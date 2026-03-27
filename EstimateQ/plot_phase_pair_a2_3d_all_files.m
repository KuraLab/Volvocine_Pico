function varargout = plot_phase_pair_a2_3d_all_files(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N, varargin)
% Overlay 3D phase-phase-a2 trajectories from all CSV files in a directory.
%
% x = phase of agent 1
% y = phase of agent 2
% z = normalized a2 of one analyzed target agent
% plotted as a simple 3D point cloud
% The function analyzes all agents listed in phase_agent_ids using the same
% Fourier fitting workflow. In addition to a2 itself, it also analyzes the
% derived quantity built from the target agent phase and normalized a2.
% The a2 normalization is done per agent using that agent's full available
% a2 series in the file: the 10th and 90th percentiles define the scale,
% their average defines the center, and their difference is normalized to 1.
% The raw a2 values are still kept separately in the output as a2_raw.
%
% Usage:
%   plot_phase_pair_a2_3d_all_files()
%   plot_phase_pair_a2_3d_all_files(dirpath)
%   plot_phase_pair_a2_3d_all_files(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices)
%   plot_phase_pair_a2_3d_all_files(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N)
%   plot_phase_pair_a2_3d_all_files(..., M, N, [m_phi2, n_phi1])
%
% Defaults:
%   dirpath = 'Spring1/255'
%   phase_agent_ids = first two non-99 agents found in the first valid CSV
%   z_agent_id = phase_agent_ids(2)  % primary target for backward-compatible top-level outputs
%   analysis_duration_sec = 15
%   analysis_start_sec = 5
%   file_indices = []  % [] means all CSV files in name-sorted order
%   M, N = 10          % default Fourier orders for automatic fitting
%   gamma_ratio = [2 1] for psi = 2*phi2 - phi1
%
% Output:
%   out: struct with plotting settings and analysis results for all target
%        agents. For backward compatibility, top-level figure/point-cloud/
%        fit fields correspond to z_agent_id, while all agent-wise results
%        are also available in out.agent_analysis.
%
% Side effect:
%   Portable gamma reconstruction data are exported to
%   fullfile(dirpath, 'gamma_exports') so that other MATLAB code can load
%   and reuse the reconstructed Gamma(psi) results.

    if nargout > 1
        error('Too many output arguments.');
    end

    % Clear any figure windows left from previous runs before starting.
    existing_figures = findall(0, 'Type', 'figure');
    if ~isempty(existing_figures)
        close(existing_figures);
        drawnow;
    end

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('Spring1', '255');
    end
    if nargin < 2
        phase_agent_ids = [];
    end
    if nargin < 3
        z_agent_id = [];
    end
    if nargin < 4 || isempty(analysis_duration_sec)
        analysis_duration_sec = 80;
    end
    if nargin < 5 || isempty(analysis_start_sec)
        analysis_start_sec = 10;
    end
    if nargin < 6
        file_indices = [];
    end
    if nargin < 7 || isempty(M)
        M = 10;
    end
    if nargin < 8 || isempty(N)
        N = 10;
    end

    % Toggle this to save only scatter/fitting figures after tuneFigure.
    enable_save_figure = false;

    gamma_ratio = [1 1];

    % Ignore legacy weighted-fit arguments if they are still passed.
    if ~isempty(varargin)
        first_extra = varargin{1};
        if isnumeric(first_extra) && numel(first_extra) == 2 && all(isfinite(first_extra(:)))
            gamma_ratio = double(first_extra(:).');
        end
    end

    % ===== Derived-signal definition: edit here =====
    derived_signal_expression = '+cos(phase_target + pi - 0.6*pi) .* a2_normalized';
    derived_signal_display_name = 'cos(phi_target + pi - 0.6*pi) * a2_norm';
    derived_signal_axis_label = 'cos(phi_target + pi - 0.6*pi) * a2_{norm}';
    derived_signal_func = @(phase_target, a2_normalized) +5*cos(phase_target + 0.4*pi) .* a2_normalized;

    sample_dt = 0.01;
    if ~isscalar(analysis_duration_sec) || analysis_duration_sec <= 0 || ~isfinite(analysis_duration_sec)
        error('analysis_duration_sec must be a positive scalar.');
    end
    if ~isscalar(analysis_start_sec) || analysis_start_sec < 0 || ~isfinite(analysis_start_sec)
        error('analysis_start_sec must be a non-negative scalar.');
    end
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'M');
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'N');
    validateattributes(gamma_ratio, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive', 'finite'}, mfilename, 'gamma_ratio');
    dirpath = resolve_analysis_dirpath(dirpath);
    if ~isfolder(dirpath)
        error('Directory not found: %s', dirpath);
    end

    csv_paths = list_csv_paths(dirpath, file_indices);
    if isempty(csv_paths)
        error('No CSV files selected in %s', dirpath);
    end

    if isempty(phase_agent_ids)
        phase_agent_ids = detect_default_phase_agents(csv_paths);
    else
        phase_agent_ids = phase_agent_ids(:).';
        if numel(phase_agent_ids) ~= 2
            error('phase_agent_ids must contain exactly two agent IDs.');
        end
        if phase_agent_ids(1) == phase_agent_ids(2)
            error('phase_agent_ids must specify two distinct agents.');
        end
    end

    if isempty(z_agent_id)
        z_agent_id = phase_agent_ids(2);
    end

    point_color = [0.0, 0.4470, 0.7410];
    marker_size = 10;
    marker_alpha = 0.18;
    gamma_settings_a2 = struct('enabled', false, 'component', 'full', 'overlay_full', false, 'show_surface_overlay', true, 'auto_save_figure', enable_save_figure);
    gamma_settings_sin_phi2_a2 = struct('enabled', true, 'component', 'full', 'overlay_full', false, 'show_surface_overlay', false, 'auto_save_figure', enable_save_figure);
    analysis_agent_ids = unique([z_agent_id, phase_agent_ids(:).'], 'stable');

    n_analysis_agents = numel(analysis_agent_ids);
    agent_analysis = struct([]);
    for agent_idx = 1:n_analysis_agents
        target_agent_id = analysis_agent_ids(agent_idx);
        agent_result = run_single_agent_analysis( ...
            csv_paths, phase_agent_ids, target_agent_id, ...
            analysis_duration_sec, analysis_start_sec, sample_dt, ...
            M, N, gamma_ratio, point_color, marker_size, marker_alpha, ...
            gamma_settings_a2, gamma_settings_sin_phi2_a2, ...
            derived_signal_func);

        if agent_idx == 1
            agent_analysis = orderfields(agent_result);
        else
            agent_analysis(agent_idx) = orderfields(agent_result, agent_analysis(1));
        end
    end

    primary_index = find([agent_analysis.agent_id] == z_agent_id, 1, 'first');
    if isempty(primary_index)
        primary_index = 1;
    end
    primary_analysis = agent_analysis(primary_index);
    true_gamma = compute_true_gamma_from_agent_analysis(agent_analysis, phase_agent_ids, gamma_ratio);
    phase_agent_omega = compute_phase_agent_mean_omega_all_files( ...
        csv_paths, phase_agent_ids, analysis_duration_sec, analysis_start_sec);

    out = struct();
    out.dirpath = dirpath;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_agent_ids = analysis_agent_ids;
    out.primary_agent_index = primary_index;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.file_indices = file_indices;
    out.fit_mode = 'unweighted';
    out.gamma_ratio = gamma_ratio;
    out.gamma_settings = struct( ...
        'a2', gamma_settings_a2, ...
        'sin_phi2_a2', gamma_settings_sin_phi2_a2);
    out.derived_signal = struct( ...
        'expression', derived_signal_expression, ...
        'display_name', derived_signal_display_name, ...
        'axis_label', derived_signal_axis_label);
    out.point_color = point_color;
    out.marker_size = marker_size;
    out.marker_alpha = marker_alpha;
    out.agent_analysis = agent_analysis;
    out.true_gamma = true_gamma;
    out.phase_agent_mean_omega = phase_agent_omega;
    out.figure_true_gamma = unwrap_scalar_field(true_gamma.figure);

    % Backward-compatible top-level fields refer to the primary z_agent_id.
    out.figure = unwrap_scalar_field(primary_analysis.figure);
    out.figure_sin_phi2_a2 = unwrap_scalar_field(primary_analysis.figure_sin_phi2_a2);
    out.used_files = unwrap_scalar_field(primary_analysis.used_files);
    out.skipped_files = unwrap_scalar_field(primary_analysis.skipped_files);
    out.per_file = unwrap_scalar_field(primary_analysis.per_file);
    out.point_cloud = unwrap_scalar_field(primary_analysis.point_cloud);
    out.fourier_fit = unwrap_scalar_field(primary_analysis.fourier_fit);
    out.fourier_fit_sin_phi2_a2 = unwrap_scalar_field(primary_analysis.fourier_fit_sin_phi2_a2);
    out.gamma_export = export_gamma_results(dirpath, out, derived_signal_func);

    if nargout == 0
        display_phase_agent_mean_omega_summary(out.phase_agent_mean_omega);
    else
        varargout{1} = out;
    end
end

function resolved_dirpath = resolve_analysis_dirpath(dirpath)
    if isfolder(dirpath)
        resolved_dirpath = dirpath;
        return;
    end

    candidate_pwd = fullfile(pwd, dirpath);
    if isfolder(candidate_pwd)
        resolved_dirpath = candidate_pwd;
        return;
    end

    base_dir = fileparts(mfilename('fullpath'));
    candidate_base = fullfile(base_dir, dirpath);
    if isfolder(candidate_base)
        resolved_dirpath = candidate_base;
        return;
    end

    resolved_dirpath = candidate_pwd;
end

function csv_paths = list_csv_paths(dirpath, file_indices)
    files = dir(fullfile(dirpath, '*.csv'));
    if isempty(files)
        csv_paths = {};
        return;
    end

    names = sort({files.name});
    if isempty(file_indices)
        selected = names;
    else
        file_indices = file_indices(:).';
        if any(file_indices < 1) || any(file_indices > numel(names))
            error('file_indices must be within 1..%d', numel(names));
        end
        selected = names(file_indices);
    end

    csv_paths = cellfun(@(name) fullfile(dirpath, name), selected, 'UniformOutput', false);
end

function phase_agent_ids = detect_default_phase_agents(csv_paths)
    for i = 1:numel(csv_paths)
        try
            T = readtable(csv_paths{i});
        catch
            continue;
        end
        if ~ismember('agent_id', T.Properties.VariableNames)
            continue;
        end
        agents = unique(T.agent_id, 'sorted').';
        agents = agents(agents ~= 99);
        if numel(agents) >= 2
            phase_agent_ids = agents(1:2);
            return;
        end
    end

    error('Could not determine default phase_agent_ids from the selected CSV files.');
end

function agent_out = run_single_agent_analysis(csv_paths, phase_agent_ids, target_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt, M, N, gamma_ratio, point_color, marker_size, marker_alpha, gamma_settings_a2, gamma_settings_sin_phi2_a2, derived_signal_func)
    fig = figure('Color', 'w', 'Name', sprintf('a2 agent %d', target_agent_id));
    ax = axes('Parent', fig);
    hold(ax, 'on');

    fig_sin_phi2_a2 = [];

    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});
    per_file = struct('file_path', {}, 'window_start_abs', {}, 'window_end_abs', {}, ...
        'n_points', {}, 'phase_agent_ids', {}, 'z_agent_id', {});
    phi1_all = [];
    phi2_all = [];
    a2_raw_all = [];
    a2_all = [];
    a2_normalized_all = [];
    sin_phi2_a2_all = [];
    time_all = [];
    file_id_all = [];
    save_after_tune = isfield(gamma_settings_a2, 'auto_save_figure') && logical(gamma_settings_a2.auto_save_figure);
    phi_idx = find(phase_agent_ids == target_agent_id, 1, 'first');
    if isempty(phi_idx)
        s_label = ['$$s_{' num2str(target_agent_id) '}(\phi_1,\phi_2)$$'];
    else
        s_label = ['$$s_' num2str(phi_idx) '(\phi_1,\phi_2)$$'];
    end

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            [point_data, meta] = compute_points_for_csv( ...
                csv_path, phase_agent_ids, target_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            warning('Skipping %s for agent %d: %s', csv_path, target_agent_id, ME.message);
            skipped_files(end+1) = struct('file_path', csv_path, 'reason', ME.message); %#ok<AGROW>
            continue;
        end

        if isempty(point_data.time)
            skipped_files(end+1) = struct('file_path', csv_path, 'reason', 'No valid points after processing.'); %#ok<AGROW>
            continue;
        end

        scatter3(ax, point_data.phase1, point_data.phase2, point_data.a2, marker_size, ...
            'filled', 'MarkerFaceColor', point_color, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha);

        phase_for_derived_signal = select_phase_for_target_agent(point_data, phase_agent_ids, target_agent_id);
        sin_phi2_a2 = derived_signal_func(phase_for_derived_signal, point_data.a2_normalized);

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        a2_raw_all = [a2_raw_all; point_data.a2_raw(:)]; %#ok<AGROW>
        a2_all = [a2_all; point_data.a2(:)]; %#ok<AGROW>
        a2_normalized_all = [a2_normalized_all; point_data.a2_normalized(:)]; %#ok<AGROW>
        sin_phi2_a2_all = [sin_phi2_a2_all; sin_phi2_a2(:)]; %#ok<AGROW>
        time_all = [time_all; point_data.time(:)]; %#ok<AGROW>
        file_id_all = [file_id_all; i * ones(numel(point_data.time), 1)]; %#ok<AGROW>

        used_files{end+1} = csv_path; %#ok<AGROW>
        per_file(end+1) = meta; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available to overlay for agent %d.', target_agent_id);
    end

    configure_phase_pair_axes(ax, s_label);

    figure(fig);
    tuneFigure;
    if save_after_tune
        saveFigure;
    end

    point_cloud = struct( ...
        'phi1', phi1_all, ...
        'phi2', phi2_all, ...
        'a2_raw', a2_raw_all, ...
        'a2', a2_all, ...
        'a2_normalized', a2_normalized_all, ...
        'sin_phi2_a2', sin_phi2_a2_all, ...
        'time', time_all, ...
        'file_index', file_id_all);

    fourier_fit = fitDoubleFourierScatter( ...
        point_cloud.phi1, point_cloud.phi2, point_cloud.a2, M, N, ...
        s_label, 'full', gamma_ratio, gamma_settings_a2);
    fourier_fit.M = M;
    fourier_fit.N = N;
    fourier_fit.fit_mode = 'unweighted';
    fourier_fit.gamma_ratio = gamma_ratio;

    fourier_fit_sin_phi2_a2 = fitDoubleFourierScatter( ...
        point_cloud.phi1, point_cloud.phi2, point_cloud.sin_phi2_a2, M, N, ...
        s_label, 'mixed-only', gamma_ratio, gamma_settings_sin_phi2_a2);
    fourier_fit_sin_phi2_a2.M = M;
    fourier_fit_sin_phi2_a2.N = N;
    fourier_fit_sin_phi2_a2.fit_mode = 'unweighted';
    fourier_fit_sin_phi2_a2.gamma_ratio = gamma_ratio;

    agent_out = struct();
    agent_out.agent_id = target_agent_id;
    agent_out.fit_mode = 'unweighted';
    agent_out.gamma_ratio = gamma_ratio;
    agent_out.gamma_settings = struct( ...
        'a2', gamma_settings_a2, ...
        'sin_phi2_a2', gamma_settings_sin_phi2_a2);
    agent_out.figure = fig;
    agent_out.figure_sin_phi2_a2 = fig_sin_phi2_a2;
    agent_out.used_files = used_files;
    agent_out.skipped_files = skipped_files;
    agent_out.per_file = per_file;
    agent_out.point_color = point_color;
    agent_out.marker_size = marker_size;
    agent_out.marker_alpha = marker_alpha;
    agent_out.point_cloud = point_cloud;
    agent_out.fourier_fit = fourier_fit;
    agent_out.fourier_fit_sin_phi2_a2 = fourier_fit_sin_phi2_a2;
end

function phase_target = select_phase_for_target_agent(point_data, phase_agent_ids, target_agent_id)
    if target_agent_id == phase_agent_ids(1)
        phase_target = point_data.phase1;
    elseif target_agent_id == phase_agent_ids(2)
        phase_target = point_data.phase2;
    else
        phase_target = point_data.phase2;
    end
end

function true_gamma = compute_true_gamma_from_agent_analysis(agent_analysis, phase_agent_ids, gamma_ratio)
    true_gamma = struct( ...
        'available', false, ...
        'reason', '', ...
        'agent_id_1', [], ...
        'agent_id_2', [], ...
        'gamma_ratio', gamma_ratio, ...
        'component', '', ...
        'psi_grid', [], ...
        'gamma_agent_1', [], ...
        'gamma_agent_2', [], ...
        'gamma_true', [], ...
        'figure', []);

    if numel(phase_agent_ids) < 2 || numel(agent_analysis) < 2
        true_gamma.reason = 'Two analyzed phase agents are required.';
        return;
    end

    idx_agent_1 = find([agent_analysis.agent_id] == phase_agent_ids(1), 1, 'first');
    idx_agent_2 = find([agent_analysis.agent_id] == phase_agent_ids(2), 1, 'first');
    if isempty(idx_agent_1) || isempty(idx_agent_2)
        true_gamma.reason = 'Could not find both phase agents in agent_analysis.';
        return;
    end

    fit_agent_1 = unwrap_scalar_field(agent_analysis(idx_agent_1).fourier_fit_sin_phi2_a2);
    fit_agent_2 = unwrap_scalar_field(agent_analysis(idx_agent_2).fourier_fit_sin_phi2_a2);
    if ~isstruct(fit_agent_1) || ~isfield(fit_agent_1, 'gamma_resonance') || ...
            ~isstruct(fit_agent_2) || ~isfield(fit_agent_2, 'gamma_resonance')
        true_gamma.reason = 'Gamma reconstruction results were not found for both agents.';
        return;
    end

    gamma_1 = fit_agent_1.gamma_resonance;
    gamma_2 = fit_agent_2.gamma_resonance;
    [psi_grid_1, gamma_values_1] = extract_gamma_curve(gamma_1, 'full');
    [psi_grid_2, gamma_values_2] = extract_gamma_curve(gamma_2, 'full');
    if isempty(psi_grid_1) || isempty(psi_grid_2) || isempty(gamma_values_1) || isempty(gamma_values_2)
        true_gamma.reason = 'Full gamma curves were empty for at least one agent.';
        return;
    end

    n_grid = max(numel(psi_grid_1), numel(psi_grid_2));
    psi_grid = linspace(-pi, pi, n_grid).';
    gamma_agent_1 = interp1(psi_grid_1, gamma_values_1, psi_grid, 'linear', 'extrap');
    gamma_agent_2 = interp1(psi_grid_2, gamma_values_2, psi_grid, 'linear', 'extrap');
    gamma_true_values = gamma_agent_2 - gamma_agent_1;

    fig_true_gamma = figure('Color', 'w', 'Name', 'true gamma');
    ax_true_gamma = axes('Parent', fig_true_gamma);
    hold(ax_true_gamma, 'on');
    plot(ax_true_gamma, psi_grid, gamma_true_values, 'LineWidth', 1.8, 'Color', [0.8500, 0.3250, 0.0980]);
    plot(ax_true_gamma, psi_grid, zeros(size(psi_grid)), ':', 'LineWidth', 1.0, 'Color', [0.5, 0.5, 0.5]);
    xlabel(ax_true_gamma, '$$\psi$$', 'Interpreter', 'latex');
    ylabel(ax_true_gamma, '$$\Gamma_{\mathrm{true}}(\psi)$$', 'Interpreter', 'latex');
    grid(ax_true_gamma, 'on');
    box(ax_true_gamma, 'on');
    xlim(ax_true_gamma, [-pi, pi]);
    xticks(ax_true_gamma, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax_true_gamma, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
    ax_true_gamma.XLabel.Interpreter = 'latex';
    ax_true_gamma.YLabel.Interpreter = 'latex';
    ax_true_gamma.TickLabelInterpreter = 'latex';
    figure(fig_true_gamma);
    tuneFigure;

    true_gamma.available = true;
    true_gamma.agent_id_1 = phase_agent_ids(1);
    true_gamma.agent_id_2 = phase_agent_ids(2);
    true_gamma.component = 'full';
    true_gamma.psi_grid = psi_grid;
    true_gamma.gamma_agent_1 = gamma_agent_1;
    true_gamma.gamma_agent_2 = gamma_agent_2;
    true_gamma.gamma_true = gamma_true_values;
    true_gamma.figure = fig_true_gamma;
end

function [psi_grid, gamma_values] = extract_gamma_curve(gamma_resonance, mode)
    psi_grid = [];
    gamma_values = [];
    if nargin < 2 || isempty(mode)
        mode = 'selected';
    end
    if ~isstruct(gamma_resonance)
        return;
    end

    if strcmp(mode, 'full')
        if isfield(gamma_resonance, 'psi_grid_centered') && isfield(gamma_resonance, 'gamma_values_full_centered') && ...
                ~isempty(gamma_resonance.psi_grid_centered) && ~isempty(gamma_resonance.gamma_values_full_centered)
            psi_grid = gamma_resonance.psi_grid_centered(:);
            gamma_values = gamma_resonance.gamma_values_full_centered(:);
            return;
        end

        if isfield(gamma_resonance, 'psi_grid') && isfield(gamma_resonance, 'gamma_values_full') && ...
                ~isempty(gamma_resonance.psi_grid) && ~isempty(gamma_resonance.gamma_values_full)
            psi_grid = gamma_resonance.psi_grid(:);
            gamma_values = gamma_resonance.gamma_values_full(:);
            return;
        end
    end

    if isfield(gamma_resonance, 'psi_grid_centered') && isfield(gamma_resonance, 'gamma_values_centered') && ...
            ~isempty(gamma_resonance.psi_grid_centered) && ~isempty(gamma_resonance.gamma_values_centered)
        psi_grid = gamma_resonance.psi_grid_centered(:);
        gamma_values = gamma_resonance.gamma_values_centered(:);
        return;
    end

    if isfield(gamma_resonance, 'psi_grid') && isfield(gamma_resonance, 'gamma_values') && ...
            ~isempty(gamma_resonance.psi_grid) && ~isempty(gamma_resonance.gamma_values)
        psi_grid = gamma_resonance.psi_grid(:);
        gamma_values = gamma_resonance.gamma_values(:);
    end
end

function display_phase_agent_mean_omega_summary(omega_summary)
    if ~isstruct(omega_summary) || ~isfield(omega_summary, 'available') || ~omega_summary.available
        if isstruct(omega_summary) && isfield(omega_summary, 'reason') && ~isempty(omega_summary.reason)
            fprintf('[INFO] Phase-agent mean angular velocity unavailable: %s\n', omega_summary.reason);
        else
            fprintf('[INFO] Phase-agent mean angular velocity unavailable.\n');
        end
        return;
    end

    agent_ids = omega_summary.phase_agent_ids(:).';
    omega_rad = omega_summary.mean_omega_rad_s(:).';
    omega_deg = omega_summary.mean_omega_deg_s(:).';
    std_rad = omega_summary.std_omega_rad_s(:).';
    std_deg = omega_summary.std_omega_deg_s(:).';
    n_valid_files = numel(omega_summary.per_file);
    n_skipped_files = numel(omega_summary.skipped_files);

    fprintf('[INFO] Phase-agent mean angular velocity summary (%d files', n_valid_files);
    if n_skipped_files > 0
        fprintf(', %d skipped', n_skipped_files);
    end
    fprintf(')\n');

    for agent_idx = 1:numel(agent_ids)
        fprintf('  Agent %d: %.12f rad/s (%.12f deg/s)', ...
            agent_ids(agent_idx), omega_rad(agent_idx), omega_deg(agent_idx));
        if agent_idx <= numel(std_rad) && isfinite(std_rad(agent_idx))
            fprintf(', std %.12f rad/s (%.12f deg/s)', std_rad(agent_idx), std_deg(agent_idx));
        end
        fprintf('\n');
    end

    if numel(agent_ids) >= 2 && isfield(omega_summary, 'mean_ratio_2_over_1') && isfinite(omega_summary.mean_ratio_2_over_1)
        fprintf('  Mean ratio omega(%d)/omega(%d): %.12f\n', ...
            agent_ids(2), agent_ids(1), omega_summary.mean_ratio_2_over_1);
    end

    if numel(agent_ids) >= 2 && isfield(omega_summary, 'median_ratio_2_over_1') && isfinite(omega_summary.median_ratio_2_over_1)
        fprintf('  Median ratio omega(%d)/omega(%d): %.12f\n', ...
            agent_ids(2), agent_ids(1), omega_summary.median_ratio_2_over_1);
    end
end

function value = unwrap_scalar_field(value)
    if iscell(value)
        if isempty(value)
            value = [];
        else
            value = value{1};
        end
    end
end

function configure_phase_pair_axes(ax, z_axis_label)
    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlabel(ax, '$$\phi_1$$', 'Interpreter', 'latex');
    ylabel(ax, '$$\phi_2$$', 'Interpreter', 'latex');
    zlabel(ax, z_axis_label, 'Interpreter', 'latex');
    grid(ax, 'on');
    view(ax, 3);
    box(ax, 'on');
end

function [point_data, meta] = compute_points_for_csv(csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt)
    requested_agents = unique([phase_agent_ids(:); z_agent_id]);
    series_by_agent = load_corrected_agent_series_from_csv( ...
        csv_path, requested_agents, {'time_pc_sec_abs', 'a0', 'a2'});
    overlap_start = -inf;
    overlap_end = inf;

    for k = 1:numel(requested_agents)
        aid = requested_agents(k);
        series = series_by_agent(aid);
        if isempty(series.time)
            error('No valid samples found for agent %d.', aid);
        end

        overlap_start = max(overlap_start, min(series.time));
        overlap_end = min(overlap_end, max(series.time));
    end

    if overlap_end <= overlap_start
        error('Selected agents do not share an overlapping time range.');
    end

    window_start_abs = overlap_start + analysis_start_sec;
    window_end_abs = min(overlap_end, window_start_abs + analysis_duration_sec);
    if window_end_abs <= window_start_abs
        error('Requested analysis window %.3f-%.3f s is outside the common overlap.', ...
            analysis_start_sec, analysis_start_sec + analysis_duration_sec);
    end

    time_abs = (window_start_abs:sample_dt:window_end_abs).';
    if numel(time_abs) < 2
        error('Analysis window is too short after applying the selected settings.');
    end
    time_rel = time_abs - window_start_abs;

    a0_1 = interp1(series_by_agent(phase_agent_ids(1)).time, series_by_agent(phase_agent_ids(1)).a0_corr, time_abs, 'linear', NaN);
    a0_2 = interp1(series_by_agent(phase_agent_ids(2)).time, series_by_agent(phase_agent_ids(2)).a0_corr, time_abs, 'linear', NaN);
    a2_z = interp1(series_by_agent(z_agent_id).time, series_by_agent(z_agent_id).a2, time_abs, 'linear', NaN);

    valid = isfinite(a0_1) & isfinite(a0_2) & isfinite(a2_z);
    time_rel = time_rel(valid);
    a0_1 = a0_1(valid);
    a0_2 = a0_2(valid);
    a2_z = a2_z(valid);

    if numel(time_rel) < 2
        error('No overlapping interpolated samples were found in the selected analysis window.');
    end

    point_data = struct();
    point_data.time = time_rel;
    point_data.phase1 = mod(a0_1, 256) * (2*pi/256);
    point_data.phase2 = mod(a0_2, 256) * (2*pi/256);
    point_data.a2_raw = a2_z;
    point_data.a2_normalized = clip_values( ...
        normalize_by_agent_percentile_span(a2_z, series_by_agent(z_agent_id).a2, 10), -0.5, 0.5);
    point_data.a2 = point_data.a2_normalized;

    meta = struct();
    meta.file_path = csv_path;
    meta.window_start_abs = window_start_abs;
    meta.window_end_abs = window_end_abs;
    meta.n_points = numel(point_data.time);
    meta.phase_agent_ids = phase_agent_ids;
    meta.z_agent_id = z_agent_id;
end

function xNorm = normalize_by_agent_percentile_span(x, reference_x, tailPct)
    x = double(x(:));
    reference_x = double(reference_x(:));
    valid = isfinite(x);
    reference_valid = isfinite(reference_x);
    xNorm = nan(size(x));

    if nnz(valid) < 1
        return;
    end

    if nnz(reference_valid) < 5
        xNorm(valid) = x(valid);
        return;
    end

    reference_values = reference_x(reference_valid);
    low_value = prctile(reference_values, tailPct);
    high_value = prctile(reference_values, 100 - tailPct);
    center_value = 0.5 * (low_value + high_value);
    span_value = high_value - low_value;

    if ~isfinite(span_value) || span_value <= 0
        xNorm(valid) = x(valid) - center_value;
        return;
    end

    xNorm(valid) = (x(valid) - center_value) / span_value;
end

function xClipped = clip_values(x, lower_bound, upper_bound)
    xClipped = min(max(x, lower_bound), upper_bound);
end

function export_info = export_gamma_results(dirpath, out, phase_sensitivity_func)
    export_info = struct( ...
        'available', false, ...
        'reason', '', ...
        'directory', '', ...
        'mat_file', '', ...
        'latest_mat_file', '', ...
        'curve_files', struct('label', {}, 'file_path', {}), ...
        'prc_snippet_files', struct('label', {}, 'file_path', {}));

    try
        export_dir = fullfile(dirpath, 'gamma_exports');
        if ~exist(export_dir, 'dir')
            mkdir(export_dir);
        end

        gamma_export = build_gamma_export_bundle(dirpath, out);
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        base_name = sprintf('gamma_export_phase%d_%d_z%d_%s', ...
            out.phase_agent_ids(1), out.phase_agent_ids(2), out.z_agent_id, timestamp);
        mat_file = fullfile(export_dir, [base_name '.mat']);
        latest_mat_file = fullfile(export_dir, 'gamma_export_latest.mat');

        save(mat_file, 'gamma_export');
        save(latest_mat_file, 'gamma_export');

        curve_files = write_gamma_curve_exports(export_dir, gamma_export);
        prc_snippet_files = write_prc_snippet_exports(export_dir, dirpath, phase_sensitivity_func, out.derived_signal.display_name, gamma_export, 10);

        export_info.available = true;
        export_info.directory = export_dir;
        export_info.mat_file = mat_file;
        export_info.latest_mat_file = latest_mat_file;
        export_info.curve_files = curve_files;
        export_info.prc_snippet_files = prc_snippet_files;

        fprintf('[INFO] Exported gamma bundle to %s\n', latest_mat_file);
    catch ME
        export_info.reason = ME.message;
        warning('Gamma export failed for %s: %s', dirpath, ME.message);
    end
end

function snippet_files = write_prc_snippet_exports(export_dir, dirpath, ~, ~, ~, max_harmonics)
    if nargin < 6 || isempty(max_harmonics)
        max_harmonics = 10;
    end

    snippet_files = struct('label', {}, 'file_path', {});

    % Reference snippets matching plot_psi_with_desined_Z first overlay source.
    delta = 0.4 * pi;
    cos_ref_file = fullfile(export_dir, 'prc_snippet_ref_cos.txt');
    if write_direct_coeff_prc_snippet(cos_ref_file, cos(delta), -sin(delta), 'cos(phi + 0.4*pi) reference', max_harmonics)
        snippet_files(end + 1) = struct('label', 'ref_cos', 'file_path', cos_ref_file); %#ok<AGROW>
    end

    w1_model = load_exported_w_model_for_snippet(dirpath, 'w1');
    w1_ref_file = fullfile(export_dir, 'prc_snippet_ref_w1.txt');
    if write_model_prc_snippet(w1_ref_file, w1_model, 'W1 exported reference', max_harmonics)
        snippet_files(end + 1) = struct('label', 'ref_w1', 'file_path', w1_ref_file); %#ok<AGROW>
    end

    w2_model = load_exported_w_model_for_snippet(dirpath, 'w2');
    w2_ref_file = fullfile(export_dir, 'prc_snippet_ref_w2.txt');
    if write_model_prc_snippet(w2_ref_file, w2_model, 'W2 exported reference', max_harmonics)
        snippet_files(end + 1) = struct('label', 'ref_w2', 'file_path', w2_ref_file); %#ok<AGROW>
    end

end

function did_write = write_direct_coeff_prc_snippet(file_path, a1, b1, source_label, max_harmonics)
    did_write = false;
    n_h = max(1, max_harmonics);
    prc_a = zeros(n_h + 1, 1);
    prc_b = zeros(n_h + 1, 1);
    prc_a(2) = a1;
    prc_b(2) = b1;

    snippet_text = build_prc_python_snippet(prc_a, prc_b, source_label);
    fid = fopen(file_path, 'w');
    if fid < 0
        return;
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', snippet_text);
    did_write = true;
end

function did_write = write_model_prc_snippet(file_path, model, source_label, max_harmonics)
    did_write = false;
    if ~isstruct(model) || ~isfield(model, 'available') || ~model.available
        return;
    end

    n_h = min(max_harmonics, model.order);
    if n_h < 1
        return;
    end

    prc_a = zeros(n_h + 1, 1);
    prc_b = zeros(n_h + 1, 1);
    prc_a(1) = model.a0;
    prc_a(2:end) = model.a(1:n_h);
    prc_b(2:end) = model.b(1:n_h);

    snippet_text = build_prc_python_snippet(prc_a, prc_b, source_label);
    fid = fopen(file_path, 'w');
    if fid < 0
        return;
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', snippet_text);
    did_write = true;
end

function model = load_exported_w_model_for_snippet(dirpath, mode_name)
    model = struct('available', false, 'source_file', '', 'order', 0, 'a0', 0, 'a', [], 'b', []);

    search_dirs = {dirpath, fullfile(dirpath, 'gamma_exports')};

    mat_path = find_latest_matching_file(search_dirs, 'W_fourier_fit_params_*.mat');
    if ~isempty(mat_path)
        model = load_w_model_from_mat_for_snippet(mat_path, mode_name);
        if model.available
            return;
        end
    end

    if strcmpi(mode_name, 'w1')
        csv_pattern = 'W1_fourier_fit_params_agent*_*.csv';
    else
        csv_pattern = 'W2_fourier_fit_params_agent*_*.csv';
    end
    csv_path = find_latest_matching_file(search_dirs, csv_pattern);
    if ~isempty(csv_path)
        model = load_w_model_from_csv_for_snippet(csv_path);
    end
end

function latest_path = find_latest_matching_file(search_dirs, pattern)
    latest_path = '';
    latest_time = -inf;
    for d = 1:numel(search_dirs)
        dir_path = search_dirs{d};
        if ~isfolder(dir_path)
            continue;
        end
        files = dir(fullfile(dir_path, pattern));
        for i = 1:numel(files)
            if files(i).datenum > latest_time
                latest_time = files(i).datenum;
                latest_path = fullfile(files(i).folder, files(i).name);
            end
        end
    end
end

function model = load_w_model_from_mat_for_snippet(mat_path, mode_name)
    model = struct('available', false, 'source_file', mat_path, 'order', 0, 'a0', 0, 'a', [], 'b', []);
    try
        S = load(mat_path);
    catch
        return;
    end
    if ~isfield(S, 'fourier_fit_export') || ~isstruct(S.fourier_fit_export)
        return;
    end

    if strcmpi(mode_name, 'w1')
        field_name = 'fit_w1';
    else
        field_name = 'fit_w2';
    end
    if ~isfield(S.fourier_fit_export, field_name) || ~isstruct(S.fourier_fit_export.(field_name))
        return;
    end

    fit_w = S.fourier_fit_export.(field_name);
    if ~all(isfield(fit_w, {'order', 'a0', 'a', 'b'}))
        return;
    end

    order = max(0, round(double(fit_w.order)));
    a = reshape(double(fit_w.a), [], 1);
    b = reshape(double(fit_w.b), [], 1);
    if numel(a) < order || numel(b) < order
        return;
    end

    model.available = true;
    model.order = order;
    model.a0 = double(fit_w.a0);
    model.a = a(1:order);
    model.b = b(1:order);
end

function model = load_w_model_from_csv_for_snippet(csv_path)
    model = struct('available', false, 'source_file', csv_path, 'order', 0, 'a0', 0, 'a', [], 'b', []);
    try
        T = readtable(csv_path);
    catch
        return;
    end
    if ~all(ismember({'k', 'a_k', 'b_k'}, T.Properties.VariableNames))
        return;
    end

    k = double(T.k(:));
    a_k = double(T.a_k(:));
    b_k = double(T.b_k(:));
    valid = isfinite(k) & isfinite(a_k) & isfinite(b_k) & (k >= 0) & (abs(k - round(k)) < eps(10));
    k = round(k(valid));
    a_k = a_k(valid);
    b_k = b_k(valid);
    if isempty(k)
        return;
    end

    order = max(k);
    a = zeros(order, 1);
    b = zeros(order, 1);
    idx0 = find(k == 0, 1, 'first');
    if isempty(idx0)
        return;
    end
    a0 = a_k(idx0);

    for n = 1:order
        idxn = find(k == n, 1, 'first');
        if ~isempty(idxn)
            a(n) = a_k(idxn);
            b(n) = b_k(idxn);
        end
    end

    model.available = true;
    model.order = order;
    model.a0 = a0;
    model.a = a;
    model.b = b;
end

function snippet_text = build_prc_python_snippet(prc_a, prc_b, source_label)
    prc_harmonics = numel(prc_a) - 1;
    lines = cell(0, 1);
    lines{end + 1} = sprintf('# Auto-generated from %s', source_label);
    lines{end + 1} = sprintf('prc_harmonics = %d', prc_harmonics);
    lines{end + 1} = 'prc_a = [0.0] * (prc_harmonics + 1)';
    lines{end + 1} = 'prc_b = [0.0] * (prc_harmonics + 1)';
    lines{end + 1} = '';
    lines{end + 1} = sprintf('prc_a[0] = %.10f', prc_a(1));
    lines{end + 1} = sprintf('prc_b[0] = %.10f', prc_b(1));
    for n = 1:prc_harmonics
        lines{end + 1} = sprintf('prc_a[%d] = %.10f', n, prc_a(n + 1));
        lines{end + 1} = sprintf('prc_b[%d] = %.10f', n, prc_b(n + 1));
    end

    snippet_text = strjoin(lines, newline);
end

function gamma_export = build_gamma_export_bundle(dirpath, out)
    gamma_export = struct();
    gamma_export.schema_version = 1;
    gamma_export.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    gamma_export.source_dirpath = dirpath;
    gamma_export.phase_agent_ids = out.phase_agent_ids;
    gamma_export.z_agent_id = out.z_agent_id;
    gamma_export.analysis_agent_ids = out.analysis_agent_ids;
    gamma_export.analysis_duration_sec = out.analysis_duration_sec;
    gamma_export.analysis_start_sec = out.analysis_start_sec;
    gamma_export.file_indices = out.file_indices;
    gamma_export.fit_mode = out.fit_mode;
    gamma_export.gamma_ratio = out.gamma_ratio;
    gamma_export.gamma_settings = out.gamma_settings;
    gamma_export.derived_signal = out.derived_signal;
    gamma_export.phase_agent_mean_omega = out.phase_agent_mean_omega;
    gamma_export.agents = build_agent_gamma_exports(out.agent_analysis);
    gamma_export.true_gamma = sanitize_true_gamma_export(out.true_gamma);
end

function agent_exports = build_agent_gamma_exports(agent_analysis)
    if isempty(agent_analysis)
        agent_exports = struct([]);
        return;
    end

    agent_exports = repmat(struct( ...
        'agent_id', [], ...
        'fit_mode', '', ...
        'gamma_ratio', [], ...
        'used_files', {{}}, ...
        'skipped_files', struct('file_path', {}, 'reason', {}), ...
        'per_file', struct('file_path', {}, 'window_start_abs', {}, 'window_end_abs', {}, 'n_points', {}, 'phase_agent_ids', {}, 'z_agent_id', {}), ...
        'n_points', 0, ...
        'a2_gamma', struct(), ...
        'derived_gamma', struct()), size(agent_analysis));

    for idx = 1:numel(agent_analysis)
        agent_exports(idx) = build_single_agent_gamma_export(agent_analysis(idx));
    end
end

function agent_export = build_single_agent_gamma_export(agent_out)
    point_cloud = unwrap_scalar_field(agent_out.point_cloud);
    n_points = 0;
    if isstruct(point_cloud) && isfield(point_cloud, 'time') && ~isempty(point_cloud.time)
        n_points = numel(point_cloud.time);
    end

    agent_export = struct();
    agent_export.agent_id = agent_out.agent_id;
    agent_export.fit_mode = agent_out.fit_mode;
    agent_export.gamma_ratio = agent_out.gamma_ratio;
    agent_export.used_files = agent_out.used_files;
    agent_export.skipped_files = agent_out.skipped_files;
    agent_export.per_file = agent_out.per_file;
    agent_export.n_points = n_points;
    agent_export.a2_gamma = build_fit_gamma_export(unwrap_scalar_field(agent_out.fourier_fit), 'a2');
    agent_export.derived_gamma = build_fit_gamma_export(unwrap_scalar_field(agent_out.fourier_fit_sin_phi2_a2), 'derived');
end

function fit_export = build_fit_gamma_export(fit_result, signal_role)
    fit_export = struct( ...
        'available', false, ...
        'signal_role', signal_role, ...
        'target_name', '', ...
        'fit_mode', '', ...
        'heatmap_mode', '', ...
        'rmse', [], ...
        'M', [], ...
        'N', [], ...
        'z_mean', [], ...
        'coeff', [], ...
        'basis_names', {{}}, ...
        'basis_groups', {{}}, ...
        'basis_types', {{}}, ...
        'basis_phi1_order', [], ...
        'basis_phi2_order', [], ...
        'basis_table', table(), ...
        'gamma_ratio', [], ...
        'gamma_ratio_reduced', [], ...
        'gamma_settings', struct(), ...
        'gamma_resonance', sanitize_gamma_resonance(struct()));

    if ~isstruct(fit_result) || isempty(fieldnames(fit_result))
        return;
    end

    fit_export.available = true;
    if isfield(fit_result, 'target_name')
        fit_export.target_name = fit_result.target_name;
    end
    if isfield(fit_result, 'fit_mode')
        fit_export.fit_mode = fit_result.fit_mode;
    end
    if isfield(fit_result, 'heatmap_mode')
        fit_export.heatmap_mode = fit_result.heatmap_mode;
    end
    if isfield(fit_result, 'rmse')
        fit_export.rmse = fit_result.rmse;
    end
    if isfield(fit_result, 'M')
        fit_export.M = fit_result.M;
    end
    if isfield(fit_result, 'N')
        fit_export.N = fit_result.N;
    end
    if isfield(fit_result, 'z_mean')
        fit_export.z_mean = fit_result.z_mean;
    end
    if isfield(fit_result, 'coeff')
        fit_export.coeff = fit_result.coeff;
    end
    if isfield(fit_result, 'basis_names')
        fit_export.basis_names = fit_result.basis_names;
    end
    if isfield(fit_result, 'basis_groups')
        fit_export.basis_groups = fit_result.basis_groups;
    end
    if isfield(fit_result, 'basis_types')
        fit_export.basis_types = fit_result.basis_types;
    end
    if isfield(fit_result, 'basis_m')
        fit_export.basis_phi1_order = fit_result.basis_m;
    end
    if isfield(fit_result, 'basis_n')
        fit_export.basis_phi2_order = fit_result.basis_n;
    end
    if isfield(fit_result, 'contribution_table') && istable(fit_result.contribution_table)
        basis_table = fit_result.contribution_table;
        keep_columns = {'term_index', 'basis_name', 'basis_group', 'phi1_order', 'phi2_order', 'coefficient'};
        keep_columns = keep_columns(ismember(keep_columns, basis_table.Properties.VariableNames));
        if ~isempty(keep_columns)
            basis_table = basis_table(:, keep_columns);
        end
        fit_export.basis_table = basis_table;

        if isempty(fit_export.basis_phi1_order) && ismember('phi1_order', basis_table.Properties.VariableNames)
            fit_export.basis_phi1_order = basis_table.phi1_order;
        end
        if isempty(fit_export.basis_phi2_order) && ismember('phi2_order', basis_table.Properties.VariableNames)
            fit_export.basis_phi2_order = basis_table.phi2_order;
        end
    end
    if isfield(fit_result, 'gamma_ratio')
        fit_export.gamma_ratio = fit_result.gamma_ratio;
    end
    if isfield(fit_result, 'gamma_ratio_reduced')
        fit_export.gamma_ratio_reduced = fit_result.gamma_ratio_reduced;
    end
    if isfield(fit_result, 'gamma_settings')
        fit_export.gamma_settings = fit_result.gamma_settings;
    end
    if isfield(fit_result, 'gamma_resonance')
        fit_export.gamma_resonance = sanitize_gamma_resonance(fit_result.gamma_resonance);
    end
end

function gamma_resonance_export = sanitize_gamma_resonance(gamma_resonance)
    gamma_resonance_export = struct( ...
        'enabled', false, ...
        'component', '', ...
        'overlay_full', false, ...
        'phi1_base', [], ...
        'phi2_base', [], ...
        'psi_label', '', ...
        'harmonic_index', [], ...
        'gamma_cos', [], ...
        'gamma_sin', [], ...
        'psi_grid', [], ...
        'gamma_values', [], ...
        'psi_grid_centered', [], ...
        'gamma_values_centered', [], ...
        'gamma_values_full', [], ...
        'gamma_values_full_centered', [], ...
        'gamma_values_symmetric', [], ...
        'gamma_values_symmetric_centered', [], ...
        'gamma_values_antisymmetric', [], ...
        'gamma_values_antisymmetric_centered', []);

    if ~isstruct(gamma_resonance) || isempty(fieldnames(gamma_resonance))
        return;
    end

    gamma_resonance_export.enabled = true;
    if isfield(gamma_resonance, 'enabled')
        gamma_resonance_export.enabled = logical(gamma_resonance.enabled);
    end

    field_names = fieldnames(gamma_resonance_export);
    for idx = 1:numel(field_names)
        field_name = field_names{idx};
        if strcmp(field_name, 'enabled')
            continue;
        end
        if isfield(gamma_resonance, field_name)
            gamma_resonance_export.(field_name) = gamma_resonance.(field_name);
        end
    end
end

function true_gamma_export = sanitize_true_gamma_export(true_gamma)
    true_gamma_export = struct( ...
        'available', false, ...
        'reason', '', ...
        'agent_id_1', [], ...
        'agent_id_2', [], ...
        'gamma_ratio', [], ...
        'component', '', ...
        'psi_grid', [], ...
        'gamma_agent_1', [], ...
        'gamma_agent_2', [], ...
        'gamma_true', []);

    if ~isstruct(true_gamma) || isempty(fieldnames(true_gamma))
        return;
    end

    field_names = fieldnames(true_gamma_export);
    for idx = 1:numel(field_names)
        field_name = field_names{idx};
        if isfield(true_gamma, field_name)
            true_gamma_export.(field_name) = true_gamma.(field_name);
        end
    end
end

function curve_files = write_gamma_curve_exports(export_dir, gamma_export)
    curve_files = struct('label', {}, 'file_path', {});

    for idx = 1:numel(gamma_export.agents)
        agent_export = gamma_export.agents(idx);

        a2_file = fullfile(export_dir, sprintf('gamma_curve_agent%d_a2.csv', agent_export.agent_id));
        if write_single_gamma_curve_csv(a2_file, agent_export.a2_gamma.gamma_resonance)
            curve_files(end + 1) = struct( ...
                'label', sprintf('agent_%d_a2', agent_export.agent_id), ...
                'file_path', a2_file); %#ok<AGROW>
        end
    end

    if isfield(gamma_export, 'true_gamma') && isstruct(gamma_export.true_gamma) && ...
            isfield(gamma_export.true_gamma, 'available') && gamma_export.true_gamma.available
        true_gamma_file = fullfile(export_dir, sprintf('gamma_true_agent%d_%d.csv', ...
            gamma_export.true_gamma.agent_id_1, gamma_export.true_gamma.agent_id_2));
        if write_single_gamma_curve_csv(true_gamma_file, gamma_export.true_gamma)
            curve_files(end + 1) = struct( ...
                'label', 'true_gamma', ...
                'file_path', true_gamma_file); %#ok<AGROW>
        end
    end
end

function did_write = write_single_gamma_curve_csv(file_path, gamma_definition)
    did_write = false;
    curve_table = build_gamma_curve_table(gamma_definition);
    if isempty(curve_table) || height(curve_table) < 1
        return;
    end

    writetable(curve_table, file_path);
    did_write = true;
end

function curve_table = build_gamma_curve_table(gamma_definition)
    curve_table = table();
    if ~isstruct(gamma_definition) || isempty(fieldnames(gamma_definition))
        return;
    end

    if isfield(gamma_definition, 'gamma_true') && isfield(gamma_definition, 'psi_grid')
        psi_grid = gamma_definition.psi_grid(:);
        if isempty(psi_grid)
            return;
        end

        curve_table = table(psi_grid, gamma_definition.gamma_agent_1(:), gamma_definition.gamma_agent_2(:), gamma_definition.gamma_true(:), ...
            'VariableNames', {'psi', 'gamma_agent_1', 'gamma_agent_2', 'gamma_true'});
        return;
    end

    psi_grid = get_export_curve_field(gamma_definition, 'psi_grid_centered', 'psi_grid');
    if isempty(psi_grid)
        return;
    end

    curve_table = table(psi_grid(:), 'VariableNames', {'psi'});
    append_gamma_curve_column('gamma_selected', get_export_curve_field(gamma_definition, 'gamma_values_centered', 'gamma_values'));
    append_gamma_curve_column('gamma_full', get_export_curve_field(gamma_definition, 'gamma_values_full_centered', 'gamma_values_full'));
    append_gamma_curve_column('gamma_symmetric', get_export_curve_field(gamma_definition, 'gamma_values_symmetric_centered', 'gamma_values_symmetric'));
    append_gamma_curve_column('gamma_antisymmetric', get_export_curve_field(gamma_definition, 'gamma_values_antisymmetric_centered', 'gamma_values_antisymmetric'));

    function append_gamma_curve_column(variable_name, values)
        if isempty(values)
            return;
        end
        values = values(:);
        if numel(values) ~= height(curve_table)
            return;
        end
        curve_table.(variable_name) = values;
    end
end

function values = get_export_curve_field(gamma_definition, primary_field, secondary_field)
    values = [];
    if isfield(gamma_definition, primary_field) && ~isempty(gamma_definition.(primary_field))
        values = gamma_definition.(primary_field);
        return;
    end
    if nargin >= 3 && isfield(gamma_definition, secondary_field) && ~isempty(gamma_definition.(secondary_field))
        values = gamma_definition.(secondary_field);
    end
end
