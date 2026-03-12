function out = plot_phase_pair_a2_3d_all_files(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N, varargin)
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
%   dirpath = 'EstimateF/Spring5/250'
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

    % Clear any figure windows left from previous runs before starting.
    existing_figures = findall(0, 'Type', 'figure');
    if ~isempty(existing_figures)
        close(existing_figures);
        drawnow;
    end

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('EstimateQ', 'Spring1', '250');
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
        analysis_start_sec = 5;
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
    derived_signal_func = @(phase_target, a2_normalized) +5*(cos(phase_target + pi - 0.6*pi) + 0 * cos(phase_target + pi + 0.0*pi)) .* a2_normalized;

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
    gamma_settings_a2 = struct('enabled', false, 'component', 'full', 'overlay_full', false, 'show_surface_overlay', true, 'resonant_only_bar_plot', false);
    gamma_settings_sin_phi2_a2 = struct('enabled', true, 'component', 'antisymmetric', 'overlay_full', true, 'show_surface_overlay', false, 'resonant_only_bar_plot', true);
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
            derived_signal_func, derived_signal_display_name);

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

    if nargout == 0
        display_phase_agent_mean_omega_summary(out.phase_agent_mean_omega);
    end
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

function agent_out = run_single_agent_analysis(csv_paths, phase_agent_ids, target_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt, M, N, gamma_ratio, point_color, marker_size, marker_alpha, gamma_settings_a2, gamma_settings_sin_phi2_a2, derived_signal_func, derived_signal_display_name) %#ok<INUSD>
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
    phase_name_for_derived_signal = select_phase_name_for_target_agent(phase_agent_ids, target_agent_id);
    derived_signal_display_name_target = strrep(derived_signal_display_name, 'phi_target', phase_name_for_derived_signal);

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

    configure_phase_pair_axes(ax, phase_agent_ids, sprintf('Agent %d a2_norm', target_agent_id));

    figure(fig);
    tuneFigure;

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
        sprintf('a2_norm(agent %d)', target_agent_id), 'full', gamma_ratio, gamma_settings_a2);
    fourier_fit.M = M;
    fourier_fit.N = N;
    fourier_fit.fit_mode = 'unweighted';
    fourier_fit.gamma_ratio = gamma_ratio;

    fourier_fit_sin_phi2_a2 = fitDoubleFourierScatter( ...
        point_cloud.phi1, point_cloud.phi2, point_cloud.sin_phi2_a2, M, N, ...
        sprintf('%s(agent %d)', derived_signal_display_name_target, target_agent_id), 'mixed-only', gamma_ratio, gamma_settings_sin_phi2_a2);
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

function phase_name = select_phase_name_for_target_agent(phase_agent_ids, target_agent_id)
    if target_agent_id == phase_agent_ids(1)
        phase_name = 'phi1';
    elseif target_agent_id == phase_agent_ids(2)
        phase_name = 'phi2';
    else
        phase_name = 'phi2';
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
    [psi_grid_1, gamma_values_1] = extract_gamma_curve(gamma_1);
    [psi_grid_2, gamma_values_2] = extract_gamma_curve(gamma_2);
    if isempty(psi_grid_1) || isempty(psi_grid_2) || isempty(gamma_values_1) || isempty(gamma_values_2)
        true_gamma.reason = 'Gamma curves were empty for at least one agent.';
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
    true_gamma.component = describe_gamma_component(gamma_1, gamma_2);
    true_gamma.psi_grid = psi_grid;
    true_gamma.gamma_agent_1 = gamma_agent_1;
    true_gamma.gamma_agent_2 = gamma_agent_2;
    true_gamma.gamma_true = gamma_true_values;
    true_gamma.figure = fig_true_gamma;
end

function [psi_grid, gamma_values] = extract_gamma_curve(gamma_resonance)
    psi_grid = [];
    gamma_values = [];
    if ~isstruct(gamma_resonance)
        return;
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

function component_label = describe_gamma_component(gamma_1, gamma_2)
    component_1 = '';
    component_2 = '';
    if isstruct(gamma_1) && isfield(gamma_1, 'component') && ~isempty(gamma_1.component)
        component_1 = char(gamma_1.component);
    end
    if isstruct(gamma_2) && isfield(gamma_2, 'component') && ~isempty(gamma_2.component)
        component_2 = char(gamma_2.component);
    end

    if strcmp(component_1, component_2) && ~isempty(component_1)
        component_label = component_1;
    else
        component_label = 'selected';
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

function configure_phase_pair_axes(ax, phase_agent_ids, z_axis_label)
    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlabel(ax, sprintf('Agent %d phase (rad)', phase_agent_ids(1)));
    ylabel(ax, sprintf('Agent %d phase (rad)', phase_agent_ids(2)));
    zlabel(ax, z_axis_label);
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
    point_data.a2_normalized = normalize_by_agent_percentile_span(a2_z, series_by_agent(z_agent_id).a2, 10);
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
