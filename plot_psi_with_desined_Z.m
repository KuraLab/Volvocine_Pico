function varargout = plot_psi_with_desined_Z(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N, varargin)
% Analyze phase-coupling functions from CSV files in a directory.
%
% The function analyzes all agents listed in phase_agent_ids using Fourier
% fitting on derived signals built from target phase and normalized a2.
% The a2 normalization is done per agent using that agent's full available
% a2 series in the file: the 10th and 90th percentiles define the scale,
% their average defines the center, and their difference is normalized to 1.
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
%   out: struct with analysis settings and results for all target agents.
%
% Side effect:
%   None. This function does not export analysis files.

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
        dirpath = fullfile('EstimateQ', 'Spring5', '255');
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

    enable_save_figure = false;

    gamma_ratio = [1 1];

    % Ignore legacy weighted-fit arguments if they are still passed.
    if ~isempty(varargin)
        first_extra = varargin{1};
        if isnumeric(first_extra) && numel(first_extra) == 2 && all(isfinite(first_extra(:)))
            gamma_ratio = double(first_extra(:).');
        end
    end

    % ===== Derived-signal definition =====
    control_gain = 1;
    derived_signal_expression = sprintf('%g*cos(phase_target + pi - 0.6*pi) .* a2_normalized', control_gain);
    derived_signal_display_name = sprintf('%g*cos(phi_target + pi - 0.6*pi) * a2_norm', control_gain);
    derived_signal_axis_label = sprintf('%g*cos(phi_target + pi - 0.6*pi) * a2_{norm}', control_gain);
    derived_signal_func = @(phase_target, a2_normalized) control_gain * cos(phase_target + pi - 0.6*pi) .* a2_normalized;

    w1_model = load_exported_w1_model(dirpath);
    use_w1_derived_signal = w1_model.available;
    if use_w1_derived_signal
        derived_signal_w1_expression = sprintf('%g*W_1_exported(phase_target) .* a2_normalized', control_gain);
        derived_signal_w1_display_name = sprintf('%g*W_1_exported(phi_target) * a2_norm', control_gain);
        derived_signal_w1_axis_label = sprintf('%g*W_1_exported(phi_target) * a2_{norm}', control_gain);
        derived_signal_w1_func = @(phase_target, a2_normalized) control_gain * evaluate_exported_w1(phase_target, w1_model) .* a2_normalized;
        fprintf('[INFO] Using exported W1 phase sensitivity from %s\n', w1_model.source_file);
    else
        derived_signal_w1_expression = '';
        derived_signal_w1_display_name = '';
        derived_signal_w1_axis_label = '';
        derived_signal_w1_func = [];
        fprintf('[INFO] Exported W1 phase sensitivity was not found in %s. Running cos-based analysis only.\n', dirpath);
    end

    w2_model = load_exported_w2_model(dirpath);
    use_w2_derived_signal = w2_model.available;
    if use_w2_derived_signal
        derived_signal_w2_expression = sprintf('%g*W_2_exported(phase_target) .* a2_normalized', control_gain);
        derived_signal_w2_display_name = sprintf('%g*W_2_exported(phi_target) * a2_norm', control_gain);
        derived_signal_w2_axis_label = sprintf('%g*W_2_exported(phi_target) * a2_{norm}', control_gain);
        derived_signal_w2_func = @(phase_target, a2_normalized) control_gain * evaluate_exported_w2(phase_target, w2_model) .* a2_normalized;
        fprintf('[INFO] Using exported W2 phase sensitivity from %s\n', w2_model.source_file);
    else
        derived_signal_w2_expression = '';
        derived_signal_w2_display_name = '';
        derived_signal_w2_axis_label = '';
        derived_signal_w2_func = [];
        fprintf('[INFO] Exported W2 phase sensitivity was not found in %s. Skipping W2-based analysis.\n', dirpath);
    end

    [figure_phase_sensitivity, phase_sensitivity_power] = plot_phase_sensitivity_overlay( ...
        control_gain, use_w1_derived_signal, w1_model, use_w2_derived_signal, w2_model);

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

    gamma_settings_derived = struct('enabled', true, 'component', 'full', 'overlay_full', false, 'show_surface_overlay', false, 'auto_save_figure', enable_save_figure);
    analysis_agent_ids = unique([z_agent_id, phase_agent_ids(:).'], 'stable');

    n_analysis_agents = numel(analysis_agent_ids);
    agent_analysis = struct([]);
    for agent_idx = 1:n_analysis_agents
        target_agent_id = analysis_agent_ids(agent_idx);
        agent_result = run_single_agent_analysis( ...
            csv_paths, phase_agent_ids, target_agent_id, ...
            analysis_duration_sec, analysis_start_sec, sample_dt, ...
            M, N, gamma_ratio, gamma_settings_derived, ...
            derived_signal_func, derived_signal_w1_func, use_w1_derived_signal, ...
            derived_signal_w2_func, use_w2_derived_signal);

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
    true_gamma = compute_true_gamma_from_agent_analysis( ...
        agent_analysis, phase_agent_ids, gamma_ratio, ...
        'fourier_fit_sin_phi2_a2', 'cos', [0.8500, 0.3250, 0.0980]);
    true_gamma_w1 = compute_true_gamma_from_agent_analysis( ...
        agent_analysis, phase_agent_ids, gamma_ratio, ...
        'fourier_fit_sin_phi2_a2_w1', 'W1', [0.0000, 0.4470, 0.7410]);
    true_gamma_w2 = compute_true_gamma_from_agent_analysis( ...
        agent_analysis, phase_agent_ids, gamma_ratio, ...
        'fourier_fit_sin_phi2_a2_w2', 'W2', [0.4660, 0.6740, 0.1880]);
    phase_agent_omega = compute_phase_agent_mean_omega_all_files( ...
        csv_paths, phase_agent_ids, analysis_duration_sec, analysis_start_sec);

    figure_true_gamma_overlay = plot_true_gamma_overlay(true_gamma, true_gamma_w1, true_gamma_w2);
    if ~isempty(figure_true_gamma_overlay)
        if true_gamma.available
            true_gamma.figure = figure_true_gamma_overlay;
        end
        if true_gamma_w1.available
            true_gamma_w1.figure = figure_true_gamma_overlay;
        end
        if true_gamma_w2.available
            true_gamma_w2.figure = figure_true_gamma_overlay;
        end
    end

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
    out.gamma_settings = struct('derived', gamma_settings_derived);
    out.derived_signal = struct( ...
        'expression', derived_signal_expression, ...
        'display_name', derived_signal_display_name, ...
        'axis_label', derived_signal_axis_label);
    out.control_gain = control_gain;
    out.w1_model = rmfield_if_exists(w1_model, 'evaluator');
    out.w2_model = rmfield_if_exists(w2_model, 'evaluator');
    out.derived_signal_w1 = struct( ...
        'available', use_w1_derived_signal, ...
        'expression', derived_signal_w1_expression, ...
        'display_name', derived_signal_w1_display_name, ...
        'axis_label', derived_signal_w1_axis_label);
    out.derived_signal_w2 = struct( ...
        'available', use_w2_derived_signal, ...
        'expression', derived_signal_w2_expression, ...
        'display_name', derived_signal_w2_display_name, ...
        'axis_label', derived_signal_w2_axis_label);
    out.agent_analysis = agent_analysis;
    out.true_gamma = true_gamma;
    out.true_gamma_w1 = true_gamma_w1;
    out.true_gamma_w2 = true_gamma_w2;
    out.phase_agent_mean_omega = phase_agent_omega;
    out.figure_phase_sensitivity = unwrap_scalar_field(figure_phase_sensitivity);
    out.phase_sensitivity_power = phase_sensitivity_power;
    out.figure_true_gamma = unwrap_scalar_field(figure_true_gamma_overlay);
    out.figure_true_gamma_w1 = unwrap_scalar_field(figure_true_gamma_overlay);
    out.figure_true_gamma_w2 = unwrap_scalar_field(figure_true_gamma_overlay);

    % Top-level fields refer to the primary z_agent_id.
    out.used_files = unwrap_scalar_field(primary_analysis.used_files);
    out.skipped_files = unwrap_scalar_field(primary_analysis.skipped_files);
    out.per_file = unwrap_scalar_field(primary_analysis.per_file);
    out.fourier_fit_sin_phi2_a2 = unwrap_scalar_field(primary_analysis.fourier_fit_sin_phi2_a2);
    out.fourier_fit_sin_phi2_a2_w1 = unwrap_scalar_field(primary_analysis.fourier_fit_sin_phi2_a2_w1);
    out.fourier_fit_sin_phi2_a2_w2 = unwrap_scalar_field(primary_analysis.fourier_fit_sin_phi2_a2_w2);
    out.gamma_export = struct('available', false, 'reason', 'Export disabled in plot_psi_with_desined_Z.');

    if nargout == 0
        display_phase_agent_mean_omega_summary(out.phase_agent_mean_omega);
    else
        varargout{1} = out;
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

function agent_out = run_single_agent_analysis(csv_paths, phase_agent_ids, target_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt, M, N, gamma_ratio, gamma_settings_derived, derived_signal_func, derived_signal_w1_func, use_w1_derived_signal, derived_signal_w2_func, use_w2_derived_signal)
    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});
    per_file = struct('file_path', {}, 'window_start_abs', {}, 'window_end_abs', {}, ...
        'n_points', {}, 'phase_agent_ids', {}, 'z_agent_id', {});
    phi1_all = [];
    phi2_all = [];
    sin_phi2_a2_all = [];
    sin_phi2_a2_w1_all = [];
    sin_phi2_a2_w2_all = [];

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

        phase_for_derived_signal = select_phase_for_target_agent(point_data, phase_agent_ids, target_agent_id);
        sin_phi2_a2 = derived_signal_func(phase_for_derived_signal, point_data.a2_normalized);
        if use_w1_derived_signal && isa(derived_signal_w1_func, 'function_handle')
            sin_phi2_a2_w1 = derived_signal_w1_func(phase_for_derived_signal, point_data.a2_normalized);
        else
            sin_phi2_a2_w1 = [];
        end
        if use_w2_derived_signal && isa(derived_signal_w2_func, 'function_handle')
            sin_phi2_a2_w2 = derived_signal_w2_func(phase_for_derived_signal, point_data.a2_normalized);
        else
            sin_phi2_a2_w2 = [];
        end

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        sin_phi2_a2_all = [sin_phi2_a2_all; sin_phi2_a2(:)]; %#ok<AGROW>
        if use_w1_derived_signal
            sin_phi2_a2_w1_all = [sin_phi2_a2_w1_all; sin_phi2_a2_w1(:)]; %#ok<AGROW>
        end
        if use_w2_derived_signal
            sin_phi2_a2_w2_all = [sin_phi2_a2_w2_all; sin_phi2_a2_w2(:)]; %#ok<AGROW>
        end

        used_files{end+1} = csv_path; %#ok<AGROW>
        per_file(end+1) = meta; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available to overlay for agent %d.', target_agent_id);
    end

    fourier_fit_sin_phi2_a2 = fitDoubleFourierScatter( ...
        phi1_all, phi2_all, sin_phi2_a2_all, M, N, ...
        s_label, 'mixed-only', gamma_ratio, gamma_settings_derived);
    fourier_fit_sin_phi2_a2.M = M;
    fourier_fit_sin_phi2_a2.N = N;
    fourier_fit_sin_phi2_a2.fit_mode = 'unweighted';
    fourier_fit_sin_phi2_a2.gamma_ratio = gamma_ratio;

    fourier_fit_sin_phi2_a2_w1 = struct();
    if use_w1_derived_signal && ~isempty(sin_phi2_a2_w1_all)
        fourier_fit_sin_phi2_a2_w1 = fitDoubleFourierScatter( ...
            phi1_all, phi2_all, sin_phi2_a2_w1_all, M, N, ...
            s_label, 'mixed-only', gamma_ratio, gamma_settings_derived);
        fourier_fit_sin_phi2_a2_w1.M = M;
        fourier_fit_sin_phi2_a2_w1.N = N;
        fourier_fit_sin_phi2_a2_w1.fit_mode = 'unweighted';
        fourier_fit_sin_phi2_a2_w1.gamma_ratio = gamma_ratio;
    end

    fourier_fit_sin_phi2_a2_w2 = struct();
    if use_w2_derived_signal && ~isempty(sin_phi2_a2_w2_all)
        fourier_fit_sin_phi2_a2_w2 = fitDoubleFourierScatter( ...
            phi1_all, phi2_all, sin_phi2_a2_w2_all, M, N, ...
            s_label, 'mixed-only', gamma_ratio, gamma_settings_derived);
        fourier_fit_sin_phi2_a2_w2.M = M;
        fourier_fit_sin_phi2_a2_w2.N = N;
        fourier_fit_sin_phi2_a2_w2.fit_mode = 'unweighted';
        fourier_fit_sin_phi2_a2_w2.gamma_ratio = gamma_ratio;
    end

    agent_out = struct();
    agent_out.agent_id = target_agent_id;
    agent_out.fit_mode = 'unweighted';
    agent_out.gamma_ratio = gamma_ratio;
    agent_out.gamma_settings = struct('derived', gamma_settings_derived);
    agent_out.used_files = used_files;
    agent_out.skipped_files = skipped_files;
    agent_out.per_file = per_file;
    agent_out.fourier_fit_sin_phi2_a2 = fourier_fit_sin_phi2_a2;
    agent_out.fourier_fit_sin_phi2_a2_w1 = fourier_fit_sin_phi2_a2_w1;
    agent_out.fourier_fit_sin_phi2_a2_w2 = fourier_fit_sin_phi2_a2_w2;
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

function true_gamma = compute_true_gamma_from_agent_analysis(agent_analysis, phase_agent_ids, gamma_ratio, fit_field_name, mode_label, curve_color)
    if nargin < 4 || isempty(fit_field_name)
        fit_field_name = 'fourier_fit_sin_phi2_a2';
    end
    if nargin < 5 || isempty(mode_label)
        mode_label = 'cos';
    end
    if nargin < 6 || isempty(curve_color)
        curve_color = [0.8500, 0.3250, 0.0980];
    end

    true_gamma = struct( ...
        'available', false, ...
        'reason', '', ...
        'mode_label', mode_label, ...
        'curve_color', curve_color, ...
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

    if ~isfield(agent_analysis, fit_field_name)
        true_gamma.reason = sprintf('Field %s is not available in agent_analysis.', fit_field_name);
        return;
    end

    fit_agent_1 = unwrap_scalar_field(agent_analysis(idx_agent_1).(fit_field_name));
    fit_agent_2 = unwrap_scalar_field(agent_analysis(idx_agent_2).(fit_field_name));
    if ~isstruct(fit_agent_1) || ~isfield(fit_agent_1, 'gamma_resonance') || ...
            ~isstruct(fit_agent_2) || ~isfield(fit_agent_2, 'gamma_resonance')
        true_gamma.reason = sprintf('Gamma reconstruction results (%s) were not found for both agents.', mode_label);
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

    true_gamma.available = true;
    true_gamma.agent_id_1 = phase_agent_ids(1);
    true_gamma.agent_id_2 = phase_agent_ids(2);
    true_gamma.component = sprintf('full_%s', mode_label);
    true_gamma.psi_grid = psi_grid;
    true_gamma.gamma_agent_1 = gamma_agent_1;
    true_gamma.gamma_agent_2 = gamma_agent_2;
    true_gamma.gamma_true = gamma_true_values;
    true_gamma.figure = [];
end

function [fig_handle, power_summary] = plot_phase_sensitivity_overlay(control_gain, use_w1, w1_model, use_w2, w2_model)
    phase_grid = linspace(-pi, pi, 801);
    z_cos = control_gain * cos(phase_grid + pi - 0.6 * pi);
    power_cos = trapz(phase_grid, z_cos .^ 2) / (2 * pi);

    power_summary = struct( ...
        'definition', '(1/(2*pi))*integral_{-pi}^{pi} Z(phi)^2 dphi', ...
        'grid_count', numel(phase_grid), ...
        'cos', struct('available', true, 'power', power_cos), ...
        'w1', struct('available', false, 'power', NaN), ...
        'w2', struct('available', false, 'power', NaN), ...
        'ratio_w1_to_cos', NaN, ...
        'ratio_w2_to_cos', NaN);

    fig_handle = figure('Color', 'w', 'Name', 'phase sensitivity (overlay)');
    ax = axes('Parent', fig_handle);
    hold(ax, 'on');

    plot(ax, phase_grid, z_cos, 'LineWidth', 1.8, 'Color', [0.8500, 0.3250, 0.0980], ...
        'DisplayName', sprintf('%g*cos(\phi+\pi-0.6\pi) | P=%.6g', control_gain, power_cos));

    if use_w1 && isstruct(w1_model) && w1_model.available
        z_w1 = control_gain * evaluate_exported_w1(phase_grid, w1_model);
        power_w1 = trapz(phase_grid, z_w1 .^ 2) / (2 * pi);
        power_summary.w1.available = true;
        power_summary.w1.power = power_w1;
        power_summary.ratio_w1_to_cos = power_w1 / max(power_cos, eps);
        plot(ax, phase_grid, z_w1, 'LineWidth', 1.8, 'Color', [0.0000, 0.4470, 0.7410], ...
            'DisplayName', sprintf('%g*W_1(\phi) | P=%.6g', control_gain, power_w1));
    end

    if use_w2 && isstruct(w2_model) && w2_model.available
        z_w2 = control_gain * evaluate_exported_w2(phase_grid, w2_model);
        power_w2 = trapz(phase_grid, z_w2 .^ 2) / (2 * pi);
        power_summary.w2.available = true;
        power_summary.w2.power = power_w2;
        power_summary.ratio_w2_to_cos = power_w2 / max(power_cos, eps);
        plot(ax, phase_grid, z_w2, 'LineWidth', 1.8, 'Color', [0.4660, 0.6740, 0.1880], ...
            'DisplayName', sprintf('%g*W_2(\phi) | P=%.6g', control_gain, power_w2));
    end

    plot(ax, phase_grid, zeros(size(phase_grid)), ':', 'LineWidth', 1.0, 'Color', [0.5, 0.5, 0.5], ...
        'HandleVisibility', 'off');
    xlabel(ax, '$$\phi$$', 'Interpreter', 'latex');
    ylabel(ax, '$$Z(\phi)$$', 'Interpreter', 'latex');
    title(ax, '$$\mathrm{Phase\ sensitivity\ functions\ (overlay)}$$', 'Interpreter', 'latex');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [-pi, pi]);
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
    legend(ax, 'Location', 'best');
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    figure(fig_handle);
    tuneFigure;

    fprintf('[INFO] Phase-sensitivity power comparison: cos=%.12f\n', power_summary.cos.power);
    if power_summary.w1.available
        fprintf('[INFO]   W1=%.12f (W1/cos=%.12f)\n', power_summary.w1.power, power_summary.ratio_w1_to_cos);
    else
        fprintf('[INFO]   W1 not available.\n');
    end
    if power_summary.w2.available
        fprintf('[INFO]   W2=%.12f (W2/cos=%.12f)\n', power_summary.w2.power, power_summary.ratio_w2_to_cos);
    else
        fprintf('[INFO]   W2 not available.\n');
    end
end

function fig_handle = plot_true_gamma_overlay(true_gamma_cos, true_gamma_w1, true_gamma_w2)
    gamma_list = {true_gamma_cos, true_gamma_w1, true_gamma_w2};
    fig_handle = [];

    active_idx = [];
    for i = 1:numel(gamma_list)
        gamma_item = gamma_list{i};
        if isstruct(gamma_item) && isfield(gamma_item, 'available') && gamma_item.available && ...
                isfield(gamma_item, 'psi_grid') && ~isempty(gamma_item.psi_grid) && ...
                isfield(gamma_item, 'gamma_true') && ~isempty(gamma_item.gamma_true)
            active_idx(end + 1) = i; %#ok<AGROW>
        end
    end

    if isempty(active_idx)
        return;
    end

    fig_handle = figure('Color', 'w', 'Name', 'true gamma (overlay)');
    ax = axes('Parent', fig_handle);
    hold(ax, 'on');

    for idx = active_idx
        gamma_item = gamma_list{idx};
        mode_label = 'unknown';
        if isfield(gamma_item, 'mode_label') && ~isempty(gamma_item.mode_label)
            mode_label = gamma_item.mode_label;
        end

        line_color = [0.0000, 0.0000, 0.0000];
        if isfield(gamma_item, 'curve_color') && isnumeric(gamma_item.curve_color) && numel(gamma_item.curve_color) == 3
            line_color = gamma_item.curve_color;
        end

        plot(ax, gamma_item.psi_grid(:), gamma_item.gamma_true(:), 'LineWidth', 1.8, 'Color', line_color, ...
            'DisplayName', sprintf('\\Gamma_{true} [%s]', mode_label));
    end

    psi_ref = linspace(-pi, pi, 201);
    plot(ax, psi_ref, zeros(size(psi_ref)), ':', 'LineWidth', 1.0, 'Color', [0.5, 0.5, 0.5], ...
        'HandleVisibility', 'off');
    xlabel(ax, '$$\psi$$', 'Interpreter', 'latex');
    ylabel(ax, '$$\Gamma_{\mathrm{true}}(\psi)$$', 'Interpreter', 'latex');
    title(ax, '$$\Gamma_{\mathrm{true}}(\psi)\ \mathrm{overlay}$$', 'Interpreter', 'latex');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [-pi, pi]);
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
    legend(ax, 'Location', 'best');
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    figure(fig_handle);
    tuneFigure;
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
    point_data.a2_normalized = clip_values( ...
        normalize_by_agent_percentile_span(a2_z, series_by_agent(z_agent_id).a2, 10), -0.5, 0.5);

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

function model = load_exported_w1_model(dirpath)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', '', ...
        'source_format', '', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    search_dirs = {dirpath, fullfile(dirpath, 'gamma_exports')};
    [mat_path, mat_found] = find_latest_file(search_dirs, 'W_fourier_fit_params_*.mat');
    if mat_found
        model = load_w1_from_mat(mat_path);
        if model.available
            return;
        end
    end

    [csv_path, csv_found] = find_latest_file(search_dirs, 'W1_fourier_fit_params_agent*_*.csv');
    if csv_found
        model = load_w1_from_csv(csv_path);
        if model.available
            return;
        end
    end

    if ~model.available && isempty(model.reason)
        model.reason = 'No exported W1 parameter file was found.';
    end
end

function y = evaluate_exported_w1(phase, model)
    if ~isstruct(model) || ~isfield(model, 'available') || ~model.available
        error('evaluate_exported_w1 requires a valid loaded model.');
    end

    y = model.a0 * ones(size(phase));
    for k = 1:model.order
        y = y + model.a(k) * cos(k * phase) + model.b(k) * sin(k * phase);
    end
end

function model = load_w1_from_mat(mat_path)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', mat_path, ...
        'source_format', 'mat', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    try
        S = load(mat_path);
    catch ME
        model.reason = sprintf('Could not load MAT file: %s', ME.message);
        return;
    end

    if ~isfield(S, 'fourier_fit_export') || ~isstruct(S.fourier_fit_export)
        model.reason = 'MAT file does not contain fourier_fit_export.';
        return;
    end

    fit_export = S.fourier_fit_export;
    if ~isfield(fit_export, 'fit_w1') || ~isstruct(fit_export.fit_w1)
        model.reason = 'MAT file does not contain fit_w1.';
        return;
    end

    fit_w1 = fit_export.fit_w1;
    required_fields = {'order', 'a0', 'a', 'b'};
    if ~all(isfield(fit_w1, required_fields))
        model.reason = 'fit_w1 does not contain required Fourier fields.';
        return;
    end

    order = max(0, round(double(fit_w1.order)));
    a = reshape(double(fit_w1.a), [], 1);
    b = reshape(double(fit_w1.b), [], 1);
    if numel(a) < order || numel(b) < order
        model.reason = 'fit_w1 Fourier vectors are shorter than the declared order.';
        return;
    end

    model.available = true;
    model.reason = '';
    if isfield(fit_export, 'signal_role')
        model.signal_role = fit_export.signal_role;
    end
    model.order = order;
    model.a0 = double(fit_w1.a0);
    model.a = a(1:order);
    model.b = b(1:order);
end

function model = load_w1_from_csv(csv_path)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', csv_path, ...
        'source_format', 'csv', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    try
        T = readtable(csv_path);
    catch ME
        model.reason = sprintf('Could not read W1 CSV file: %s', ME.message);
        return;
    end

    required_columns = {'k', 'a_k', 'b_k'};
    if ~all(ismember(required_columns, T.Properties.VariableNames))
        model.reason = 'W1 CSV must contain columns: k, a_k, b_k.';
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
        model.reason = 'W1 CSV has no valid Fourier rows.';
        return;
    end

    order = max(k);
    a = zeros(order, 1);
    b = zeros(order, 1);

    idx0 = find(k == 0, 1, 'first');
    if isempty(idx0)
        model.reason = 'W1 CSV is missing k=0 (a0) row.';
        return;
    end
    a0 = a_k(idx0);

    for n = 1:order
        idxn = find(k == n, 1, 'first');
        if isempty(idxn)
            continue;
        end
        a(n) = a_k(idxn);
        b(n) = b_k(idxn);
    end

    model.available = true;
    model.reason = '';
    model.order = order;
    model.a0 = a0;
    model.a = a;
    model.b = b;
end

function model = load_exported_w2_model(dirpath)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', '', ...
        'source_format', '', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    search_dirs = {dirpath, fullfile(dirpath, 'gamma_exports')};
    [mat_path, mat_found] = find_latest_file(search_dirs, 'W_fourier_fit_params_*.mat');
    if mat_found
        model = load_w2_from_mat(mat_path);
        if model.available
            return;
        end
    end

    [csv_path, csv_found] = find_latest_file(search_dirs, 'W2_fourier_fit_params_agent*_*.csv');
    if csv_found
        model = load_w2_from_csv(csv_path);
        if model.available
            return;
        end
    end

    if ~model.available && isempty(model.reason)
        model.reason = 'No exported W2 parameter file was found.';
    end
end

function y = evaluate_exported_w2(phase, model)
    if ~isstruct(model) || ~isfield(model, 'available') || ~model.available
        error('evaluate_exported_w2 requires a valid loaded model.');
    end

    y = model.a0 * ones(size(phase));
    for k = 1:model.order
        y = y + model.a(k) * cos(k * phase) + model.b(k) * sin(k * phase);
    end
end

function model = load_w2_from_mat(mat_path)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', mat_path, ...
        'source_format', 'mat', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    try
        S = load(mat_path);
    catch ME
        model.reason = sprintf('Could not load MAT file: %s', ME.message);
        return;
    end

    if ~isfield(S, 'fourier_fit_export') || ~isstruct(S.fourier_fit_export)
        model.reason = 'MAT file does not contain fourier_fit_export.';
        return;
    end

    fit_export = S.fourier_fit_export;
    if ~isfield(fit_export, 'fit_w2') || ~isstruct(fit_export.fit_w2)
        model.reason = 'MAT file does not contain fit_w2.';
        return;
    end

    fit_w2 = fit_export.fit_w2;
    required_fields = {'order', 'a0', 'a', 'b'};
    if ~all(isfield(fit_w2, required_fields))
        model.reason = 'fit_w2 does not contain required Fourier fields.';
        return;
    end

    order = max(0, round(double(fit_w2.order)));
    a = reshape(double(fit_w2.a), [], 1);
    b = reshape(double(fit_w2.b), [], 1);
    if numel(a) < order || numel(b) < order
        model.reason = 'fit_w2 Fourier vectors are shorter than the declared order.';
        return;
    end

    model.available = true;
    model.reason = '';
    if isfield(fit_export, 'signal_role')
        model.signal_role = fit_export.signal_role;
    end
    model.order = order;
    model.a0 = double(fit_w2.a0);
    model.a = a(1:order);
    model.b = b(1:order);
end

function model = load_w2_from_csv(csv_path)
    model = struct( ...
        'available', false, ...
        'reason', '', ...
        'source_file', csv_path, ...
        'source_format', 'csv', ...
        'signal_role', '', ...
        'order', 0, ...
        'a0', 0, ...
        'a', [], ...
        'b', []);

    try
        T = readtable(csv_path);
    catch ME
        model.reason = sprintf('Could not read W2 CSV file: %s', ME.message);
        return;
    end

    required_columns = {'k', 'a_k', 'b_k'};
    if ~all(ismember(required_columns, T.Properties.VariableNames))
        model.reason = 'W2 CSV must contain columns: k, a_k, b_k.';
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
        model.reason = 'W2 CSV has no valid Fourier rows.';
        return;
    end

    order = max(k);
    a = zeros(order, 1);
    b = zeros(order, 1);

    idx0 = find(k == 0, 1, 'first');
    if isempty(idx0)
        model.reason = 'W2 CSV is missing k=0 (a0) row.';
        return;
    end
    a0 = a_k(idx0);

    for n = 1:order
        idxn = find(k == n, 1, 'first');
        if isempty(idxn)
            continue;
        end
        a(n) = a_k(idxn);
        b(n) = b_k(idxn);
    end

    model.available = true;
    model.reason = '';
    model.order = order;
    model.a0 = a0;
    model.a = a;
    model.b = b;
end

function [file_path, is_found] = find_latest_file(search_dirs, pattern)
    file_path = '';
    is_found = false;
    latest_datenum = -inf;

    for idx = 1:numel(search_dirs)
        search_dir = search_dirs{idx};
        if ~isfolder(search_dir)
            continue;
        end
        candidates = dir(fullfile(search_dir, pattern));
        for k = 1:numel(candidates)
            candidate = candidates(k);
            if candidate.datenum > latest_datenum
                latest_datenum = candidate.datenum;
                file_path = fullfile(candidate.folder, candidate.name);
                is_found = true;
            end
        end
    end
end

function S = rmfield_if_exists(S, field_name)
    if isstruct(S) && isfield(S, field_name)
        S = rmfield(S, field_name);
    end
end
