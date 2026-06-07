function out = plot_spring3_a2_reference_fit_error(varargin)
% plot_spring3_a2_reference_fit_error Compare Spring3 folders ending in
% _0531 against a Fourier fit learned from Spring3/255.
%
% The reference folder Spring3/255 is used to fit an unweighted double
% Fourier model for a2 on the usual (phi1, phi2) phase-pair features.
% The learned coefficients are then evaluated on the other immediate
% subfolders under Spring3 whose names end with _0531, and the prediction
% error is summarized and plotted.
%
% Usage:
%   out = plot_spring3_a2_reference_fit_error()
%   out = plot_spring3_a2_reference_fit_error('ShowFigures', true)
%   out = plot_spring3_a2_reference_fit_error('M', 12, 'N', 12)
%
% Name-value options:
%   'ParentDir'             : parent Spring3 directory (default: local EstimateQ/Spring3)
%   'ReferenceFolder'       : reference subfolder name (default: '255')
%   'AnalysisDurationSec'   : analysis window length (default: 80)
%   'AnalysisStartSec'      : analysis window start offset (default: 10)
%   'SampleDt'              : interpolation step in seconds (default: 0.01)
%   'M'                     : Fourier order in phi1 (default: 10)
%   'N'                     : Fourier order in phi2 (default: 10)
%   'PhaseAgentIds'         : optional two phase-agent IDs to force (default: auto detect)
%   'ZAgentId'              : optional a2 target agent ID to force (default: second phase agent)
%   'IncludeReferenceFolder': include the reference folder in the comparison (default: true)
%   'ShowFigures'           : draw summary figures (default: true)

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'ParentDir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ReferenceFolder', '240_0531', @(x) ischar(x) || isstring(x));
    addParameter(p, 'AnalysisDurationSec', 80, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'finite'}));
    addParameter(p, 'AnalysisStartSec', 10, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'finite'}));
    addParameter(p, 'SampleDt', 0.01, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'finite'}));
    addParameter(p, 'M', 10, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}));
    addParameter(p, 'N', 10, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}));
    addParameter(p, 'PhaseAgentIds', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    addParameter(p, 'ZAgentId', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'IncludeReferenceFolder', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowFigures', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    parent_dir = char(p.Results.ParentDir);
    if isempty(parent_dir)
        parent_dir = fullfile(fileparts(mfilename('fullpath')), 'Spring3');
    end
    reference_folder = char(p.Results.ReferenceFolder);
    analysis_duration_sec = double(p.Results.AnalysisDurationSec);
    analysis_start_sec = double(p.Results.AnalysisStartSec);
    sample_dt = double(p.Results.SampleDt);
    M = double(p.Results.M);
    N = double(p.Results.N);
    phase_agent_ids = p.Results.PhaseAgentIds;
    z_agent_id = p.Results.ZAgentId;
    include_reference_folder = p.Results.IncludeReferenceFolder;
    show_figures = p.Results.ShowFigures;

    ensure_local_function_folder_on_path();

    parent_dir = resolve_parent_dirpath(parent_dir);
    if ~isfolder(parent_dir)
        error('Directory not found: %s', parent_dir);
    end

    reference_dir = fullfile(parent_dir, reference_folder);
    if ~isfolder(reference_dir)
        error('Reference directory not found: %s', reference_dir);
    end

    reference_csv_paths = list_csv_paths(reference_dir);
    if isempty(reference_csv_paths)
        error('No CSV files found in reference directory: %s', reference_dir);
    end

    if isempty(phase_agent_ids)
        phase_agent_ids = detect_default_phase_agents(reference_csv_paths);
    else
        phase_agent_ids = phase_agent_ids(:).';
        if numel(phase_agent_ids) ~= 2 || phase_agent_ids(1) == phase_agent_ids(2)
            error('PhaseAgentIds must contain two distinct agent IDs.');
        end
    end

    if isempty(z_agent_id)
        z_agent_id = phase_agent_ids(2);
    end

    reference_point_cloud = build_point_cloud( ...
        reference_csv_paths, phase_agent_ids, z_agent_id, ...
        analysis_duration_sec, analysis_start_sec, sample_dt);

    reference_gamma_settings = struct( ...
        'enabled', false, ...
        'component', 'full', ...
        'overlay_full', false, ...
        'show_surface_overlay', false, ...
        'auto_save_figure', false);
    reference_fit = fitDoubleFourierScatter( ...
        reference_point_cloud.phi1, reference_point_cloud.phi2, reference_point_cloud.a2, ...
        M, N, sprintf('reference a2 (%s)', reference_folder), 'full', [1 1], reference_gamma_settings);

    folder_infos = list_spring3_subfolders(parent_dir);
    if isempty(folder_infos)
        error('No comparable subfolders with CSV files were found under %s.', parent_dir);
    end

    if include_reference_folder
        selected_mask = true(size(folder_infos));
    else
        selected_mask = ~[folder_infos.is_reference];
    end
    folder_infos = folder_infos(selected_mask);

    if isempty(folder_infos)
        error('No folders selected for comparison after applying IncludeReferenceFolder.');
    end

    comparison_results = repmat(struct( ...
        'folder_name', '', ...
        'folder_path', '', ...
        'csv_paths', {{}}, ...
        'used_files', {{}}, ...
        'skipped_files', struct('file_path', {}, 'reason', {}), ...
        'point_cloud', struct(), ...
        'prediction', struct(), ...
        'rmse', NaN, ...
        'mae', NaN, ...
        'bias', NaN, ...
        'corrcoef', NaN, ...
        'omega1_rad_s', NaN, ...
        'omega2_rad_s', NaN, ...
        'omega1_deg_s', NaN, ...
        'omega2_deg_s', NaN, ...
        'omega_delta_rad_s', NaN, ...
        'omega_delta_deg_s', NaN, ...
        'omega_ratio_2_over_1', NaN, ...
        'n_points', 0), numel(folder_infos), 1);

    for idx = 1:numel(folder_infos)
        folder_info = folder_infos(idx);
        csv_paths = list_csv_paths(folder_info.folder_path);
        if isempty(csv_paths)
            warning('Skipping %s: no CSV files were found.', folder_info.folder_path);
            continue;
        end

        try
            point_cloud = build_point_cloud( ...
                csv_paths, phase_agent_ids, z_agent_id, ...
                analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            warning('Skipping %s: %s', folder_info.folder_path, ME.message);
            comparison_results(idx).folder_name = folder_info.folder_name;
            comparison_results(idx).folder_path = folder_info.folder_path;
            comparison_results(idx).csv_paths = csv_paths;
            comparison_results(idx).skipped_files = struct('file_path', folder_info.folder_path, 'reason', ME.message);
            continue;
        end

        [prediction, metrics] = evaluate_reference_fit(point_cloud, reference_fit);

        comparison_results(idx).folder_name = folder_info.folder_name;
        comparison_results(idx).folder_path = folder_info.folder_path;
        comparison_results(idx).csv_paths = csv_paths;
        comparison_results(idx).used_files = point_cloud.used_files;
        comparison_results(idx).skipped_files = point_cloud.skipped_files;
        comparison_results(idx).point_cloud = point_cloud;
        comparison_results(idx).prediction = prediction;
        comparison_results(idx).rmse = metrics.rmse;
        comparison_results(idx).mae = metrics.mae;
        comparison_results(idx).bias = metrics.bias;
        comparison_results(idx).corrcoef = metrics.corrcoef;
        comparison_results(idx).omega1_rad_s = metrics.omega1_rad_s;
        comparison_results(idx).omega2_rad_s = metrics.omega2_rad_s;
        comparison_results(idx).omega1_deg_s = metrics.omega1_deg_s;
        comparison_results(idx).omega2_deg_s = metrics.omega2_deg_s;
        comparison_results(idx).omega_delta_rad_s = metrics.omega_delta_rad_s;
        comparison_results(idx).omega_delta_deg_s = metrics.omega_delta_deg_s;
        comparison_results(idx).omega_ratio_2_over_1 = metrics.omega_ratio_2_over_1;
        comparison_results(idx).n_points = numel(point_cloud.a2);
    end

    valid_mask = arrayfun(@(s) isfinite(s.rmse) && s.n_points > 0, comparison_results);
    comparison_results = comparison_results(valid_mask);
    if isempty(comparison_results)
        error('No valid folders produced comparison results.');
    end

    summary_table = build_summary_table(comparison_results);

    figure_rmse = [];
    figure_mae = [];
    figure_residual = [];
    figure_omega = [];
    if show_figures
        figure_rmse = plot_single_error_summary(summary_table, reference_folder, 'rmse', 'RMSE', 'Spring3 a2 reference RMSE');
        figure_mae = plot_single_error_summary(summary_table, reference_folder, 'mae', 'MAE', 'Spring3 a2 reference MAE');
        figure_residual = plot_residual_summary(comparison_results, reference_folder);
        figure_omega = plot_omega_summary(summary_table, reference_folder, phase_agent_ids);
    end

    out = struct();
    out.parent_dir = parent_dir;
    out.reference_folder = reference_folder;
    out.reference_dir = reference_dir;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.sample_dt = sample_dt;
    out.M = M;
    out.N = N;
    out.reference_point_cloud = reference_point_cloud;
    out.reference_fit = reference_fit;
    out.comparison_results = comparison_results;
    out.summary_table = summary_table;
    out.figure_error = figure_rmse;
    out.figure_rmse = figure_rmse;
    out.figure_mae = figure_mae;
    out.figure_residual = figure_residual;
    out.figure_omega = figure_omega;

    fprintf('[INFO] Reference fit built from %s (%d samples).\n', reference_dir, numel(reference_point_cloud.a2));
    fprintf('[INFO] Compared %d folder(s) under %s.\n', numel(comparison_results), parent_dir);
    disp(summary_table);
end

function parent_dir = resolve_parent_dirpath(parent_dir)
    if isfolder(parent_dir)
        return;
    end

    candidate_pwd = fullfile(pwd, parent_dir);
    if isfolder(candidate_pwd)
        parent_dir = candidate_pwd;
        return;
    end

    base_dir = fileparts(mfilename('fullpath'));
    candidate_base = fullfile(base_dir, parent_dir);
    if isfolder(candidate_base)
        parent_dir = candidate_base;
        return;
    end

    parent_dir = candidate_pwd;
end

function ensure_local_function_folder_on_path()
    local_dir = fileparts(mfilename('fullpath'));
    parent_dir = fileparts(local_dir);

    if isempty(which('load_corrected_agent_series_from_csv')) && ~contains(path, parent_dir)
        addpath(parent_dir);
    end
    if isempty(which('fitDoubleFourierScatter')) && ~contains(path, local_dir)
        addpath(local_dir);
    end
end

function csv_paths = list_csv_paths(dirpath)
    files = dir(fullfile(dirpath, '*.csv'));
    if isempty(files)
        csv_paths = {};
        return;
    end

    names = sort({files.name});
    csv_paths = cellfun(@(name) fullfile(dirpath, name), names, 'UniformOutput', false);
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

    error('Could not determine default phase_agent_ids from the reference CSV files.');
end

function folder_infos = list_spring3_subfolders(parent_dir)
    listing = dir(parent_dir);
    listing = listing([listing.isdir]);
    listing = listing(~ismember({listing.name}, {'.', '..'}));

    folder_infos = struct('folder_name', {}, 'folder_path', {}, 'is_reference', {});
    for i = 1:numel(listing)
        folder_name = listing(i).name;
        if isempty(regexp(folder_name, '_0531$', 'once'))
            continue;
        end
        folder_path = fullfile(parent_dir, folder_name);
        if isempty(dir(fullfile(folder_path, '*.csv')))
            continue;
        end

        folder_infos(end + 1) = struct( ...
            'folder_name', folder_name, ...
            'folder_path', folder_path, ...
            'is_reference', strcmp(folder_name, '255')); %#ok<AGROW>
    end

    if isempty(folder_infos)
        return;
    end

    folder_infos = sort_folder_infos(folder_infos);
end

function folder_infos = sort_folder_infos(folder_infos)
    folder_names = {folder_infos.folder_name};
    numeric_keys = nan(size(folder_names));
    for i = 1:numel(folder_names)
        numeric_keys(i) = str2double(folder_names{i});
    end

    if all(isfinite(numeric_keys))
        [~, order] = sort(numeric_keys);
    else
        [~, order] = sort(lower(folder_names));
    end

    folder_infos = folder_infos(order);
end

function point_cloud = build_point_cloud(csv_paths, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt)
    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});
    phi1_all = [];
    phi2_all = [];
    a2_raw_all = [];
    a2_all = [];
    time_all = [];
    file_id_all = [];
    omega1_all = [];
    omega2_all = [];

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            [point_data, ~] = compute_points_for_csv( ...
                csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            warning('Skipping %s: %s', csv_path, ME.message);
            skipped_files(end + 1) = struct('file_path', csv_path, 'reason', ME.message); %#ok<AGROW>
            continue;
        end

        if isempty(point_data.time)
            skipped_files(end + 1) = struct('file_path', csv_path, 'reason', 'No valid points after processing.'); %#ok<AGROW>
            continue;
        end

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        a2_raw_all = [a2_raw_all; point_data.a2_raw(:)]; %#ok<AGROW>
        a2_all = [a2_all; point_data.a2(:)]; %#ok<AGROW>
        time_all = [time_all; point_data.time(:)]; %#ok<AGROW>
        file_id_all = [file_id_all; i * ones(numel(point_data.time), 1)]; %#ok<AGROW>
        [omega1_rad_s, omega2_rad_s] = estimate_file_omegas(point_data.time, point_data.phase1, point_data.phase2);
        omega1_all = [omega1_all; omega1_rad_s]; %#ok<AGROW>
        omega2_all = [omega2_all; omega2_rad_s]; %#ok<AGROW>
        used_files{end + 1} = csv_path; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available to build the point cloud.');
    end

    point_cloud = struct();
    point_cloud.phi1 = phi1_all;
    point_cloud.phi2 = phi2_all;
    point_cloud.a2_raw = a2_raw_all;
    point_cloud.a2 = a2_all;
    point_cloud.time = time_all;
    point_cloud.file_index = file_id_all;
    point_cloud.omega1_rad_s = omega1_all;
    point_cloud.omega2_rad_s = omega2_all;
    point_cloud.omega1_deg_s = omega1_all * (180 / pi);
    point_cloud.omega2_deg_s = omega2_all * (180 / pi);
    point_cloud.used_files = used_files;
    point_cloud.skipped_files = skipped_files;
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

function x_norm = normalize_by_agent_percentile_span(x, reference_x, tail_pct)
    x = double(x(:));
    reference_x = double(reference_x(:));

    valid = isfinite(x);
    reference_valid = isfinite(reference_x);
    x_norm = nan(size(x));

    if nnz(valid) < 1
        return;
    end

    if nnz(reference_valid) < 5
        x_norm(valid) = x(valid);
        return;
    end

    reference_values = reference_x(reference_valid);
    low_value = prctile(reference_values, tail_pct);
    high_value = prctile(reference_values, 100 - tail_pct);
    center_value = 0.5 * (low_value + high_value);
    span_value = high_value - low_value;

    if ~isfinite(span_value) || span_value <= 0
        x_norm(valid) = x(valid) - center_value;
        return;
    end

    x_norm(valid) = (x(valid) - center_value) / span_value;
end

function x_clipped = clip_values(x, lower_bound, upper_bound)
    x_clipped = min(max(x, lower_bound), upper_bound);
end

function [prediction, metrics] = evaluate_reference_fit(point_cloud, reference_fit)
    A = build_double_fourier_design_matrix(point_cloud.phi1, point_cloud.phi2, reference_fit.M, reference_fit.N);
    z_hat = A * reference_fit.coeff + reference_fit.z_mean;
    residual = point_cloud.a2 - z_hat;

    prediction = struct();
    prediction.z_hat = z_hat;
    prediction.residual = residual;

    metrics = struct();
    metrics.rmse = sqrt(mean(residual .^ 2));
    metrics.mae = mean(abs(residual));
    metrics.bias = mean(residual);
    if numel(point_cloud.a2) >= 2
        c = corrcoef(point_cloud.a2, z_hat);
        metrics.corrcoef = c(1, 2);
    else
        metrics.corrcoef = NaN;
    end

    metrics.omega1_rad_s = mean(point_cloud.omega1_rad_s, 'omitnan');
    metrics.omega2_rad_s = mean(point_cloud.omega2_rad_s, 'omitnan');
    metrics.omega1_deg_s = mean(point_cloud.omega1_deg_s, 'omitnan');
    metrics.omega2_deg_s = mean(point_cloud.omega2_deg_s, 'omitnan');
    metrics.omega_delta_rad_s = metrics.omega2_rad_s - metrics.omega1_rad_s;
    metrics.omega_delta_deg_s = metrics.omega2_deg_s - metrics.omega1_deg_s;
    if isfinite(metrics.omega1_rad_s) && metrics.omega1_rad_s ~= 0
        metrics.omega_ratio_2_over_1 = metrics.omega2_rad_s / metrics.omega1_rad_s;
    else
        metrics.omega_ratio_2_over_1 = NaN;
    end
end

function summary_table = build_summary_table(comparison_results)
    folder_name = string({comparison_results.folder_name}).';
    n_points = [comparison_results.n_points].';
    rmse = [comparison_results.rmse].';
    mae = [comparison_results.mae].';
    bias = [comparison_results.bias].';
    corrcoef_values = [comparison_results.corrcoef].';
    omega1_rad_s = [comparison_results.omega1_rad_s].';
    omega2_rad_s = [comparison_results.omega2_rad_s].';
    omega1_deg_s = [comparison_results.omega1_deg_s].';
    omega2_deg_s = [comparison_results.omega2_deg_s].';
    omega_delta_rad_s = [comparison_results.omega_delta_rad_s].';
    omega_delta_deg_s = [comparison_results.omega_delta_deg_s].';
    omega_ratio_2_over_1 = [comparison_results.omega_ratio_2_over_1].';

    summary_table = table(folder_name, n_points, rmse, mae, bias, corrcoef_values, ...
        omega1_rad_s, omega2_rad_s, omega1_deg_s, omega2_deg_s, omega_delta_rad_s, omega_delta_deg_s, omega_ratio_2_over_1, ...
        'VariableNames', {'folder_name', 'n_points', 'rmse', 'mae', 'bias', 'corrcoef', ...
        'omega1_rad_s', 'omega2_rad_s', 'omega1_deg_s', 'omega2_deg_s', 'omega_delta_rad_s', 'omega_delta_deg_s', 'omega_ratio_2_over_1'});
end

function fig = plot_single_error_summary(summary_table, reference_folder, metric_field, metric_label, fig_name)
    fig = figure('Color', 'w', 'Name', fig_name);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    x = categorical(summary_table.folder_name);
    x = reordercats(x, cellstr(summary_table.folder_name));
    bar(ax, x, summary_table.(metric_field));
    ylabel(ax, metric_label);
    xlabel(ax, '$$\Delta\omega$$','Interpreter', 'latex');
    %title(ax, sprintf('Reference fit from Spring3/%s', reference_folder));
    grid(ax, 'on');
    xtickangle(ax, 30);
    labels = format_folder_xticklabels(cellstr(summary_table.folder_name));
    xticklabels(ax, labels);
    set(ax, 'TickLabelInterpreter', 'latex');
    tuneFigure;
end

function fig = plot_residual_summary(comparison_results, reference_folder)
    fig = figure('Color', 'w', 'Name', 'Spring3 a2 residuals');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    colors = lines(numel(comparison_results));

    for i = 1:numel(comparison_results)
        residual = comparison_results(i).prediction.residual(:);
        if isempty(residual)
            continue;
        end

        jitter = linspace(-0.16, 0.16, numel(residual)).';
        scatter(ax, i + jitter, residual, 10, 'filled', ...
            'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.12, ...
            'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'none');
        plot(ax, i, mean(residual), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    end

    yline(ax, 0, ':', 'Color', [0.4, 0.4, 0.4]);
    grid(ax, 'on');
    xlim(ax, [0.5, numel(comparison_results) + 0.5]);
    xticks(ax, 1:numel(comparison_results));
    names = {comparison_results.folder_name};
    labels = format_folder_xticklabels(names);
    xticklabels(ax, labels);
    xtickangle(ax, 30);
    xlabel(ax, '$$\Delta\omega$$','Interpreter', 'latex');
    ylabel(ax, 'Residual (a2 - prediction)');
    %title(ax, sprintf('Residuals vs reference fit from Spring3/%s', reference_folder));
end

function fig = plot_omega_summary(summary_table, reference_folder, phase_agent_ids)
    fig = figure('Color', 'w', 'Name', 'Spring3 angular velocity summary');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    x = categorical(summary_table.folder_name);
    x = reordercats(x, cellstr(summary_table.folder_name));
    bar(ax, x, [summary_table.omega1_rad_s, summary_table.omega2_rad_s], 'grouped');
    ylabel(ax, 'Mean angular velocity (rad/s)');
    xlabel(ax, '$$\Delta\omega$$','Interpreter', 'latex');
    legend(ax, {sprintf('Agent %d', phase_agent_ids(1)), sprintf('Agent %d', phase_agent_ids(2))}, 'Location', 'northwest');
    title(ax, sprintf('Per-folder mean angular velocity using reference window from Spring3/%s', reference_folder));
    grid(ax, 'on');
    xtickangle(ax, 30);
    labels = format_folder_xticklabels(cellstr(summary_table.folder_name));
    xticklabels(ax, labels);
    set(ax, 'TickLabelInterpreter', 'latex');
end

function [omega1_rad_s, omega2_rad_s] = estimate_file_omegas(time_sec, phase1_rad, phase2_rad)
    omega1_rad_s = estimate_mean_omega_from_phase(time_sec, phase1_rad);
    omega2_rad_s = estimate_mean_omega_from_phase(time_sec, phase2_rad);
end

function labels = format_folder_xticklabels(folder_names)
    % Convert folder names like '260_0531' to LaTeX labels of the form
    % '$<value>/\pi$' where value = 2.5 - (numeric_part / 100).
    if isstring(folder_names)
        folder_names = cellstr(folder_names);
    end
    labels = cell(size(folder_names));
    for k = 1:numel(folder_names)
        name = folder_names{k};
        % remove trailing _0531 if present
        numpart = regexprep(name, '_0531$', '');
        val = str2double(numpart);
        if isnan(val)
            % fallback: escape underscores for LaTeX
            labels{k} = strrep(name, '_', '\_');
            continue;
        end
        scaled = val / 100;
        result = -2.5 + scaled;
        labels{k} = sprintf('$%0.3g\\pi$', result);
    end
end

function omega_rad_s = estimate_mean_omega_from_phase(time_sec, phase_rad)
    time_sec = time_sec(:);
    phase_rad = phase_rad(:);
    valid = isfinite(time_sec) & isfinite(phase_rad);
    time_sec = time_sec(valid);
    phase_rad = phase_rad(valid);

    if numel(time_sec) < 2
        omega_rad_s = NaN;
        return;
    end

    phase_unwrapped = unwrap(phase_rad);
    p = polyfit(time_sec - time_sec(1), phase_unwrapped, 1);
    omega_rad_s = p(1);
end

function [A, basis_names, basis_groups, basis_m, basis_n, basis_types] = build_double_fourier_design_matrix(phi1, phi2, M, N)
% build_double_fourier_design_matrix Construct the real trigonometric basis matrix.

    phi1 = phi1(:);
    phi2 = phi2(:);

    n_samples = numel(phi1);
    n_basis = 1 + 2 * M + 2 * N + 4 * M * N;
    A = zeros(n_samples, n_basis);
    basis_names = cell(n_basis, 1);
    basis_groups = cell(n_basis, 1);
    basis_m = zeros(n_basis, 1);
    basis_n = zeros(n_basis, 1);
    basis_types = cell(n_basis, 1);

    col = 1;

    A(:, col) = 1;
    basis_names{col} = '1';
    basis_groups{col} = 'constant';
    basis_m(col) = 0;
    basis_n(col) = 0;
    basis_types{col} = 'constant';
    col = col + 1;

    for m = 1:M
        A(:, col) = cos(m * phi1);
        basis_names{col} = sprintf('cos(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        basis_n(col) = 0;
        basis_types{col} = 'phi1_cos';
        col = col + 1;

        A(:, col) = sin(m * phi1);
        basis_names{col} = sprintf('sin(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        basis_n(col) = 0;
        basis_types{col} = 'phi1_sin';
        col = col + 1;
    end

    for n = 1:N
        A(:, col) = cos(n * phi2);
        basis_names{col} = sprintf('cos(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_m(col) = 0;
        basis_n(col) = n;
        basis_types{col} = 'phi2_cos';
        col = col + 1;

        A(:, col) = sin(n * phi2);
        basis_names{col} = sprintf('sin(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_m(col) = 0;
        basis_n(col) = n;
        basis_types{col} = 'phi2_sin';
        col = col + 1;
    end

    for m = 1:M
        c1 = cos(m * phi1);
        s1 = sin(m * phi1);
        for n = 1:N
            c2 = cos(n * phi2);
            s2 = sin(n * phi2);

            A(:, col) = c1 .* c2;
            basis_names{col} = sprintf('cos(%d*phi1)cos(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_cc';
            col = col + 1;

            A(:, col) = c1 .* s2;
            basis_names{col} = sprintf('cos(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_cs';
            col = col + 1;

            A(:, col) = s1 .* c2;
            basis_names{col} = sprintf('sin(%d*phi1)cos(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_sc';
            col = col + 1;

            A(:, col) = s1 .* s2;
            basis_names{col} = sprintf('sin(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_ss';
            col = col + 1;
        end
    end
end