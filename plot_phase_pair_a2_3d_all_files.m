function out = plot_phase_pair_a2_3d_all_files(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N, varargin)
% Overlay 3D phase-phase-a2 trajectories from all CSV files in a directory.
%
% x = phase of agent 1
% y = phase of agent 2
% z = a2 of z_agent_id
% plotted as a simple 3D point cloud
% In addition to a2 itself, the function also analyzes sin(phi2) * a2_norm
% using the same Fourier fitting workflow, where a2_norm is normalized by
% a percentile-range scaling.
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
%   z_agent_id = phase_agent_ids(2)
%   analysis_duration_sec = 15
%   analysis_start_sec = 5
%   file_indices = []  % [] means all CSV files in name-sorted order
%   M, N = 10          % default Fourier orders for automatic fitting
%   gamma_ratio = [2 1] for psi = 2*phi2 - phi1
%
% Output:
%   out: struct with figure handle, used/skipped files, plotting settings,
%        aggregated point cloud, and Fourier-fit results.

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('EstimateF', 'Spring5', '250');
    end
    if nargin < 2
        phase_agent_ids = [];
    end
    if nargin < 3
        z_agent_id = [];
    end
    if nargin < 4 || isempty(analysis_duration_sec)
        analysis_duration_sec = 15;
    end
    if nargin < 5 || isempty(analysis_start_sec)
        analysis_start_sec = 5;
    end
    if nargin < 6
        file_indices = [];
    end
    if nargin < 7 || isempty(M)
        M = 5;
    end
    if nargin < 8 || isempty(N)
        N = 5;
    end

    gamma_ratio = [1 1];

    % Ignore legacy weighted-fit arguments if they are still passed.
    if ~isempty(varargin)
        first_extra = varargin{1};
        if isnumeric(first_extra) && numel(first_extra) == 2 && all(isfinite(first_extra(:)))
            gamma_ratio = double(first_extra(:).');
        end
    end

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

    fig = figure('Color', 'w');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    fig_sin_phi2_a2 = figure('Color', 'w');
    ax_sin_phi2_a2 = axes('Parent', fig_sin_phi2_a2);
    hold(ax_sin_phi2_a2, 'on');
    point_color = [0.0, 0.4470, 0.7410];
    marker_size = 10;
    marker_alpha = 0.18;

    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});
    per_file = struct('file_path', {}, 'window_start_abs', {}, 'window_end_abs', {}, ...
        'n_points', {}, 'phase_agent_ids', {}, 'z_agent_id', {});
    phi1_all = [];
    phi2_all = [];
    a2_all = [];
    a2_normalized_all = [];
    sin_phi2_a2_all = [];
    time_all = [];
    file_id_all = [];

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            [point_data, meta] = compute_points_for_csv( ...
                csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            warning('Skipping %s: %s', csv_path, ME.message);
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

        sin_phi2_a2 = -sin(point_data.phase2+0.6*pi) .* point_data.a2_normalized;
        scatter3(ax_sin_phi2_a2, point_data.phase1, point_data.phase2, sin_phi2_a2, marker_size, ...
            'filled', 'MarkerFaceColor', point_color, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha);

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        a2_all = [a2_all; point_data.a2(:)]; %#ok<AGROW>
        a2_normalized_all = [a2_normalized_all; point_data.a2_normalized(:)]; %#ok<AGROW>
        sin_phi2_a2_all = [sin_phi2_a2_all; sin_phi2_a2(:)]; %#ok<AGROW>
        time_all = [time_all; point_data.time(:)]; %#ok<AGROW>
        file_id_all = [file_id_all; i * ones(numel(point_data.time), 1)]; %#ok<AGROW>

        used_files{end+1} = csv_path; %#ok<AGROW>
        per_file(end+1) = meta; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available to overlay.');
    end

    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlabel(ax, sprintf('Agent %d phase (rad)', phase_agent_ids(1)));
    ylabel(ax, sprintf('Agent %d phase (rad)', phase_agent_ids(2)));
    zlabel(ax, sprintf('Agent %d a2', z_agent_id));
    title(ax, sprintf('Overlayed 3D point cloud: agents [%d, %d], z=a2(agent %d)', ...
        phase_agent_ids(1), phase_agent_ids(2), z_agent_id));
    grid(ax, 'on');
    view(ax, 3);
    box(ax, 'on');

    xlim(ax_sin_phi2_a2, [0, 2*pi]);
    ylim(ax_sin_phi2_a2, [0, 2*pi]);
    xticks(ax_sin_phi2_a2, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax_sin_phi2_a2, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax_sin_phi2_a2, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax_sin_phi2_a2, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlabel(ax_sin_phi2_a2, sprintf('Agent %d phase (rad)', phase_agent_ids(1)));
    ylabel(ax_sin_phi2_a2, sprintf('Agent %d phase (rad)', phase_agent_ids(2)));
    zlabel(ax_sin_phi2_a2, sprintf('sin(phi2) * a2_{norm}(agent %d)', z_agent_id));
    title(ax_sin_phi2_a2, sprintf('Overlayed 3D point cloud: agents [%d, %d], z=sin(phi2) * a2_{norm}(agent %d)', ...
        phase_agent_ids(1), phase_agent_ids(2), z_agent_id));
    grid(ax_sin_phi2_a2, 'on');
    view(ax_sin_phi2_a2, 3);
    box(ax_sin_phi2_a2, 'on');

    sgtitle(fig, sprintf('%s: %d/%d files overlaid, window %.0f-%.0fs', ...
        dirpath, numel(used_files), numel(csv_paths), analysis_start_sec, analysis_start_sec + analysis_duration_sec));
    sgtitle(fig_sin_phi2_a2, sprintf('%s: %d/%d files overlaid, window %.0f-%.0fs', ...
        dirpath, numel(used_files), numel(csv_paths), analysis_start_sec, analysis_start_sec + analysis_duration_sec));

    out = struct();
    out.dirpath = dirpath;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.file_indices = file_indices;
    out.fit_mode = 'unweighted';
    out.gamma_ratio = gamma_ratio;
    out.figure = fig;
    out.figure_sin_phi2_a2 = fig_sin_phi2_a2;
    out.used_files = used_files;
    out.skipped_files = skipped_files;
    out.per_file = per_file;
    out.point_color = point_color;
    out.marker_size = marker_size;
    out.point_cloud = struct( ...
        'phi1', phi1_all, ...
        'phi2', phi2_all, ...
        'a2', a2_all, ...
        'a2_normalized', a2_normalized_all, ...
        'sin_phi2_a2', sin_phi2_a2_all, ...
        'time', time_all, ...
        'file_index', file_id_all);

    out.fourier_fit = fitDoubleFourierScatter( ...
        out.point_cloud.phi1, out.point_cloud.phi2, out.point_cloud.a2, M, N, ...
        sprintf('a2(agent %d)', z_agent_id), 'full', gamma_ratio);
    out.fourier_fit.M = M;
    out.fourier_fit.N = N;
    out.fourier_fit.fit_mode = 'unweighted';
    out.fourier_fit.gamma_ratio = gamma_ratio;

    out.fourier_fit_sin_phi2_a2 = fitDoubleFourierScatter( ...
        out.point_cloud.phi1, out.point_cloud.phi2, out.point_cloud.sin_phi2_a2, M, N, ...
        sprintf('sin(phi2) * a2_norm(agent %d)', z_agent_id), 'mixed-only', gamma_ratio);
    out.fourier_fit_sin_phi2_a2.M = M;
    out.fourier_fit_sin_phi2_a2.N = N;
    out.fourier_fit_sin_phi2_a2.fit_mode = 'unweighted';
    out.fourier_fit_sin_phi2_a2.gamma_ratio = gamma_ratio;
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

function [point_data, meta] = compute_points_for_csv(csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt)
    T = readtable(csv_path);

    required_cols = {'time_pc_sec_abs', 'a0', 'a2'};
    if ~all(ismember(required_cols, T.Properties.VariableNames))
        error('CSV missing required columns. Required: %s', strjoin(required_cols, ', '));
    end

    T.time_pc_sec_abs = double(T.time_pc_sec_abs);

    has_agent = ismember('agent_id', T.Properties.VariableNames);
    has_chunk = ismember('chunk_id', T.Properties.VariableNames);

    if ~has_agent
        T.agent_id = ones(height(T), 1);
    end

    if has_agent && has_chunk
        T_OVERFLOW = 2^32 / 1e6;
        T_TOL = 5.0;
        threshold_sec = T_OVERFLOW - T_TOL;
        jump_sec = T_OVERFLOW;

        T = correct_large_jump_matlab(T, threshold_sec, jump_sec);
        T = correct_chunk_start_times_matlab(T, 4000.0, T_OVERFLOW);
    end

    all_agents = unique(T.agent_id, 'sorted').';
    requested_agents = unique([phase_agent_ids(:); z_agent_id]);
    if ~all(ismember(requested_agents, all_agents))
        error('Requested agents %s are not all present. Available agents: %s', ...
            mat2str(requested_agents.'), mat2str(all_agents));
    end

    series_by_agent = struct();
    overlap_start = -inf;
    overlap_end = inf;

    for k = 1:numel(requested_agents)
        aid = requested_agents(k);
        series = get_agent_series(T, aid);
        if isempty(series.time)
            error('No valid samples found for agent %d.', aid);
        end

        series_by_agent(aid).time = series.time;
        series_by_agent(aid).a0_corr = series.a0_corr;
        series_by_agent(aid).a2 = series.a2;

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
    point_data.a2 = a2_z;
    point_data.a2_normalized = normalize_by_percentile_range(a2_z, 1, 99);

    meta = struct();
    meta.file_path = csv_path;
    meta.window_start_abs = window_start_abs;
    meta.window_end_abs = window_end_abs;
    meta.n_points = numel(point_data.time);
    meta.phase_agent_ids = phase_agent_ids;
    meta.z_agent_id = z_agent_id;
end

function xNorm = normalize_by_percentile_range(x, lowPct, highPct)
    x = double(x(:));
    valid = isfinite(x);
    xNorm = nan(size(x));

    if nnz(valid) < 5
        xNorm(valid) = x(valid);
        return;
    end

    xValid = x(valid);
    pLow = prctile(xValid, lowPct);
    pHigh = prctile(xValid, highPct);
    denom = pHigh - pLow;
    center = 0.5 * (pLow + pHigh);

    if ~isfinite(denom) || denom <= 0
        xNorm(valid) = xValid - center;
        return;
    end

    xNorm(valid) = (xValid - center) / denom;
end

function series = get_agent_series(T, agent_id)
    sub = T(T.agent_id == agent_id, :);
    sub = sortrows(sub, 'time_pc_sec_abs');

    if isempty(sub)
        series = struct('time', [], 'a0_corr', [], 'a2', []);
        return;
    end

    time_vals = double(sub.time_pc_sec_abs(:));
    a0_vals = correct_phase_discontinuity(double(sub.a0(:)));
    a2_vals = double(sub.a2(:));

    valid = isfinite(time_vals) & isfinite(a0_vals) & isfinite(a2_vals);
    time_vals = time_vals(valid);
    a0_vals = a0_vals(valid);
    a2_vals = a2_vals(valid);

    [time_vals, ia] = unique(time_vals, 'stable');
    a0_vals = a0_vals(ia);
    a2_vals = a2_vals(ia);

    series = struct('time', time_vals, 'a0_corr', a0_vals, 'a2', a2_vals);
end

function corrected_phase = correct_phase_discontinuity(phase_data)
    corrected_phase = phase_data(:);
    for i = 2:length(corrected_phase)
        diffv = corrected_phase(i) - corrected_phase(i - 1);
        if diffv < -128
            corrected_phase(i:end) = corrected_phase(i:end) + 256;
        elseif diffv > 128
            corrected_phase(i:end) = corrected_phase(i:end) - 256;
        end
    end
end

function T = correct_large_jump_matlab(T, threshold_sec, jump_sec)
    [G, ~] = findgroups(T.agent_id, T.chunk_id);

    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end

        [~, rel] = sort(T.time_pc_sec_abs(idx));
        idx = idx(rel);

        time_series = T.time_pc_sec_abs(idx);
        time_diff = [0; diff(time_series)];
        jump_idx = find(time_diff > threshold_sec);

        for j = 1:numel(jump_idx)
            fix_range = jump_idx(j):numel(time_series);
            T.time_pc_sec_abs(idx(fix_range)) = T.time_pc_sec_abs(idx(fix_range)) - jump_sec;
        end
    end
end

function T = correct_chunk_start_times_matlab(T, threshold_sec, jump_sec)
    [G, ~] = findgroups(T.agent_id, T.chunk_id);
    chunk_start = splitapply(@(x) min(x), T.time_pc_sec_abs, G);
    median_start = median(chunk_start, 'omitnan');

    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end

        start_time = min(T.time_pc_sec_abs(idx));
        if (start_time - median_start) > threshold_sec
            T.time_pc_sec_abs(idx) = T.time_pc_sec_abs(idx) - jump_sec;
        end
    end
end
