function out = plot_phase_pair_a2_3d_nth_file(dirpath, file_index, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec)
% Plot a 3D trajectory with
%   x = phase of agent 1,
%   y = phase of agent 2,
%   z = a2 of z_agent_id,
% colored by time.
%
% Usage:
%   plot_phase_pair_a2_3d_nth_file()
%   plot_phase_pair_a2_3d_nth_file(dirpath)
%   plot_phase_pair_a2_3d_nth_file(dirpath, file_index, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec)
%
% Defaults:
%   dirpath = 'EstimateF/Spring3/250'
%   file_index = 5
%   phase_agent_ids = first two non-99 agents in the file
%   z_agent_id = phase_agent_ids(2)
%   analysis_duration_sec = 15
%   analysis_start_sec = 5
%
% Output:
%   out: struct with selected file path, used agents, time window, and figure handle.

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('EstimateF', 'Spring3/250');
    end
    if nargin < 2 || isempty(file_index)
        file_index = 5;
    end
    if nargin < 3
        phase_agent_ids = [];
    end
    if nargin < 4
        z_agent_id = [];
    end
    if nargin < 5 || isempty(analysis_duration_sec)
        analysis_duration_sec = 15;
    end
    if nargin < 6 || isempty(analysis_start_sec)
        analysis_start_sec = 5;
    end

    sample_dt = 0.01;
    jump_threshold = pi;

    if ~isscalar(analysis_duration_sec) || analysis_duration_sec <= 0 || ~isfinite(analysis_duration_sec)
        error('analysis_duration_sec must be a positive scalar.');
    end
    if ~isscalar(analysis_start_sec) || analysis_start_sec < 0 || ~isfinite(analysis_start_sec)
        error('analysis_start_sec must be a non-negative scalar.');
    end

    if ~isfolder(dirpath)
        error('Directory not found: %s', dirpath);
    end

    csv_path = pick_nth_csv(dirpath, file_index);
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
    non99_agents = all_agents(all_agents ~= 99);
    if isempty(non99_agents)
        non99_agents = all_agents;
    end

    if isempty(phase_agent_ids)
        if numel(non99_agents) < 2
            error('At least two agents are required. Found: %s', mat2str(non99_agents));
        end
        phase_agent_ids = non99_agents(1:2);
    else
        phase_agent_ids = phase_agent_ids(:).';
        if numel(phase_agent_ids) ~= 2
            error('phase_agent_ids must contain exactly two agent IDs.');
        end
        if phase_agent_ids(1) == phase_agent_ids(2)
            error('phase_agent_ids must specify two distinct agents.');
        end
        if ~all(ismember(phase_agent_ids, all_agents))
            error('Some phase_agent_ids are not present in the file. Available agents: %s', mat2str(all_agents));
        end
    end

    if isempty(z_agent_id)
        z_agent_id = phase_agent_ids(2);
    end
    if ~ismember(z_agent_id, all_agents)
        error('z_agent_id %d is not present in the file. Available agents: %s', z_agent_id, mat2str(all_agents));
    end

    needed_agents = unique([phase_agent_ids(:); z_agent_id]);
    series_by_agent = struct();
    overlap_start = -inf;
    overlap_end = inf;

    for k = 1:numel(needed_agents)
        aid = needed_agents(k);
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

    phase_1 = mod(a0_1, 256) * (2*pi/256);
    phase_2 = mod(a0_2, 256) * (2*pi/256);

    [time_plot, phase_1_plot, phase_2_plot, a2_plot] = break_wrapped_phase_jumps_3d( ...
        time_rel, phase_1, phase_2, a2_z, jump_threshold);

    segments = split_nan_segments_3d(time_plot, phase_1_plot, phase_2_plot, a2_plot);
    if isempty(segments)
        error('No valid 3D segments were available for plotting.');
    end

    tmin_color = min(time_rel);
    tmax_color = max(time_rel);
    if tmax_color <= tmin_color
        tmax_color = tmin_color + eps;
    end

    fig = figure('Color', 'w');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    cmap = turbo(256);
    colormap(ax, cmap);
    caxis(ax, [tmin_color, tmax_color]);

    for s = 1:numel(segments)
        seg = segments(s);
        plot_time_gradient_line3(ax, seg.phase1, seg.phase2, seg.a2, seg.time, 1.4);
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
    title(ax, sprintf('3D phase-phase-a2 trajectory: agents [%d, %d], z=a2(agent %d)', ...
        phase_agent_ids(1), phase_agent_ids(2), z_agent_id));
    grid(ax, 'on');
    view(ax, 3);
    box(ax, 'on');

    cb = colorbar(ax);
    cb.Label.String = 'Time (s)';

    [~, filename, ext] = fileparts(csv_path);
    sgtitle(fig, sprintf('%s%s (index %d): window %.0f-%.0fs', ...
        filename, ext, file_index, analysis_start_sec, analysis_start_sec + analysis_duration_sec));

    out = struct();
    out.dirpath = dirpath;
    out.file_index = file_index;
    out.csv_path = csv_path;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.window_start_abs = window_start_abs;
    out.window_end_abs = window_end_abs;
    out.figure = fig;
    out.n_segments = numel(segments);
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

function [time_out, phase1_out, phase2_out, z_out] = break_wrapped_phase_jumps_3d(time_in, phase1_in, phase2_in, z_in, jump_threshold)
    time_out = time_in(:);
    phase1_out = phase1_in(:);
    phase2_out = phase2_in(:);
    z_out = z_in(:);

    if numel(time_out) < 2
        return;
    end

    jump_idx = find((abs(diff(phase1_out)) > jump_threshold) | (abs(diff(phase2_out)) > jump_threshold)) + 1;
    time_out(jump_idx) = NaN;
    phase1_out(jump_idx) = NaN;
    phase2_out(jump_idx) = NaN;
    z_out(jump_idx) = NaN;
end

function segments = split_nan_segments_3d(time_in, x_in, y_in, z_in)
    t = time_in(:);
    x = x_in(:);
    y = y_in(:);
    z = z_in(:);

    finite_mask = isfinite(t) & isfinite(x) & isfinite(y) & isfinite(z);
    if ~any(finite_mask)
        segments = struct('time', {}, 'phase1', {}, 'phase2', {}, 'a2', {});
        return;
    end

    seg_start = find(finite_mask & [true; ~finite_mask(1:end-1)]);
    seg_end = find(finite_mask & [~finite_mask(2:end); true]);

    n_seg = numel(seg_start);
    segments = repmat(struct('time', [], 'phase1', [], 'phase2', [], 'a2', []), 1, n_seg);
    for k = 1:n_seg
        idx = seg_start(k):seg_end(k);
        segments(k).time = t(idx);
        segments(k).phase1 = x(idx);
        segments(k).phase2 = y(idx);
        segments(k).a2 = z(idx);
    end
end

function plot_time_gradient_line3(ax, x_in, y_in, z_in, t_in, line_width)
    x = x_in(:);
    y = y_in(:);
    z = z_in(:);
    t = t_in(:);

    if numel(x) < 2
        scatter3(ax, x, y, z, 14, t, 'filled');
        return;
    end

    surface(ax, [x x], [y y], [z z], [t t], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', line_width);
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

function csvPath = pick_nth_csv(parentDir, indexN)
    files = dir(fullfile(parentDir, '*.csv'));
    if isempty(files)
        error('No CSV files found in %s', parentDir);
    end

    names = sort({files.name});
    if indexN < 1 || indexN > numel(names)
        error('file_index out of range for %s: 1..%d', parentDir, numel(names));
    end

    csvPath = fullfile(parentDir, names{indexN});
    fprintf('[INFO] Selected CSV #%d: %s\n', indexN, csvPath);
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
