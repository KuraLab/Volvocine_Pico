function out = plot_phase_a2_early_late_nth_file(dirpath, file_index, segment_ratio, target_agents, thin_step, analysis_duration_sec, analysis_start_sec)
% Plot phase-a2 relationship over full time in the Nth CSV file.
%
% Usage:
%   plot_phase_a2_early_late_nth_file()
%   plot_phase_a2_early_late_nth_file(dirpath)
%   plot_phase_a2_early_late_nth_file(dirpath, file_index, segment_ratio, target_agents, thin_step, analysis_duration_sec, analysis_start_sec)
%
% Defaults:
%   dirpath = 'EstimateF/omega245'
%   file_index = 9
%   segment_ratio = []      % unused (kept for backward compatibility)
%   target_agents = []     % [] means all agents except 99
%   thin_step = 1.0        % keep about 100%% of NaN-separated segments
%                           % 0<thin_step<1: keep ratio, thin_step>=1: legacy step
%   analysis_duration_sec = 30  % process N seconds
%   analysis_start_sec = 0      % start time (seconds from beginning)
%
% Output:
%   out: struct with selected file path and plotting metadata.

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('EstimateF', 'Spring3/250');
    end
    if nargin < 2 || isempty(file_index)
        file_index = 5;
    end
    if nargin < 3 || isempty(segment_ratio)
        segment_ratio = [];
    end
    if nargin < 4
        target_agents = [];
    end
    if nargin < 5 || isempty(thin_step)
        thin_step = 0.7;
    end
    if nargin < 6 || isempty(analysis_duration_sec)
        analysis_duration_sec = 15;
    end
    if nargin < 7 || isempty(analysis_start_sec)
        analysis_start_sec = 5;
    end

    if ~isscalar(thin_step) || thin_step <= 0 || ~isfinite(thin_step)
        error('thin_step must be a positive scalar. Use 0.1 (~10%%) or 10 (legacy step).');
    end
    if ~isscalar(analysis_duration_sec) || analysis_duration_sec <= 0 || ~isfinite(analysis_duration_sec)
        error('analysis_duration_sec must be a positive scalar.');
    end
    if ~isscalar(analysis_start_sec) || analysis_start_sec < 0 || ~isfinite(analysis_start_sec)
        error('analysis_start_sec must be a non-negative scalar.');
    end
    thin_step_input = double(thin_step);
    if thin_step_input < 1
        segment_keep_ratio = min(1, thin_step_input);
        segment_thin_step = max(1, round(1 / segment_keep_ratio));
    else
        segment_thin_step = max(1, round(thin_step_input));
        segment_keep_ratio = min(1, 1 / segment_thin_step);
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
        has_agent = true;
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
    if isempty(target_agents)
        agents = all_agents(all_agents ~= 99);
        if isempty(agents)
            agents = all_agents;
        end
    else
        target_agents = target_agents(:).';
        agents = target_agents(ismember(target_agents, all_agents));
        if isempty(agents)
            error('No target_agents found in file. Available agents: %s', mat2str(all_agents));
        end
    end

    non99_agents = all_agents(all_agents ~= 99);
    if isempty(non99_agents)
        non99_agents = all_agents;
    end
    base_agent_id = min(non99_agents);

    base_sub = T(T.agent_id == base_agent_id, :);
    base_sub = sortrows(base_sub, 'time_pc_sec_abs');
    base_time = double(base_sub.time_pc_sec_abs(:));
    base_a0_corr = correct_phase_discontinuity(double(base_sub.a0(:)));
    [base_time, ia_base] = unique(base_time, 'stable');
    base_a0_corr = base_a0_corr(ia_base);

    n_agents = numel(agents);
    n_cols = min(3, n_agents);
    n_rows = ceil(n_agents / n_cols);

    fig_phase_sensor = figure('Color', 'w');
    tl_phase_sensor = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_phase_diff_time = figure('Color', 'w');
    tl_phase_diff_time = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_phase_sensor_delta = figure('Color', 'w');
    tl_phase_sensor_delta = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    stats = repmat(struct('agent_id', [], 'mean_phase', [], 'mean_phase_diff', [], ...
        'n_segments', [], 'segment_mean_phase_diff', [], 'segment_mean_time_sec', []), 1, n_agents);

    for i = 1:n_agents
        ax = nexttile(tl_phase_sensor, i); %#ok<LAXES>
        ax_diff = nexttile(tl_phase_diff_time, i); %#ok<LAXES>
        ax_delta = nexttile(tl_phase_sensor_delta, i); %#ok<LAXES>
        aid = agents(i);

        sub = T(T.agent_id == aid, :);
        sub = sortrows(sub, 'time_pc_sec_abs');

        t = double(sub.time_pc_sec_abs(:));
        a0_raw = double(sub.a0(:));
        a0_corr = correct_phase_discontinuity(a0_raw);
        phase = mod(a0_raw, 256) * (2*pi/256);
        sensor = double(sub.a2(:));

        valid = isfinite(t) & isfinite(phase) & isfinite(sensor) & isfinite(a0_corr);
        t = t(valid);
        a0_corr = a0_corr(valid);
        phase = phase(valid);
        sensor = sensor(valid);

        if isempty(t)
            text(ax, 0.5, 0.5, sprintf('agent %g: no valid data', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            text(ax_diff, 0.5, 0.5, sprintf('agent %g: no valid data', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_diff, 'off');
            text(ax_delta, 0.5, 0.5, sprintf('agent %g: no valid data', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_delta, 'off');
            continue;
        end

        t_zero = min(t);
        t_start_proc = t_zero + analysis_start_sec;
        t_end_proc = t_start_proc + analysis_duration_sec;
        proc_mask = (t >= t_start_proc) & (t <= t_end_proc);
        t = t(proc_mask);
        a0_corr = a0_corr(proc_mask);
        phase = phase(proc_mask);
        sensor = sensor(proc_mask);

        if isempty(t)
            text(ax, 0.5, 0.5, sprintf('agent %g: no data in window %.1f-%.1fs', ...
                aid, analysis_start_sec, analysis_start_sec + analysis_duration_sec), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            text(ax_diff, 0.5, 0.5, sprintf('agent %g: no data in window %.1f-%.1fs', ...
                aid, analysis_start_sec, analysis_start_sec + analysis_duration_sec), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_diff, 'off');
            text(ax_delta, 0.5, 0.5, sprintf('agent %g: no data in window %.1f-%.1fs', ...
                aid, analysis_start_sec, analysis_start_sec + analysis_duration_sec), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_delta, 'off');
            continue;
        end

        if aid == base_agent_id
            phase_diff = zeros(size(a0_corr));
        else
            if numel(base_time) < 2
                phase_diff = nan(size(a0_corr));
            else
                base_interp = interp1(base_time, base_a0_corr, t, 'linear', NaN);
                phase_raw = a0_corr - base_interp;
                phase_diff = mod(phase_raw + 128, 256) - 128;
                phase_diff = phase_diff * (2*pi/256);
            end
        end

        valid_diff = isfinite(phase_diff);
        t = t(valid_diff);
        phase = phase(valid_diff);
        sensor = sensor(valid_diff);
        phase_diff = phase_diff(valid_diff);

        if isempty(t)
            text(ax, 0.5, 0.5, sprintf('agent %g: no overlap with base agent %g', aid, base_agent_id), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            text(ax_diff, 0.5, 0.5, sprintf('agent %g: no overlap with base agent %g', aid, base_agent_id), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_diff, 'off');
            text(ax_delta, 0.5, 0.5, sprintf('agent %g: no overlap with base agent %g', aid, base_agent_id), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_delta, 'off');
            continue;
        end

        phase_jump_threshold = pi;
        t_rel = t - t(1);

        [time_all, phase_all, sensor_all, phase_diff_all] = break_phase_jumps_for_plot( ...
            t_rel, phase, sensor, phase_diff, phase_jump_threshold);

        [time_all, phase_all, sensor_all, phase_diff_all] = thin_nan_segments( ...
            time_all, phase_all, sensor_all, phase_diff_all, segment_keep_ratio);

        segments = split_nan_segments(time_all, phase_all, sensor_all, phase_diff_all);
        if isempty(segments)
            text(ax, 0.5, 0.5, sprintf('agent %g: no segments after thinning', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            text(ax_diff, 0.5, 0.5, sprintf('agent %g: no segments after thinning', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_diff, 'off');
            text(ax_delta, 0.5, 0.5, sprintf('agent %g: no segments after thinning', aid), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ax_delta, 'off');
            continue;
        end

        mean_phase = circular_mean_0_2pi(phase_all);
        mean_phase_diff = circular_mean_pi(phase_diff_all);

        n_segments = numel(segments);
        segment_mean_phase_diff = nan(1, n_segments);
        segment_mean_time = nan(1, n_segments);

        valid_time_all = time_all(isfinite(time_all));
        tmin_color = min(valid_time_all);
        tmax_color = max(valid_time_all);
        if ~isfinite(tmin_color) || ~isfinite(tmax_color)
            tmin_color = 0;
            tmax_color = 1;
        elseif tmax_color <= tmin_color
            tmax_color = tmin_color + eps;
        end

        time_cmap = turbo(256);
        colormap(ax, time_cmap);
        caxis(ax, [tmin_color, tmax_color]);
        colormap(ax_diff, time_cmap);
        caxis(ax_diff, [tmin_color, tmax_color]);
        colormap(ax_delta, time_cmap);
        caxis(ax_delta, [tmin_color, tmax_color]);

        hold(ax, 'on');
        hold(ax_diff, 'on');
        hold(ax_delta, 'on');
        for s = 1:n_segments
            seg = segments(s);
            segment_mean_phase_diff(s) = circular_mean_pi(seg.phase_diff);
            segment_mean_time(s) = mean(seg.time, 'omitnan');

            plot_time_gradient_line(ax, seg.phase, seg.sensor, seg.time, 1.2);
            plot_time_gradient_line(ax_diff, seg.time, seg.phase_diff, seg.time, 1.2);
        end

        ref_seg = select_latest_complete_segment(segments);
        yline(ax_delta, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
        for s = 1:n_segments
            seg = segments(s);
            [phase_delta, sensor_delta, time_delta] = compute_segment_sensor_difference(seg, ref_seg);
            if isempty(phase_delta)
                continue;
            end
            plot_time_gradient_line(ax_delta, phase_delta, sensor_delta, time_delta, 1.2);
        end
        hold(ax, 'off');
        hold(ax_diff, 'off');
        hold(ax_delta, 'off');

        xlim(ax, [0, 2*pi]);
        xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
        xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        xlabel(ax, 'Phase (rad)');
        ylabel(ax, 'a2');
        title(ax, sprintf('agent %g vs base %g (Nseg=%d, mean Δ\phi=%.2f rad)', ...
            aid, base_agent_id, n_segments, mean_phase_diff));
        grid(ax, 'on');
        cb = colorbar(ax);
        cb.Label.String = 'Time (s)';

        xlabel(ax_diff, 'Time (s)');
        ylabel(ax_diff, '\Delta\phi (rad)');
        xlim(ax_diff, [tmin_color, tmax_color]);
        ylim(ax_diff, [-pi, pi]);
        yticks(ax_diff, [-pi, 0, pi]);
        yticklabels(ax_diff, {'-\pi', '0', '\pi'});
        title(ax_diff, sprintf('agent %g: \Delta\phi_{%d-%d}(t)', aid, aid, base_agent_id));
        grid(ax_diff, 'on');
        cb_diff = colorbar(ax_diff);
        cb_diff.Label.String = 'Time (s)';

        xlim(ax_delta, [0, 2*pi]);
        xticks(ax_delta, [0, pi/2, pi, 3*pi/2, 2*pi]);
        xticklabels(ax_delta, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        xlabel(ax_delta, 'Phase (rad)');
        ylabel(ax_delta, '\Delta a_2 vs latest');
        title(ax_delta, sprintf('agent %g: a2 - latest complete segment', aid));
        grid(ax_delta, 'on');
        cb_delta = colorbar(ax_delta);
        cb_delta.Label.String = 'Time (s)';

        stats(i).agent_id = aid;
        stats(i).mean_phase = mean_phase;
        stats(i).mean_phase_diff = mean_phase_diff;
        stats(i).n_segments = n_segments;
        stats(i).segment_mean_phase_diff = segment_mean_phase_diff;
        stats(i).segment_mean_time_sec = segment_mean_time;
    end

    [~, filename, ext] = fileparts(csv_path);
    sgtitle(tl_phase_sensor, sprintf('%s%s (index %d): phase vs a2, base=%d, window %.0f-%.0fs, segment-keep~%.0f%%', ...
        filename, ext, file_index, base_agent_id, analysis_start_sec, analysis_start_sec + analysis_duration_sec, segment_keep_ratio * 100));
    sgtitle(tl_phase_diff_time, sprintf('%s%s (index %d): phase-diff time series, base=%d, window %.0f-%.0fs, segment-keep~%.0f%%', ...
        filename, ext, file_index, base_agent_id, analysis_start_sec, analysis_start_sec + analysis_duration_sec, segment_keep_ratio * 100));
    sgtitle(tl_phase_sensor_delta, sprintf('%s%s (index %d): a2 difference vs latest segment, window %.0f-%.0fs, segment-keep~%.0f%%', ...
        filename, ext, file_index, analysis_start_sec, analysis_start_sec + analysis_duration_sec, segment_keep_ratio * 100));

    out = struct();
    out.dirpath = dirpath;
    out.file_index = file_index;
    out.csv_path = csv_path;
    out.segment_ratio = segment_ratio;
    out.use_all_time = true;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.base_agent_id = base_agent_id;
    out.figure_phase_sensor = fig_phase_sensor;
    out.figure_phase_diff_time = fig_phase_diff_time;
    out.figure_phase_sensor_delta = fig_phase_sensor_delta;
    out.thin_step = thin_step_input;
    out.segment_thin_step = segment_thin_step;
    out.segment_keep_ratio = segment_keep_ratio;
    out.agents = agents;
    out.stats = stats;
end

function [time_out, phase_out, sensor_out, phase_diff_out] = break_phase_jumps_for_plot(time_in, phase_in, sensor_in, phase_diff_in, jump_threshold)
    time_out = time_in(:);
    phase_out = phase_in(:);
    sensor_out = sensor_in(:);
    phase_diff_out = phase_diff_in(:);

    if numel(phase_out) < 2
        return;
    end

    jump_idx = find(abs(diff(phase_out)) > jump_threshold) + 1;
    time_out(jump_idx) = NaN;
    phase_out(jump_idx) = NaN;
    sensor_out(jump_idx) = NaN;
    phase_diff_out(jump_idx) = NaN;
end

function [t_out, x_out, y_out, d_out] = thin_nan_segments(t_in, x_in, y_in, d_in, keep_ratio)
    t = t_in(:);
    x = x_in(:);
    y = y_in(:);
    d = d_in(:);

    if keep_ratio >= 1
        t_out = t;
        x_out = x;
        y_out = y;
        d_out = d;
        return;
    end

    finite_mask = isfinite(t) & isfinite(x) & isfinite(y) & isfinite(d);
    if ~any(finite_mask)
        t_out = t;
        x_out = x;
        y_out = y;
        d_out = d;
        return;
    end

    seg_start = find(finite_mask & [true; ~finite_mask(1:end-1)]);
    seg_end = find(finite_mask & [~finite_mask(2:end); true]);

    n_seg = numel(seg_start);
    n_keep = max(1, round(n_seg * keep_ratio));
    keep_idx = unique(round(linspace(1, n_seg, n_keep)));

    keep_segment = false(n_seg, 1);
    keep_segment(keep_idx) = true;

    keep_points = false(size(finite_mask));
    for k = 1:numel(seg_start)
        if keep_segment(k)
            keep_points(seg_start(k):seg_end(k)) = true;
        end
    end

    drop_points = finite_mask & ~keep_points;
    t(drop_points) = NaN;
    x(drop_points) = NaN;
    y(drop_points) = NaN;
    d(drop_points) = NaN;

    t_out = t;
    x_out = x;
    y_out = y;
    d_out = d;
end

function segments = split_nan_segments(time_in, phase_in, sensor_in, phase_diff_in)
    t = time_in(:);
    p = phase_in(:);
    s = sensor_in(:);
    d = phase_diff_in(:);

    finite_mask = isfinite(t) & isfinite(p) & isfinite(s) & isfinite(d);
    if ~any(finite_mask)
        segments = struct('time', {}, 'phase', {}, 'sensor', {}, 'phase_diff', {}, 'ends_with_nan', {});
        return;
    end

    seg_start = find(finite_mask & [true; ~finite_mask(1:end-1)]);
    seg_end = find(finite_mask & [~finite_mask(2:end); true]);

    n_seg = numel(seg_start);
    segments = repmat(struct('time', [], 'phase', [], 'sensor', [], 'phase_diff', [], 'ends_with_nan', false), 1, n_seg);
    for k = 1:n_seg
        idx = seg_start(k):seg_end(k);
        segments(k).time = t(idx);
        segments(k).phase = p(idx);
        segments(k).sensor = s(idx);
        segments(k).phase_diff = d(idx);
        segments(k).ends_with_nan = seg_end(k) < numel(finite_mask) && ~finite_mask(seg_end(k) + 1);
    end
end

function ref_seg = select_latest_complete_segment(segments)
    if isempty(segments)
        error('No segments available for reference selection.');
    end

    complete_mask = [segments.ends_with_nan];
    if any(complete_mask)
        ref_seg = segments(find(complete_mask, 1, 'last'));
    else
        ref_seg = segments(end);
    end
end

function plot_time_gradient_line(ax, x_in, y_in, t_in, line_width)
    x = x_in(:);
    y = y_in(:);
    t = t_in(:);

    if numel(x) < 2
        scatter(ax, x, y, 10, t, 'filled');
        return;
    end

    surface(ax, [x x], [y y], zeros(numel(x), 2), [t t], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', line_width);
end

function [phase_common, sensor_delta, time_common] = compute_segment_sensor_difference(seg, ref_seg)
    [phase_seg, sensor_seg, time_seg] = prepare_segment_curve(seg.phase, seg.sensor, seg.time);
    [phase_ref, sensor_ref, ~] = prepare_segment_curve(ref_seg.phase, ref_seg.sensor, ref_seg.time);

    if isempty(phase_seg) || isempty(phase_ref)
        phase_common = [];
        sensor_delta = [];
        time_common = [];
        return;
    end

    phase_min = max(min(phase_seg), min(phase_ref));
    phase_max = min(max(phase_seg), max(phase_ref));
    if ~isfinite(phase_min) || ~isfinite(phase_max) || phase_max <= phase_min
        phase_common = [];
        sensor_delta = [];
        time_common = [];
        return;
    end

    n_grid = max(300, 4 * max(numel(phase_seg), numel(phase_ref)));
    n_grid = min(n_grid, 4000);
    phase_common = linspace(phase_min, phase_max, n_grid).';

    sensor_seg_common = interp1(phase_seg, sensor_seg, phase_common, 'linear');
    sensor_ref_common = interp1(phase_ref, sensor_ref, phase_common, 'linear');
    time_common = interp1(phase_seg, time_seg, phase_common, 'linear');
    valid = isfinite(sensor_seg_common) & isfinite(sensor_ref_common) & isfinite(time_common);
    phase_common = phase_common(valid);
    sensor_delta = sensor_seg_common(valid) - sensor_ref_common(valid);
    time_common = time_common(valid);

    if isempty(phase_common)
        sensor_delta = [];
        time_common = [];
    end
end

function [phase_u, sensor_u, time_u] = prepare_segment_curve(phase_in, sensor_in, time_in)
    phase_vals = phase_in(:);
    sensor_vals = sensor_in(:);
    time_vals = time_in(:);

    valid = isfinite(phase_vals) & isfinite(sensor_vals) & isfinite(time_vals);
    phase_vals = phase_vals(valid);
    sensor_vals = sensor_vals(valid);
    time_vals = time_vals(valid);

    if isempty(phase_vals)
        phase_u = [];
        sensor_u = [];
        time_u = [];
        return;
    end

    [phase_vals, order] = sort(phase_vals);
    sensor_vals = sensor_vals(order);
    time_vals = time_vals(order);

    [G, phase_u] = findgroups(phase_vals);
    sensor_u = splitapply(@mean, sensor_vals, G);
    time_u = splitapply(@mean, time_vals, G);
end

function mu = circular_mean_0_2pi(phase_vals)
    vals = phase_vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        mu = NaN;
        return;
    end

    mu = atan2(mean(sin(vals)), mean(cos(vals)));
    mu = mod(mu, 2*pi);
end

function mu = circular_mean_pi(phase_vals)
    vals = phase_vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        mu = NaN;
        return;
    end

    mu = atan2(mean(sin(vals)), mean(cos(vals)));
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
