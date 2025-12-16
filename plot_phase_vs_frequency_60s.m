function plot_phase_vs_frequency_60s(rootdir, t_sample, plot_duration, apply_filter, filter_window_size, sample_window)
% Plot mean phase-difference at a fixed time across frequency folders
%
% Usage:
%   plot_phase_vs_frequency_60s()
%   plot_phase_vs_frequency_60s(rootdir, t_sample)
%
% Defaults:
%   rootdir = 'EstimateF'
%   t_sample = 60 (seconds after start)
%   sample_window = [t_sample - 10, t_sample]
%   plot_duration = 120
%   apply_filter = true
%   filter_window_size = 1

    if nargin < 1 || isempty(rootdir)
        rootdir = 'EstimateF';
    end
    if nargin < 2 || isempty(t_sample)
        t_sample = 60;
    end
    if nargin < 6 || isempty(sample_window)
        sample_window = [max(0, t_sample - 10), t_sample];
    end
    if nargin < 3 || isempty(plot_duration)
        plot_duration = 120;
    end
    if nargin < 4 || isempty(apply_filter)
        apply_filter = true;
    end
    if nargin < 5 || isempty(filter_window_size)
        filter_window_size = 1;
    end

    if ~isfolder(rootdir)
        error('Root directory not found: %s', rootdir);
    end

    root_dirs = resolve_root_dirs(rootdir);
    if isempty(root_dirs)
        error('No Spring directories found under %s', rootdir);
    end

    % arrays for plotting (may contain multiple points per frequency when clustered)
    freq_plot = [];
    mean_plot = [];
    stderr_plot = [];
    dataset_idx = [];
    dataset_labels = cell(numel(root_dirs),1);

    % constants for overflow correction (same as other scripts)
    T_OVERFLOW = 2^32 / 1e6;
    T_TOL = 5.0;
    threshold_sec = T_OVERFLOW - T_TOL;
    jump_sec = T_OVERFLOW;
    CLUSTER_VAR_THRESHOLD = 0.35;

    for r = 1:numel(root_dirs)
        root_path = root_dirs{r};
        [~, dataset_labels{r}] = fileparts(root_path);

        [freq_r, mean_r, stderr_r] = collect_phase_points_for_root( ...
            root_path, sample_window, plot_duration, apply_filter, filter_window_size, ...
            threshold_sec, jump_sec, T_OVERFLOW, CLUSTER_VAR_THRESHOLD);

        if isempty(freq_r)
            warning('No valid frequency data under %s', root_path);
            continue;
        end

        freq_plot = [freq_plot, freq_r]; %#ok<AGROW>
        mean_plot = [mean_plot, mean_r]; %#ok<AGROW>
        stderr_plot = [stderr_plot, stderr_r]; %#ok<AGROW>
        dataset_idx = [dataset_idx, r * ones(1, numel(freq_r))]; %#ok<AGROW>
    end

    if isempty(freq_plot)
        error('No valid frequency folders with CSV files found under %s', rootdir);
    end

    % sort by frequency for nicer plotting
    [freq_plot, I] = sort(freq_plot);
    mean_plot = mean_plot(I);
    stderr_plot = stderr_plot(I);
    dataset_idx = dataset_idx(I);

    % Plot (phase on x-axis, frequency on y-axis)
    figure;
    hold on;
    colors = lines(max(1, numel(root_dirs)));
    marker_size = 40;
    legend_handles = gobjects(0);
    legend_labels_drawn = {};
    for d = 1:numel(root_dirs)
        mask = dataset_idx == d;
        if ~any(mask)
            continue;
        end
        xvals = mean_plot(mask);
        yvals = freq_plot(mask);
        h = scatter(xvals, yvals, marker_size, 'MarkerFaceColor', colors(d,:), ...
            'MarkerEdgeColor', colors(d,:), 'DisplayName', dataset_labels{d});
        legend_handles(end+1) = h; %#ok<AGROW>
        legend_labels_drawn{end+1} = dataset_labels{d}; %#ok<AGROW>
    end
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels_drawn, 'Location', 'best');
    end
    ylabel('Frequency (folder name)');
    xlabel('Mean phase difference (rad)');
    grid on;
    xlim([-pi, pi]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    hold off;
    tuneFigure;
end

function root_dirs = resolve_root_dirs(root_input)
    if iscell(root_input)
        candidates = root_input;
    else
        candidates = {root_input};
    end

    root_dirs = {};
    for i = 1:numel(candidates)
        cand = candidates{i};
        if ~(ischar(cand) || isstring(cand))
            continue;
        end
        cand = char(cand);
        if ~isfolder(cand)
            warning('Skipping %s: not a folder', cand);
            continue;
        end
        entries = dir(cand);
        entries = entries([entries.isdir]);
        entries = entries(~ismember({entries.name},{'.','..'}));
        names = {entries.name};
        numeric_mask = ~isnan(str2double(names));
        if any(numeric_mask)
            root_dirs{end+1} = cand; %#ok<AGROW>
        else
            spring_mask = strncmpi(names, 'Spring', numel('Spring'));
            spring_entries = entries(spring_mask);
            for s = 1:numel(spring_entries)
                root_dirs{end+1} = fullfile(cand, spring_entries(s).name); %#ok<AGROW>
            end
        end
    end
    root_dirs = unique(root_dirs, 'stable');
end

function [freq_vals, mean_vals, stderr_vals] = collect_phase_points_for_root(root_folder, sample_window, plot_duration, apply_filter, filter_window_size, threshold_sec, jump_sec, T_OVERFLOW, cluster_var_threshold)
    freq_vals = [];
    mean_vals = [];
    stderr_vals = [];

    if numel(sample_window) ~= 2
        error('sample_window must contain exactly two elements [t_start, t_end].');
    end
    sample_window = sort(sample_window(:));
    window_start = sample_window(1);
    window_end = sample_window(2);
    if window_end <= window_start
        error('sample_window end must be greater than start.');
    end

    D = dir(root_folder);
    D = D([D.isdir]);
    D = D(~ismember({D.name},{'.','..'}));
    freq_mask = ~isnan(str2double({D.name}));
    freq_dirs = D(freq_mask);

    for s = 1:numel(freq_dirs)
        freq = str2double(freq_dirs(s).name);
        folder = fullfile(root_folder, freq_dirs(s).name);
        csvs = dir(fullfile(folder, '*.csv'));
        if isempty(csvs)
            continue;
        end

        per_file_means = nan(numel(csvs),1);
        all_samples = [];
        file_tables = cell(numel(csvs),1);
        agent_sets = cell(numel(csvs),1);

        for f = 1:numel(csvs)
            fp = fullfile(folder, csvs(f).name);
            try
                T = readtable(fp);
            catch
                warning('Failed to read %s', fp);
                continue;
            end
            if ~all(ismember({'time_pc_sec_abs','a0','agent_id'}, T.Properties.VariableNames))
                warning('Skipping %s: missing required columns', fp);
                continue;
            end
            if ~ismember('chunk_id', T.Properties.VariableNames)
                T.chunk_id = ones(height(T),1);
            end

            T = correct_large_jump_matlab(T, threshold_sec, jump_sec);
            T = correct_chunk_start_times_matlab(T, 4000.0, T_OVERFLOW);

            file_tables{f} = T;
            agent_sets{f} = unique(T.agent_id);
        end

        valid_idx = ~cellfun(@isempty, file_tables);
        file_tables = file_tables(valid_idx);
        agent_sets = agent_sets(valid_idx);
        per_file_means = per_file_means(valid_idx);

        if isempty(file_tables)
            continue;
        end

        common_agents = agent_sets{1};
        for f = 2:numel(agent_sets)
            common_agents = intersect(common_agents, agent_sets{f});
        end
        if isempty(common_agents)
            common_agents = agent_sets{1};
            warning('No common agents in %s; using agents from first file.', folder);
        end
        base_agent = min(common_agents);

        for f = 1:numel(file_tables)
            T = file_tables{f};
            series_struct = compute_phase_series_for_file(T, base_agent, 0, plot_duration, 1, apply_filter, filter_window_size, [], 1, 1);

            all_agents = find(~cellfun(@isempty,{series_struct.time}));
            if isempty(all_agents)
                continue;
            end

            sampled_vals = [];
            for ai = all_agents
                if ai == base_agent
                    continue; % skip base agent (always zero phase)
                end
                times = series_struct(ai).time;
                phases = series_struct(ai).phase;
                if isempty(times) || isempty(phases)
                    continue;
                end
                mask = (times >= window_start) & (times <= window_end);
                if ~any(mask)
                    continue;
                end
                window_phases = phases(mask);
                window_phases = window_phases(~isnan(window_phases));
                if isempty(window_phases)
                    continue;
                end
                agent_mean = atan2(mean(sin(window_phases)), mean(cos(window_phases)));
                sampled_vals(end+1) = agent_mean; %#ok<AGROW>
            end

            if isempty(sampled_vals)
                per_file_means(f) = NaN;
            else
                per_file_means(f) = mean(sampled_vals);
                all_samples = [all_samples, sampled_vals]; %#ok<AGROW>
            end
        end

        valid = ~isnan(per_file_means);
        if ~any(valid)
            continue;
        end

        if isempty(all_samples)
            continue;
        end

        Nsamples = numel(all_samples);
        R = abs(mean(exp(1i * all_samples(:))));
        circ_var = 1 - R;

        if circ_var > cluster_var_threshold && Nsamples >= 6
            k = 2;
            pts = [cos(all_samples(:)), sin(all_samples(:))];
            try
                opts = statset('MaxIter',500);
                [idx, ~] = kmeans(pts, k, 'Replicates', 5, 'Options', opts);
            catch
                idx = ones(size(all_samples(:)));
            end

            offsets = linspace(-0.12, 0.12, k);
            for c = 1:k
                ths = all_samples(idx == c);
                if isempty(ths)
                    continue;
                end
                mmean = atan2(mean(sin(ths)), mean(cos(ths)));
                Rcl = abs(mean(exp(1i * ths)));
                circ_std = sqrt(max(0, -2 * log(max(Rcl, eps))));
                stderr = circ_std / sqrt(max(1, numel(ths)));

                freq_vals(end+1) = freq + offsets(c); %#ok<AGROW>
                mean_vals(end+1) = mmean; %#ok<AGROW>
                stderr_vals(end+1) = stderr; %#ok<AGROW>
            end
        else
            all_vals = all_samples(:);
            mmean = atan2(mean(sin(all_vals)), mean(cos(all_vals)));
            Rall = abs(mean(exp(1i * all_vals)));
            circ_std = sqrt(max(0, -2 * log(max(Rall, eps))));
            stderr = circ_std / sqrt(max(1, numel(all_vals)));

            freq_vals(end+1) = freq; %#ok<AGROW>
            mean_vals(end+1) = mmean; %#ok<AGROW>
            stderr_vals(end+1) = stderr; %#ok<AGROW>
        end
    end

    freq_vals = freq_vals(:).';
    mean_vals = mean_vals(:).';
    stderr_vals = stderr_vals(:).';
end

% ----------------------- helper functions (copied/adapted) -----------------
function series_struct = compute_phase_series_for_file(df_all, base_agent_id, n_seconds_to_cut, plot_duration, allow_missing_agents, apply_filter, filter_window_size, max_agent_id, n_sync, m_sync)
    if nargin < 5 || isempty(allow_missing_agents)
        allow_missing_agents = 1;
    end
    if nargin < 6 || isempty(apply_filter)
        apply_filter = true;
    end
    if nargin < 7 || isempty(filter_window_size)
        filter_window_size = 1;
    end
    if nargin < 8 || isempty(max_agent_id)
        max_agent_id = max(df_all.agent_id);
    end
    if nargin < 9 || isempty(n_sync)
        n_sync = 2;
    end
    if nargin < 10 || isempty(m_sync)
        m_sync = 1;
    end

    series_struct = repmat(struct('time', [], 'phase', []), 1, max_agent_id);

    if isempty(df_all)
        return;
    end

    df_main = df_all(df_all.agent_id ~= 99, :);
    if isempty(df_main)
        return;
    end

    min_time = min(df_main.time_pc_sec_abs);
    max_time = max(df_main.time_pc_sec_abs);
    agents = unique(df_main.agent_id, 'sorted').';

    agent_ranges = zeros(length(agents), 2);
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_main(df_main.agent_id == agent_id, :);
        agent_ranges(i,1) = min(sub.time_pc_sec_abs);
        agent_ranges(i,2) = max(sub.time_pc_sec_abs);
    end

    if isempty(base_agent_id) || ~ismember(base_agent_id, agents)
        base_agent_id = min(agents);
    end

    start_time_abs = min_time + n_seconds_to_cut;
    new_time_series = (start_time_abs:0.01:max_time) - start_time_abs;
    if isempty(new_time_series)
        return;
    end

    valid_counts = zeros(size(new_time_series));
    for t_idx = 1:length(new_time_series)
        t_abs = new_time_series(t_idx) + start_time_abs;
        valid_counts(t_idx) = sum(agent_ranges(:,1) <= t_abs & agent_ranges(:,2) >= t_abs);
    end

    min_valid = length(agents) - allow_missing_agents;
    valid_idx = find(valid_counts >= min_valid);
    if isempty(valid_idx)
        return;
    end

    new_time_series = new_time_series(valid_idx);

    max_allowed_time = plot_duration - n_seconds_to_cut;
    time_mask = new_time_series <= max_allowed_time;
    new_time_series = new_time_series(time_mask);
    if isempty(new_time_series)
        return;
    end

    interpolated_data = struct();
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_main(df_main.agent_id == agent_id, :);
        sub = sortrows(sub, 'time_pc_sec_abs');
        sub.a0 = correct_phase_discontinuity(sub.a0);
        [~, ia] = unique(sub.time_pc_sec_abs);
        sub = sub(ia, :);
        t_min_agent = min(sub.time_pc_sec_abs) - start_time_abs;
        t_max_agent = max(sub.time_pc_sec_abs) - start_time_abs;
        valid_mask = (new_time_series >= t_min_agent) & (new_time_series <= t_max_agent);
        interp_a0 = nan(size(new_time_series));
        interp_a0(valid_mask) = interp1(sub.time_pc_sec_abs - start_time_abs, sub.a0, new_time_series(valid_mask), 'linear', 'extrap');
        if apply_filter
            interp_a0 = movmean(interp_a0, filter_window_size, 'omitnan');
        end
        interpolated_data(agent_id).a0 = interp_a0;
    end

    base_agent_a0 = interpolated_data(base_agent_id).a0;

    for i = 1:length(agents)
        agent_id = agents(i);
        phase_raw = n_sync * interpolated_data(agent_id).a0 - m_sync * base_agent_a0;
        phase_diff = mod(phase_raw + 128, 256) - 128;
        phase_diff = phase_diff * (2*pi/256);
        phase_diff_with_nan = phase_diff;
        for j = 2:length(phase_diff_with_nan)
            if isnan(phase_diff_with_nan(j)) || isnan(phase_diff_with_nan(j-1))
                continue;
            end
            if abs(phase_diff_with_nan(j) - phase_diff_with_nan(j-1)) > pi
                phase_diff_with_nan(j) = NaN;
            end
        end
        series_struct(agent_id).time = new_time_series;
        if agent_id == base_agent_id
            series_struct(agent_id).phase = zeros(size(new_time_series));
        else
            series_struct(agent_id).phase = phase_diff_with_nan;
        end
    end
end

function corrected_phase = correct_phase_discontinuity(phase_data)
    corrected_phase = phase_data;
    for i = 2:length(corrected_phase)
        diffv = corrected_phase(i) - corrected_phase(i - 1);
        if diffv < -128
            corrected_phase(i:end) = corrected_phase(i:end) + 256;
        elseif diffv > 128
            corrected_phase(i:end) = corrected_phase(i:end) - 256;
        end
    end
end

function df_all = correct_large_jump_matlab(df_all, threshold_sec, jump_sec)
    [G, ~] = findgroups(df_all.agent_id, df_all.chunk_id);
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end
        [~, sorted_idx_rel] = sort(df_all.time_pc_sec_abs(idx));
        idx = idx(sorted_idx_rel);
        time_series = df_all.time_pc_sec_abs(idx);
        time_diff = [0; diff(time_series)];
        jump_idx = find(time_diff > threshold_sec);
        for j = 1:length(jump_idx)
            fix_range = jump_idx(j):length(time_series);
            df_all.time_pc_sec_abs(idx(fix_range)) = df_all.time_pc_sec_abs(idx(fix_range)) - jump_sec;
        end
    end
end

function df_all = correct_chunk_start_times_matlab(df_all, threshold_sec, jump_sec)
    [G, chunk_keys] = findgroups(df_all.agent_id, df_all.chunk_id);
    chunk_start = splitapply(@(x) min(x), df_all.time_pc_sec_abs, G);
    median_start = median(chunk_start);
    has_chunk_id = size(chunk_keys, 2) >= 2;
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end
        start_time = df_all.time_pc_sec_abs(idx(1));
        if start_time - median_start > threshold_sec
            df_all.time_pc_sec_abs(idx) = df_all.time_pc_sec_abs(idx) - jump_sec;
        end
    end
end
