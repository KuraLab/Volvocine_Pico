function varargout = plot_phase_evolution_omega245(dirpath, n_seconds_to_cut, plot_duration, apply_filter, filter_window_size, do_save_figure, n_sync, m_sync, sample_window)
% Overlay phase-relationship time evolutions from all files in a directory
%
% Usage:
%   plot_phase_evolution_omega250()
%   plot_phase_evolution_omega250(dirpath, n_seconds_to_cut, plot_duration)
%
% Defaults:
%   dirpath = 'EstimateF/omega250'
%   n_seconds_to_cut = 0
%   plot_duration = 105
save%   apply_filter = true
%   filter_window_size = 10
%   do_save_figure = false
%   sample_window = [50, 60]

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('EstimateF','Spring5/265');
    end
    if nargin < 2 || isempty(n_seconds_to_cut)
        n_seconds_to_cut = 0;0
    end
    if nargin < 3 || isempty(plot_duration)
        plot_duration = 105;
    end
    if nargin < 4 || isempty(apply_filter)
        apply_filter = true;
    end
    if nargin < 5 || isempty(filter_window_size)
        filter_window_size = 1;
    end
    if nargin < 6 || isempty(do_save_figure)
        do_save_figure = false;
    end
    if nargin < 7 || isempty(n_sync)
        n_sync = 1; % default: check 2:1 synchronization
    end
    if nargin < 8 || isempty(m_sync)
        m_sync = 1;
    end
    if nargin < 9 || isempty(sample_window)
        sample_window = [50, 60];
    end

    if numel(sample_window) ~= 2
        error('sample_window must contain exactly two elements [t_start, t_end].');
    end
    sample_window = sort(sample_window(:)).';
    if sample_window(2) <= sample_window(1)
        error('sample_window end must be greater than start.');
    end

    if ~isfolder(dirpath)
        error('Directory not found: %s', dirpath);
    end

    csvs = dir(fullfile(dirpath, '*.csv'));
    if isempty(csvs)
        error('No CSV files found in %s', dirpath);
    end

    file_list = arrayfun(@(d) fullfile(dirpath,d.name), csvs, 'UniformOutput', false);

    % --- overflow / chunk jump correction constants (same as original) ---
    T_OVERFLOW = 2^32 / 1e6; % 約4294.967296秒
    T_TOL = 5.0;             % 許容誤差（秒）
    threshold_sec = T_OVERFLOW - T_TOL;
    jump_sec = T_OVERFLOW;
    CLUSTER_VAR_THRESHOLD = 0.35;

    % Read all files and collect agent sets
    file_tables = cell(size(file_list));
    agent_sets = cell(size(file_list));

    for i = 1:numel(file_list)
        try
            T = readtable(file_list{i});
        catch ME
            warning('Failed to read %s: %s', file_list{i}, ME.message);
            continue;
        end
        % Keep expected columns if present
        if ~all(ismember({'time_pc_sec_abs','a0','agent_id'}, T.Properties.VariableNames))
            warning('Skipping %s: missing required columns', file_list{i});
            continue;
        end
        % Ensure chunk_id exists (if not, create a dummy single chunk)
        if ~ismember('chunk_id', T.Properties.VariableNames)
            T.chunk_id = ones(height(T),1);
        end
        file_tables{i} = T(:,{'time_pc_sec_abs','a0','agent_id','chunk_id'});
        % Apply overflow / chunk-start corrections per-file (same logic as original)
        file_tables{i} = correct_large_jump_matlab(file_tables{i}, threshold_sec, jump_sec);
        file_tables{i} = correct_chunk_start_times_matlab(file_tables{i}, 4000.0, T_OVERFLOW);
        FT = file_tables{i};
        agent_sets{i} = unique(FT.agent_id);
    end

    valid_idx = ~cellfun(@isempty, file_tables);
    file_tables = file_tables(valid_idx);
    file_list = file_list(valid_idx);
    agent_sets = agent_sets(valid_idx);

    if isempty(file_tables)
        error('No valid data files to plot.');
    end

    % Determine agents common to all files (prefer intersection)
    common_agents = agent_sets{1};
    for i = 2:numel(agent_sets)
        common_agents = intersect(common_agents, agent_sets{i});
    end
    if isempty(common_agents)
        % fallback: use agents from first file
        common_agents = agent_sets{1};
        warning('No common agents across files; using agents from first file.');
    end

    % Decide base agent
    base_agent = min(common_agents);
    other_agents = setdiff(common_agents, base_agent);
    if isempty(other_agents)
        error('Not enough agents to compute relative phases (need at least 2).');
    end

    % Prepare per-file phase series using original pipeline
    allow_missing_agents = 1;
    max_agent_id = max(common_agents);
    phase_series_by_file = cell(numel(file_tables), 1);
    for f = 1:numel(file_tables)
        phase_series_by_file{f} = compute_phase_series_for_file( ...
            file_tables{f}, base_agent, n_seconds_to_cut, plot_duration, ...
            allow_missing_agents, apply_filter, filter_window_size, max_agent_id, n_sync, m_sync);
    end

    % Prepare figure with one subplot per other agent
    n_plots = numel(other_agents);
    %figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
    figure('Visible','on');
    colors = lines(numel(file_list));
    max_plot_time = min(plot_duration - n_seconds_to_cut, 60);

    % Display label for modified phase-combination using n:m notation
    y_label_str = sprintf('$$%d\\phi_{j} - %d\\phi_{\\mathrm{1}}$$', n_sync, m_sync);

    for p = 1:n_plots
        ag = other_agents(p);
        subplot(n_plots,1,p); hold on;
%        title(sprintf('Agent %d - Agent %d', ag, base_agent));
        ylabel(y_label_str,'Interpreter','latex');
        ylim([-pi, pi]);
        yticks([-pi,0,pi]);
        yticklabels({'-\pi','0','\pi'});
        set(gca,'TickLabelInterpreter','latex');
        yline(0, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 0.8);

        for f = 1:numel(file_list)
            series_struct = phase_series_by_file{f};
            if ag > numel(series_struct)
                continue;
            end
            if isempty(series_struct(ag).time) || isempty(series_struct(ag).phase)
                continue;
            end
            plot(series_struct(ag).time, series_struct(ag).phase, ...
                'Color', colors(f,:), 'LineWidth', 0.8);
        end

        xlim([0, max_plot_time]);
        if p == n_plots
            xlabel('Time (s)','Interpreter','latex');
        end
        grid on;
        hold off;
    end

    % no legend requested — overlays only

    tuneFigure();

    cluster_info = cluster_phase_window_means(phase_series_by_file, other_agents, sample_window, CLUSTER_VAR_THRESHOLD);

    if do_save_figure
        saveFigure();
    end

    if nargout >= 1
        varargout{1} = cluster_info;
    end
    if nargout >= 2
        varargout{2} = phase_series_by_file;
    end
end

function cluster_info = cluster_phase_window_means(phase_series_by_file, other_agents, sample_window, cluster_var_threshold)
    window_start = sample_window(1);
    window_end = sample_window(2);
    all_samples = [];

    for f = 1:numel(phase_series_by_file)
        series_struct = phase_series_by_file{f};
        if isempty(series_struct)
            continue;
        end
        for ag = other_agents(:).'
            if ag > numel(series_struct)
                continue;
            end
            times = series_struct(ag).time;
            phases = series_struct(ag).phase;
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
            all_samples(end+1) = agent_mean; %#ok<AGROW>
        end
    end

    cluster_info = struct('cluster_id', {}, 'n_clusters', {}, 'mean', {}, 'stderr', {}, 'count', {}, 'circ_var', {}, 'nsamples_total', {});

    if isempty(all_samples)
        fprintf('[INFO] No valid phase samples found between %.2f and %.2f seconds.\n', window_start, window_end);
        return;
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
            k = 1;
        end
    else
        k = 1;
        idx = ones(size(all_samples(:)));
    end

    fprintf('[INFO] Phase clusters over %.2f-%.2f s (Nsamples=%d, circVar=%.3f):\n', ...
        window_start, window_end, Nsamples, circ_var);

    for c = 1:k
        mask = idx == c;
        ths = all_samples(mask);
        if isempty(ths)
            continue;
        end
        mmean = atan2(mean(sin(ths)), mean(cos(ths)));
        Rcl = abs(mean(exp(1i * ths)));
        circ_std = sqrt(max(0, -2 * log(max(Rcl, eps))));
        stderr = circ_std / sqrt(max(1, numel(ths)));
        fprintf('  Cluster %d/%d: mean = %.3f rad (%.2f deg), stderr = %.3f (n=%d)\n', ...
            c, k, mmean, mmean * 180/pi, stderr, numel(ths));
        cluster_info(end+1) = struct( ...
            'cluster_id', c, ...
            'n_clusters', k, ...
            'mean', mmean, ...
            'stderr', stderr, ...
            'count', numel(ths), ...
            'circ_var', circ_var, ...
            'nsamples_total', Nsamples); %#ok<AGROW>
    end
end

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
        % Use n_sync*phi_agent - m_sync*phi_base for n:m synchronization checks
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
    % same logic as in plot_relative_phase_matlab_2modules.m
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
    % グループ化（agent_id, chunk_id 単位）
    [G, ~] = findgroups(df_all.agent_id, df_all.chunk_id);
    fprintf('[INFO] Found %d unique chunks.\n', max(G));

    % 該当ブロックを修正
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end

        % 時系列を並び替え
        [~, sorted_idx_rel] = sort(df_all.time_pc_sec_abs(idx));
        idx = idx(sorted_idx_rel);
        
        time_series = df_all.time_pc_sec_abs(idx);
        time_diff = [0; diff(time_series)];

        jump_idx = find(time_diff > threshold_sec);
        for j = 1:length(jump_idx)
            fix_range = jump_idx(j):length(time_series);
            df_all.time_pc_sec_abs(idx(fix_range)) = df_all.time_pc_sec_abs(idx(fix_range)) - jump_sec;
            fprintf('[FIX] Corrected overflow at index %d, subtracted %.6f sec.\n', idx(jump_idx(j)), jump_sec);
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
            continue;  % 空グループスキップ
        end
        start_time = df_all.time_pc_sec_abs(idx(1));
        if start_time - median_start > threshold_sec
            df_all.time_pc_sec_abs(idx) = df_all.time_pc_sec_abs(idx) - jump_sec;
            % agent_id, chunk_idの表示（chunk_idが無い場合も対応）
            if has_chunk_id
                aid = chunk_keys(i,1);
                cid = chunk_keys(i,2);
            else
                aid = chunk_keys(i);
                cid = -1;  % または NaN
            end
            fprintf('[FIX] Corrected chunk time for agent %d, chunk %d: %.3f → %.3f\n', ...
                aid, cid, start_time, start_time - jump_sec);
        end
    end
end