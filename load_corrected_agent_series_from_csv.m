function [series_by_agent, all_agents] = load_corrected_agent_series_from_csv(csv_path, requested_agents, required_cols)
% Load one CSV, apply shared preprocessing, and return per-agent time series.

    if nargin < 2
        requested_agents = [];
    end
    if nargin < 3 || isempty(required_cols)
        required_cols = {'time_pc_sec_abs', 'a0'};
    end

    required_cols = [{'time_pc_sec_abs'}, required_cols(:).'];
    required_cols = unique(required_cols, 'stable');

    T = readtable(csv_path);
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
        t_overflow = 2^32 / 1e6;
        t_tol = 5.0;
        threshold_sec = t_overflow - t_tol;
        jump_sec = t_overflow;

        T = correct_large_jump_matlab(T, threshold_sec, jump_sec);
        T = correct_chunk_start_times_matlab(T, 4000.0, t_overflow);
    end

    all_agents = unique(T.agent_id, 'sorted').';

    if isempty(requested_agents)
        requested_agents = all_agents;
    else
        requested_agents = requested_agents(:).';
    end

    if isempty(requested_agents)
        series_by_agent = struct('time', {}, 'a0_corr', {}, 'a1', {}, 'a2', {});
        return;
    end

    if any(~isfinite(requested_agents)) || any(requested_agents < 1) || any(requested_agents ~= round(requested_agents))
        error('requested_agents must contain positive integer agent IDs.');
    end

    if ~all(ismember(requested_agents, all_agents))
        error('Requested agents %s are not all present. Available agents: %s', ...
            mat2str(requested_agents), mat2str(all_agents));
    end

    series_template = struct('time', [], 'a0_corr', [], 'a1', [], 'a2', []);
    series_by_agent = repmat(series_template, 1, max(requested_agents));
    for agent_idx = 1:numel(requested_agents)
        aid = requested_agents(agent_idx);
        series_by_agent(aid) = get_agent_series(T, aid);
    end
end

function series = get_agent_series(T, agent_id)
    sub = T(T.agent_id == agent_id, :);
    sub = sortrows(sub, 'time_pc_sec_abs');

    if isempty(sub)
        series = struct('time', [], 'a0_corr', [], 'a2', []);
        return;
    end

    time_vals = double(sub.time_pc_sec_abs(:));
    valid = isfinite(time_vals);

    a0_vals = [];
    if ismember('a0', sub.Properties.VariableNames)
        a0_vals = correct_phase_discontinuity(double(sub.a0(:)));
        valid = valid & isfinite(a0_vals);
    end

    a1_vals = [];
    if ismember('a1', sub.Properties.VariableNames)
        a1_vals = double(sub.a1(:));
        valid = valid & isfinite(a1_vals);
    end

    a2_vals = [];
    if ismember('a2', sub.Properties.VariableNames)
        a2_vals = double(sub.a2(:));
        valid = valid & isfinite(a2_vals);
    end

    time_vals = time_vals(valid);
    if ~isempty(a0_vals)
        a0_vals = a0_vals(valid);
    end
    if ~isempty(a2_vals)
        a2_vals = a2_vals(valid);
    end

    [time_vals, ia] = unique(time_vals, 'stable');
    if ~isempty(a0_vals)
        a0_vals = a0_vals(ia);
    end
    if ~isempty(a1_vals)
        a1_vals = a1_vals(ia);
    end
    if ~isempty(a2_vals)
        a2_vals = a2_vals(ia);
    end

    series = struct('time', time_vals, 'a0_corr', [], 'a1', [], 'a2', []);
    if ~isempty(a0_vals)
        series.a0_corr = a0_vals;
    end
    if ~isempty(a1_vals)
        series.a1 = a1_vals;
    end
    if ~isempty(a2_vals)
        series.a2 = a2_vals;
    end
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