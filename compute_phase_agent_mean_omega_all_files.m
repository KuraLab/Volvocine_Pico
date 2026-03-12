function omega_summary = compute_phase_agent_mean_omega_all_files(csv_paths, phase_agent_ids, analysis_duration_sec, analysis_start_sec)
% Compute mean angular velocity for the two analyzed phase agents across files.

    omega_summary = struct( ...
        'available', false, ...
        'reason', '', ...
        'phase_agent_ids', phase_agent_ids(:).', ...
        'analysis_duration_sec', analysis_duration_sec, ...
        'analysis_start_sec', analysis_start_sec, ...
        'per_file', struct('file_path', {}, 'window_start_abs', {}, 'window_end_abs', {}, ...
            'duration_used_s', {}, 'mean_omega_rad_s', {}, 'mean_omega_deg_s', {}, 'omega_ratio_2_over_1', {}), ...
        'skipped_files', struct('file_path', {}, 'reason', {}), ...
        'mean_omega_rad_s', [NaN, NaN], ...
        'mean_omega_deg_s', [NaN, NaN], ...
        'std_omega_rad_s', [NaN, NaN], ...
        'std_omega_deg_s', [NaN, NaN], ...
        'mean_ratio_2_over_1', NaN, ...
        'median_ratio_2_over_1', NaN);

    per_file = omega_summary.per_file;
    skipped_files = omega_summary.skipped_files;
    omega_values_rad = [];
    ratio_values = [];

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            file_result = compute_phase_agent_mean_omega_for_csv( ...
                csv_path, phase_agent_ids, analysis_duration_sec, analysis_start_sec);
        catch ME
            skipped_files(end+1) = struct('file_path', csv_path, 'reason', ME.message); %#ok<AGROW>
            continue;
        end

        per_file(end+1) = file_result; %#ok<AGROW>
        omega_values_rad = [omega_values_rad; file_result.mean_omega_rad_s(:).']; %#ok<AGROW>
        ratio_values(end+1, 1) = file_result.omega_ratio_2_over_1; %#ok<AGROW>
    end

    if isempty(per_file)
        omega_summary.reason = 'No valid files were available to compute phase-agent mean angular velocity.';
        omega_summary.skipped_files = skipped_files;
        return;
    end

    omega_summary.available = true;
    omega_summary.per_file = per_file;
    omega_summary.skipped_files = skipped_files;
    omega_summary.mean_omega_rad_s = mean(omega_values_rad, 1, 'omitnan');
    omega_summary.mean_omega_deg_s = omega_summary.mean_omega_rad_s * (180 / pi);
    omega_summary.std_omega_rad_s = std(omega_values_rad, 0, 1, 'omitnan');
    omega_summary.std_omega_deg_s = omega_summary.std_omega_rad_s * (180 / pi);
    omega_summary.mean_ratio_2_over_1 = mean(ratio_values, 'omitnan');
    omega_summary.median_ratio_2_over_1 = median(ratio_values, 'omitnan');
end

function file_result = compute_phase_agent_mean_omega_for_csv(csv_path, phase_agent_ids, analysis_duration_sec, analysis_start_sec)
    requested_agents = phase_agent_ids(:).';
    series_by_agent = load_corrected_agent_series_from_csv( ...
        csv_path, requested_agents, {'time_pc_sec_abs', 'a0'});
    overlap_start = -inf;
    overlap_end = inf;
    for agent_idx = 1:numel(requested_agents)
        aid = requested_agents(agent_idx);
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

    duration_used_s = window_end_abs - window_start_abs;
    if duration_used_s <= 0
        error('No usable duration remained after applying the analysis window.');
    end

    mean_omega_rad_s = nan(1, numel(requested_agents));
    for agent_idx = 1:numel(requested_agents)
        aid = requested_agents(agent_idx);
        time_vals = series_by_agent(aid).time;
        theta_vals = double(series_by_agent(aid).a0_corr) * (2 * pi / 256);

        theta_start = interp1(time_vals, theta_vals, window_start_abs, 'linear', 'extrap');
        theta_end = interp1(time_vals, theta_vals, window_end_abs, 'linear', 'extrap');
        mean_omega_rad_s(agent_idx) = (theta_end - theta_start) / duration_used_s;
    end

    omega_ratio_2_over_1 = NaN;
    if abs(mean_omega_rad_s(1)) >= 1e-12
        omega_ratio_2_over_1 = mean_omega_rad_s(2) / mean_omega_rad_s(1);
    end

    file_result = struct( ...
        'file_path', csv_path, ...
        'window_start_abs', window_start_abs, ...
        'window_end_abs', window_end_abs, ...
        'duration_used_s', duration_used_s, ...
        'mean_omega_rad_s', mean_omega_rad_s, ...
        'mean_omega_deg_s', mean_omega_rad_s * (180 / pi), ...
        'omega_ratio_2_over_1', omega_ratio_2_over_1);
end