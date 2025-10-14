function resultsTable = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, agent_id, export_outputs, make_plot)
% Compute average angular velocity for a given agent (default id=6)
% between [t_start, t_end] seconds for all CSV logs in folder_path.
% Time is set relative to the first timestamp in each file (experiment time 0).
% Uses the same correction logic as plot_relative_phase_matlab_2modules.m.
%
% Usage:
%   resultsTable = compute_avg_omega_id6_kappa5();
%   resultsTable = compute_avg_omega_id6_kappa5('ArnoldPlot/kappa5', 10, 50, 6);

    if nargin < 1 || isempty(folder_path)
        folder_path = fullfile(pwd, 'ArnoldPlot', 'kappa5');
    end
    if nargin < 2 || isempty(t_start)
        t_start = 10;  % [s]
    end
    if nargin < 3 || isempty(t_end)
        t_end = 50;    % [s]
    end
    if nargin < 4 || isempty(agent_id)
        agent_id = 6;
    end
    if nargin < 5 || isempty(export_outputs)
        export_outputs = false; % only save when true
    end
    if nargin < 6 || isempty(make_plot)
        make_plot = true; % draw figure by default
    end

    if ~isfolder(folder_path)
        error('Folder not found: %s', folder_path);
    end

    files = dir(fullfile(folder_path, '*.csv'));
    if isempty(files)
        error('No CSV files found in %s', folder_path);
    end

    % Constants (same as reference)
    T_OVERFLOW = 2^32 / 1e6; % ~4294.967296 sec
    T_TOL = 5.0;             % tolerance [s]
    threshold_sec = T_OVERFLOW - T_TOL;
    jump_sec = T_OVERFLOW;

    results = struct('file', {}, 'agent_id', {}, 'mean_omega_rad_s', {}, 'mean_omega_deg_s', {}, ...
                     'base_agent_id', {}, 'base_mean_omega_rad_s', {}, 'ratio_to_base', {}, ...
                     'duration_used_s', {}, 'status', {});

    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        try
            T = readtable(fpath);
        catch ME
            warning('Failed to read %s: %s', fpath, ME.message);
            results(end+1) = struct('file', files(i).name, ...
                                    'agent_id', agent_id, ...
                                    'mean_omega_rad_s', NaN, ...
                                    'mean_omega_deg_s', NaN, ...
                                    'base_agent_id', NaN, ...
                                    'base_mean_omega_rad_s', NaN, ...
                                    'ratio_to_base', NaN, ...
                                    'duration_used_s', 0, ...
                                    'status', 'read_error'); %#ok<AGROW>
            continue;
        end

        required = {'agent_id','chunk_id','time_pc_sec_abs','a0'};
        if ~all(ismember(required, T.Properties.VariableNames))
            warning('Skipping %s: missing required columns %s', files(i).name, strjoin(required, ', '));
            results(end+1) = struct('file', files(i).name, ...
                                    'agent_id', agent_id, ...
                                    'mean_omega_rad_s', NaN, ...
                                    'mean_omega_deg_s', NaN, ...
                                    'base_agent_id', NaN, ...
                                    'base_mean_omega_rad_s', NaN, ...
                                    'ratio_to_base', NaN, ...
                                    'duration_used_s', 0, ...
                                    'status', 'missing_columns'); %#ok<AGROW>
            continue;
        end

        % Compute target agent mean omega
        [ok_t, mean_omega, mean_omega_deg, dur_used, status_t] = ...
            compute_agent_mean_omega(T, agent_id, t_start, t_end, threshold_sec, jump_sec, T_OVERFLOW);
        if ~ok_t
            if ~strcmp(status_t, 'no_agent_data')
                warning('Agent %d in %s: %s', agent_id, files(i).name, status_t);
            end
            results(end+1) = struct('file', files(i).name, ...
                                    'agent_id', agent_id, ...
                                    'mean_omega_rad_s', NaN, ...
                                    'mean_omega_deg_s', NaN, ...
                                    'base_agent_id', NaN, ...
                                    'base_mean_omega_rad_s', NaN, ...
                                    'ratio_to_base', NaN, ...
                                    'duration_used_s', dur_used, ...
                                    'status', status_t); %#ok<AGROW>
            continue;
        end

        % Determine baseline as smallest agent id in file (excluding 99)
        agents_in_file = unique(T.agent_id);
        agents_in_file = agents_in_file(agents_in_file ~= 99);
        if isempty(agents_in_file)
            base_id = agent_id; % fallback
        else
            base_id = min(agents_in_file);
        end

        [ok_b, base_mean_omega, ~, ~, status_b] = ...
            compute_agent_mean_omega(T, base_id, t_start, t_end, threshold_sec, jump_sec, T_OVERFLOW);
        if ~ok_b
            warning('Base agent %d in %s: %s', base_id, files(i).name, status_b);
            results(end+1) = struct('file', files(i).name, ...
                                    'agent_id', agent_id, ...
                                    'mean_omega_rad_s', mean_omega, ...
                                    'mean_omega_deg_s', mean_omega_deg, ...
                                    'base_agent_id', base_id, ...
                                    'base_mean_omega_rad_s', NaN, ...
                                    'ratio_to_base', NaN, ...
                                    'duration_used_s', dur_used, ...
                                    'status', ['base_' status_b]); %#ok<AGROW>
            continue;
        end

        if abs(base_mean_omega) < 1e-9
            ratio_val = NaN;
        else
            ratio_val = mean_omega / base_mean_omega;
        end

        results(end+1) = struct('file', files(i).name, ...
                                'agent_id', agent_id, ...
                                'mean_omega_rad_s', mean_omega, ...
                                'mean_omega_deg_s', mean_omega_deg, ...
                                'base_agent_id', base_id, ...
                                'base_mean_omega_rad_s', base_mean_omega, ...
                                'ratio_to_base', ratio_val, ...
                                'duration_used_s', t_end - t_start, ...
                                'status', 'ok'); %#ok<AGROW>
    end

    % Build table
    if isempty(results)
        resultsTable = table();
    else
        resultsTable = struct2table(results);
    end

    % Sort by filename naturally
    try
        resultsTable = sortrows(resultsTable, 'file');
    catch
    end

    % Plot line of ratio (omega_target / omega_base)
    if make_plot
        ratio_all = resultsTable.ratio_to_base;
        figure('Color','w');
        if ~isempty(ratio_all)
            % X-axis: 1.0 .. 4.0 by 0.1 when count is 31; else fallback to linspace
            if numel(ratio_all) == 31
                x = 1.0:0.1:4.0;
            else
                x = linspace(1.0, 4.0, numel(ratio_all));
            end
            plot(x, ratio_all, '-o', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.4 0.8]);
            grid on;
            xlabel('X');
            ylabel('Ratio: \\omega_{agent} / \\omega_{base}');
            title(sprintf('Agent %d ratio to base (min agent id) in [%g, %g] s across %d files', agent_id, t_start, t_end, numel(ratio_all)));
            xlim([min(x) max(x)]);
            xticks(x);
        else
            text(0.5, 0.5, 'No valid results', 'HorizontalAlignment','center');
            axis off;
        end
    end

    % Save outputs (only when export_outputs is true)
    if export_outputs && make_plot
        outDir = fullfile(pwd, 'exports');
        if ~isfolder(outDir)
            try
                mkdir(outDir);
            catch
            end
        end
        ts = datestr(now, 'yyyymmdd_HHMMSS');
        outFig = fullfile(outDir, sprintf('ratio_omega_agent%d_to_base_%s.png', agent_id, ts));
        outCsv = fullfile(outDir, sprintf('ratio_omega_agent%d_to_base_%s.csv', agent_id, ts));
        try
            saveas(gcf, outFig);
        catch
        end
        try
            writetable(resultsTable, outCsv);
        catch
        end
    end
end

% --- helpers copied from plot_relative_phase_matlab_2modules.m ---
function corrected_phase = correct_phase_discontinuity(phase_data)
    % Unwrap-like correction for 8-bit phase with wrap at 256
    corrected_phase = double(phase_data);
    for i = 2:length(corrected_phase)
        d = corrected_phase(i) - corrected_phase(i - 1);
        if d < -128
            corrected_phase(i:end) = corrected_phase(i:end) + 256;
        elseif d > 128
            corrected_phase(i:end) = corrected_phase(i:end) - 256;
        end
    end
end

function df_all = correct_large_jump_matlab(df_all, threshold_sec, jump_sec)
    % Correct large time jumps (overflow) per (agent_id, chunk_id)
    [G, ~] = findgroups(df_all.agent_id, df_all.chunk_id);
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx), continue; end
        [~, sidx_rel] = sort(df_all.time_pc_sec_abs(idx));
        idx = idx(sidx_rel);
        time_series = df_all.time_pc_sec_abs(idx);
        time_diff = [0; diff(time_series)];
        jump_idx = find(time_diff > threshold_sec);
        for j = 1:numel(jump_idx)
            fix_range = jump_idx(j):numel(time_series);
            df_all.time_pc_sec_abs(idx(fix_range)) = df_all.time_pc_sec_abs(idx(fix_range)) - jump_sec;
        end
    end
end

function df_all = correct_chunk_start_times_matlab(df_all, threshold_sec, jump_sec)
    [G, ~] = findgroups(df_all.agent_id, df_all.chunk_id);
    if isempty(G)
        return;
    end
    chunk_start = splitapply(@(x) min(x), df_all.time_pc_sec_abs, G);
    median_start = median(chunk_start);
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx), continue; end
        start_time = df_all.time_pc_sec_abs(idx(1));
        if start_time - median_start > threshold_sec
            df_all.time_pc_sec_abs(idx) = df_all.time_pc_sec_abs(idx) - jump_sec;
        end
    end
end

function [ok, mean_omega, mean_omega_deg, duration_used, status] = compute_agent_mean_omega(T, aid, t_start, t_end, threshold_sec, jump_sec, T_OVERFLOW)
    ok = false; mean_omega = NaN; mean_omega_deg = NaN; duration_used = 0; status = 'unknown';
    try
        sub = T(T.agent_id == aid, {'agent_id','chunk_id','time_pc_sec_abs','a0'});
        if isempty(sub)
            status = 'no_agent_data'; return; end

        sub = correct_large_jump_matlab(sub, threshold_sec, jump_sec);
        sub = correct_chunk_start_times_matlab(sub, 4000.0, T_OVERFLOW);
        sub = sortrows(sub, 'time_pc_sec_abs');
        [~, ia] = unique(sub.time_pc_sec_abs);
        sub = sub(ia, :);
        t_rel = sub.time_pc_sec_abs - sub.time_pc_sec_abs(1);
        if t_rel(end) < t_end
            duration_used = t_rel(end);
            status = 'short_duration';
            return;
        end
        a0_corr = correct_phase_discontinuity(sub.a0);
        theta = double(a0_corr) * (2*pi/256);
        theta_start = interp1(t_rel, theta, t_start, 'linear', 'extrap');
        theta_end   = interp1(t_rel, theta, t_end,   'linear', 'extrap');
        mean_omega = (theta_end - theta_start) / (t_end - t_start);
        mean_omega_deg = mean_omega * (180/pi);
        duration_used = t_end - t_start;
        ok = true; status = 'ok';
    catch ME
        status = ['error_' ME.identifier];
    end
end
