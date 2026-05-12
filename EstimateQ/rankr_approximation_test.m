function rank_analysis = rankr_approximation_test(dirpath, phase_agent_ids, varargin)
% Minimal rank-R approximation test for s_j(phi1, phi2)
% Uses extract_phase_pairs_for_rank_analysis to get data from CSV
% Then performs SVD-based low-rank decomposition analysis
%
% Usage:
%   rank_analysis = rankr_approximation_test();
%   rank_analysis = rankr_approximation_test('Spring3/255', [2, 4]);
%   rank_analysis = rankr_approximation_test(..., 'M', 10, 'N', 10, 'SignalField', 'a2_normalized_agent');
%   rank_analysis = rankr_approximation_test(..., 'SignalField', 'sin_phi2_a2_agent');  % Phase-modulated version
%   rank_analysis = rankr_approximation_test(..., 'SignalField', 'a2_normalized_all');  % Common signal for both agents
%
%   % Remove self-phase-only terms and plot self-only profile
%   rank_analysis = rankr_approximation_test(..., ...
%       'RemoveSelfOnly', true, 'RemoveConstant', false, 'RemoveOtherOnly', false, ...
%       'PlotSelfOnlyProfile', true);
%
%   % Analyze pure interaction components only (remove all marginal terms and constant)
%   rank_analysis = rankr_approximation_test(..., ...
%       'RemoveSelfOnly', true, 'RemoveConstant', true, 'RemoveOtherOnly', true, ...
%       'PlotSelfOnlyProfile', true);
%
%   % Check saved results
%   numel(rank_analysis)
%   unique([rank_analysis.agent_id])
%   [rank_analysis([rank_analysis.agent_id] == 2).rank]
%   [rank_analysis([rank_analysis.agent_id] == 4).rank]
%
% ========== SETTINGS (Edit here to change defaults) ==========
    DEFAULT_M = 20;                              % Fourier order M
    DEFAULT_N = 20;                              % Fourier order N
    DEFAULT_SIGNAL_FIELD = 'a2_normalized_agent'; % Signal field: 'a2_normalized_agent' for agent-specific raw signal, 'sin_phi2_a2_agent' for phase-modulated, 'a2_normalized_all' for common signal
    DEFAULT_ANALYSIS_DURATION = 80;             % Analysis duration in seconds
    DEFAULT_ANALYSIS_START = 10;                % Analysis start time in seconds
    PLOT_SHOW_POINTS = false;                    % Show original data points on 3D plots (true/false)
    DEFAULT_COLORMAP = 'jet';                   % Default colormap: 'jet', 'parula', 'cool', 'hot', etc.
    DEFAULT_PROFILE_RANK = 3;                   % Number of rank components to profile (a_r, b_r)
    DEFAULT_REMOVE_SELF_ONLY = true;            % Remove self-phase-only terms before rank analysis
    DEFAULT_REMOVE_CONSTANT = true;            % Remove constant term (m=0, n=0)
    DEFAULT_REMOVE_OTHER_ONLY = false;          % Remove other-phase-only terms
    DEFAULT_PLOT_SELF_ONLY_PROFILE = true;      % Plot self-phase-only profile
% =============================================================="

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('Spring3', '255');
    end

    if nargin < 2 || isempty(phase_agent_ids)
        phase_agent_ids = [];
    end

    parser = inputParser;
    addParameter(parser, 'M', DEFAULT_M, @isnumeric);
    addParameter(parser, 'N', DEFAULT_N, @isnumeric);
    addParameter(parser, 'SignalField', DEFAULT_SIGNAL_FIELD, @ischar);
    addParameter(parser, 'AnalysisDuration', DEFAULT_ANALYSIS_DURATION, @isnumeric);
    addParameter(parser, 'AnalysisStart', DEFAULT_ANALYSIS_START, @isnumeric);
    addParameter(parser, 'ProfileRank', DEFAULT_PROFILE_RANK, @isnumeric);
    addParameter(parser, 'RemoveSelfOnly', DEFAULT_REMOVE_SELF_ONLY, @islogical);
    addParameter(parser, 'RemoveConstant', DEFAULT_REMOVE_CONSTANT, @islogical);
    addParameter(parser, 'RemoveOtherOnly', DEFAULT_REMOVE_OTHER_ONLY, @islogical);
    addParameter(parser, 'PlotSelfOnlyProfile', DEFAULT_PLOT_SELF_ONLY_PROFILE, @islogical);

    parse(parser, varargin{:});

    M = round(parser.Results.M);
    N = round(parser.Results.N);
    signal_field = parser.Results.SignalField;
    profile_rank = round(parser.Results.ProfileRank);
    remove_self_only = parser.Results.RemoveSelfOnly;
    remove_constant = parser.Results.RemoveConstant;
    remove_other_only = parser.Results.RemoveOtherOnly;
    flag_plot_self_only_profile = parser.Results.PlotSelfOnlyProfile;
    
    % Extract data inline - for all agents
    out_full = extract_phase_pairs_for_rank_analysis_inline(dirpath, phase_agent_ids, [], ...
        parser.Results.AnalysisDuration, parser.Results.AnalysisStart);

    phi1_all = out_full.phi1_all;
    phi2_all = out_full.phi2_all;
    
    % Analyze each target agent
    target_agents = out_full.phase_agent_ids;
    fprintf('[Rank-R Approximation Analysis]\n');
    fprintf('  Analyzing target agents: %s\n', sprintf('%d ', target_agents));
    fprintf('  Samples: %d\n', numel(phi1_all));
    fprintf('  Fourier order: M=%d, N=%d\n', M, N);
    fprintf('  Remove self-only terms: %s\n', mat2str(remove_self_only));
    fprintf('  Remove constant term: %s\n', mat2str(remove_constant));
    fprintf('  Remove other-only terms: %s\n', mat2str(remove_other_only));
    fprintf('  Plot self-only profile: %s\n\n', mat2str(flag_plot_self_only_profile));
    
    rank_analysis = struct([]);
    expected_result_count = 0;
    
    for agent_idx = 1:numel(target_agents)
        target_agent_id = target_agents(agent_idx);
        
        % Select signal for this agent
        if strcmp(signal_field, 'sin_phi2_a2_agent')
            % Use agent-specific derived signal sin(phi2 + tau) * a2_normalized
            field_name = sprintf('sin_phi2_a2_agent%d', target_agent_id);
            if isfield(out_full, field_name)
                y_all = out_full.(field_name);
            else
                error('Agent-specific signal field %s not found for agent %d', field_name, target_agent_id);
            end
        elseif strcmp(signal_field, 'a2_normalized_agent')
            % Use agent-specific raw a2_normalized (no phase modulation)
            field_name = sprintf('a2_normalized_agent%d', target_agent_id);
            if isfield(out_full, field_name)
                y_all = out_full.(field_name);
            else
                error('Agent-specific raw signal field %s not found for agent %d', field_name, target_agent_id);
            end
        elseif strcmp(signal_field, 'a2_normalized_all')
            % Use raw a2_normalized for both agents
            y_all = out_full.a2_normalized_all;
        elseif isfield(out_full, signal_field)
            y_all = out_full.(signal_field);
        else
            error('Signal field %s for agent %d not found', signal_field, target_agent_id);
        end

        if isempty(y_all)
            fprintf('  [Agent %d] No valid data.\n\n', target_agent_id);
            continue;
        end
        
        % Determine signal index (s_1 or s_2)
        if target_agent_id == target_agents(1)
            signal_index = 1;
        else
            signal_index = 2;
        end
        
        % Determine self phase
        if target_agent_id == out_full.phase_agent_ids(1)
            self_phase_index = 1;  % phi1
            self_phase_name = 'phi1';
        else
            self_phase_index = 2;  % phi2
            self_phase_name = 'phi2';
        end

        fprintf('[Agent %d (s_%d)] - Self phase: %s\n', target_agent_id, signal_index, self_phase_name);
        fprintf('  Rank | Energy Ratio | NRMSE\n');
        fprintf('  -----+--------------+----------\n');

        % Estimate full Fourier coefficients via regularized LS
        [C_full, m_values, n_values] = estimate_fourier_coeff_matrix(phi1_all, phi2_all, y_all, M, N);
        
        % Remove phase marginal terms if requested
        [C_analysis, y_analysis, removed_info] = remove_phase_marginal_terms( ...
            C_full, phi1_all, phi2_all, y_all, m_values, n_values, ...
            target_agent_id, out_full.phase_agent_ids, ...
            remove_self_only, remove_constant, remove_other_only);
        
        % Print removal information
        fprintf('  Removed self-only energy ratio:  %.6f\n', removed_info.removed_self_energy_ratio);
        fprintf('  Removed constant energy ratio:   %.6f\n', removed_info.removed_constant_energy_ratio);
        fprintf('  Removed other-only energy ratio: %.6f\n', removed_info.removed_other_energy_ratio);
        fprintf('  Remaining coefficient energy ratio: %.6f\n', removed_info.remaining_energy_ratio);
        fprintf('  Reconstruction mismatch RMSE: %.6e\n', removed_info.reconstruction_mismatch_rmse);

        % Compute rank-R analysis using SVD
        [U, S, V] = svd(C_analysis, 'econ');
        sigma = diag(S);
        energy = sum(sigma.^2);

        term_m = repelem(m_values(:), numel(n_values));
        term_n = repmat(n_values(:), numel(m_values), 1);

        n_rank_to_test = min(5, numel(sigma));
        expected_result_count = expected_result_count + n_rank_to_test;
        
        % Build separable components once for all ranks
        components_all = build_separable_components(U, S, V, m_values, n_values, min(5, numel(sigma)));
        
        % Compute self-only profile once for this agent
        self_profile = compute_self_only_profile(C_full, m_values, n_values, target_agent_id, out_full.phase_agent_ids);
        
        % Extract energy metrics from removed_info
        total_coeff_energy = removed_info.total_coeff_energy;
        removed_total_coeff_energy = removed_info.removed_total_coeff_energy;
        remaining_coeff_energy = removed_info.remaining_coeff_energy;
        removed_total_energy_ratio = removed_info.removed_total_energy_ratio;
        remaining_energy_ratio = removed_info.remaining_energy_ratio;
        
        fprintf('  Removed total energy ratio: %.6f\n', removed_total_energy_ratio);
        fprintf('  Remaining energy ratio: %.6f\n\n', remaining_energy_ratio);
        fprintf('  Rank | Residual Energy | Total Energy | Residual NRMSE | Total NRMSE\n');
        fprintf('  -----+-----------------+--------------+---------------+------------\n');

        rank_results = struct([]);  % Store results for plotting
        
        for R = 1:n_rank_to_test
            C_rankr = U(:, 1:R) * S(1:R, 1:R) * V(:, 1:R)';
            coef_rankr = reshape(C_rankr.', [], 1);

            energy_r = sum(sigma(1:R).^2);
            
            % Residual metrics (C_analysis based)
            residual_energy_ratio = energy_r / energy;
            residual_nrmse = compute_nrmse_fast(phi1_all, phi2_all, y_analysis, coef_rankr, term_m, term_n);
            
            % Total metrics (C_full based, with removed components added back)
            C_total_rankr = removed_info.C_removed_total + C_rankr;
            coef_total_rankr = reshape(C_total_rankr.', [], 1);
            total_energy_ratio = (removed_total_coeff_energy + energy_r) / total_coeff_energy;
            total_nrmse = compute_nrmse_fast(phi1_all, phi2_all, y_all, coef_total_rankr, term_m, term_n);

            fprintf('  %2d   |     %.6f     |    %.6f    |     %.6f    |   %.6f\n', ...
                R, residual_energy_ratio, total_energy_ratio, residual_nrmse, total_nrmse);

            result = struct();
            result.agent_id = target_agent_id;
            result.signal_index = signal_index;
            result.self_phase = self_phase_name;
            result.rank = R;
            result.residual_energy_ratio = residual_energy_ratio;
            result.total_energy_ratio = total_energy_ratio;
            result.residual_nrmse = residual_nrmse;
            result.total_nrmse = total_nrmse;
            result.energy_ratio = residual_energy_ratio;  % Backward compatibility
            result.nrmse = residual_nrmse;  % Backward compatibility
            result.removed_total_energy_ratio = removed_total_energy_ratio;
            result.remaining_energy_ratio = remaining_energy_ratio;
            result.coef = coef_rankr;
            result.C_rankr = C_rankr;
            result.C_full = C_full;
            result.C_analysis = C_analysis;
            result.removed_info = removed_info;
            result.self_profile = self_profile;
            result.components = components_all(1:min(R, numel(components_all)));
            
            % Append to rank_results
            if isempty(rank_results)
                rank_results = orderfields(result);
            else
                rank_results(end + 1) = orderfields(result, rank_results(1)); %#ok<AGROW>
            end
            
            % Append to rank_analysis
            if isempty(rank_analysis)
                rank_analysis = orderfields(result);
            else
                rank_analysis(end + 1) = orderfields(result, rank_analysis(1)); %#ok<AGROW>
            end
        end
        fprintf('\n');
        
        % Plot colormap for ranks 1, 2, 3
        if numel(rank_results) >= 3
            title_suffix = '';
            if remove_self_only || remove_constant || remove_other_only
                title_suffix = ' - marginal removed';
            end
            
            plot_rank_approximations(phi1_all, phi2_all, y_analysis, rank_results(1:3), m_values, n_values, target_agent_id, title_suffix);
            plot_rank_approximations_3d(phi1_all, phi2_all, y_analysis, rank_results(1:3), m_values, n_values, target_agent_id, PLOT_SHOW_POINTS, DEFAULT_COLORMAP, title_suffix);
            
            % Plot separable profiles for ranks 1 to profile_rank
            plot_separable_profiles(components_all(1:min(profile_rank, numel(components_all))), m_values, n_values, target_agent_id, profile_rank, title_suffix);
            
            % Plot self-only profile if requested
            if flag_plot_self_only_profile
                plot_self_only_profile(self_profile);
            end
        end
    end
    
    % Display summary
    fprintf('[INFO] Total saved rank results: %d\n', numel(rank_analysis));
    fprintf('[INFO] Expected rank results: %d\n', expected_result_count);
    
    % Check saved results
    if numel(rank_analysis) > 0
        fprintf('[CHECK] Unique agents: %s\n', sprintf('%d ', unique([rank_analysis.agent_id])));
        agent_ids = [rank_analysis.agent_id];
        for aid = unique(agent_ids)
            n_for_agent = sum(agent_ids == aid);
            ranks_for_agent = [rank_analysis(agent_ids == aid).rank];
            fprintf('[CHECK] Agent %d: %d results, ranks = %s\n', aid, n_for_agent, sprintf('%d ', ranks_for_agent));
        end
    end
end

function [C, m_values, n_values] = estimate_fourier_coeff_matrix(phi1, phi2, y, M, N)
    m_values = (-M:M).';
    n_values = (-N:N).';

    term_m = repelem(m_values, numel(n_values));
    term_n = repmat(n_values, numel(m_values), 1);

    n_terms = numel(term_m);
    n_samples = numel(y);

    G = complex(zeros(n_terms, n_terms));
    h = complex(zeros(n_terms, 1));
    chunk_size = 20000;

    for k0 = 1:chunk_size:n_samples
        k1 = min(k0 + chunk_size - 1, n_samples);
        idx = k0:k1;

        A = exp(1i * (phi1(idx) * term_m.' + phi2(idx) * term_n.'));
        G = G + A' * A;
        h = h + A' * y(idx(:));
    end

    % Regularization: lambda = 1e-12 * trace(G) / n_terms
    lambda = 1e-12 * trace(G) / n_terms;
    Greg = G + lambda * eye(n_terms);

    % Solve regularized LS
    coef = Greg \ h;
    C = reshape(coef, numel(n_values), numel(m_values)).';
end

function nrmse = compute_nrmse_fast(phi1, phi2, y, coef, term_m, term_n)
    y = y(:);
    n_samples = numel(y);

    A = exp(1i * (phi1(:) * term_m.' + phi2(:) * term_n.'));
    y_pred = real(A * coef);
    residual = y - y_pred;
    sse = sum(residual.^2);
    y_std = std(y);

    if y_std > 0
        nrmse = sqrt(sse / n_samples) / y_std;
    else
        nrmse = sqrt(sse / n_samples);
    end
end

% ==================== Integration of extract_phase_pairs_for_rank_analysis ====================

function out = extract_phase_pairs_for_rank_analysis_inline(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, M, N)
    if nargin < 6 || isempty(file_indices)
        file_indices = [];
    end
    if nargin < 7 || isempty(M)
        M = 20;
    end
    if nargin < 8 || isempty(N)
        N = 20;
    end

    ensure_local_function_folder_on_path();

    dirpath = resolve_analysis_dirpath(dirpath);
    if ~isfolder(dirpath)
        error('Directory not found: %s', dirpath);
    end

    csv_paths = list_csv_paths(dirpath, file_indices);
    if isempty(csv_paths)
        error('No CSV files found in %s', dirpath);
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
        z_agent_id = phase_agent_ids;  % Use both phase agents as target agents
    end

    fprintf('[Extract Phase Pairs for Rank Analysis]\n');
    fprintf('  Directory: %s\n', dirpath);
    fprintf('  Phase agents: %s\n', sprintf('%d ', phase_agent_ids));
    fprintf('  Target signal agents: %s\n', sprintf('%d ', z_agent_id));
    fprintf('  Analysis window: %.1f-%.1f sec\n', analysis_start_sec, analysis_start_sec + analysis_duration_sec);
    fprintf('  Fourier order: M=%d, N=%d\n\n', M, N);

    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});
    
    phi1_all = [];
    phi2_all = [];
    a2_raw_all = [];
    a2_normalized_all = [];
    a2_normalized_agent2_all = [];
    a2_normalized_agent4_all = [];
    sin_phi2_a2_all = [];
    sin_phi2_a2_agent2_all = [];
    sin_phi2_a2_agent4_all = [];
    time_all = [];
    file_id_all = [];

    sample_dt = 0.01;
    tau = 3.246312408709;
    derived_signal_func = @(phase_target, a2_normalized) +sin(phase_target + tau) .* a2_normalized;

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            % Get data for agent 2
            [point_data_2, meta_2] = compute_points_for_csv( ...
                csv_path, phase_agent_ids, 2, analysis_duration_sec, analysis_start_sec, sample_dt);
            % Get data for agent 4
            [point_data_4, meta_4] = compute_points_for_csv( ...
                csv_path, phase_agent_ids, 4, analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            fprintf('  Skipped %s: %s\n', csv_path, ME.message);
            skipped_files(end+1) = struct('file_path', csv_path, 'reason', ME.message); %#ok<AGROW>
            continue;
        end

        if isempty(point_data_2.time) && isempty(point_data_4.time)
            fprintf('  Skipped %s: No valid points after processing.\n', csv_path);
            skipped_files(end+1) = struct('file_path', csv_path, 'reason', 'No valid points after processing.'); %#ok<AGROW>
            continue;
        end
        
        % Use agent 2 data as reference for phase
        point_data = point_data_2;

        % Calculate derived signals for both agents
        phase_for_derived_signal_2 = select_phase_for_target_agent(point_data_2, phase_agent_ids, 2);
        sin_phi2_a2_2 = derived_signal_func(phase_for_derived_signal_2, point_data_2.a2_normalized);
        
        phase_for_derived_signal_4 = select_phase_for_target_agent(point_data_4, phase_agent_ids, 4);
        sin_phi2_a2_4 = derived_signal_func(phase_for_derived_signal_4, point_data_4.a2_normalized);

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        a2_raw_all = [a2_raw_all; point_data.a2_raw(:)]; %#ok<AGROW>
        a2_normalized_all = [a2_normalized_all; point_data.a2_normalized(:)]; %#ok<AGROW>
        a2_normalized_agent2_all = [a2_normalized_agent2_all; point_data_2.a2_normalized(:)]; %#ok<AGROW>
        a2_normalized_agent4_all = [a2_normalized_agent4_all; point_data_4.a2_normalized(:)]; %#ok<AGROW>
        sin_phi2_a2_all = [sin_phi2_a2_all; sin_phi2_a2_2(:)]; %#ok<AGROW>
        sin_phi2_a2_agent2_all = [sin_phi2_a2_agent2_all; sin_phi2_a2_2(:)]; %#ok<AGROW>
        sin_phi2_a2_agent4_all = [sin_phi2_a2_agent4_all; sin_phi2_a2_4(:)]; %#ok<AGROW>
        time_all = [time_all; point_data.time(:)]; %#ok<AGROW>
        file_id_all = [file_id_all; i * ones(numel(point_data.time), 1)]; %#ok<AGROW>

        used_files{end+1} = csv_path; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available.');
    end

    fprintf('  Loaded: %d samples from %d files\n\n', numel(phi1_all), numel(used_files));

    out = struct();
    out.dirpath = dirpath;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.M = M;
    out.N = N;
    out.used_files = used_files;
    out.skipped_files = skipped_files;
    out.phi1_all = phi1_all;
    out.phi2_all = phi2_all;
    out.a2_raw_all = a2_raw_all;
    out.a2_normalized_all = a2_normalized_all;
    out.a2_normalized_agent2 = a2_normalized_agent2_all;
    out.a2_normalized_agent4 = a2_normalized_agent4_all;
    out.sin_phi2_a2_all = sin_phi2_a2_all;
    out.sin_phi2_a2_agent2 = sin_phi2_a2_agent2_all;
    out.sin_phi2_a2_agent4 = sin_phi2_a2_agent4_all;
    out.signal_agent_1 = sin_phi2_a2_agent2_all;  % Alias for agent 2
    out.signal_agent_2 = sin_phi2_a2_agent4_all;  % Alias for agent 4
    out.time_all = time_all;
    out.file_id_all = file_id_all;
end

function resolved_dirpath = resolve_analysis_dirpath(dirpath)
    if isfolder(dirpath)
        resolved_dirpath = dirpath;
        return;
    end

    candidate_pwd = fullfile(pwd, dirpath);
    if isfolder(candidate_pwd)
        resolved_dirpath = candidate_pwd;
        return;
    end

    base_dir = fileparts(mfilename('fullpath'));
    candidate_base = fullfile(base_dir, dirpath);
    if isfolder(candidate_base)
        resolved_dirpath = candidate_base;
        return;
    end

    resolved_dirpath = candidate_pwd;
end

function ensure_local_function_folder_on_path()
    local_dir = fileparts(mfilename('fullpath'));
    if isempty(which('fitDoubleFourierScatter')) && ~contains(path, local_dir)
        addpath(local_dir);
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

function phase_target = select_phase_for_target_agent(point_data, phase_agent_ids, target_agent_id)
    if target_agent_id == phase_agent_ids(1)
        phase_target = point_data.phase1;
    elseif target_agent_id == phase_agent_ids(2)
        phase_target = point_data.phase2;
    else
        phase_target = point_data.phase2;
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

% ==================== Plotting Functions ====================

function plot_rank_approximations(phi1_all, phi2_all, y_all, rank_results, m_values, n_values, agent_id, title_suffix)
    if nargin < 8
        title_suffix = '';
    end
    % Plot colormap of s_j(phi1, phi2) for ranks 1, 2, 3
    
    % Create grid for evaluation
    n_grid = 64;
    phi1_grid = linspace(0, 2*pi, n_grid);
    phi2_grid = linspace(0, 2*pi, n_grid);
    [PHI1_grid, PHI2_grid] = meshgrid(phi1_grid, phi2_grid);
    
    % Create figure with 3 subplots
    fig = figure('Color', 'w', 'Position', [100, 100, 1200, 400], ...
        'Name', sprintf('Rank Approximations - Agent %d%s', agent_id, title_suffix));
    
    for rank_idx = 1:3
        if rank_idx > numel(rank_results)
            break;
        end
        
        result = rank_results(rank_idx);
        rank = result.rank;
        coef = result.coef;
        
        % Evaluate on grid
        term_m = repelem(m_values(:), numel(n_values));
        term_n = repmat(n_values(:), numel(m_values), 1);
        
        n_phi1 = numel(phi1_grid);
        n_phi2 = numel(phi2_grid);
        Z = zeros(n_phi2, n_phi1);
        
        for i = 1:n_phi1
            for j = 1:n_phi2
                A = exp(1i * (phi1_grid(i) * term_m.' + phi2_grid(j) * term_n.'));
                Z(j, i) = real(A * coef);
            end
        end
        
        % Plot
        ax = subplot(1, 3, rank_idx);
        imagesc(phi1_grid, phi2_grid, Z);
        set(gca, 'YDir', 'normal');
        colorbar;
        
        % Create symmetric red-blue colormap
        n_colors = 256;
        cmap_blue = [linspace(0, 1, n_colors/2).' linspace(0, 0, n_colors/2).' linspace(1, 0.5, n_colors/2).'];
        cmap_red = [linspace(0.5, 1, n_colors/2).' linspace(0, 0, n_colors/2).' linspace(0, 0, n_colors/2).'];
        cmap = [flipud(cmap_blue); cmap_red];
        colormap(ax, cmap);
        
        % Set limits symmetrically
        z_max = max(abs(Z(:)));
        caxis([-z_max, z_max]);
        
        xlabel(ax, '$$\phi_1$$', 'Interpreter', 'latex');
        ylabel(ax, '$$\phi_2$$', 'Interpreter', 'latex');
        title(ax, sprintf('Rank %d (NRMSE=%.4f)', rank, result.nrmse), 'Interpreter', 'latex');
        
        % Set axis ticks
        set(ax, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax, 'YTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax, 'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        set(ax, 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        
        axis(ax, 'square');
    end
end

function plot_rank_approximations_3d(phi1_all, phi2_all, y_all, rank_results, m_values, n_values, agent_id, show_points, colormap_name, title_suffix)
    if nargin < 10
        title_suffix = '';
    end
    % Plot 3D surface of s_j(phi1, phi2) for ranks 1, 2, 3 with optional original data points
    %   show_points: logical flag to display original data (default: true)
    %   colormap_name: string name of colormap to use (default: 'jet')
    
    % Create grid for evaluation
    n_grid = 64;
    phi1_grid = linspace(0, 2*pi, n_grid);
    phi2_grid = linspace(0, 2*pi, n_grid);
    [PHI1_grid, PHI2_grid] = meshgrid(phi1_grid, phi2_grid);
    
    % Create figure with 3 subplots
    fig = figure('Color', 'w', 'Position', [100, 550, 1200, 400], ...
        'Name', sprintf('Rank Approximations 3D - Agent %d%s', agent_id, title_suffix));
    
    for rank_idx = 1:3
        if rank_idx > numel(rank_results)
            break;
        end
        
        result = rank_results(rank_idx);
        rank = result.rank;
        coef = result.coef;
        
        % Evaluate on grid
        term_m = repelem(m_values(:), numel(n_values));
        term_n = repmat(n_values(:), numel(m_values), 1);
        
        n_phi1 = numel(phi1_grid);
        n_phi2 = numel(phi2_grid);
        Z = zeros(n_phi2, n_phi1);
        
        for i = 1:n_phi1
            for j = 1:n_phi2
                A = exp(1i * (phi1_grid(i) * term_m.' + phi2_grid(j) * term_n.'));
                Z(j, i) = real(A * coef);
            end
        end
        
        % Plot 3D surface
        ax = subplot(1, 3, rank_idx);
        surf(ax, PHI1_grid, PHI2_grid, Z, 'EdgeColor', 'none');
        hold(ax, 'on');
        
        % Overlay original data points as scatter (optional)
        if show_points
            scatter3(ax, phi1_all, phi2_all, y_all, 15, y_all, 'filled', 'MarkerEdgeColor', 'k', ...
                'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.5);
        end
        
        % Apply colormap
        colormap(ax, colormap_name);
        colorbar;
        
        xlabel(ax, '$$\phi_1$$', 'Interpreter', 'latex');
        ylabel(ax, '$$\phi_2$$', 'Interpreter', 'latex');
        zlabel(ax, '$$s_j(\phi_1,\phi_2)$$', 'Interpreter', 'latex');
        title(ax, sprintf('Rank %d (NRMSE=%.4f)', rank, result.nrmse), 'Interpreter', 'latex');
        
        % Set axis ticks
        set(ax, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax, 'YTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax, 'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        set(ax, 'YTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
        
        view(ax, 45, 30);
        grid(ax, 'on');
        hold(ax, 'off');
    end
end

function components = build_separable_components(U, S, V, m_values, n_values, profile_rank)
    % Build separable components from SVD
    % Output: structure array with fields .alpha, .beta, .sigma for each rank component
    
    components = struct('alpha', {}, 'beta', {}, 'sigma', {});
    n_ranks = min(profile_rank, size(U, 2));
    
    for r = 1:n_ranks
        sigma_r = S(r, r);
        alpha_r = sqrt(sigma_r) * U(:, r);
        beta_r = sqrt(sigma_r) * conj(V(:, r));
        
        % Phase adjustment: make maximum of |alpha_r| be positive real
        phi_grid_temp = linspace(0, 2*pi, 512);
        a_r_temp = exp(1i * phi_grid_temp(:) * m_values(:).') * alpha_r;
        [~, idx_max] = max(abs(a_r_temp));
        phase_shift = exp(-1i * angle(a_r_temp(idx_max)));
        alpha_r = phase_shift * alpha_r;
        beta_r = conj(phase_shift) * beta_r;
        
        components(r).alpha = alpha_r;
        components(r).beta = beta_r;
        components(r).sigma = sigma_r;
    end
end

function plot_separable_profiles(components, m_values, n_values, agent_id, profile_rank, title_suffix)
    if nargin < 6
        title_suffix = '';
    end
    % Plot separable components a_r(phi1) and b_r(phi2)
    % Each rank gets a row: left column = a_r(phi1), right column = b_r(phi2)
    
    n_components = numel(components);
    n_display = min(n_components, profile_rank);
    
    % Create figure with tiled layout
    fig = figure('Color', 'w', 'Position', [100, 100, 1000, 300*n_display], ...
        'Name', sprintf('Separable Profiles - Agent %d%s', agent_id, title_suffix));
    tiledlayout(n_display, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Grid for evaluation
    phi_grid = linspace(0, 2*pi, 512);
    
    for r = 1:n_display
        comp = components(r);
        sigma_r = comp.sigma;
        
        % Evaluate a_r(phi1) and b_r(phi2)
        a_r = exp(1i * phi_grid(:) * m_values(:).') * comp.alpha;
        b_r = exp(1i * phi_grid(:) * n_values(:).') * comp.beta;
        
        % Left subplot: a_r(phi1)
        ax_a = nexttile;
        hold(ax_a, 'on');
        plot(ax_a, phi_grid, real(a_r), '-', 'DisplayName', 'Real', 'LineWidth', 1.5);
        plot(ax_a, phi_grid, imag(a_r), '--', 'DisplayName', 'Imag', 'LineWidth', 1.5);
        plot(ax_a, phi_grid, abs(a_r), ':', 'DisplayName', 'Abs', 'LineWidth', 2);
        hold(ax_a, 'off');
        
        xlabel(ax_a, '$$\phi_1$$', 'Interpreter', 'latex');
        ylabel(ax_a, '$$a_r(\phi_1)$$', 'Interpreter', 'latex');
        title(ax_a, sprintf('Rank %d: $a_{%d}(\\phi_1)$ ($\\sigma_{%d}=%.4e$)', r, r, r, sigma_r), 'Interpreter', 'latex');
        set(ax_a, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax_a, 'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, 'TickLabelInterpreter', 'latex');
        grid(ax_a, 'on');
        legend(ax_a, 'Location', 'best');
        
        % Right subplot: b_r(phi2)
        ax_b = nexttile;
        hold(ax_b, 'on');
        plot(ax_b, phi_grid, real(b_r), '-', 'DisplayName', 'Real', 'LineWidth', 1.5);
        plot(ax_b, phi_grid, imag(b_r), '--', 'DisplayName', 'Imag', 'LineWidth', 1.5);
        plot(ax_b, phi_grid, abs(b_r), ':', 'DisplayName', 'Abs', 'LineWidth', 2);
        hold(ax_b, 'off');
        
        xlabel(ax_b, '$$\phi_2$$', 'Interpreter', 'latex');
        ylabel(ax_b, '$$b_r(\phi_2)$$', 'Interpreter', 'latex');
        title(ax_b, sprintf('Rank %d: $b_{%d}(\\phi_2)$', r, r), 'Interpreter', 'latex');
        set(ax_b, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
        set(ax_b, 'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, 'TickLabelInterpreter', 'latex');
        grid(ax_b, 'on');
        legend(ax_b, 'Location', 'best');
    end
end

function [C_analysis, y_analysis, removed_info] = remove_phase_marginal_terms( ...
    C_full, phi1, phi2, y, m_values, n_values, ...
    target_agent_id, phase_agent_ids, ...
    remove_self_only, remove_constant, remove_other_only)
    
    % Initialize
    C_analysis = C_full;
    C_removed_self = zeros(size(C_full));
    C_removed_constant = zeros(size(C_full));
    C_removed_other = zeros(size(C_full));
    
    idx_m0 = find(m_values == 0);
    idx_n0 = find(n_values == 0);
    
    % Determine which signal this is (s_1 or s_2)
    if target_agent_id == phase_agent_ids(1)
        signal_index = 1;
        self_phase_index = 1;  % phi1
        self_phase_name = 'phi1';
    else
        signal_index = 2;
        self_phase_index = 2;  % phi2
        self_phase_name = 'phi2';
    end
    
    % Remove terms based on signal index and options
    if signal_index == 1
        % s_1: self phase is phi1 (m index)
        % Self-only: m ~= 0, n = 0
        % Other-only: m = 0, n ~= 0
        
        if remove_self_only
            C_removed_self(m_values ~= 0, idx_n0) = C_full(m_values ~= 0, idx_n0);
            C_analysis(m_values ~= 0, idx_n0) = 0;
        end
        
        if remove_other_only
            C_removed_other(idx_m0, n_values ~= 0) = C_full(idx_m0, n_values ~= 0);
            C_analysis(idx_m0, n_values ~= 0) = 0;
        end
    else
        % s_2: self phase is phi2 (n index)
        % Self-only: m = 0, n ~= 0
        % Other-only: m ~= 0, n = 0
        
        if remove_self_only
            C_removed_self(idx_m0, n_values ~= 0) = C_full(idx_m0, n_values ~= 0);
            C_analysis(idx_m0, n_values ~= 0) = 0;
        end
        
        if remove_other_only
            C_removed_other(m_values ~= 0, idx_n0) = C_full(m_values ~= 0, idx_n0);
            C_analysis(m_values ~= 0, idx_n0) = 0;
        end
    end
    
    % Remove constant term
    if remove_constant
        C_removed_constant(idx_m0, idx_n0) = C_full(idx_m0, idx_n0);
        C_analysis(idx_m0, idx_n0) = 0;
    end
    
    % Combine removed components
    C_removed_total = C_removed_self + C_removed_constant + C_removed_other;
    
    % Evaluate removed and residual signals
    y_removed_total = evaluate_fourier_from_C(phi1, phi2, C_removed_total, m_values, n_values);
    y_analysis = y(:) - y_removed_total(:);
    
    % Verify reconstruction
    y_analysis_from_C = evaluate_fourier_from_C(phi1, phi2, C_analysis, m_values, n_values);
    reconstruction_mismatch_rmse = sqrt(mean((y_analysis(:) - y_analysis_from_C(:)).^2));
    
    % Calculate individual removed signals
    y_removed_self = evaluate_fourier_from_C(phi1, phi2, C_removed_self, m_values, n_values);
    y_removed_constant = evaluate_fourier_from_C(phi1, phi2, C_removed_constant, m_values, n_values);
    y_removed_other = evaluate_fourier_from_C(phi1, phi2, C_removed_other, m_values, n_values);
    
    % Calculate energy metrics
    total_coeff_energy = sum(abs(C_full(:)).^2);
    removed_self_coeff_energy = sum(abs(C_removed_self(:)).^2);
    removed_constant_energy = sum(abs(C_removed_constant(:)).^2);
    removed_other_coeff_energy = sum(abs(C_removed_other(:)).^2);
    removed_total_coeff_energy = sum(abs(C_removed_total(:)).^2);
    remaining_coeff_energy = sum(abs(C_analysis(:)).^2);
    
    % Store results
    removed_info = struct();
    removed_info.remove_self_only = remove_self_only;
    removed_info.remove_constant = remove_constant;
    removed_info.remove_other_only = remove_other_only;
    removed_info.self_phase = self_phase_name;
    removed_info.self_phase_index = self_phase_index;
    removed_info.target_agent_id = target_agent_id;
    removed_info.phase_agent_ids = phase_agent_ids;
    removed_info.signal_index = signal_index;
    removed_info.C_removed_self = C_removed_self;
    removed_info.C_removed_constant = C_removed_constant;
    removed_info.C_removed_other = C_removed_other;
    removed_info.C_removed_total = C_removed_total;
    removed_info.y_removed_self = y_removed_self;
    removed_info.y_removed_constant = y_removed_constant;
    removed_info.y_removed_other = y_removed_other;
    removed_info.y_removed_total = y_removed_total;
    removed_info.y_analysis = y_analysis;
    removed_info.removed_self_coeff_energy = removed_self_coeff_energy;
    removed_info.removed_constant_energy = removed_constant_energy;
    removed_info.removed_other_coeff_energy = removed_other_coeff_energy;
    removed_info.removed_total_coeff_energy = removed_total_coeff_energy;
    removed_info.remaining_coeff_energy = remaining_coeff_energy;
    removed_info.total_coeff_energy = total_coeff_energy;
    removed_info.removed_self_energy_ratio = removed_self_coeff_energy / total_coeff_energy;
    removed_info.removed_constant_energy_ratio = removed_constant_energy / total_coeff_energy;
    removed_info.removed_other_energy_ratio = removed_other_coeff_energy / total_coeff_energy;
    removed_info.removed_total_energy_ratio = removed_total_coeff_energy / total_coeff_energy;
    removed_info.remaining_energy_ratio = remaining_coeff_energy / total_coeff_energy;
    removed_info.reconstruction_mismatch_rmse = reconstruction_mismatch_rmse;
end

function y_hat = evaluate_fourier_from_C(phi1, phi2, C, m_values, n_values)
    % Evaluate Fourier reconstruction from coefficient matrix C
    % Uses chunk processing for memory efficiency
    
    coef = reshape(C.', [], 1);
    term_m = repelem(m_values(:), numel(n_values));
    term_n = repmat(n_values(:), numel(m_values), 1);
    term_n = term_n(:);
    
    phi1 = phi1(:);
    phi2 = phi2(:);
    n_samples = numel(phi1);
    chunk_size = 20000;
    y_hat = zeros(n_samples, 1);
    
    for k0 = 1:chunk_size:n_samples
        k_end = min(k0 + chunk_size - 1, n_samples);
        k_idx = k0:k_end;
        
        % Create design matrix for chunk
        A_chunk = exp(1i * (phi1(k_idx) * term_m.' + phi2(k_idx) * term_n.'));
        y_hat(k_idx) = real(A_chunk * coef);
    end
end

function self_profile = compute_self_only_profile(C_full, m_values, n_values, target_agent_id, phase_agent_ids)
    % Compute self-phase-only profile
    
    phi_grid = linspace(0, 2*pi, 512).';
    
    % Determine which signal and self phase
    if target_agent_id == phase_agent_ids(1)
        % Signal s_1: self phase is phi1 (m index)
        % Self-only coefficients: m ~= 0, n = 0
        signal_index = 1;
        self_phase_index = 1;  % phi1
        self_phase_name = 'phi1';
        
        idx_n0 = find(n_values == 0);
        idx_m_nonzero = find(m_values ~= 0);
        
        if ~isempty(idx_n0) && ~isempty(idx_m_nonzero)
            % Extract coefficients for m~=0, n=0
            coeff_vec = C_full(idx_m_nonzero, idx_n0(:));  % [n_m_nonzero, n_n0]
            coeff_vec = coeff_vec(:);  % Flatten to column vector
            m_basis = m_values(idx_m_nonzero);
            E = exp(1i * phi_grid * m_basis(:).');
        else
            coeff_vec = 0;
            E = zeros(numel(phi_grid), 1);
        end
    else
        % Signal s_2: self phase is phi2 (n index)
        % Self-only coefficients: m = 0, n ~= 0
        signal_index = 2;
        self_phase_index = 2;  % phi2
        self_phase_name = 'phi2';
        
        idx_m0 = find(m_values == 0);
        idx_n_nonzero = find(n_values ~= 0);
        
        if ~isempty(idx_m0) && ~isempty(idx_n_nonzero)
            % Extract coefficients for m=0, n~=0
            coeff_vec = C_full(idx_m0(:), idx_n_nonzero);  % [n_m0, n_n_nonzero]
            coeff_vec = coeff_vec(:);  % Flatten to column vector
            n_basis = n_values(idx_n_nonzero);
            E = exp(1i * phi_grid * n_basis(:).');
        else
            coeff_vec = 0;
            E = zeros(numel(phi_grid), 1);
        end
    end
    
    % Evaluate profile
    if numel(E) == 1
        profile_values = zeros(size(phi_grid));
    else
        profile_values = E * coeff_vec;
    end
    
    % Store results
    self_profile = struct();
    self_profile.phi_grid = phi_grid;
    self_profile.values = profile_values;
    self_profile.real = real(profile_values);
    self_profile.imag = imag(profile_values);
    self_profile.abs = abs(profile_values);
    self_profile.self_phase = self_phase_name;
    self_profile.self_phase_index = self_phase_index;
    self_profile.target_agent_id = target_agent_id;
    self_profile.signal_index = signal_index;
    self_profile.max_abs_imag = max(abs(imag(profile_values)));
end

function fig = plot_self_only_profile(self_profile)
    % Plot self-phase-only profile
    
    phi_grid = self_profile.phi_grid;
    signal_index = self_profile.signal_index;
    target_agent_id = self_profile.target_agent_id;
    self_phase = self_profile.self_phase;
    
    % Create figure
    fig = figure('Color', 'w', 'Position', [100, 50, 800, 400], ...
        'Name', sprintf('Self-only Profile - s_%d (Agent %d)', signal_index, target_agent_id));
    
    ax = axes(fig);
    hold(ax, 'on');
    plot(ax, phi_grid, self_profile.real, '-', 'DisplayName', 'Real', 'LineWidth', 1.5);
    plot(ax, phi_grid, self_profile.imag, '--', 'DisplayName', 'Imag', 'LineWidth', 1.5);
    plot(ax, phi_grid, self_profile.abs, ':', 'DisplayName', 'Abs', 'LineWidth', 2);
    hold(ax, 'off');
    
    % Labels and title
    xlabel(ax, sprintf('$$%s$$', self_phase), 'Interpreter', 'latex');
    ylabel(ax, '$$s_{j,\mathrm{self}}$$', 'Interpreter', 'latex');
    title(ax, sprintf('Self-only profile removed from $s_{%d}$ (Agent %d, self phase = %s, max |imag| = %.2e)', ...
        signal_index, target_agent_id, self_phase, self_profile.max_abs_imag), ...
        'Interpreter', 'latex');
    
    % X-axis ticks
    set(ax, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
    set(ax, 'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, 'TickLabelInterpreter', 'latex');
    
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end
