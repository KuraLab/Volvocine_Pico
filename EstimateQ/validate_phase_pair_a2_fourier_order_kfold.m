function out = validate_phase_pair_a2_fourier_order_kfold(dirpath, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, file_indices, varargin)
% Validate Fourier truncation order for phase-pair -> a2 fitting via K-fold CV.
%
% This function extracts only the first fitting pipeline used in
% plot_phase_pair_a2_3d_all_files.m:
%   1) build (phi1, phi2, a2_norm) point cloud from CSVs
%   2) fit unweighted double-Fourier model Q(phi1,phi2)
%
% Then it performs model-order validation:
%   - K-fold CV RMSE for orders k = 0..Mmax with M = N = k
%   - per-order shell RMS contribution at Mmax
%   - comparison of shell RMS against residual noise floor
%
% Usage:
%   out = validate_phase_pair_a2_fourier_order_kfold()
%   out = validate_phase_pair_a2_fourier_order_kfold(dirpath)
%   out = validate_phase_pair_a2_fourier_order_kfold(..., 'Mmax', 20, 'Kfold', 5)
%
% Name-value options:
%   'Mmax'        : maximum order (default 20)
%   'Kfold'       : number of folds (default 5)
%   'SampleDt'    : interpolation step [sec] (default 0.01)
%   'RandomSeed'  : RNG seed for fold split (default 42)
%   'ShowFigure'  : draw summary figure (default true)
%
% Output:
%   out: struct with configuration, CV statistics, shell-vs-noise summary,
%        selected-order suggestions, and collected point cloud.

    if nargin < 1 || isempty(dirpath)
        dirpath = fullfile('Spring3', '255');
    end
    if nargin < 2
        phase_agent_ids = [];
    end
    if nargin < 3
        z_agent_id = [];
    end
    if nargin < 4 || isempty(analysis_duration_sec)
        analysis_duration_sec = 80;
    end
    if nargin < 5 || isempty(analysis_start_sec)
        analysis_start_sec = 10;
    end
    if nargin < 6
        file_indices = [];
    end

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'Mmax', 20, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}));
    addParameter(p, 'Kfold', 5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 2, 'finite'}));
    addParameter(p, 'SampleDt', 0.01, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, 'finite'}));
    addParameter(p, 'RandomSeed', 42, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'finite'}));
    addParameter(p, 'ShowFigure', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    Mmax = double(p.Results.Mmax);
    Kfold = double(p.Results.Kfold);
    sample_dt = double(p.Results.SampleDt);
    random_seed = double(p.Results.RandomSeed);
    show_figure = p.Results.ShowFigure;

    ensure_required_paths_on_path();

    if ~isscalar(analysis_duration_sec) || analysis_duration_sec <= 0 || ~isfinite(analysis_duration_sec)
        error('analysis_duration_sec must be a positive scalar.');
    end
    if ~isscalar(analysis_start_sec) || analysis_start_sec < 0 || ~isfinite(analysis_start_sec)
        error('analysis_start_sec must be a non-negative scalar.');
    end

    dirpath = resolve_analysis_dirpath(dirpath);
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

    [point_cloud, used_files, skipped_files] = build_point_cloud( ...
        csv_paths, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt);

    n_samples = numel(point_cloud.a2);
    if n_samples < Kfold
        error('Not enough samples (%d) for Kfold=%d.', n_samples, Kfold);
    end

    fold_id = create_kfold_partition(n_samples, Kfold, random_seed);
    orders = (0:Mmax).';
    order_col_end = arrayfun(@(k) num_double_fourier_terms(k, k), orders);

    fold_rmse = nan(numel(orders), Kfold);
    n_train = zeros(numel(orders), Kfold);
    n_valid = zeros(numel(orders), Kfold);
    n_basis = order_col_end;

    A_train_max_by_fold = cell(Kfold, 1);
    A_valid_max_by_fold = cell(Kfold, 1);
    z_train_by_fold = cell(Kfold, 1);
    z_valid_by_fold = cell(Kfold, 1);
    n_train_by_fold = zeros(Kfold, 1);
    n_valid_by_fold = zeros(Kfold, 1);

    for f = 1:Kfold
        valid_mask = (fold_id == f);
        train_mask = ~valid_mask;

        phi1_train = point_cloud.phi1(train_mask);
        phi2_train = point_cloud.phi2(train_mask);
        z_train = point_cloud.a2(train_mask);

        phi1_valid = point_cloud.phi1(valid_mask);
        phi2_valid = point_cloud.phi2(valid_mask);
        z_valid = point_cloud.a2(valid_mask);

        A_train_max_by_fold{f} = build_double_fourier_design_matrix_equal_order(phi1_train, phi2_train, Mmax);
        A_valid_max_by_fold{f} = build_double_fourier_design_matrix_equal_order(phi1_valid, phi2_valid, Mmax);
        z_train_by_fold{f} = z_train;
        z_valid_by_fold{f} = z_valid;
        n_train_by_fold(f) = nnz(train_mask);
        n_valid_by_fold(f) = nnz(valid_mask);
    end

    for oi = 1:numel(orders)
        n_cols = order_col_end(oi);

        for f = 1:Kfold
            z_train = z_train_by_fold{f};
            z_valid = z_valid_by_fold{f};

            n_train(oi, f) = n_train_by_fold(f);
            n_valid(oi, f) = n_valid_by_fold(f);

            if n_train(oi, f) <= n_cols
                continue;
            end

            z_mean_train = mean(z_train);
            z_train_centered = z_train - z_mean_train;

            A_train = A_train_max_by_fold{f}(:, 1:n_cols);
            coeff = A_train \ z_train_centered;

            A_valid = A_valid_max_by_fold{f}(:, 1:n_cols);
            z_valid_hat = A_valid * coeff + z_mean_train;

            fold_rmse(oi, f) = sqrt(mean((z_valid - z_valid_hat) .^ 2));
        end
    end

    cv_rmse_mean = mean(fold_rmse, 2, 'omitnan');
    cv_rmse_std = std(fold_rmse, 0, 2, 'omitnan');
    cv_n_valid_folds = sum(isfinite(fold_rmse), 2);

    [best_cv_rmse, best_idx] = min(cv_rmse_mean);
    best_order = orders(best_idx);

    one_se_threshold = best_cv_rmse + cv_rmse_std(best_idx);
    one_se_candidates = find(cv_rmse_mean <= one_se_threshold);
    if isempty(one_se_candidates)
        one_se_order = best_order;
    else
        one_se_order = orders(one_se_candidates(1));
    end

    % Full-data fit at Mmax for shell-vs-noise-floor analysis.
    phi1_all = point_cloud.phi1;
    phi2_all = point_cloud.phi2;
    z_all = point_cloud.a2;

    [A_full, basis_meta] = build_double_fourier_design_matrix_equal_order(phi1_all, phi2_all, Mmax);
    z_mean_full = mean(z_all);
    z_centered_full = z_all - z_mean_full;
    coeff_full = A_full \ z_centered_full;
    residual_full = z_centered_full - A_full * coeff_full;

    noise_floor_rms = sqrt(mean(residual_full .^ 2));
    noise_floor_std = std(residual_full, 0);

    shell_orders = (0:Mmax).';
    shell_n_terms = zeros(numel(shell_orders), 1);
    shell_rms = zeros(numel(shell_orders), 1);

    for si = 1:numel(shell_orders)
        shell_order = shell_orders(si);
        shell_mask = (basis_meta.shell_order == shell_order);
        shell_n_terms(si) = nnz(shell_mask);
        if shell_n_terms(si) == 0
            shell_rms(si) = 0;
            continue;
        end
        shell_signal = A_full(:, shell_mask) * coeff_full(shell_mask);
        shell_rms(si) = sqrt(mean(shell_signal .^ 2));
    end

    shell_to_noise_ratio = shell_rms / max(noise_floor_rms, eps);
    above_noise_idx = find(shell_to_noise_ratio > 1, 1, 'last');
    if isempty(above_noise_idx)
        noise_floor_order = 0;
    else
        noise_floor_order = shell_orders(above_noise_idx);
    end

    cv_table = table( ...
        orders, n_basis, cv_n_valid_folds, cv_rmse_mean, cv_rmse_std, ...
        'VariableNames', {'order', 'n_basis', 'n_valid_folds', 'cv_rmse_mean', 'cv_rmse_std'});

    shell_table = table( ...
        shell_orders, shell_n_terms, shell_rms, shell_to_noise_ratio, ...
        'VariableNames', {'order', 'n_terms', 'shell_rms', 'shell_to_noise_ratio'});

    if show_figure
        fig_cv = figure('Color', 'w', 'Name', 'Fourier order validation: CV');
        ax_cv = axes('Parent', fig_cv);
        hold(ax_cv, 'on');
        errorbar(ax_cv, orders, cv_rmse_mean, cv_rmse_std, 'o-', 'LineWidth', 3, 'MarkerSize', 4, ...
            'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
        h_best = xline(ax_cv, best_order, '--', 'Best', 'LabelOrientation', 'horizontal', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth', 3);
        h_one_se = xline(ax_cv, one_se_order, ':', '1SE', 'LabelOrientation', 'horizontal', 'Color', [0.4940, 0.1840, 0.5560],'LineWidth', 3);
        h_one_se_th = yline(ax_cv, one_se_threshold, ':', '1SE threshold', 'Color', [0.3, 0.3, 0.3],'LineWidth', 3);
        set([h_best, h_one_se, h_one_se_th], 'FontSize', 20);
        grid(ax_cv, 'on');
        xlabel(ax_cv, 'Order k (M=N=k)');
        ylabel(ax_cv, 'Validation RMSE');
        %title(ax_cv, sprintf('K-fold CV (K=%d)', Kfold));
        xlim(ax_cv, [0, Mmax]);
        tuneFigure;

        fig_shell = figure('Color', 'w', 'Name', 'Fourier order validation: shell vs noise');
        ax_shell = axes('Parent', fig_shell);
        hold(ax_shell, 'on');
        stem(ax_shell, shell_orders, shell_rms, 'filled', 'LineWidth', 3, 'Color', [0.4660, 0.6740, 0.1880]);
        h_noise_floor = yline(ax_shell, noise_floor_rms, '--', 'noise floor (residual RMS at Mmax)', 'Color', [0.8500, 0.3250, 0.0980],'LineWidth', 3);
        h_last_above = xline(ax_shell, noise_floor_order, ':', 'last > noise', 'LabelOrientation', 'horizontal', 'Color', [0.2, 0.2, 0.2],'LineWidth', 3);
        set([h_noise_floor, h_last_above], 'FontSize', 20);
        grid(ax_shell, 'on');
        xlabel(ax_shell, 'Shell order h = max(m,n)');
        ylabel(ax_shell, 'Shell RMS contribution');
        %title(ax_shell, 'Per-order shell contribution vs noise floor');
        xlim(ax_shell, [0, Mmax]);
        tuneFigure;
    else
        fig_cv = [];
        fig_shell = [];
    end

    out = struct();
    out.dirpath = dirpath;
    out.phase_agent_ids = phase_agent_ids;
    out.z_agent_id = z_agent_id;
    out.analysis_duration_sec = analysis_duration_sec;
    out.analysis_start_sec = analysis_start_sec;
    out.file_indices = file_indices;
    out.sample_dt = sample_dt;
    out.used_files = used_files;
    out.skipped_files = skipped_files;
    out.point_cloud = point_cloud;

    out.Mmax = Mmax;
    out.Kfold = Kfold;
    out.random_seed = random_seed;
    out.cv = struct();
    out.cv.fold_rmse = fold_rmse;
    out.cv.table = cv_table;
    out.cv.best_order = best_order;
    out.cv.best_rmse = best_cv_rmse;
    out.cv.one_se_order = one_se_order;
    out.cv.one_se_threshold = one_se_threshold;

    out.noise_floor = struct();
    out.noise_floor.rms = noise_floor_rms;
    out.noise_floor.std = noise_floor_std;
    out.noise_floor.order_last_above_noise = noise_floor_order;
    out.noise_floor.shell_table = shell_table;

    out.recommendation = struct();
    out.recommendation.from_cv_best = best_order;
    out.recommendation.from_cv_one_se = one_se_order;
    out.recommendation.from_noise_floor = noise_floor_order;
    out.figure_cv = fig_cv;
    out.figure_shell = fig_shell;
    out.figure = fig_cv;

    fprintf('[INFO] K-fold CV completed: best order=%d, 1SE order=%d, noise-floor order=%d\n', ...
        best_order, one_se_order, noise_floor_order);
    fprintf('[INFO] CV RMSE(best)=%.6g, CV threshold(1SE)=%.6g, residual noise RMS(Mmax) = %.6g\n', ...
        best_cv_rmse, one_se_threshold, noise_floor_rms);
end

function ensure_required_paths_on_path()
    local_dir = fileparts(mfilename('fullpath'));
    parent_dir = fileparts(local_dir);

    if isempty(which('load_corrected_agent_series_from_csv'))
        addpath(parent_dir);
    end
    if isempty(which('fitDoubleFourierScatter'))
        addpath(local_dir);
    end
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

    error('Could not determine default phase_agent_ids from selected CSV files.');
end

function [point_cloud, used_files, skipped_files] = build_point_cloud(csv_paths, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt)
    used_files = {};
    skipped_files = struct('file_path', {}, 'reason', {});

    phi1_all = [];
    phi2_all = [];
    a2_all = [];
    time_all = [];
    file_id_all = [];

    for i = 1:numel(csv_paths)
        csv_path = csv_paths{i};
        try
            point_data = compute_points_for_csv( ...
                csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt);
        catch ME
            skipped_files(end + 1) = struct('file_path', csv_path, 'reason', ME.message); %#ok<AGROW>
            continue;
        end

        if isempty(point_data.time)
            skipped_files(end + 1) = struct('file_path', csv_path, 'reason', 'No valid points after processing.'); %#ok<AGROW>
            continue;
        end

        phi1_all = [phi1_all; point_data.phase1(:)]; %#ok<AGROW>
        phi2_all = [phi2_all; point_data.phase2(:)]; %#ok<AGROW>
        a2_all = [a2_all; point_data.a2(:)]; %#ok<AGROW>
        time_all = [time_all; point_data.time(:)]; %#ok<AGROW>
        file_id_all = [file_id_all; i * ones(numel(point_data.time), 1)]; %#ok<AGROW>

        used_files{end + 1} = csv_path; %#ok<AGROW>
    end

    if isempty(used_files)
        error('No valid files were available to build point cloud.');
    end

    point_cloud = struct();
    point_cloud.phi1 = phi1_all;
    point_cloud.phi2 = phi2_all;
    point_cloud.a2 = a2_all;
    point_cloud.time = time_all;
    point_cloud.file_index = file_id_all;
end

function point_data = compute_points_for_csv(csv_path, phase_agent_ids, z_agent_id, analysis_duration_sec, analysis_start_sec, sample_dt)
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
        error('Requested analysis window is outside the common overlap.');
    end

    time_abs = (window_start_abs:sample_dt:window_end_abs).';
    if numel(time_abs) < 2
        error('Analysis window is too short after applying selected settings.');
    end

    a0_1 = interp1(series_by_agent(phase_agent_ids(1)).time, series_by_agent(phase_agent_ids(1)).a0_corr, time_abs, 'linear', NaN);
    a0_2 = interp1(series_by_agent(phase_agent_ids(2)).time, series_by_agent(phase_agent_ids(2)).a0_corr, time_abs, 'linear', NaN);
    a2_z = interp1(series_by_agent(z_agent_id).time, series_by_agent(z_agent_id).a2, time_abs, 'linear', NaN);

    valid = isfinite(a0_1) & isfinite(a0_2) & isfinite(a2_z);
    time_abs = time_abs(valid);
    a0_1 = a0_1(valid);
    a0_2 = a0_2(valid);
    a2_z = a2_z(valid);

    if numel(time_abs) < 2
        error('No overlapping interpolated samples in selected analysis window.');
    end

    point_data = struct();
    point_data.time = time_abs - window_start_abs;
    point_data.phase1 = mod(a0_1, 256) * (2 * pi / 256);
    point_data.phase2 = mod(a0_2, 256) * (2 * pi / 256);
    point_data.a2 = clip_values( ...
        normalize_by_agent_percentile_span(a2_z, series_by_agent(z_agent_id).a2, 10), -0.5, 0.5);
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

function fold_id = create_kfold_partition(n_samples, kfold, random_seed)
    rng(random_seed, 'twister');
    perm = randperm(n_samples);

    fold_id = zeros(n_samples, 1);
    for i = 1:n_samples
        fold_id(perm(i)) = mod(i - 1, kfold) + 1;
    end
end

function n_terms = num_double_fourier_terms(M, N)
    n_terms = 1 + 2 * M + 2 * N + 4 * M * N;
end

function [A, basis_meta] = build_double_fourier_design_matrix_equal_order(phi1, phi2, Kmax)
    phi1 = mod(phi1(:), 2 * pi);
    phi2 = mod(phi2(:), 2 * pi);
    n = numel(phi1);

    n_cols = num_double_fourier_terms(Kmax, Kmax);
    A = zeros(n, n_cols);

    if nargout > 1
        basis_meta = struct();
        basis_meta.m = zeros(n_cols, 1);
        basis_meta.n = zeros(n_cols, 1);
        basis_meta.shell_order = zeros(n_cols, 1);
        basis_meta.group = strings(n_cols, 1);
    else
        basis_meta = struct();
    end

    col = 0;

    col = col + 1;
    A(:, col) = 1;
    if nargout > 1
        basis_meta.group(col) = "constant";
    end

    % Shell-ordered basis so columns 1:num_double_fourier_terms(k,k)
    % correspond exactly to the order-k model (M=N=k).
    for h = 1:Kmax
        col = col + 1;
        A(:, col) = cos(h * phi1);
        if nargout > 1
            basis_meta.m(col) = h;
            basis_meta.group(col) = "phi1_only";
            basis_meta.shell_order(col) = h;
        end

        col = col + 1;
        A(:, col) = sin(h * phi1);
        if nargout > 1
            basis_meta.m(col) = h;
            basis_meta.group(col) = "phi1_only";
            basis_meta.shell_order(col) = h;
        end

        col = col + 1;
        A(:, col) = cos(h * phi2);
        if nargout > 1
            basis_meta.n(col) = h;
            basis_meta.group(col) = "phi2_only";
            basis_meta.shell_order(col) = h;
        end

        col = col + 1;
        A(:, col) = sin(h * phi2);
        if nargout > 1
            basis_meta.n(col) = h;
            basis_meta.group(col) = "phi2_only";
            basis_meta.shell_order(col) = h;
        end

        for m = 1:h
            n_idx = h;
            c1 = cos(m * phi1);
            s1 = sin(m * phi1);
            c2 = cos(n_idx * phi2);
            s2 = sin(n_idx * phi2);

            col = col + 1;
            A(:, col) = c1 .* c2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = c1 .* s2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = s1 .* c2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = s1 .* s2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end
        end

        for n_idx = 1:(h - 1)
            m = h;
            c1 = cos(m * phi1);
            s1 = sin(m * phi1);
            c2 = cos(n_idx * phi2);
            s2 = sin(n_idx * phi2);

            col = col + 1;
            A(:, col) = c1 .* c2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = c1 .* s2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = s1 .* c2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end

            col = col + 1;
            A(:, col) = s1 .* s2;
            if nargout > 1
                basis_meta.m(col) = m;
                basis_meta.n(col) = n_idx;
                basis_meta.group(col) = "mixed";
                basis_meta.shell_order(col) = h;
            end
        end
    end

    if col ~= n_cols
        error('Internal column-count mismatch in design matrix builder.');
    end
end
