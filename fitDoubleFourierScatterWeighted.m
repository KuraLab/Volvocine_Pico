function result = fitDoubleFourierScatterWeighted(phi1, phi2, z, M, N, beta, r, clip_percentiles)
% fitDoubleFourierScatterWeighted
% Fit a real double Fourier series to scattered periodic data using both
% ordinary least squares and k-nearest-neighbor-distance-based weighted
% least squares.
%
%   result = fitDoubleFourierScatterWeighted(phi1, phi2, z, M, N, beta, r)
%   result = fitDoubleFourierScatterWeighted(phi1, phi2, z, M, N, beta, r, clip_percentiles)
%
% Inputs:
%   phi1, phi2        : angle variables (radians). They are treated as
%                       2*pi-periodic and normalized to [0, 2*pi).
%   z                 : observed values on scattered points Q(phi1, phi2)
%   M, N              : Fourier truncation orders in phi1 and phi2
%   beta              : exponent in the density-compensation weight
%                       w_k = d_k^(2*beta)
%   r                 : neighbor order used to define d_k, the r-th nearest
%                       neighbor distance under the periodic 2D metric
%   clip_percentiles  : two-element vector [p_low, p_high] used to clip the
%                       raw weights before normalization. Default: [1 99]
%
% Periodic distance definition:
%   dphi = atan2(sin(a-b), cos(a-b))
%   dist = sqrt(dphi1.^2 + dphi2.^2)
%
% Weight definition:
%   d_k   = r-th nearest-neighbor distance for sample k
%   w_k   = d_k^(2*beta)
%   w_k is then percentile-clipped and normalized so that mean(w)=1.
%
% The function computes both:
%   coeff_unweighted = A \ z
%   coeff_weighted   = (A .* sqrt(w)) \ (z .* sqrt(w))
%
% The output struct contains coefficients, fitted values, residuals,
% unweighted/weighted RMSE values, weighted contribution summaries,
% weight diagnostics, basis metadata, evaluated fitted surfaces, and
% figure handles.
%
% Example 1:
%   out = plot_phase_pair_a2_3d_all_files(fullfile('EstimateF','Spring3','250'), [2 3], 3, 15, 5);
%   result = fitDoubleFourierScatterWeighted( ...
%       out.point_cloud.phi1, out.point_cloud.phi2, out.point_cloud.a2, ...
%       6, 6, 1.0, 10, [1 99]);
%
% Example 2:
%   result = fitDoubleFourierScatterWeighted(phi1, phi2, z, 4, 4, 0.5, 8);
%
% Notes:
% - The periodic k-nearest-neighbor computation is implemented in base
%   MATLAB using block processing, so no toolbox-specific kNN routine is
%   required.
% - If duplicated points are present and beta < 0, the raw weight formula
%   can diverge. To avoid numerical overflow, d_k is lower-bounded by a
%   small positive number when weights are formed.

    if nargin < 7
        error('Usage: fitDoubleFourierScatterWeighted(phi1, phi2, z, M, N, beta, r, clip_percentiles)');
    end
    if nargin < 8 || isempty(clip_percentiles)
        clip_percentiles = [1 99];
    end

    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'M');
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'N');
    validateattributes(beta, {'numeric'}, {'scalar', 'finite', 'real'}, mfilename, 'beta');
    validateattributes(r, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'r');
    validateattributes(clip_percentiles, {'numeric'}, {'vector', 'numel', 2, 'nondecreasing', '>=', 0, '<=', 100}, mfilename, 'clip_percentiles');

    phi1 = phi1(:);
    phi2 = phi2(:);
    z = z(:);

    if ~(numel(phi1) == numel(phi2) && numel(phi1) == numel(z))
        error('phi1, phi2, and z must have the same number of elements.');
    end

    % Normalize the angular variables to the 2*pi-periodic domain.
    phi1 = mod(phi1, 2*pi);
    phi2 = mod(phi2, 2*pi);

    % Remove invalid samples before fitting.
    valid_mask = isfinite(phi1) & isfinite(phi2) & isfinite(z);
    phi1 = phi1(valid_mask);
    phi2 = phi2(valid_mask);
    z = z(valid_mask);

    n_samples = numel(z);
    if n_samples == 0
        error('No valid samples remain after removing NaN/Inf values.');
    end
    if n_samples <= r
        error('r must satisfy r <= number_of_valid_samples - 1. Current values: r=%d, N=%d.', r, n_samples);
    end

    % Mean-center the observed values before fitting so the constant term
    % does not dominate the contribution analysis.
    z_original = z;
    z_mean = mean(z_original);
    z = z_original - z_mean;

    % Build the real trigonometric design matrix.
    [A, basis_names, basis_groups, basis_m, basis_n] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N);

    % Compute k-nearest-neighbor-distance-based weights on the periodic
    % phase torus.
    weight_info = computePeriodicKnnWeights(phi1, phi2, r, beta, clip_percentiles);
    w = weight_info.weights;
    sqrt_w = sqrt(w);

    % Ordinary least squares fit.
    coeff_unweighted = A \ z;
    z_hat_unweighted = A * coeff_unweighted;
    residual_unweighted = z - z_hat_unweighted;
    rmse_unweighted = sqrt(mean(residual_unweighted .^ 2));
    wrmse_unweighted = computeWeightedRMSE(residual_unweighted, w);

    % Weighted least squares fit using Aw = A .* sqrt(w), zw = z .* sqrt(w).
    Aw = A .* sqrt_w;
    zw = z .* sqrt_w;
    coeff_weighted = Aw \ zw;
    z_hat_weighted = A * coeff_weighted;
    residual_weighted = z - z_hat_weighted;
    rmse_weighted = sqrt(mean(residual_weighted .^ 2));
    wrmse_weighted = computeWeightedRMSE(residual_weighted, w);
    z_hat_unweighted_original_scale = z_hat_unweighted + z_mean;
    z_hat_weighted_original_scale = z_hat_weighted + z_mean;

    % Diagnostics for both systems.
    rankA = rank(A);
    condA = cond(A);
    rankAw = rank(Aw);
    condAw = cond(Aw);

    if rankA < size(A, 2)
        warning('fitDoubleFourierScatterWeighted:RankDeficientUnweighted', ...
            'Unweighted design matrix may be rank deficient: rank(A) = %d < %d.', rankA, size(A, 2));
    end
    if rankAw < size(Aw, 2)
        warning('fitDoubleFourierScatterWeighted:RankDeficientWeighted', ...
            'Weighted design matrix may be rank deficient: rank(Aw) = %d < %d.', rankAw, size(Aw, 2));
    end
    if isfinite(condA) && condA > 1e10
        warning('fitDoubleFourierScatterWeighted:IllConditionedUnweighted', ...
            'Unweighted design matrix may be ill-conditioned: cond(A) = %.3e.', condA);
    end
    if isfinite(condAw) && condAw > 1e10
        warning('fitDoubleFourierScatterWeighted:IllConditionedWeighted', ...
            'Weighted design matrix may be ill-conditioned: cond(Aw) = %.3e.', condAw);
    end

    % Weighted contribution analysis for the weighted fit.
    % Each fitted basis term is A(:,k) * coeff_weighted(k). Its contribution
    % is summarized by the weighted mean-square amplitude over the samples.
    term_matrix_weighted = bsxfun(@times, A, coeff_weighted.');
    term_energy_contribution_weighted = sum(w .* (term_matrix_weighted .^ 2), 1).' / sum(w);
    term_rms_contribution_weighted = sqrt(term_energy_contribution_weighted);

    total_term_energy_weighted = sum(term_energy_contribution_weighted);
    if total_term_energy_weighted > 0
        term_contribution_ratio_weighted = term_energy_contribution_weighted / total_term_energy_weighted;
    else
        term_contribution_ratio_weighted = zeros(size(term_energy_contribution_weighted));
    end

    contribution_table_weighted = table((1:size(A, 2)).', basis_names, basis_groups, basis_m, basis_n, coeff_weighted, ...
        term_rms_contribution_weighted, term_energy_contribution_weighted, term_contribution_ratio_weighted, ...
        'VariableNames', {'term_index', 'basis_name', 'basis_group', 'phi1_order', 'phi2_order', 'coefficient', ...
        'rms_contribution', 'energy_contribution', 'contribution_ratio'});

    group_names = {'constant'; 'phi1_only'; 'phi2_only'; 'mixed'};
    group_display_names = {'constant'; 'phi1 only'; 'phi2 only'; 'phi1-phi2 mixed'};
    group_term_count = zeros(numel(group_names), 1);
    group_energy_contribution_weighted = zeros(numel(group_names), 1);
    group_contribution_ratio_weighted = zeros(numel(group_names), 1);
    for g = 1:numel(group_names)
        group_mask = strcmp(basis_groups, group_names{g});
        group_term_count(g) = nnz(group_mask);
        group_energy_contribution_weighted(g) = sum(term_energy_contribution_weighted(group_mask));
        if total_term_energy_weighted > 0
            group_contribution_ratio_weighted(g) = group_energy_contribution_weighted(g) / total_term_energy_weighted;
        end
    end

    group_contribution_summary_weighted = table(group_names, group_display_names, group_term_count, ...
        group_energy_contribution_weighted, group_contribution_ratio_weighted, ...
        'VariableNames', {'group_name', 'display_name', 'n_terms', ...
        'energy_contribution', 'contribution_ratio'});

    phi1_harmonic_summary_weighted = summarizeSingleAxisContributions( ...
        basis_groups, basis_m, term_energy_contribution_weighted, total_term_energy_weighted, 'phi1_only');
    phi2_harmonic_summary_weighted = summarizeSingleAxisContributions( ...
        basis_groups, basis_n, term_energy_contribution_weighted, total_term_energy_weighted, 'phi2_only');
    mixed_pair_summary_weighted = summarizeMixedPairContributions( ...
        basis_groups, basis_m, basis_n, term_energy_contribution_weighted, total_term_energy_weighted);

    % Evaluate the fitted surfaces on a regular periodic grid for
    % visualization.
    surface_grid_size = 81;
    surface_phi1_values = linspace(0, 2*pi, surface_grid_size);
    surface_phi2_values = linspace(0, 2*pi, surface_grid_size);
    [surface_phi1, surface_phi2] = meshgrid(surface_phi1_values, surface_phi2_values);
    A_surface = buildDoubleFourierDesignMatrix(surface_phi1(:), surface_phi2(:), M, N);
    surface_z_unweighted = reshape(A_surface * coeff_unweighted, size(surface_phi1)) + z_mean;
    surface_z_weighted = reshape(A_surface * coeff_weighted, size(surface_phi1)) + z_mean;

    % Visualize the normalized weights and their distribution.
    fig_weights = figure('Color', 'w');
    tiled_weights = tiledlayout(fig_weights, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_weight_map = nexttile(tiled_weights);
    scatter(ax_weight_map, phi1, phi2, 18, w, 'filled');
    formatPeriodicPhasePlane(ax_weight_map);
    title(ax_weight_map, sprintf('Normalized weights on periodic phase plane (r=%d, beta=%.3g)', r, beta), 'Interpreter', 'none');
    cb = colorbar(ax_weight_map);
    ylabel(cb, 'normalized weight');

    ax_weight_hist = nexttile(tiled_weights);
    histogram(ax_weight_hist, w, 40, 'FaceColor', [0.0, 0.4470, 0.7410], 'EdgeColor', 'none');
    xlabel(ax_weight_hist, 'normalized weight');
    ylabel(ax_weight_hist, 'count');
    title(ax_weight_hist, sprintf('Weight distribution after clipping [%g, %g]%%', clip_percentiles(1), clip_percentiles(2)), 'Interpreter', 'none');
    grid(ax_weight_hist, 'on');
    box(ax_weight_hist, 'on');
    xline(ax_weight_hist, 1.0, '--k', 'mean(w)=1', 'LabelVerticalAlignment', 'bottom');

    % Visualize the original scattered data with the unweighted and weighted
    % fitted surfaces.
    fig_fit_compare = figure('Color', 'w');
    tiled_fit = tiledlayout(fig_fit_compare, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    z_limits = [min([z_original; surface_z_unweighted(:); surface_z_weighted(:)]), ...
                max([z_original; surface_z_unweighted(:); surface_z_weighted(:)])];

    ax_unweighted = nexttile(tiled_fit);
    plotScatterAndSurfaceOverlay(ax_unweighted, phi1, phi2, z_original, surface_phi1, surface_phi2, surface_z_unweighted, ...
        'Unweighted fit', sprintf('OLS: RMSE=%.4g, weighted RMSE=%.4g', rmse_unweighted, wrmse_unweighted), z_limits);

    ax_weighted = nexttile(tiled_fit);
    plotScatterAndSurfaceOverlay(ax_weighted, phi1, phi2, z_original, surface_phi1, surface_phi2, surface_z_weighted, ...
        'Weighted fit', sprintf('WLS: RMSE=%.4g, weighted RMSE=%.4g', rmse_weighted, wrmse_weighted), z_limits);

    sgtitle(tiled_fit, sprintf('Double Fourier comparison (M=%d, N=%d)', M, N), 'Interpreter', 'none');

    fig_contribution = figure('Color', 'w');
    tiled_contribution = tiledlayout(fig_contribution, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    max_display_items = 10;
    plot_group_colors = [
        0.0, 0.4470, 0.7410
        0.4660, 0.6740, 0.1880
        0.8500, 0.3250, 0.0980];

    ax_group = nexttile(tiled_contribution);
    plotContributionSummaryBars(ax_group, phi1_harmonic_summary_weighted, max_display_items, ...
        '\phi_1-only contribution by harmonic order', 'Harmonic order m', plot_group_colors(1, :));

    ax_group = nexttile(tiled_contribution);
    plotContributionSummaryBars(ax_group, phi2_harmonic_summary_weighted, max_display_items, ...
        '\phi_2-only contribution by harmonic order', 'Harmonic order n', plot_group_colors(2, :));

    ax_group = nexttile(tiled_contribution);
    plotContributionSummaryBars(ax_group, mixed_pair_summary_weighted, max_display_items, ...
        '\phi_1-\phi_2 mixed contribution by harmonic pair', 'Harmonic pair (m,n)', plot_group_colors(3, :));

    constant_row = strcmp(group_contribution_summary_weighted.group_name, 'constant');
    constant_ratio_percent = 100 * group_contribution_summary_weighted.contribution_ratio(constant_row);
    sgtitle(tiled_contribution, sprintf([ ...
        'Weighted mean-centered contribution ratio (M=%d, N=%d, beta=%.3g, r=%d) | ', ...
        'constant term = %.2f%%'], M, N, beta, r, constant_ratio_percent), 'Interpreter', 'tex');

    fprintf(['[INFO] Weighted double Fourier fit completed: Nsamples=%d, Nbasis=%d, ', ...
        'rank(A)=%d, cond(A)=%.3e, rank(Aw)=%d, cond(Aw)=%.3e, ', ...
        'RMSE(OLS)=%.6g, RMSE(WLS)=%.6g, WRMSE(OLS)=%.6g, WRMSE(WLS)=%.6g\n'], ...
        n_samples, size(A, 2), rankA, condA, rankAw, condAw, ...
        rmse_unweighted, rmse_weighted, wrmse_unweighted, wrmse_weighted);

    result = struct();
    result.phi1 = phi1;
    result.phi2 = phi2;
    result.z = z;
    result.z_original = z_original;
    result.z_centered = z;
    result.z_mean = z_mean;
    result.valid_mask = valid_mask;

    result.M = M;
    result.N = N;
    result.beta = beta;
    result.r = r;
    result.clip_percentiles = clip_percentiles;

    result.A = A;
    result.basis_names = basis_names;
    result.basis_groups = basis_groups;
    result.basis_m = basis_m;
    result.basis_n = basis_n;

    result.weight_info = weight_info;
    result.dk = weight_info.dk;
    result.weights_raw = weight_info.weights_raw;
    result.weights_clipped = weight_info.weights_clipped;
    result.weights = weight_info.weights;

    result.coeff_unweighted = coeff_unweighted;
    result.coeff_weighted = coeff_weighted;
    result.z_hat_unweighted = z_hat_unweighted;
    result.z_hat_weighted = z_hat_weighted;
    result.z_hat_unweighted_original_scale = z_hat_unweighted_original_scale;
    result.z_hat_weighted_original_scale = z_hat_weighted_original_scale;
    result.residual_unweighted = residual_unweighted;
    result.residual_weighted = residual_weighted;
    result.rmse_unweighted = rmse_unweighted;
    result.rmse_weighted = rmse_weighted;
    result.weighted_rmse_unweighted = wrmse_unweighted;
    result.weighted_rmse_weighted = wrmse_weighted;
    result.term_rms_contribution = term_rms_contribution_weighted;
    result.term_energy_contribution = term_energy_contribution_weighted;
    result.term_contribution_ratio = term_contribution_ratio_weighted;
    result.contribution_table = contribution_table_weighted;
    result.group_contribution_summary = group_contribution_summary_weighted;
    result.phi1_harmonic_summary = phi1_harmonic_summary_weighted;
    result.phi2_harmonic_summary = phi2_harmonic_summary_weighted;
    result.mixed_pair_summary = mixed_pair_summary_weighted;

    result.rankA = rankA;
    result.condA = condA;
    result.rankAw = rankAw;
    result.condAw = condAw;

    result.surface_phi1 = surface_phi1;
    result.surface_phi2 = surface_phi2;
    result.surface_z_unweighted = surface_z_unweighted;
    result.surface_z_weighted = surface_z_weighted;

    result.fig_weights = fig_weights;
    result.fig_fit_compare = fig_fit_compare;
    result.fig_contribution = fig_contribution;
end

function [A, basis_names, basis_groups, basis_m, basis_n] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N)
% buildDoubleFourierDesignMatrix
% Construct the real trigonometric basis matrix for a double Fourier model.
%
% Basis order:
%   1
%   cos(m*phi1), sin(m*phi1)                            m = 1..M
%   cos(n*phi2), sin(n*phi2)                            n = 1..N
%   cos(m*phi1)cos(n*phi2), cos(m*phi1)sin(n*phi2)
%   sin(m*phi1)cos(n*phi2), sin(m*phi1)sin(n*phi2)      m = 1..M, n = 1..N

    phi1 = phi1(:);
    phi2 = phi2(:);
    n_samples = numel(phi1);

    n_basis = 1 + 2 * M + 2 * N + 4 * M * N;
    A = zeros(n_samples, n_basis);
    basis_names = cell(n_basis, 1);
    basis_groups = cell(n_basis, 1);
    basis_m = zeros(n_basis, 1);
    basis_n = zeros(n_basis, 1);

    col = 1;

    A(:, col) = 1;
    basis_names{col} = '1';
    basis_groups{col} = 'constant';
    col = col + 1;

    for m = 1:M
        A(:, col) = cos(m * phi1);
        basis_names{col} = sprintf('cos(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        col = col + 1;

        A(:, col) = sin(m * phi1);
        basis_names{col} = sprintf('sin(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        col = col + 1;
    end

    for n = 1:N
        A(:, col) = cos(n * phi2);
        basis_names{col} = sprintf('cos(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_n(col) = n;
        col = col + 1;

        A(:, col) = sin(n * phi2);
        basis_names{col} = sprintf('sin(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_n(col) = n;
        col = col + 1;
    end

    for m = 1:M
        c1 = cos(m * phi1);
        s1 = sin(m * phi1);
        for n = 1:N
            c2 = cos(n * phi2);
            s2 = sin(n * phi2);

            A(:, col) = c1 .* c2;
            basis_names{col} = sprintf('cos(%d*phi1)cos(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            col = col + 1;

            A(:, col) = c1 .* s2;
            basis_names{col} = sprintf('cos(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            col = col + 1;

            A(:, col) = s1 .* c2;
            basis_names{col} = sprintf('sin(%d*phi1)cos(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            col = col + 1;

            A(:, col) = s1 .* s2;
            basis_names{col} = sprintf('sin(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            col = col + 1;
        end
    end
end

function weight_info = computePeriodicKnnWeights(phi1, phi2, r, beta, clip_percentiles)
% computePeriodicKnnWeights
% Compute periodic kNN-distance-based weights.
%
% Steps:
%   1. Compute the r-th nearest-neighbor distance d_k under the periodic
%      2D metric.
%   2. Form raw weights w_k = d_k^(2*beta).
%   3. Clip raw weights by percentile.
%   4. Normalize so mean(w)=1.

    dk = computePeriodicKthNeighborDistance(phi1, phi2, r);

    % Protect against duplicated points when beta < 0.
    dk_safe = max(dk, sqrt(eps(class(dk))));
    weights_raw = dk_safe .^ (2 * beta);

    clip_bounds = computePercentileBounds(weights_raw, clip_percentiles);
    low_bound = clip_bounds(1);
    high_bound = clip_bounds(2);

    weights_clipped = min(max(weights_raw, low_bound), high_bound);
    mean_weight = mean(weights_clipped);
    if ~isfinite(mean_weight) || mean_weight <= 0
        error('Normalized weight mean is invalid. Check beta and the input data.');
    end
    weights = weights_clipped / mean_weight;

    weight_info = struct();
    weight_info.r = r;
    weight_info.beta = beta;
    weight_info.dk = dk;
    weight_info.dk_safe = dk_safe;
    weight_info.weights_raw = weights_raw;
    weight_info.weights_clipped = weights_clipped;
    weight_info.weights = weights;
    weight_info.clip_percentiles = clip_percentiles;
    weight_info.clip_bounds = [low_bound, high_bound];
end

function dk = computePeriodicKthNeighborDistance(phi1, phi2, r)
% computePeriodicKthNeighborDistance
% Compute the r-th nearest-neighbor distance under the periodic metric
%   dist = sqrt(dphi1.^2 + dphi2.^2)
% where
%   dphi = atan2(sin(a-b), cos(a-b)).
%
% The implementation uses block processing to avoid forming the full
% N-by-N distance matrix at once.

    n_samples = numel(phi1);
    dk = zeros(n_samples, 1);

    target_block_entries = 2e6;
    block_size = max(1, min(n_samples, floor(target_block_entries / max(n_samples, 1))));

    phi1_row = reshape(phi1, 1, []);
    phi2_row = reshape(phi2, 1, []);

    for start_idx = 1:block_size:n_samples
        stop_idx = min(start_idx + block_size - 1, n_samples);
        row_idx = start_idx:stop_idx;

        dphi1 = periodicAngleDifference(phi1(row_idx), phi1_row);
        dphi2 = periodicAngleDifference(phi2(row_idx), phi2_row);
        dist_block = sqrt(dphi1 .^ 2 + dphi2 .^ 2);

        local_rows = 1:numel(row_idx);
        self_linear_idx = sub2ind(size(dist_block), local_rows, row_idx);
        dist_block(self_linear_idx) = inf;

        nearest_distances = mink(dist_block, r, 2);
        dk(row_idx) = nearest_distances(:, end);
    end
end

function dphi = periodicAngleDifference(a, b)
% periodicAngleDifference
% Return the shortest signed periodic angle difference using
%   atan2(sin(a-b), cos(a-b)).
%
% The inputs can be scalars, vectors, or arrays with compatible sizes.

    dphi = atan2(sin(a - b), cos(a - b));
end

function percentile_values = computePercentileBounds(x, percentiles)
% computePercentileBounds
% Base-MATLAB percentile evaluation by linear interpolation on sorted data.
% This avoids a toolbox dependency for percentile clipping.

    x = sort(x(:));
    n = numel(x);
    percentile_values = zeros(size(percentiles));

    if n == 0
        error('Cannot compute percentiles of an empty array.');
    end
    if n == 1
        percentile_values(:) = x;
        return;
    end

    for i = 1:numel(percentiles)
        p = percentiles(i);
        if p <= 0
            percentile_values(i) = x(1);
        elseif p >= 100
            percentile_values(i) = x(end);
        else
            position = 1 + (n - 1) * (p / 100);
            idx_low = floor(position);
            idx_high = ceil(position);
            if idx_low == idx_high
                percentile_values(i) = x(idx_low);
            else
                alpha = position - idx_low;
                percentile_values(i) = (1 - alpha) * x(idx_low) + alpha * x(idx_high);
            end
        end
    end
end

function wrmse = computeWeightedRMSE(residual, weights)
% computeWeightedRMSE
% Weighted root-mean-square error.

    wrmse = sqrt(sum(weights .* (residual .^ 2)) / sum(weights));
end

function summary_table = summarizeSingleAxisContributions(basis_groups, order_values, term_energy_contribution, total_term_energy, target_group)
    group_mask = strcmp(basis_groups, target_group) & order_values > 0;
    unique_orders = unique(order_values(group_mask));

    if isempty(unique_orders)
        summary_table = table([], cell(0, 1), [], [], 'VariableNames', ...
            {'order_value', 'label', 'energy_contribution', 'contribution_ratio'});
        return;
    end

    n_orders = numel(unique_orders);
    labels = cell(n_orders, 1);
    energy_contribution = zeros(n_orders, 1);
    contribution_ratio = zeros(n_orders, 1);
    for i = 1:n_orders
        order_value = unique_orders(i);
        order_mask = group_mask & order_values == order_value;
        energy_contribution(i) = sum(term_energy_contribution(order_mask));
        if total_term_energy > 0
            contribution_ratio(i) = energy_contribution(i) / total_term_energy;
        end
        labels{i} = sprintf('%d', order_value);
    end

    summary_table = table(unique_orders, labels, energy_contribution, contribution_ratio, ...
        'VariableNames', {'order_value', 'label', 'energy_contribution', 'contribution_ratio'});
    summary_table = sortrows(summary_table, 'contribution_ratio', 'descend');
end

function summary_table = summarizeMixedPairContributions(basis_groups, basis_m, basis_n, term_energy_contribution, total_term_energy)
    group_mask = strcmp(basis_groups, 'mixed') & basis_m > 0 & basis_n > 0;
    pair_values = unique([basis_m(group_mask), basis_n(group_mask)], 'rows');

    if isempty(pair_values)
        summary_table = table([], [], cell(0, 1), [], [], 'VariableNames', ...
            {'m_order', 'n_order', 'label', 'energy_contribution', 'contribution_ratio'});
        return;
    end

    n_pairs = size(pair_values, 1);
    labels = cell(n_pairs, 1);
    energy_contribution = zeros(n_pairs, 1);
    contribution_ratio = zeros(n_pairs, 1);
    for i = 1:n_pairs
        m_order = pair_values(i, 1);
        n_order = pair_values(i, 2);
        pair_mask = group_mask & basis_m == m_order & basis_n == n_order;
        energy_contribution(i) = sum(term_energy_contribution(pair_mask));
        if total_term_energy > 0
            contribution_ratio(i) = energy_contribution(i) / total_term_energy;
        end
        labels{i} = sprintf('(%d,%d)', m_order, n_order);
    end

    summary_table = table(pair_values(:, 1), pair_values(:, 2), labels, ...
        energy_contribution, contribution_ratio, 'VariableNames', ...
        {'m_order', 'n_order', 'label', 'energy_contribution', 'contribution_ratio'});
    summary_table = sortrows(summary_table, 'contribution_ratio', 'descend');
end

function plotContributionSummaryBars(ax, summary_table, max_display_items, title_prefix, xlabel_text, bar_color)
    if isempty(summary_table) || height(summary_table) == 0
        axis(ax, 'off');
        title(ax, title_prefix, 'Interpreter', 'tex');
        text(ax, 0.5, 0.5, 'No terms in this group', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        return;
    end

    n_display = min(max_display_items, height(summary_table));
    summary_top = summary_table(1:n_display, :);
    x_positions = 1:n_display;

    bar(ax, x_positions, 100 * summary_top.contribution_ratio, ...
        'FaceColor', bar_color, 'EdgeColor', 'none');
    xlabel(ax, xlabel_text);
    ylabel(ax, 'Contribution (%)');
    title(ax, sprintf('%s (top %d, total %.2f%%)', ...
        title_prefix, n_display, 100 * sum(summary_table.contribution_ratio)), 'Interpreter', 'tex');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'XTick', x_positions, 'XTickLabel', summary_top.label, 'TickLabelInterpreter', 'none');
    xtickangle(ax, 25);
end

function plotScatterAndSurfaceOverlay(ax, phi1, phi2, z_scatter, surface_phi1, surface_phi2, surface_z, panel_title, subtitle_text, z_limits)
% plotScatterAndSurfaceOverlay
% Overlay the original 3D scattered data with a fitted Fourier surface.

    hold(ax, 'on');
    surf(ax, surface_phi1, surface_phi2, surface_z, surface_z, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.78, 'FaceLighting', 'gouraud');
    scatter3(ax, phi1, phi2, z_scatter, 12, ...
        'filled', 'MarkerFaceColor', [0.08, 0.08, 0.08], 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

    xlabel(ax, '\phi_1 (rad)', 'Interpreter', 'tex');
    ylabel(ax, '\phi_2 (rad)', 'Interpreter', 'tex');
    zlabel(ax, 'z', 'Interpreter', 'none');
    title(ax, {panel_title, subtitle_text}, 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');
    view(ax, [40, 28]);

    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    if nargin >= 9 && ~isempty(z_limits) && all(isfinite(z_limits)) && z_limits(1) < z_limits(2)
        zlim(ax, z_limits);
    end
    colormap(ax, parula);
    colorbar(ax);
end

function formatPeriodicPhasePlane(ax)
% formatPeriodicPhasePlane
% Standard formatting for the 2*pi-periodic phi1-phi2 plane.

    xlabel(ax, '\phi_1 (rad)', 'Interpreter', 'tex');
    ylabel(ax, '\phi_2 (rad)', 'Interpreter', 'tex');
    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    grid(ax, 'on');
    box(ax, 'on');
end
