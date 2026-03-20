function result = fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name, heatmap_mode, gamma_ratio, gamma_settings)
% fitDoubleFourierScatter Fit a real double Fourier series to scattered data.
%
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N)
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name)
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name, heatmap_mode)
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name, heatmap_mode, gamma_ratio)
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name, heatmap_mode, gamma_ratio, gamma_settings)
%
% Inputs:
%   phi1, phi2 : angle variables (radians). They are normalized to [0, 2*pi).
%   z          : observed values Q(phi1, phi2)
%   M, N       : truncation orders for phi1 and phi2 directions
%   target_name: optional label used in figure titles and z-axis labels
%   heatmap_mode:
%       'full'       -> show (0,0), (m,0), (0,n), (m,n)
%       'mixed-only' -> show only the mixed-term block on m=1:M, n=1:N
%   gamma_ratio: two-element integer vector [m_phi2, n_phi1] defining
%       psi = m_phi2 * phi2 - n_phi1 * phi1 for the reconstructed
%       one-variable resonant interaction Gamma(psi). The ratio is reduced
%       to lowest terms internally, so [1 1] and [2 2] are treated as the
%       same resonance. Default: [2 1]
%   gamma_settings: optional struct with fields
%       enabled   : whether Gamma(psi) is reconstructed/plotted (default true)
%       component : 'full', 'symmetric', or 'antisymmetric' (default 'full')
%       overlay_full : when true, overlay the full reconstructed Gamma(psi)
%           together with the selected component on the same axes
%           (default false)
%       show_surface_overlay : whether to create the 3D scatter + fitted
%           surface overlay figure (default true)
%       auto_save_figure : whether to call saveFigure after tuneFigure for
%           the 3D scatter + fitted surface overlay figure (default false)
%
% The fitted model uses the real trigonometric basis
%   1
%   cos(m*phi1), sin(m*phi1)                    for m = 1..M
%   cos(n*phi2), sin(n*phi2)                    for n = 1..N
%   cos(m*phi1)cos(n*phi2)                      for m = 1..M, n = 1..N
%   cos(m*phi1)sin(n*phi2)                      for m = 1..M, n = 1..N
%   sin(m*phi1)cos(n*phi2)                      for m = 1..M, n = 1..N
%   sin(m*phi1)sin(n*phi2)                      for m = 1..M, n = 1..N
%
% The coefficient vector is estimated by linear least squares:
%   coeff = A \ z
%
% Output result contains:
%   result.coeff
%   result.z_hat
%   result.z_hat_original_scale
%   result.z_mean
%   result.z_original
%   result.z_centered
%   result.rmse
%   result.A
%   result.basis_names
%   result.basis_groups
%   result.phi1_harmonic_summary
%   result.phi2_harmonic_summary
%   result.mixed_pair_summary
%   result.mixed_pair_summary_resonant
%   result.term_rms_contribution
%   result.term_energy_contribution
%   result.term_contribution_ratio
%   result.group_contribution_summary
%   result.contribution_map
%   result.contribution_map_percent
%   result.contribution_map_mixed
%   result.contribution_map_mixed_percent
%   result.surface_phi1
%   result.surface_phi2
%   result.surface_z
%   result.fig_contribution_bars
%   result.target_name
%   result.heatmap_mode
%   result.gamma_ratio
%   result.gamma_ratio_reduced
%   result.gamma_settings
%   result.gamma_resonance
%   result.fig_gamma_resonance
% and additional diagnostic fields.
%
% Example:
%   out = plot_phase_pair_a2_3d_all_files(fullfile('EstimateF','Spring3','250'), [2 3], 3, 15, 5);
%   result = fitDoubleFourierScatter(out.point_cloud.phi1, out.point_cloud.phi2, out.point_cloud.a2, 3, 3);
%
%   % If you already have scattered data directly:
%   % result = fitDoubleFourierScatter(phi1, phi2, z, 4, 4);

    if nargin < 5
        error('Usage: fitDoubleFourierScatter(phi1, phi2, z, M, N, target_name, heatmap_mode, gamma_ratio, gamma_settings)');
    end
    if nargin < 6 || isempty(target_name)
        target_name = 'z';
    end
    if nargin < 7 || isempty(heatmap_mode)
        heatmap_mode = 'full';
    end
    if nargin < 8 || isempty(gamma_ratio)
        gamma_ratio = [2 1];
    end
    if nargin < 9 || isempty(gamma_settings)
        gamma_settings = struct();
    end

    if isstring(target_name)
        target_name = char(target_name);
    end
    if isstring(heatmap_mode)
        heatmap_mode = char(heatmap_mode);
    end
    heatmap_mode = lower(strtrim(heatmap_mode));
    gamma_ratio = double(gamma_ratio(:).');
    gamma_ratio_reduced = reducePositiveIntegerRatio(gamma_ratio);
    gamma_settings = normalizeGammaSettings(gamma_settings);

    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'M');
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'N');
    validateattributes(gamma_ratio, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive', 'finite'}, mfilename, 'gamma_ratio');

    phi1 = phi1(:);
    phi2 = phi2(:);
    z = z(:);

    if ~(numel(phi1) == numel(phi2) && numel(phi1) == numel(z))
        error('phi1, phi2, and z must have the same number of elements.');
    end

    % Normalize angle variables to the periodic domain [0, 2*pi)
    phi1 = mod(phi1, 2*pi);
    phi2 = mod(phi2, 2*pi);

    % Remove NaN / Inf samples before fitting
    valid_mask = isfinite(phi1) & isfinite(phi2) & isfinite(z);
    phi1 = phi1(valid_mask);
    phi2 = phi2(valid_mask);
    z = z(valid_mask);

    if isempty(z)
        error('No valid samples remain after removing NaN/Inf values.');
    end

    % Mean-center the observed values before fitting so the constant term
    % does not dominate the contribution analysis.
    z_original = z;
    z_mean = mean(z_original);
    z = z_original - z_mean;

    [A, basis_names, basis_groups, basis_m, basis_n, basis_types] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N);

    % Solve the linear least-squares problem A * coeff ≈ z
    coeff = A \ z;
    z_hat = A * coeff;
    residual = z - z_hat;
    rmse = sqrt(mean(residual .^ 2));
    z_hat_original_scale = z_hat + z_mean;

    rankA = rank(A);
    n_cols = size(A, 2);
    condA = cond(A);

    if rankA < n_cols
        warning('fitDoubleFourierScatter:RankDeficient', ...
            'Design matrix may be rank deficient: rank(A) = %d < %d.', rankA, n_cols);
    end
    if isfinite(condA) && condA > 1e10
        warning('fitDoubleFourierScatter:IllConditioned', ...
            'Design matrix may be ill-conditioned: cond(A) = %.3e.', condA);
    end

    % Per-term contribution on the observed samples.
    % Each fitted basis term is A(:,k) * coeff(k). We summarize its
    % contribution by the sample RMS / mean-square amplitude.
    term_matrix = bsxfun(@times, A, coeff.');
    term_energy_contribution = mean(term_matrix .^ 2, 1).';
    term_rms_contribution = sqrt(term_energy_contribution);

    total_term_energy = sum(term_energy_contribution);
    if total_term_energy > 0
        term_contribution_ratio = term_energy_contribution / total_term_energy;
    else
        term_contribution_ratio = zeros(size(term_energy_contribution));
    end

    contribution_table = table((1:n_cols).', basis_names, basis_groups, basis_m, basis_n, coeff, ...
        term_rms_contribution, term_energy_contribution, term_contribution_ratio, ...
        'VariableNames', {'term_index', 'basis_name', 'basis_group', 'phi1_order', 'phi2_order', 'coefficient', ...
        'rms_contribution', 'energy_contribution', 'contribution_ratio'});

    group_names = {'constant'; 'phi1_only'; 'phi2_only'; 'mixed'};
    group_display_names = {'constant'; 'phi1 only'; 'phi2 only'; 'phi1-phi2 mixed'};
    group_term_count = zeros(numel(group_names), 1);
    group_energy_contribution = zeros(numel(group_names), 1);
    group_contribution_ratio = zeros(numel(group_names), 1);
    for g = 1:numel(group_names)
        group_mask = strcmp(basis_groups, group_names{g});
        group_term_count(g) = nnz(group_mask);
        group_energy_contribution(g) = sum(term_energy_contribution(group_mask));
        if total_term_energy > 0
            group_contribution_ratio(g) = group_energy_contribution(g) / total_term_energy;
        end
    end

    group_contribution_summary = table(group_names, group_display_names, group_term_count, ...
        group_energy_contribution, group_contribution_ratio, ...
        'VariableNames', {'group_name', 'display_name', 'n_terms', ...
        'energy_contribution', 'contribution_ratio'});

    phi1_harmonic_summary = summarizeSingleAxisContributions( ...
        basis_groups, basis_m, term_energy_contribution, total_term_energy, 'phi1_only');
    phi2_harmonic_summary = summarizeSingleAxisContributions( ...
        basis_groups, basis_n, term_energy_contribution, total_term_energy, 'phi2_only');
    mixed_pair_summary = summarizeMixedPairContributions( ...
        basis_groups, basis_m, basis_n, term_energy_contribution, total_term_energy);
    mixed_pair_summary_resonant = filterMixedPairSummaryByResonance( ...
        mixed_pair_summary, gamma_ratio_reduced(2), gamma_ratio_reduced(1));
    contribution_map = buildContributionHeatmap( ...
        M, N, group_contribution_summary, phi1_harmonic_summary, phi2_harmonic_summary, mixed_pair_summary);
    contribution_map_percent = 100 * contribution_map;
    contribution_map_mixed = contribution_map(2:end, 2:end);
    contribution_map_mixed_percent = 100 * contribution_map_mixed;

    surface_grid_size = 81;
    surface_phi1_values = linspace(0, 2*pi, surface_grid_size);
    surface_phi2_values = linspace(0, 2*pi, surface_grid_size);
    [surface_phi1, surface_phi2] = meshgrid(surface_phi1_values, surface_phi2_values);
    A_surface = buildDoubleFourierDesignMatrix(surface_phi1(:), surface_phi2(:), M, N);
    surface_z_centered = reshape(A_surface * coeff, size(surface_phi1));
    surface_z = surface_z_centered + z_mean;

    if gamma_settings.show_surface_overlay
        fig_original = figure('Color', 'w');
        ax_original = axes('Parent', fig_original);
        plotScatterAndSurfaceOverlay(ax_original, phi1, phi2, z_original, surface_phi1, surface_phi2, surface_z, ...
            target_name, [-37.5, 30]);
        colorbar(ax_original);
        figure(fig_original);
        tuneFigure;
        if gamma_settings.auto_save_figure
            saveFigure;
        end
    else
        fig_original = [];
    end

    fig_fit = [];
    fig_residual = [];

    fig_contribution = [];
    fig_contribution_bars = [];

    if gamma_settings.enabled
        gamma_resonance = reconstructResonantGamma( ...
            coeff, basis_groups, basis_types, basis_m, basis_n, gamma_ratio_reduced(2), gamma_ratio_reduced(1));
        [gamma_values_full, gamma_values_symmetric, gamma_values_antisymmetric] = ...
            evaluateGammaComponents(gamma_resonance.psi_grid, gamma_resonance.harmonic_index, ...
                gamma_resonance.gamma_cos, gamma_resonance.gamma_sin);
        psi_grid_centered = linspace(-pi, pi, numel(gamma_resonance.psi_grid)).';
        [gamma_values_full_centered, gamma_values_symmetric_centered, gamma_values_antisymmetric_centered] = ...
            evaluateGammaComponents(psi_grid_centered, gamma_resonance.harmonic_index, ...
                gamma_resonance.gamma_cos, gamma_resonance.gamma_sin);

        gamma_resonance.gamma_values_full = gamma_values_full;
        gamma_resonance.gamma_values_symmetric = gamma_values_symmetric;
        gamma_resonance.gamma_values_antisymmetric = gamma_values_antisymmetric;
        gamma_resonance.psi_grid_centered = psi_grid_centered;
        gamma_resonance.gamma_values_full_centered = gamma_values_full_centered;
        gamma_resonance.gamma_values_symmetric_centered = gamma_values_symmetric_centered;
        gamma_resonance.gamma_values_antisymmetric_centered = gamma_values_antisymmetric_centered;
        gamma_resonance.component = gamma_settings.component;
        gamma_resonance.overlay_full = gamma_settings.overlay_full;

        switch gamma_settings.component
            case 'full'
                gamma_resonance.gamma_values = gamma_values_full;
                gamma_resonance.gamma_values_centered = gamma_values_full_centered;
                gamma_component_label = 'full';
            case 'symmetric'
                gamma_resonance.gamma_values = gamma_values_symmetric;
                gamma_resonance.gamma_values_centered = gamma_values_symmetric_centered;
                gamma_component_label = 'symmetric part';
            case 'antisymmetric'
                gamma_resonance.gamma_values = gamma_values_antisymmetric;
                gamma_resonance.gamma_values_centered = gamma_values_antisymmetric_centered;
                gamma_component_label = 'antisymmetric part';
            otherwise
                error('Unsupported gamma component: %s', gamma_settings.component);
        end

        fig_gamma_resonance = figure('Color', 'w');
        ax_gamma = axes('Parent', fig_gamma_resonance);
        if isempty(gamma_resonance.harmonic_index)
            axis(ax_gamma, 'off');
            text(ax_gamma, 0.5, 0.5, sprintf('No resonant mixed terms found for \\psi = %g\\phi_2 - %g\\phi_1', ...
                gamma_ratio_reduced(1), gamma_ratio_reduced(2)), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Interpreter', 'tex');
        else
            hold(ax_gamma, 'on');
            overlay_full_curve = gamma_settings.overlay_full && ~strcmp(gamma_settings.component, 'full');
            if overlay_full_curve
                plot(ax_gamma, gamma_resonance.psi_grid_centered, gamma_resonance.gamma_values_full_centered, '--', ...
                    'LineWidth', 1.4, 'Color', [0.45, 0.45, 0.45], 'DisplayName', 'full reconstructed Gamma(psi)');
            end
            plot(ax_gamma, gamma_resonance.psi_grid_centered, gamma_resonance.gamma_values_centered, 'LineWidth', 1.8, ...
                'Color', [0.0, 0.4470, 0.7410], 'DisplayName', sprintf('%s of Gamma(psi)', gamma_component_label));
            xlabel(ax_gamma, '$$\psi$$', 'Interpreter', 'latex');
            ylabel(ax_gamma, '$$\Gamma(\psi)$$', 'Interpreter', 'latex');
            if overlay_full_curve
                legend(ax_gamma, 'Location', 'best');
            end
            grid(ax_gamma, 'on');
            box(ax_gamma, 'on');
            xlim(ax_gamma, [-pi, pi]);
            xticks(ax_gamma, [-pi, -pi/2, 0, pi/2, pi]);
            xticklabels(ax_gamma, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
            ax_gamma.XLabel.Interpreter = 'latex';
            ax_gamma.YLabel.Interpreter = 'latex';
            ax_gamma.TickLabelInterpreter = 'latex';
            text(ax_gamma, 0.02, 0.98, sprintf('Harmonics used: %s', mat2str(gamma_resonance.harmonic_index.')), ...
                'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'Interpreter', 'none', 'BackgroundColor', 'w', 'Margin', 3);
        end
        figure(fig_gamma_resonance);
        tuneFigure;
    else
        gamma_resonance = struct( ...
            'enabled', false, ...
            'component', gamma_settings.component, ...
            'overlay_full', gamma_settings.overlay_full, ...
            'phi1_base', gamma_ratio_reduced(2), ...
            'phi2_base', gamma_ratio_reduced(1), ...
            'psi_label', sprintf('%g*phi2 - %g*phi1', gamma_ratio_reduced(1), gamma_ratio_reduced(2)), ...
            'harmonic_index', [], ...
            'gamma_cos', [], ...
            'gamma_sin', [], ...
            'psi_grid', [], ...
            'gamma_values', [], ...
            'psi_grid_centered', [], ...
            'gamma_values_centered', [], ...
            'gamma_values_full', [], ...
            'gamma_values_full_centered', [], ...
            'gamma_values_symmetric', [], ...
            'gamma_values_symmetric_centered', [], ...
            'gamma_values_antisymmetric', [], ...
            'gamma_values_antisymmetric_centered', []);
        fig_gamma_resonance = [];
    end

    fprintf('[INFO] Double Fourier fit completed after mean-centering: mean(z)=%.6g, Nsamples=%d, Nbasis=%d, rank(A)=%d, cond(A)=%.3e, RMSE=%.6g\n', ...
        z_mean, numel(z), n_cols, rankA, condA, rmse);

    result = struct();
    result.coeff = coeff;
    result.z_hat = z_hat;
    result.z_hat_original_scale = z_hat_original_scale;
    result.z_mean = z_mean;
    result.z_original = z_original;
    result.z_centered = z;
    result.rmse = rmse;
    result.A = A;
    result.basis_names = basis_names;
    result.basis_groups = basis_groups;
    result.basis_types = basis_types;
    result.phi1_harmonic_summary = phi1_harmonic_summary;
    result.phi2_harmonic_summary = phi2_harmonic_summary;
    result.mixed_pair_summary = mixed_pair_summary;
    result.mixed_pair_summary_resonant = mixed_pair_summary_resonant;
    result.term_rms_contribution = term_rms_contribution;
    result.term_energy_contribution = term_energy_contribution;
    result.term_contribution_ratio = term_contribution_ratio;
    result.contribution_table = contribution_table;
    result.group_contribution_summary = group_contribution_summary;
    result.contribution_map = contribution_map;
    result.contribution_map_percent = contribution_map_percent;
    result.contribution_map_mixed = contribution_map_mixed;
    result.contribution_map_mixed_percent = contribution_map_mixed_percent;
    result.surface_phi1 = surface_phi1;
    result.surface_phi2 = surface_phi2;
    result.surface_z = surface_z;
    result.residual = residual;
    result.phi1 = phi1;
    result.phi2 = phi2;
    result.z = z;
    result.M = M;
    result.N = N;
    result.target_name = target_name;
    result.heatmap_mode = heatmap_mode;
    result.gamma_ratio = gamma_ratio;
    result.gamma_ratio_reduced = gamma_ratio_reduced;
    result.gamma_settings = gamma_settings;
    result.gamma_resonance = gamma_resonance;
    result.rankA = rankA;
    result.condA = condA;
    result.valid_mask = valid_mask;
    result.fig_original = fig_original;
    result.fig_fit = fig_fit;
    result.fig_residual = fig_residual;
    result.fig_contribution = fig_contribution;
    result.fig_contribution_bars = fig_contribution_bars;
    result.fig_gamma_resonance = fig_gamma_resonance;
end

function [A, basis_names, basis_groups, basis_m, basis_n, basis_types] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N)
% buildDoubleFourierDesignMatrix Construct the real trigonometric basis matrix.
%
% Each column of A corresponds to one named basis function in basis_names.

    n_samples = numel(phi1);
    n_basis = 1 + 2 * M + 2 * N + 4 * M * N;
    A = zeros(n_samples, n_basis);
    basis_names = cell(n_basis, 1);
    basis_groups = cell(n_basis, 1);
    basis_m = zeros(n_basis, 1);
    basis_n = zeros(n_basis, 1);
    basis_types = cell(n_basis, 1);

    col = 1;

    % Constant term: 1
    A(:, col) = 1;
    basis_names{col} = '1';
    basis_groups{col} = 'constant';
    basis_m(col) = 0;
    basis_n(col) = 0;
    basis_types{col} = 'constant';
    col = col + 1;

    % phi1-only terms: cos(m*phi1), sin(m*phi1)
    for m = 1:M
        A(:, col) = cos(m * phi1);
        basis_names{col} = sprintf('cos(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        basis_n(col) = 0;
        basis_types{col} = 'phi1_cos';
        col = col + 1;

        A(:, col) = sin(m * phi1);
        basis_names{col} = sprintf('sin(%d*phi1)', m);
        basis_groups{col} = 'phi1_only';
        basis_m(col) = m;
        basis_n(col) = 0;
        basis_types{col} = 'phi1_sin';
        col = col + 1;
    end

    % phi2-only terms: cos(n*phi2), sin(n*phi2)
    for n = 1:N
        A(:, col) = cos(n * phi2);
        basis_names{col} = sprintf('cos(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_m(col) = 0;
        basis_n(col) = n;
        basis_types{col} = 'phi2_cos';
        col = col + 1;

        A(:, col) = sin(n * phi2);
        basis_names{col} = sprintf('sin(%d*phi2)', n);
        basis_groups{col} = 'phi2_only';
        basis_m(col) = 0;
        basis_n(col) = n;
        basis_types{col} = 'phi2_sin';
        col = col + 1;
    end

    % Mixed terms:
    %   cos(m*phi1)cos(n*phi2)
    %   cos(m*phi1)sin(n*phi2)
    %   sin(m*phi1)cos(n*phi2)
    %   sin(m*phi1)sin(n*phi2)
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
            basis_types{col} = 'mixed_cc';
            col = col + 1;

            A(:, col) = c1 .* s2;
            basis_names{col} = sprintf('cos(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_cs';
            col = col + 1;

            A(:, col) = s1 .* c2;
            basis_names{col} = sprintf('sin(%d*phi1)cos(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_sc';
            col = col + 1;

            A(:, col) = s1 .* s2;
            basis_names{col} = sprintf('sin(%d*phi1)sin(%d*phi2)', m, n);
            basis_groups{col} = 'mixed';
            basis_m(col) = m;
            basis_n(col) = n;
            basis_types{col} = 'mixed_ss';
            col = col + 1;
        end
    end
end

function gamma_resonance = reconstructResonantGamma(coeff, basis_groups, basis_types, basis_m, basis_n, phi1_base, phi2_base)
    max_harmonic = min(floor(max(basis_m) / max(phi1_base, 1)), floor(max(basis_n) / max(phi2_base, 1)));
    harmonic_index = [];
    gamma_cos = [];
    gamma_sin = [];

    for k = 1:max_harmonic
        m_order = k * phi1_base;
        n_order = k * phi2_base;
        if m_order <= 0 || n_order <= 0
            continue;
        end

        cc = getMixedCoefficient(coeff, basis_groups, basis_types, basis_m, basis_n, m_order, n_order, 'mixed_cc');
        cs = getMixedCoefficient(coeff, basis_groups, basis_types, basis_m, basis_n, m_order, n_order, 'mixed_cs');
        sc = getMixedCoefficient(coeff, basis_groups, basis_types, basis_m, basis_n, m_order, n_order, 'mixed_sc');
        ss = getMixedCoefficient(coeff, basis_groups, basis_types, basis_m, basis_n, m_order, n_order, 'mixed_ss');

        diff_cos = 0.5 * (cc + ss);
        diff_sin = 0.5 * (cs - sc);

        if abs(diff_cos) > 0 || abs(diff_sin) > 0
            harmonic_index(end + 1, 1) = k; %#ok<AGROW>
            gamma_cos(end + 1, 1) = diff_cos; %#ok<AGROW>
            gamma_sin(end + 1, 1) = diff_sin; %#ok<AGROW>
        end
    end

    psi_grid = linspace(0, 2*pi, 512).';
    gamma_values = zeros(size(psi_grid));
    for i = 1:numel(harmonic_index)
        k = harmonic_index(i);
        gamma_values = gamma_values + gamma_cos(i) * cos(k * psi_grid) + gamma_sin(i) * sin(k * psi_grid);
    end

    gamma_resonance = struct();
    gamma_resonance.phi1_base = phi1_base;
    gamma_resonance.phi2_base = phi2_base;
    gamma_resonance.psi_label = sprintf('%d*phi2 - %d*phi1', phi2_base, phi1_base);
    gamma_resonance.harmonic_index = harmonic_index;
    gamma_resonance.gamma_cos = gamma_cos;
    gamma_resonance.gamma_sin = gamma_sin;
    gamma_resonance.psi_grid = psi_grid;
    gamma_resonance.gamma_values = gamma_values;
end

function [gamma_values_full, gamma_values_symmetric, gamma_values_antisymmetric] = evaluateGammaComponents(psi_grid, harmonic_index, gamma_cos, gamma_sin)
    gamma_values_symmetric = zeros(size(psi_grid));
    gamma_values_antisymmetric = zeros(size(psi_grid));

    for i = 1:numel(harmonic_index)
        k = harmonic_index(i);
        gamma_values_symmetric = gamma_values_symmetric + gamma_cos(i) * cos(k * psi_grid);
        gamma_values_antisymmetric = gamma_values_antisymmetric + gamma_sin(i) * sin(k * psi_grid);
    end

    gamma_values_full = gamma_values_symmetric + gamma_values_antisymmetric;
end

function coeff_value = getMixedCoefficient(coeff, basis_groups, basis_types, basis_m, basis_n, m_order, n_order, basis_type)
    mask = strcmp(basis_groups, 'mixed') & strcmp(basis_types, basis_type) & ...
        basis_m == m_order & basis_n == n_order;
    if any(mask)
        coeff_value = coeff(find(mask, 1, 'first'));
    else
        coeff_value = 0;
    end
end

function ratio_reduced = reducePositiveIntegerRatio(ratio)
    a = round(ratio(1));
    b = round(ratio(2));

    while b ~= 0
        t = mod(a, b);
        a = b;
        b = t;
    end

    ratio_gcd = max(a, 1);
    ratio_reduced = ratio / ratio_gcd;
end

function gamma_settings = normalizeGammaSettings(gamma_settings)
    if ~isstruct(gamma_settings)
        error('gamma_settings must be a struct when provided.');
    end

    if ~isfield(gamma_settings, 'enabled') || isempty(gamma_settings.enabled)
        gamma_settings.enabled = true;
    end
    gamma_settings.enabled = logical(gamma_settings.enabled);

    if ~isfield(gamma_settings, 'component') || isempty(gamma_settings.component)
        gamma_settings.component = 'full';
    end
    if isstring(gamma_settings.component)
        gamma_settings.component = char(gamma_settings.component);
    end
    gamma_settings.component = lower(strtrim(gamma_settings.component));

    allowed_components = {'full', 'symmetric', 'antisymmetric'};
    if ~ismember(gamma_settings.component, allowed_components)
        error('Unsupported gamma_settings.component: %s', gamma_settings.component);
    end

    if ~isfield(gamma_settings, 'overlay_full') || isempty(gamma_settings.overlay_full)
        gamma_settings.overlay_full = false;
    end
    gamma_settings.overlay_full = logical(gamma_settings.overlay_full);

    if ~isfield(gamma_settings, 'show_surface_overlay') || isempty(gamma_settings.show_surface_overlay)
        gamma_settings.show_surface_overlay = true;
    end
    gamma_settings.show_surface_overlay = logical(gamma_settings.show_surface_overlay);

    if ~isfield(gamma_settings, 'auto_save_figure') || isempty(gamma_settings.auto_save_figure)
        gamma_settings.auto_save_figure = false;
    end
    gamma_settings.auto_save_figure = logical(gamma_settings.auto_save_figure);
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

function summary_table = filterMixedPairSummaryByResonance(mixed_pair_summary, phi1_base, phi2_base)
    if isempty(mixed_pair_summary) || height(mixed_pair_summary) == 0
        summary_table = mixed_pair_summary;
        return;
    end

    phi1_base = double(phi1_base);
    phi2_base = double(phi2_base);
    k_phi1 = mixed_pair_summary.m_order / phi1_base;
    k_phi2 = mixed_pair_summary.n_order / phi2_base;
    resonant_mask = abs(k_phi1 - round(k_phi1)) < 1e-12 & ...
        abs(k_phi2 - round(k_phi2)) < 1e-12 & ...
        round(k_phi1) == round(k_phi2) & round(k_phi1) >= 1;

    summary_table = mixed_pair_summary(resonant_mask, :);
    if ~isempty(summary_table)
        summary_table = sortrows(summary_table, 'contribution_ratio', 'descend');
    end
end

function contribution_map = buildContributionHeatmap(M, N, group_contribution_summary, phi1_harmonic_summary, phi2_harmonic_summary, mixed_pair_summary)
    contribution_map = zeros(M + 1, N + 1);

    constant_row = strcmp(group_contribution_summary.group_name, 'constant');
    if any(constant_row)
        contribution_map(1, 1) = group_contribution_summary.contribution_ratio(constant_row);
    end

    for i = 1:height(phi1_harmonic_summary)
        m_order = phi1_harmonic_summary.order_value(i);
        if m_order >= 1 && m_order <= M
            contribution_map(m_order + 1, 1) = phi1_harmonic_summary.contribution_ratio(i);
        end
    end

    for i = 1:height(phi2_harmonic_summary)
        n_order = phi2_harmonic_summary.order_value(i);
        if n_order >= 1 && n_order <= N
            contribution_map(1, n_order + 1) = phi2_harmonic_summary.contribution_ratio(i);
        end
    end

    for i = 1:height(mixed_pair_summary)
        m_order = mixed_pair_summary.m_order(i);
        n_order = mixed_pair_summary.n_order(i);
        if m_order >= 1 && m_order <= M && n_order >= 1 && n_order <= N
            contribution_map(m_order + 1, n_order + 1) = mixed_pair_summary.contribution_ratio(i);
        end
    end
end

function plotScatterAndSurfaceOverlay(ax, phi1, phi2, z_scatter, surface_phi1, surface_phi2, surface_z, zlabel_text, view_angles)
    hold(ax, 'on');
    surf(ax, surface_phi1, surface_phi2, surface_z, surface_z, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.72, 'FaceLighting', 'gouraud');
    scatter3(ax, phi1, phi2, z_scatter, 10, ...
        'filled', 'MarkerFaceColor', [0.1, 0.1, 0.1], 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22);
    formatPhaseScatterAxes(ax, zlabel_text);
    view(ax, view_angles);
end

function formatPhaseScatterAxes(ax, zlabel_text) %#ok<INUSD>
    xlabel(ax, '$$\phi_1$$', 'Interpreter', 'latex');
    ylabel(ax, '$$\phi_2$$', 'Interpreter', 'latex');
    zlabel(ax, zlabel_text, 'Interpreter', 'latex');
    grid(ax, 'on');
    view(ax, 3);
    box(ax, 'on');

    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, 2*pi]);
    xticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    yticks(ax, [0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticklabels(ax, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
end
