function result = fitDoubleFourierScatter(phi1, phi2, z, M, N)
% fitDoubleFourierScatter Fit a real double Fourier series to scattered data.
%
%   result = fitDoubleFourierScatter(phi1, phi2, z, M, N)
%
% Inputs:
%   phi1, phi2 : angle variables (radians). They are normalized to [0, 2*pi).
%   z          : observed values Q(phi1, phi2)
%   M, N       : truncation orders for phi1 and phi2 directions
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
%   result.rmse
%   result.A
%   result.basis_names
% and additional diagnostic fields.
%
% Example:
%   out = plot_phase_pair_a2_3d_all_files(fullfile('EstimateF','Spring3','250'), [2 3], 3, 15, 5);
%   result = fitDoubleFourierScatter(out.point_cloud.phi1, out.point_cloud.phi2, out.point_cloud.a2, 3, 3);
%
%   % If you already have scattered data directly:
%   % result = fitDoubleFourierScatter(phi1, phi2, z, 4, 4);

    if nargin < 5
        error('Usage: fitDoubleFourierScatter(phi1, phi2, z, M, N)');
    end

    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'M');
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'N');

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

    [A, basis_names] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N);

    % Solve the linear least-squares problem A * coeff ≈ z
    coeff = A \ z;
    z_hat = A * coeff;
    residual = z - z_hat;
    rmse = sqrt(mean(residual .^ 2));

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

    fig_original = figure('Color', 'w');
    ax_original = axes('Parent', fig_original);
    scatter3(ax_original, phi1, phi2, z, 12, z, 'filled');
    formatPhaseScatterAxes(ax_original, '$Q(\phi_1,\phi_2)$', 'Original scattered data');
    colorbar(ax_original);

    fig_fit = figure('Color', 'w');
    ax_fit = axes('Parent', fig_fit);
    scatter3(ax_fit, phi1, phi2, z_hat, 12, z_hat, 'filled');
    formatPhaseScatterAxes(ax_fit, '$\hat{Q}(\phi_1,\phi_2)$', ...
        sprintf('Double Fourier fit (M=%d, N=%d), RMSE=%.4g', M, N, rmse));
    colorbar(ax_fit);

    fig_residual = figure('Color', 'w');
    ax_residual = axes('Parent', fig_residual);
    scatter3(ax_residual, phi1, phi2, residual, 12, residual, 'filled');
    formatPhaseScatterAxes(ax_residual, '$Q - \hat{Q}$', 'Residual scatter');
    colorbar(ax_residual);

    fprintf('[INFO] Double Fourier fit completed: Nsamples=%d, Nbasis=%d, rank(A)=%d, cond(A)=%.3e, RMSE=%.6g\n', ...
        numel(z), n_cols, rankA, condA, rmse);

    result = struct();
    result.coeff = coeff;
    result.z_hat = z_hat;
    result.rmse = rmse;
    result.A = A;
    result.basis_names = basis_names;
    result.residual = residual;
    result.phi1 = phi1;
    result.phi2 = phi2;
    result.z = z;
    result.M = M;
    result.N = N;
    result.rankA = rankA;
    result.condA = condA;
    result.valid_mask = valid_mask;
    result.fig_original = fig_original;
    result.fig_fit = fig_fit;
    result.fig_residual = fig_residual;
end

function [A, basis_names] = buildDoubleFourierDesignMatrix(phi1, phi2, M, N)
% buildDoubleFourierDesignMatrix Construct the real trigonometric basis matrix.
%
% Each column of A corresponds to one named basis function in basis_names.

    n_samples = numel(phi1);
    n_basis = 1 + 2 * M + 2 * N + 4 * M * N;
    A = zeros(n_samples, n_basis);
    basis_names = cell(n_basis, 1);

    col = 1;

    % Constant term: 1
    A(:, col) = 1;
    basis_names{col} = '1';
    col = col + 1;

    % phi1-only terms: cos(m*phi1), sin(m*phi1)
    for m = 1:M
        A(:, col) = cos(m * phi1);
        basis_names{col} = sprintf('cos(%d*phi1)', m);
        col = col + 1;

        A(:, col) = sin(m * phi1);
        basis_names{col} = sprintf('sin(%d*phi1)', m);
        col = col + 1;
    end

    % phi2-only terms: cos(n*phi2), sin(n*phi2)
    for n = 1:N
        A(:, col) = cos(n * phi2);
        basis_names{col} = sprintf('cos(%d*phi2)', n);
        col = col + 1;

        A(:, col) = sin(n * phi2);
        basis_names{col} = sprintf('sin(%d*phi2)', n);
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
            col = col + 1;

            A(:, col) = c1 .* s2;
            basis_names{col} = sprintf('cos(%d*phi1)sin(%d*phi2)', m, n);
            col = col + 1;

            A(:, col) = s1 .* c2;
            basis_names{col} = sprintf('sin(%d*phi1)cos(%d*phi2)', m, n);
            col = col + 1;

            A(:, col) = s1 .* s2;
            basis_names{col} = sprintf('sin(%d*phi1)sin(%d*phi2)', m, n);
            col = col + 1;
        end
    end
end

function formatPhaseScatterAxes(ax, zlabel_text, title_text)
    xlabel(ax, '$\phi_1$ (rad)', 'Interpreter', 'latex');
    ylabel(ax, '$\phi_2$ (rad)', 'Interpreter', 'latex');
    zlabel(ax, zlabel_text, 'Interpreter', 'latex');
    title(ax, title_text, 'Interpreter', 'none');
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
