function out = plot_estimated_equilibria_with_true_gamma(spring_ids, psi_grid)
% plot_estimated_equilibria_with_true_gamma Overlay estimated equilibria and true gamma by spring.
%
% Usage:
%   plot_estimated_equilibria_with_true_gamma
%   plot_estimated_equilibria_with_true_gamma([1 2 3 5])
%   plot_estimated_equilibria_with_true_gamma([1 2], linspace(-pi, pi, 512))
%
% Inputs:
%   spring_ids : optional vector of spring IDs. Default is [1 2 3 5].
%   psi_grid   : optional psi grid used for true gamma reconstruction.
%
% Output:
%   out : struct array with figure handles and source details.

    % Toggle this to control whether figures are saved after tuneFigure.
    enable_save_figure = true;

    existing_figures = findall(groot, 'Type', 'figure');
    if ~isempty(existing_figures)
        close(existing_figures);
        drawnow;
    end

    if nargin < 1 || isempty(spring_ids)
        spring_ids = [1 2 3 5];
    else
        validateattributes(spring_ids, {'numeric'}, {'vector', 'integer', 'positive'}, mfilename, 'spring_ids');
        spring_ids = unique(spring_ids(:).', 'stable');
    end

    if nargin < 2 || isempty(psi_grid)
        psi_grid = [];
    else
        validateattributes(psi_grid, {'numeric'}, {'vector', 'finite'}, mfilename, 'psi_grid');
        psi_grid = psi_grid(:);
    end

    freq_center = 250;
    freq_scale = 100;
    y_scale = pi;
    scale_freq = @(freq_vals) y_scale * (freq_vals - freq_center) / freq_scale;

    spring_data_all = build_estimated_spring_data(scale_freq);
    available_ids = [spring_data_all.id];

    missing_ids = spring_ids(~ismember(spring_ids, available_ids));
    if ~isempty(missing_ids)
        warning('plot_estimated_equilibria_with_true_gamma:missingSpringData', ...
            'No estimated equilibrium data found for spring IDs: %s', mat2str(missing_ids));
    end

    spring_ids = spring_ids(ismember(spring_ids, available_ids));
    if isempty(spring_ids)
        error('No plottable spring IDs remain after filtering.');
    end

    out = repmat(struct( ...
        'spring_id', [], ...
        'figure', gobjects(1), ...
        'true_gamma_path', '', ...
        'true_gamma_loaded', false), 1, numel(spring_ids));

    for idx = 1:numel(spring_ids)
        spring_id = spring_ids(idx);
        spring_data = spring_data_all([spring_data_all.id] == spring_id);

        true_gamma_path = fullfile('EstimateQ', sprintf('Spring%d', spring_id), '255', 'gamma_exports', 'gamma_export_latest.mat');
        [psi_true, gamma_true, true_loaded] = load_true_gamma_curve(true_gamma_path, psi_grid);

        fig = figure('Color', 'w', 'Name', sprintf('Spring %d: estimated equilibria + true gamma', spring_id));
        ax = axes('Parent', fig);
        hold(ax, 'on');

        h_stable = scatter(ax, spring_data.phi_stable, spring_data.freq_stable_norm, 72, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerFaceColor', [0.00, 0.45, 0.74], 'MarkerEdgeColor', [0.00, 0.45, 0.74], ...
            'DisplayName', 'Stable EQP');
        h_unstable = scatter(ax, spring_data.phi_unstable, spring_data.freq_unstable_norm, 72, 's', 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerFaceColor', [0.85, 0.33, 0.10], 'MarkerEdgeColor', [0.85, 0.33, 0.10], ...
            'DisplayName', 'Unstable EQP');
        if true_loaded
            h_true = plot(ax, psi_true, gamma_true, '-', 'LineWidth', 2.0, ...
                'Color', [0.20, 0.20, 0.20], 'DisplayName', '$$f(\psi)$$');
            legend(ax, [h_stable, h_unstable, h_true], ...
                {'Stable EQP', 'Unstable EQP', '$$f(\psi)$$'}, 'Location', 'best');
            apply_padded_ylim(ax, [spring_data.freq_stable_norm(:); spring_data.freq_unstable_norm(:); gamma_true(:)]);
        else
            legend(ax, [h_stable, h_unstable], ...
                {'Stable EQP', 'Unstable EQP'}, 'Location', 'best');
            apply_padded_ylim(ax, [spring_data.freq_stable_norm(:); spring_data.freq_unstable_norm(:)]);
        end
        ylabel(ax, '$$f(\psi)$$', 'Interpreter', 'latex');

        xlabel(ax, '$$\psi$$', 'Interpreter', 'latex');
        title(ax, sprintf('Spring %d', spring_id));
        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, [-pi, pi]);
        xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
        xticklabels(ax, {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
        ax.TickLabelInterpreter = 'latex';
        set(ax, 'FontSize', 13, 'LineWidth', 1.2);

        if exist('tuneFigure', 'file') == 2
            tuneFigure;
            if enable_save_figure && exist('saveFigure', 'file') == 2
                saveFigure;
            end
        elseif enable_save_figure && exist('saveFigure', 'file') == 2
            saveFigure;
        end

        out(idx).spring_id = spring_id;
        out(idx).figure = fig;
        out(idx).true_gamma_path = true_gamma_path;
        out(idx).true_gamma_loaded = true_loaded;
    end
end

function spring_data = build_estimated_spring_data(scale_freq)
    spring_data = struct('id', {}, 'phi_stable', {}, 'phi_unstable', {}, 'freq_stable_norm', {}, 'freq_unstable_norm', {});

    % Spring 1
    freq1_all = [ ...
        240 240, ...
        245 245 245 245, ...
        250 250 250 250, ...
        255 255 255 255, ...
        260 260];

    phi1_all = [ ...
       -0.51722,  -1.032, ...
       -0.349878, 2.69335, 1.35462, -0.9891, ...
       -0.195281, 3.0993,  1.18133, -1.402885, ...
        0.015736, -2.91941, 0.976656, -1.793484, ...
        0.266414, 0.463211];
    phi1_all = -phi1_all;

    [freq1_stable, freq1_unstable, phi1_stable, phi1_unstable] = split_stable_unstable_by_count(...
        freq1_all, phi1_all, [2, 4, 4, 4, 2]);

    spring_data(end + 1) = struct( ...
        'id', 1, ...
        'phi_stable', phi1_stable, ...
        'phi_unstable', phi1_unstable, ...
        'freq_stable_norm', scale_freq(freq1_stable), ...
        'freq_unstable_norm', scale_freq(freq1_unstable)); %#ok<AGROW>

    % Spring 2
    freq2 = [230 235 240 245 250 255 260 265];
    phi2_stable = [-0.790055, -0.666673, -0.453907, -0.283704, 0.049815, 0.125068, 0.411129, 0.641414];
    phi2_unstable = [-1.36336, -1.58014, -2.01999, -2.26026, -3.07159, 2.32331, 1.86263, 1.56455];
    phi2_stable = -phi2_stable;
    phi2_unstable = -phi2_unstable;

    spring_data(end + 1) = struct( ...
        'id', 2, ...
        'phi_stable', phi2_stable, ...
        'phi_unstable', phi2_unstable, ...
        'freq_stable_norm', scale_freq(freq2), ...
        'freq_unstable_norm', scale_freq(freq2)); %#ok<AGROW>

    % Spring 3
    freq3 = [230 235 240 245 250 255 260 265];
    phi3_stable = [-0.861, -0.637, -0.475, -0.282, 0.000, 0.246, 0.474, 0.669];
    phi3_unstable = [-1.226, -2.133, -2.612, -2.787, -2.946, 2.254, 1.356, 1.114];
    phi3_stable = -phi3_stable;
    phi3_unstable = -phi3_unstable;

    spring_data(end + 1) = struct( ...
        'id', 3, ...
        'phi_stable', phi3_stable, ...
        'phi_unstable', phi3_unstable, ...
        'freq_stable_norm', scale_freq(freq3), ...
        'freq_unstable_norm', scale_freq(freq3)); %#ok<AGROW>

    % Spring 5
    freq5 = [235 240 245 250 255 260 265 270];
    phi5_stable = [-0.969, -0.687, -0.389, -0.090, 0.439, 0.672, 0.871, 1.257];
    phi5_unstable = [-1.934, -2.789, -3.038, 3.035, 2.688, 2.226, 2.178, 1.957];
    phi5_stable = -phi5_stable;
    phi5_unstable = -phi5_unstable;

    spring_data(end + 1) = struct( ...
        'id', 5, ...
        'phi_stable', phi5_stable, ...
        'phi_unstable', phi5_unstable, ...
        'freq_stable_norm', scale_freq(freq5), ...
        'freq_unstable_norm', scale_freq(freq5)); %#ok<AGROW>
end

function [freq_stable, freq_unstable, phi_stable, phi_unstable] = split_stable_unstable_by_count(freq_all, phi_all, per_freq_counts)
    validateattributes(freq_all, {'numeric'}, {'vector', 'finite'}, mfilename, 'freq_all');
    validateattributes(phi_all, {'numeric'}, {'vector', 'finite'}, mfilename, 'phi_all');
    validateattributes(per_freq_counts, {'numeric'}, {'vector', 'integer', 'positive'}, mfilename, 'per_freq_counts');

    if numel(freq_all) ~= numel(phi_all)
        error('freq_all and phi_all must have identical lengths.');
    end
    if sum(per_freq_counts) ~= numel(freq_all)
        error('sum(per_freq_counts) must match numel(freq_all).');
    end

    freq_stable = [];
    freq_unstable = [];
    phi_stable = [];
    phi_unstable = [];

    cursor = 1;
    for idx = 1:numel(per_freq_counts)
        n = per_freq_counts(idx);
        if mod(n, 2) ~= 0
            error('Each count in per_freq_counts must be even.');
        end

        ids = cursor:(cursor + n - 1);
        split_idx = n / 2;
        stable_ids = ids(1:split_idx);
        unstable_ids = ids(split_idx + 1:end);

        freq_stable = [freq_stable, freq_all(stable_ids)]; %#ok<AGROW>
        freq_unstable = [freq_unstable, freq_all(unstable_ids)]; %#ok<AGROW>
        phi_stable = [phi_stable, phi_all(stable_ids)]; %#ok<AGROW>
        phi_unstable = [phi_unstable, phi_all(unstable_ids)]; %#ok<AGROW>

        cursor = cursor + n;
    end
end

function [psi_true, gamma_true, loaded] = load_true_gamma_curve(true_gamma_mat_path, psi_grid)
    psi_true = [];
    gamma_true = [];
    loaded = false;

    if exist(true_gamma_mat_path, 'file') ~= 2
        warning('plot_estimated_equilibria_with_true_gamma:fileMissing', ...
            'true gamma file not found: %s', true_gamma_mat_path);
        return;
    end

    loaded_mat = load(true_gamma_mat_path);
    if ~isfield(loaded_mat, 'gamma_export') || ~isstruct(loaded_mat.gamma_export)
        warning('plot_estimated_equilibria_with_true_gamma:noGammaExport', ...
            'gamma_export not found in: %s', true_gamma_mat_path);
        return;
    end
    if ~isfield(loaded_mat.gamma_export, 'true_gamma')
        warning('plot_estimated_equilibria_with_true_gamma:noTrueGamma', ...
            'true_gamma field not found in: %s', true_gamma_mat_path);
        return;
    end

    true_gamma = loaded_mat.gamma_export.true_gamma;
    if ~isstruct(true_gamma) || ~isfield(true_gamma, 'available') || ~true_gamma.available
        warning('plot_estimated_equilibria_with_true_gamma:trueGammaUnavailable', ...
            'true_gamma is not available in: %s', true_gamma_mat_path);
        return;
    end

    try
        if isempty(psi_grid)
            [psi_true, gamma_true] = reconstruct_exported_gamma(true_gamma, [], 'true');
        else
            [psi_true, gamma_true] = reconstruct_exported_gamma(true_gamma, psi_grid, 'true');
        end
    catch err
        warning('plot_estimated_equilibria_with_true_gamma:reconstructFailed', ...
            'Failed to reconstruct true gamma from %s (%s)', true_gamma_mat_path, err.message);
        psi_true = [];
        gamma_true = [];
        loaded = false;
        return;
    end

    psi_true = psi_true(:);
    gamma_true = gamma_true(:);
    loaded = ~isempty(psi_true) && ~isempty(gamma_true);
end

function apply_padded_ylim(ax, y_values)
    y_values = y_values(:);
    if isempty(y_values)
        return;
    end
    y_min = min(y_values);
    y_max = max(y_values);
    if ~isfinite(y_min) || ~isfinite(y_max)
        return;
    end
    pad = max(1e-6, 0.10 * (y_max - y_min));
    ylim(ax, [y_min - pad, y_max + pad]);
end