function figs = plot_arnold_tongue_from_saved_result(mat_path)
%PLOT_ARNOLD_TONGUE_FROM_SAVED_RESULT Redraw plots from a saved result MAT file.
%
% Usage:
%   plot_arnold_tongue_from_saved_result()
%   plot_arnold_tongue_from_saved_result('path/to/arnold_tongue_result_*.mat')
%
% If no path is given, the latest arnold_tongue_result_*.mat in this folder
% is used.

% --- User settings for the 1:1 boundary ---
% lock_mode  = 'proportional' -> |ratio - 1| <= lock_value * sigma
% lock_mode  = 'fixed'        -> |ratio - 1| <= lock_value
%lock_mode = 'proportional';
lock_mode = 'fixed';
lock_value = 0.01;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

if nargin < 1 || isempty(mat_path)
    mat_path = find_latest_result_file(script_dir);
end

assert(isfile(mat_path), 'Result MAT file not found: %s', mat_path);
S = load(mat_path);

required_fields = {'ratio_map', 'deltaomega_values', 'sigma_values'};
for k = 1:numel(required_fields)
    assert(isfield(S, required_fields{k}), 'Missing field ''%s'' in MAT file.', required_fields{k});
end

ratio_map = S.ratio_map;
deltaomega_values = S.deltaomega_values;
sigma_values = S.sigma_values;

figs = struct();

%% Heatmap and boundary
figs.heatmap = figure('Color', 'w', 'Name', 'Arnold tongue from saved result');
imagesc(deltaomega_values, sigma_values, ratio_map);
axis xy;
colormap(turbo);
cb = colorbar;
cb.Label.String = '$$\bar{\dot{\phi}}_2/\bar{\dot{\phi}}_1$$';
cb.Label.Interpreter = 'latex';
set(gca, 'Layer', 'top');

xlabel('$$\Delta\omega$$', 'Interpreter', 'latex');
ylabel('$$\sigma$$', 'Interpreter', 'latex');

hold on;
deltaomega_line = linspace(min(deltaomega_values), max(deltaomega_values), 400);
sigma_line = 1/0.184 * abs(deltaomega_line);
h_guide = plot(deltaomega_line, sigma_line, 'b--', 'LineWidth', 1.8, ...
    'DisplayName', '$$\sigma = 0.184|\Delta\omega|$$');
lock_tol_map = compute_lock_tol_map(lock_mode, lock_value, sigma_values);
lock_mask = isfinite(ratio_map) & abs(ratio_map - 1) <= lock_tol_map;
first_hit_deltaomega = nan(size(sigma_values));
last_hit_deltaomega = nan(size(sigma_values));
for i_sigma = 1:numel(sigma_values)
    hit_idx = find(lock_mask(i_sigma, :), 1, 'first');
    if ~isempty(hit_idx)
        first_hit_deltaomega(i_sigma) = deltaomega_values(hit_idx);
    end
    hit_idx = find(lock_mask(i_sigma, :), 1, 'last');
    if ~isempty(hit_idx)
        last_hit_deltaomega(i_sigma) = deltaomega_values(hit_idx);
    end
end

% --- Print 1:1 boundary for requested sigma values (value-based)
% Specify sigma values you want to query here
requested_sigma_values = [3, 5];
ts = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
for ii = 1:numel(requested_sigma_values)
    req = requested_sigma_values(ii);
    [~, idx] = min(abs(sigma_values - req));
    if isempty(idx) || idx < 1 || idx > numel(sigma_values)
        continue;
    end
    matched_sigma = sigma_values(idx);
    f_delta = first_hit_deltaomega(idx);
    l_delta = last_hit_deltaomega(idx);
    if isfinite(f_delta)
        f_str = sprintf('%.6g', f_delta);
    else
        f_str = 'NaN';
    end
    if isfinite(l_delta)
        l_str = sprintf('%.6g', l_delta);
    else
        l_str = 'NaN';
    end
    fprintf('%s\t%s\trequested_sigma=%.6g\tmatched_index=%d\tmatched_sigma=%.6g\tfirst_deltaomega=%s\tlast_deltaomega=%s\n', ...
        ts, mat_path, req, idx, matched_sigma, f_str, l_str);
end
valid_hit = isfinite(first_hit_deltaomega);
first_label = sprintf('1:1 first-hit curve (%s, %.4g)', lock_mode, lock_value);
last_label = sprintf('1:1 last-hit curve (%s, %.4g)', lock_mode, lock_value);
h_first = plot(first_hit_deltaomega(valid_hit), sigma_values(valid_hit), 'k-', 'LineWidth', 2.0, ...
    'DisplayName', first_label);
valid_hit = isfinite(last_hit_deltaomega);
h_last = plot(last_hit_deltaomega(valid_hit), sigma_values(valid_hit), 'k-', 'LineWidth', 2.0, ...
    'DisplayName', last_label);
set([h_first, h_last], 'HandleVisibility', 'off');
h_curve_legend = plot(nan, nan, 'k-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('1:1 curve (%s, %.4g)', lock_mode, lock_value));
set(h_curve_legend, 'Visible', 'off');
legend([h_guide, h_curve_legend], ...
    {'$$\sigma = 0.184|\Delta\omega|$$', ...
    sprintf('1:1 curve (%s, %.4g)', lock_mode, lock_value)}, ...
    'Location', 'best');
tuneFigure;

%% z_opt plot if available
if isfield(S, 'recon') && isstruct(S.recon) && isfield(S.recon, 'z_opt_grid') ...
        && isfield(S.recon.z_opt_grid, 'theta') && isfield(S.recon.z_opt_grid, 'value')
    figs.z_opt = figure('Color', 'w', 'Name', 'z_opt from saved result');
    plot(S.recon.z_opt_grid.theta, S.recon.z_opt_grid.value, '-k', 'LineWidth', 1.5);
    grid on;
    xlabel('\theta');
    ylabel('z_{opt}(\theta)');
    title('reconstructed z_{opt}');
end

%% s2 plot: reuse saved function handle if present; otherwise reconstruct from export path in meta.
have_s2 = false;
if isfield(S, 'recon') && isstruct(S.recon) && isfield(S.recon, 's2_eval') ...
        && isa(S.recon.s2_eval, 'function_handle')
    s2_eval = S.recon.s2_eval;
    have_s2 = true;
elseif isfield(S, 'recon') && isstruct(S.recon) && isfield(S.recon, 'meta') ...
        && isfield(S.recon.meta, 'mat_path') && isfile(S.recon.meta.mat_path)
    try
        opts = struct();
        opts.mat_path = S.recon.meta.mat_path;
        opts.signal_role = S.recon.meta.signal_role;
        if isfield(S.recon.meta, 'agent_id')
            opts.agent_mode = 'agent_id';
            opts.agent_id = S.recon.meta.agent_id;
        end
        recon2 = simulate_winfree_from_export(opts);
        s2_eval = recon2.s2_eval;
        have_s2 = true;
    catch ME
        warning('ArnoldTongueReload:S2ReconstructFailed', ...
            'Could not reconstruct s2 for plotting: %s', ME.message);
    end
end

if have_s2
    n_plot = 121;
    phi_plot = linspace(0, 2*pi, n_plot + 1);
    phi_plot = phi_plot(1:end-1);
    [P1p, P2p] = meshgrid(phi_plot, phi_plot);
    Splot = s2_eval(P1p, P2p);

    figs.s2 = figure('Color', 'w', 'Name', 's2 from saved result');
    imagesc(phi_plot, phi_plot, Splot);
    axis xy;
    axis tight;
    colormap(turbo);
    colorbar;
    xlabel('\phi_1');
    ylabel('\phi_2');
    title('reconstructed s_2(\phi_1,\phi_2)');
else
    warning('ArnoldTongueReload:NoS2', 's2 plot was skipped because s2_eval could not be obtained.');
end

end

function mat_path = find_latest_result_file(script_dir)
files = dir(fullfile(script_dir, 'arnold_tongue_result_*.mat'));
assert(~isempty(files), 'No arnold_tongue_result_*.mat file found in %s', script_dir);
[~, idx] = max([files.datenum]);
mat_path = fullfile(files(idx).folder, files(idx).name);
end

function lock_tol_map = compute_lock_tol_map(lock_mode, lock_value, sigma_values)
switch lower(string(lock_mode))
    case "proportional"
        lock_tol_map = lock_value .* sigma_values(:);
    case "fixed"
        lock_tol_map = lock_value .* ones(size(sigma_values(:)));
    otherwise
        error('Unsupported lock_mode: %s. Use ''proportional'' or ''fixed''.', lock_mode);
end
end