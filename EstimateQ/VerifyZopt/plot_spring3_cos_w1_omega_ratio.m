function plot_spring3_cos_w1_omega_ratio(enable_save_figure)
% Plot omega ratio curves for both:
%   EstimateQ/VerifyZopt/Spring3/cos/**/merged_*.csv
%   EstimateQ/VerifyZopt/Spring3/w1/**/merged_*.csv
% Reuses compute_avg_omega_id6_kappa5 in per-file max-vs-min mode.

t_start = 60;
t_end = 120;
t_start_series = 0;  % start time for time-series overlay only
mean_phase_band_half_width_rad = 0.10*pi;
if nargin < 1 || isempty(enable_save_figure)
    enable_save_figure = true;
end

% Clear any existing figures before generating new plots.
close all;

rootDirCos = fullfile(pwd, 'EstimateQ', 'VerifyZopt', 'Spring3', 'k3w2');
rootDirW1  = fullfile(pwd, 'EstimateQ', 'VerifyZopt', 'Spring3', 'k5w2');

[xCos, yCos, phaseStatsCos, folderLabelsCos] = local_collect_ratio_curve(rootDirCos, t_start, t_end);
[xW1, yW1, phaseStatsW1, folderLabelsW1]  = local_collect_ratio_curve(rootDirW1, t_start, t_end);

isConvergedCos = local_build_folder_convergence_mask(folderLabelsCos, phaseStatsCos, mean_phase_band_half_width_rad);
isConvergedW1 = local_build_folder_convergence_mask(folderLabelsW1, phaseStatsW1, mean_phase_band_half_width_rad);

fprintf('\n[INFO] Per-file phase-difference stats for Spring3/cos (window %.2f-%.2f s)\n', t_start, t_end);
if isempty(phaseStatsCos)
    fprintf('[INFO] No valid phase-difference samples found for Spring3/cos.\n');
else
    disp(phaseStatsCos);
end

fprintf('\n[INFO] Per-file phase-difference stats for Spring3/w1 (window %.2f-%.2f s)\n', t_start, t_end);
if isempty(phaseStatsW1)
    fprintf('[INFO] No valid phase-difference samples found for Spring3/w1.\n');
else
    disp(phaseStatsW1);
end

figure('Color', 'w');
hW1Trend = plot(xW1, yW1, '-', 'LineWidth', 1.2, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '$z_{\mathrm{opt}}(\theta), q^*$');
hold on;
hCosTrend = plot(xCos, yCos, '-', 'LineWidth', 1.2, 'Color', [0.0, 0.4470, 0.7410], 'DisplayName', '$z_{\mathrm{sin}}(\theta), \tau_1^*$');

hW1Conv = plot(xW1(isConvergedW1), yW1(isConvergedW1), 'p', ...
    'LineStyle', 'none', 'MarkerSize', 10, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'DisplayName', '$z_{\mathrm{opt}}(\theta), q^*$ converged');
plot(xW1(~isConvergedW1), yW1(~isConvergedW1), 's', ...
    'LineStyle', 'none', 'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]);

hCosConv = plot(xCos(isConvergedCos), yCos(isConvergedCos), 'd', ...
    'LineStyle', 'none', 'MarkerSize', 6, 'MarkerFaceColor', [0.0, 0.4470, 0.7410], ...
    'MarkerEdgeColor', [0.0, 0.4470, 0.7410], 'DisplayName', '$z_{\mathrm{sin}}(\theta), \tau_1^*$ converged');
plot(xCos(~isConvergedCos), yCos(~isConvergedCos), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.0, 0.4470, 0.7410]);
grid on;
x_limits = [min([xCos(:); xW1(:)]), max([xCos(:); xW1(:)])];
xlim(x_limits);
xtick_min = floor(x_limits(1) * 10) / 10;
xtick_max = ceil(x_limits(2) * 10) / 10;
xticks(-2:1:2);
xlabel('$$\Delta\omega$$');
ylabel('$$\bar{\dot{\phi_2}} / \bar{\dot{\phi_1}}$$');
legend([hCosTrend, hW1Trend, hCosConv, hW1Conv], 'Location', 'best', 'Interpreter', 'latex');
lgd = legend(gca);
if ~isempty(lgd) && isgraphics(lgd)
    lgd.Location = 'eastoutside';
end
hold off;

if exist('tuneFigure', 'file') == 2
    tuneFigure;
end
if enable_save_figure
    if exist('saveFigure', 'file') == 2
        saveFigure;
    else
        warning('enable_save_figure is true, but saveFigure.m was not found.');
    end
end

local_plot_mean_phase_vs_delta_omega(phaseStatsCos, phaseStatsW1, mean_phase_band_half_width_rad);
local_plot_phase_evolution_in_mean_band(phaseStatsCos, phaseStatsW1, t_start_series, t_end, mean_phase_band_half_width_rad);

if enable_save_figure
    if exist('saveFigure', 'file') == 2
        saveFigure;
    else
        warning('enable_save_figure is true, but saveFigure.m was not found.');
    end
end
end


function [x, y, phase_stats, folder_labels] = local_collect_ratio_curve(rootDir, t_start, t_end)
if ~isfolder(rootDir)
    error('Folder not found: %s', rootDir);
end

dirs = dir(rootDir);
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
if isempty(dirs)
    error('No subdirectories found under %s', rootDir);
end

% Keep only numeric folder names (e.g., 240..260) and sort ascending.
labelsNum = nan(1, numel(dirs));
for i = 1:numel(dirs)
    v = str2double(dirs(i).name);
    if ~isnan(v)
        labelsNum(i) = v;
    end
end
valid = ~isnan(labelsNum);
dirs = dirs(valid);
labelsNum = labelsNum(valid);
[labelsNum, ord] = sort(labelsNum);
dirs = dirs(ord);

if isempty(dirs)
    error('No numeric subdirectories found under %s', rootDir);
end

omega2 = (labelsNum ./ 100) * pi;
x = omega2 - 2.5 * pi;  % DeltaOmega
y = nan(size(x));
folder_labels = labelsNum(:);
phase_rows = struct('folder_label', {}, 'file', {}, 'base_agent_id', {}, 'target_agent_id', {}, ...
    'mean_phase_diff_rad', {}, 'mean_phase_diff_deg', {}, ...
    'max_abs_phase_deviation_rad', {}, 'n_samples', {}, 'status', {});

for i = 1:numel(dirs)
    folderPath = fullfile(rootDir, dirs(i).name);
    try
        T = compute_avg_omega_id6_kappa5(folderPath, t_start, t_end, NaN, false, false, 'per_file_max_vs_min');
    catch ME
        warning('Failed on %s: %s', folderPath, ME.message);
        continue;
    end

    if isempty(T) || ~ismember('ratio_to_base', T.Properties.VariableNames)
        continue;
    end

    for j = 1:height(T)
        if ~ismember('file', T.Properties.VariableNames)
            continue;
        end
        csvFileName = local_table_file_name(T.file, j);
        if isempty(csvFileName)
            continue;
        end
        csvPath = fullfile(folderPath, csvFileName);
        stats = local_compute_phase_diff_stats(csvPath, t_start, t_end);
        stats.folder_label = labelsNum(i);
        phase_rows(end + 1) = stats; %#ok<AGROW>
    end

    r = T.ratio_to_base;
    r = r(~isnan(r));
    if isempty(r)
        continue;
    end

    % Use mean in case multiple CSV files exist in a folder.
    y(i) = mean(r);
end

if isempty(phase_rows)
    phase_stats = table();
else
    phase_stats = struct2table(phase_rows);
    phase_stats = sortrows(phase_stats, {'folder_label', 'file'});
end
end

function is_converged = local_build_folder_convergence_mask(folder_labels, phase_stats, mean_phase_band_half_width_rad)
folder_labels = folder_labels(:);
is_converged = false(size(folder_labels));

if isempty(folder_labels) || isempty(phase_stats) || ~istable(phase_stats)
    return;
end

required = {'folder_label', 'max_abs_phase_deviation_rad', 'status'};
if ~all(ismember(required, phase_stats.Properties.VariableNames))
    return;
end

status_col = phase_stats.status;
if iscell(status_col)
    ok_mask = strcmp(status_col, 'ok');
elseif isstring(status_col)
    ok_mask = strcmp(cellstr(status_col), 'ok');
else
    ok_mask = false(height(phase_stats), 1);
end

stats_folder = double(phase_stats.folder_label);
max_dev = double(phase_stats.max_abs_phase_deviation_rad);
valid = ok_mask & isfinite(stats_folder) & isfinite(max_dev);
if ~any(valid)
    return;
end

stats_folder = stats_folder(valid);
max_dev = max_dev(valid);

for ii = 1:numel(folder_labels)
    mask_i = (stats_folder == folder_labels(ii));
    if ~any(mask_i)
        continue;
    end

    % A folder is treated as converged only when all valid files in it
    % satisfy the phase-band condition.
    is_converged(ii) = all(max_dev(mask_i) <= mean_phase_band_half_width_rad);
end
end

function stats = local_compute_phase_diff_stats(csvPath, t_start, t_end)
stats = struct( ...
    'folder_label', NaN, ...
    'file', csvPath, ...
    'base_agent_id', NaN, ...
    'target_agent_id', NaN, ...
    'mean_phase_diff_rad', NaN, ...
    'mean_phase_diff_deg', NaN, ...
    'max_abs_phase_deviation_rad', NaN, ...
    'n_samples', 0, ...
    'status', 'unknown');

try
    [series_by_agent, all_agents] = load_corrected_agent_series_from_csv(csvPath, [], {'time_pc_sec_abs', 'a0'});
catch ME
    stats.status = 'read_error';
    warning('Failed to load %s: %s', csvPath, ME.message);
    return;
end

agents = all_agents(:).';
agents = agents(agents ~= 99);
if numel(agents) < 2
    stats.status = 'not_enough_agents';
    return;
end

base_id = min(agents);
target_id = max(agents);
stats.base_agent_id = base_id;
stats.target_agent_id = target_id;

series_base = series_by_agent(base_id);
series_target = series_by_agent(target_id);
if isempty(series_base.time) || isempty(series_target.time)
    stats.status = 'missing_series';
    return;
end

overlap_start = max(min(series_base.time), min(series_target.time));
overlap_end = min(max(series_base.time), max(series_target.time));
if overlap_end <= overlap_start
    stats.status = 'no_overlap';
    return;
end

window_start_abs = overlap_start + t_start;
window_end_abs = min(overlap_end, overlap_start + t_end);
if window_end_abs <= window_start_abs
    stats.status = 'window_outside_overlap';
    return;
end

tb = series_base.time(:);
tt = series_target.time(:);
base_mask = tb >= window_start_abs & tb <= window_end_abs;
target_mask = tt >= window_start_abs & tt <= window_end_abs;

t_grid = unique([window_start_abs; window_end_abs; tb(base_mask); tt(target_mask)]);
if numel(t_grid) < 2
    stats.status = 'insufficient_points';
    return;
end

theta_base = interp1(tb, double(series_base.a0_corr(:)) * (2*pi/256), t_grid, 'linear', NaN);
theta_target = interp1(tt, double(series_target.a0_corr(:)) * (2*pi/256), t_grid, 'linear', NaN);
valid = isfinite(theta_base) & isfinite(theta_target);
if nnz(valid) < 2
    stats.status = 'insufficient_interpolated_points';
    return;
end

delta_phase = atan2(sin(theta_target(valid) - theta_base(valid)), cos(theta_target(valid) - theta_base(valid)));
if isempty(delta_phase)
    stats.status = 'empty_delta';
    return;
end

mean_phase = atan2(mean(sin(delta_phase)), mean(cos(delta_phase)));
phase_dev = atan2(sin(delta_phase - mean_phase), cos(delta_phase - mean_phase));

stats.mean_phase_diff_rad = mean_phase;
stats.mean_phase_diff_deg = mean_phase * (180 / pi);
stats.max_abs_phase_deviation_rad = max(abs(phase_dev));
stats.n_samples = numel(delta_phase);
stats.status = 'ok';
end

function file_name = local_table_file_name(file_col, idx)
file_name = '';

if iscell(file_col)
    value = file_col{idx};
elseif isstring(file_col)
    value = file_col(idx);
elseif ischar(file_col)
    if size(file_col, 1) < idx
        return;
    end
    value = file_col(idx, :);
else
    value = file_col(idx);
end

if isstring(value)
    value = char(value);
elseif isnumeric(value)
    value = num2str(value);
end

if ischar(value)
    file_name = strtrim(value);
end
end

function local_plot_mean_phase_vs_delta_omega(phaseStatsCos, phaseStatsW1, mean_phase_band_half_width_rad)
if nargin < 3 || isempty(mean_phase_band_half_width_rad)
    mean_phase_band_half_width_rad = 0.2;
end

[filesCos, ~] = local_select_mean_band_files(phaseStatsCos, mean_phase_band_half_width_rad, 'cos');
[filesW1, ~] = local_select_mean_band_files(phaseStatsW1, mean_phase_band_half_width_rad, 'w1');
phaseStatsCosFiltered = local_filter_phase_stats_by_files(phaseStatsCos, filesCos);
phaseStatsW1Filtered = local_filter_phase_stats_by_files(phaseStatsW1, filesW1);

figure('Color', 'w');
hold on;
grid on;

[xCosRaw, yCosRaw, xCosMean, yCosMean] = local_extract_mean_phase_series(phaseStatsCosFiltered);
[xW1Raw, yW1Raw, xW1Mean, yW1Mean] = local_extract_mean_phase_series(phaseStatsW1Filtered);

if ~isempty(xCosRaw)
    scatter(xCosRaw, yCosRaw, 20, 'o', 'MarkerEdgeColor', [0.0, 0.4470, 0.7410], ...
        'MarkerFaceColor', 'none', 'DisplayName', 'Spring3/cos per-file');
end
if ~isempty(xW1Raw)
    scatter(xW1Raw, yW1Raw, 20, 's', 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
        'MarkerFaceColor', 'none', 'DisplayName', 'Spring3/w1 per-file');
end
if ~isempty(xCosMean)
    plot(xCosMean, yCosMean, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', 'auto', 'DisplayName', 'Spring3/cos mean');
end
if ~isempty(xW1Mean)
    plot(xW1Mean, yW1Mean, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', 'auto', 'DisplayName', 'Spring3/w1 mean');
end

overlay_data_phase_mean = struct();
overlay_data_phase_mean.xCosMean = xCosMean(:);
overlay_data_phase_mean.yCosMean = yCosMean(:);
overlay_data_phase_mean.xW1Mean = xW1Mean(:);
overlay_data_phase_mean.yW1Mean = yW1Mean(:);
assignin('base', 'overlay_data_spring3_phase_mean', overlay_data_phase_mean);

x_all = [xCosRaw(:); xW1Raw(:); xCosMean(:); xW1Mean(:)];
if ~isempty(x_all)
    x_limits = [min(x_all), max(x_all)];
    xlim(x_limits);
    xtick_min = floor(x_limits(1) * 10) / 10;
    xtick_max = ceil(x_limits(2) * 10) / 10;
    xticks(xtick_min:0.1:xtick_max);
end

ylim([-pi, pi]);
yticks([-pi, 0, pi]);
yticklabels({'-\pi', '0', '\pi'});
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('$$\Delta\omega$$');
ylabel('$$\langle\phi_{\mathrm{target}}-\phi_{\mathrm{base}}\rangle$$');
title(sprintf('Mean phase difference for files within mean phase band (|\Delta\phi-\mu| \le %.3f rad)', mean_phase_band_half_width_rad), 'Interpreter', 'none');
legend('Location', 'best');
hold off;

if exist('tuneFigure', 'file') == 2
    tuneFigure;
end
end

function phase_stats_filtered = local_filter_phase_stats_by_files(phase_stats, file_paths)
if isempty(phase_stats) || ~istable(phase_stats)
    phase_stats_filtered = table();
    return;
end
if ~ismember('file', phase_stats.Properties.VariableNames)
    phase_stats_filtered = table();
    return;
end
if isempty(file_paths)
    phase_stats_filtered = table();
    return;
end

file_list = cellfun(@char, file_paths, 'UniformOutput', false);
mask = false(height(phase_stats), 1);
for i = 1:height(phase_stats)
    file_path = local_table_file_name(phase_stats.file, i);
    if isempty(file_path)
        continue;
    end
    mask(i) = any(strcmp(file_list, file_path));
end

phase_stats_filtered = phase_stats(mask, :);
end

function [x_raw, y_raw, x_mean, y_mean] = local_extract_mean_phase_series(phase_stats)
x_raw = [];
y_raw = [];
x_mean = [];
y_mean = [];

if isempty(phase_stats) || ~istable(phase_stats)
    return;
end

required = {'folder_label', 'mean_phase_diff_rad', 'status'};
if ~all(ismember(required, phase_stats.Properties.VariableNames))
    return;
end

status_col = phase_stats.status;
if iscell(status_col)
    ok_mask = strcmp(status_col, 'ok');
elseif isstring(status_col)
    ok_mask = strcmp(cellstr(status_col), 'ok');
else
    ok_mask = false(height(phase_stats), 1);
end

mu = double(phase_stats.mean_phase_diff_rad);
folder_label = double(phase_stats.folder_label);
ok_mask = ok_mask & isfinite(mu) & isfinite(folder_label);
if ~any(ok_mask)
    return;
end

folder_ok = folder_label(ok_mask);
mu_ok = mu(ok_mask);

x_raw = (folder_ok ./ 100) * pi - 2.5 * pi;
y_raw = mu_ok;

u = unique(folder_ok(:).');
x_mean = nan(size(u));
y_mean = nan(size(u));
for i = 1:numel(u)
    mask_i = folder_ok == u(i);
    if ~any(mask_i)
        continue;
    end
    x_mean(i) = (u(i) / 100) * pi - 2.5 * pi;
    y_mean(i) = atan2(mean(sin(mu_ok(mask_i))), mean(cos(mu_ok(mask_i))));
end
end

function local_plot_phase_evolution_in_mean_band(phaseStatsCos, phaseStatsW1, t_start, t_end, mean_phase_band_half_width_rad)
[filesCos, labelsCos] = local_select_mean_band_files(phaseStatsCos, mean_phase_band_half_width_rad, 'cos');
[filesW1, labelsW1] = local_select_mean_band_files(phaseStatsW1, mean_phase_band_half_width_rad, 'w1');

if isempty(filesCos) && isempty(filesW1)
    fprintf('[INFO] No files found with max phase deviation from mean <= %.3f rad.\n', mean_phase_band_half_width_rad);
    return;
end

figure('Color', 'w');
ax = axes('Parent', gcf);
hold(ax, 'on');

line_handles = gobjects(0, 1);
legend_labels = {};

for i = 1:numel(filesCos)
    [t_rel, delta_phase, ok] = local_compute_phase_diff_timeseries(filesCos{i}, t_start, t_end);
    if ~ok
        continue;
    end
    h = plot(ax, t_rel, delta_phase, 'Color', [0.0, 0.4470, 0.7410], 'LineWidth', 0.8, 'LineStyle', '-');
    if i == 1
        line_handles(end + 1) = h; %#ok<AGROW>
        legend_labels{end + 1} = sprintf('cos (N=%d)', numel(filesCos)); %#ok<AGROW>
    end
    if mod(i, 20) == 0
        fprintf('[INFO] cos overlay progress: %d/%d\n', i, numel(filesCos));
    end
end

for i = 1:numel(filesW1)
    [t_rel, delta_phase, ok] = local_compute_phase_diff_timeseries(filesW1{i}, t_start, t_end);
    if ~ok
        continue;
    end
    h = plot(ax, t_rel, delta_phase, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 0.8, 'LineStyle', '--');
    if i == 1
        line_handles(end + 1) = h; %#ok<AGROW>
        legend_labels{end + 1} = sprintf('w1 (N=%d)', numel(filesW1)); %#ok<AGROW>
    end
    if mod(i, 20) == 0
        fprintf('[INFO] w1 overlay progress: %d/%d\n', i, numel(filesW1));
    end
end

yline(ax, 0, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 0.8);
ylim(ax, [-pi, pi]);
yticks(ax, [-pi, 0, pi]);
yticklabels(ax, {'-\pi', '0', '\pi'});
set(ax, 'TickLabelInterpreter', 'latex');
xlabel(ax, 'Time since t_{start} (s)', 'Interpreter', 'latex');
ylabel(ax, '$$\phi_{\mathrm{target}}-\phi_{\mathrm{base}}$$', 'Interpreter', 'latex');
title(ax, sprintf('Phase-difference overlays (|\Delta\phi-\mu| \le %.2f rad)', mean_phase_band_half_width_rad), 'Interpreter', 'none');
grid(ax, 'on');

if ~isempty(line_handles)
    legend(ax, line_handles, legend_labels, 'Location', 'best');
end

hold(ax, 'off');

if exist('tuneFigure', 'file') == 2
    tuneFigure;
end

fprintf('[INFO] Overlayed mean-band files: cos=%d, w1=%d\n', numel(filesCos), numel(filesW1));

if ~isempty(labelsCos)
    fprintf('[INFO] cos files (max |phase-mean| <= %.3f rad):\n', mean_phase_band_half_width_rad);
    disp(labelsCos(:));
end
if ~isempty(labelsW1)
    fprintf('[INFO] w1 files (max |phase-mean| <= %.3f rad):\n', mean_phase_band_half_width_rad);
    disp(labelsW1(:));
end
end

function [file_paths, labels] = local_select_mean_band_files(phase_stats, mean_phase_band_half_width_rad, prefix_label)
file_paths = {};
labels = {};

if isempty(phase_stats) || ~istable(phase_stats)
    return;
end

required = {'file', 'max_abs_phase_deviation_rad', 'status'};
if ~all(ismember(required, phase_stats.Properties.VariableNames))
    return;
end

status_col = phase_stats.status;
if iscell(status_col)
    ok_mask = strcmp(status_col, 'ok');
elseif isstring(status_col)
    ok_mask = strcmp(cellstr(status_col), 'ok');
else
    ok_mask = false(height(phase_stats), 1);
end

max_dev = double(phase_stats.max_abs_phase_deviation_rad);
valid_mask = ok_mask & isfinite(max_dev) & (max_dev <= mean_phase_band_half_width_rad);
if ~any(valid_mask)
    return;
end

rows = find(valid_mask);
for i = 1:numel(rows)
    idx = rows(i);
    file_path = local_table_file_name(phase_stats.file, idx);
    if isempty(file_path)
        continue;
    end
    if ~isfile(file_path)
        continue;
    end
    file_paths{end + 1} = file_path; %#ok<AGROW>
    labels{end + 1} = sprintf('%s | max|phase-mean|=%.4f rad | %s', prefix_label, max_dev(idx), file_path); %#ok<AGROW>
end

if isempty(file_paths)
    return;
end

[file_paths, ia] = unique(file_paths, 'stable');
labels = labels(ia);
end

function [t_rel, delta_phase, ok] = local_compute_phase_diff_timeseries(csvPath, t_start, t_end)
t_rel = [];
delta_phase = [];
ok = false;

try
    [series_by_agent, all_agents] = load_corrected_agent_series_from_csv(csvPath, [], {'time_pc_sec_abs', 'a0'});
catch
    return;
end

agents = all_agents(:).';
agents = agents(agents ~= 99);
if numel(agents) < 2
    return;
end

base_id = min(agents);
target_id = max(agents);

series_base = series_by_agent(base_id);
series_target = series_by_agent(target_id);
if isempty(series_base.time) || isempty(series_target.time)
    return;
end

overlap_start = max(min(series_base.time), min(series_target.time));
overlap_end = min(max(series_base.time), max(series_target.time));
if overlap_end <= overlap_start
    return;
end

window_start_abs = overlap_start + t_start;
window_end_abs = min(overlap_end, overlap_start + t_end);
if window_end_abs <= window_start_abs
    return;
end

tb = series_base.time(:);
tt = series_target.time(:);
base_mask = tb >= window_start_abs & tb <= window_end_abs;
target_mask = tt >= window_start_abs & tt <= window_end_abs;

t_grid = unique([window_start_abs; window_end_abs; tb(base_mask); tt(target_mask)]);
if numel(t_grid) < 2
    return;
end

theta_base = interp1(tb, double(series_base.a0_corr(:)) * (2*pi/256), t_grid, 'linear', NaN);
theta_target = interp1(tt, double(series_target.a0_corr(:)) * (2*pi/256), t_grid, 'linear', NaN);
valid = isfinite(theta_base) & isfinite(theta_target);
if nnz(valid) < 2
    return;
end

t_valid = t_grid(valid);
delta_phase = atan2(sin(theta_target(valid) - theta_base(valid)), cos(theta_target(valid) - theta_base(valid)));
t_rel = t_valid - window_start_abs;
t_rel = max(t_rel, 0);
ok = true;
end
