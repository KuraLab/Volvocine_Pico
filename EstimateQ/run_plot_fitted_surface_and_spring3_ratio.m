function run_plot_fitted_surface_and_spring3_ratio(run_surface_plot, run_ratio_plot, run_overlay_plot, run_overlay_s1, run_overlay_s2)
% Run both analysis scripts in one entry point:
%   1) EstimateQ/plot_fitted_surface_from_export.m
%   2) EstimateQ/VerifyZopt/plot_spring3_cos_w1_omega_ratio.m
%
% Usage:
%   run_plot_fitted_surface_and_spring3_ratio
%   run_plot_fitted_surface_and_spring3_ratio(true, true)
%   run_plot_fitted_surface_and_spring3_ratio(true, false)
%   run_plot_fitted_surface_and_spring3_ratio(false, true)
%   run_plot_fitted_surface_and_spring3_ratio(true, true, true)
%   run_plot_fitted_surface_and_spring3_ratio(true, true, true, true, false)

if nargin < 1 || isempty(run_surface_plot)
    run_surface_plot = true;
end
if nargin < 2 || isempty(run_ratio_plot)
    run_ratio_plot = true;
end
if nargin < 3 || isempty(run_overlay_plot)
    run_overlay_plot = true;
end
if nargin < 4 || isempty(run_overlay_s1)
    run_overlay_s1 = false;
end
if nargin < 5 || isempty(run_overlay_s2)
    run_overlay_s2 = true;
end

if ~run_surface_plot && ~run_ratio_plot
    warning('Nothing to run: both run_surface_plot and run_ratio_plot are false.');
    return;
end

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end

surface_script = fullfile(this_dir, 'plot_fitted_surface_from_export.m');
ratio_func_dir = fullfile(this_dir, 'VerifyZopt');
ratio_func_name = 'plot_spring3_cos_w1_omega_ratio';

if run_surface_plot
    if ~isfile(surface_script)
        error('Surface script not found: %s', surface_script);
    end
    fprintf('[RUN] %s\n', surface_script);

    % Keep this first because the script internally calls close all.
    % Run in a helper function workspace so clearvars in the script does not
    % clear this function's local variables.
    local_run_script_in_isolated_workspace(surface_script);
end

if run_ratio_plot
    if ~isfolder(ratio_func_dir)
        error('Function directory not found: %s', ratio_func_dir);
    end

    old_path = path;
    cleanup_obj = onCleanup(@() path(old_path)); %#ok<NASGU>
    addpath(ratio_func_dir);

    if exist(ratio_func_name, 'file') ~= 2
        error('Function not found on path: %s', ratio_func_name);
    end

    fprintf('[RUN] %s (%s)\n', ratio_func_name, ratio_func_dir);
    plot_spring3_cos_w1_omega_ratio();
end

if run_overlay_plot
    local_plot_requested_overlay(run_overlay_s1, run_overlay_s2);
end

fprintf('[DONE] run_plot_fitted_surface_and_spring3_ratio finished.\n');
end

function local_run_script_in_isolated_workspace(script_path)
run(script_path);
end

function local_plot_requested_overlay(run_overlay_s1, run_overlay_s2)
have_s2 = evalin('base', 'exist(''overlay_data_s2_gamma2_odd'', ''var'') == 1');
have_phase_mean = evalin('base', 'exist(''overlay_data_spring3_phase_mean'', ''var'') == 1');

if ~(have_s2 && have_phase_mean)
    warning(['Overlay plot skipped: required data were not found in base workspace. ', ...
        'Run both source analyses in this session before overlay.']);
    return;
end

s2 = evalin('base', 'overlay_data_s2_gamma2_odd');
pm = evalin('base', 'overlay_data_spring3_phase_mean');

if ~isstruct(s2) || ~isstruct(pm)
    warning('Overlay plot skipped: exported overlay data format is invalid.');
    return;
end

if nargin < 1 || isempty(run_overlay_s1)
    run_overlay_s1 = true;
end
if nargin < 2 || isempty(run_overlay_s2)
    run_overlay_s2 = true;
end

required_s2 = {'psi_scan_values', 'gamma1_odd_from_w', 'gamma2_odd_from_w', ...
    'gamma1_odd_from_sine', 'gamma2_odd_from_sine_pair'};
required_pm = {'xCosMean', 'yCosMean', 'xW1Mean', 'yW1Mean'};
if ~all(isfield(s2, required_s2)) || ~all(isfield(pm, required_pm))
    warning('Overlay plot skipped: exported overlay data is missing required fields.');
    return;
end

figure('Color', 'w');
hold on;
grid on;

overlay_color_sin_cos = [0.0000, 0.4470, 0.7410];
overlay_color_w_opt = [0.8500, 0.3250, 0.0980];

psi_scan = s2.psi_scan_values(:);
gamma1_opt = s2.gamma1_odd_from_w(:);
gamma2_opt = s2.gamma2_odd_from_w(:);
gamma1_sine = s2.gamma1_odd_from_sine(:);
gamma2_sine_pair = s2.gamma2_odd_from_sine_pair(:);

% Requested preprocessing for 6th-plot lines:
% cut data from max(y)-index to min(y)-index, then swap x/y.
if run_overlay_s1
    [psi_opt_cut_s1, gamma_opt_cut_s1] = local_cut_series_from_max_to_min(psi_scan, gamma1_opt);
    [psi_sine_cut_s1, gamma_sine_cut_s1] = local_cut_series_from_max_to_min(psi_scan, gamma1_sine);
end
if run_overlay_s2
    [psi_opt_cut_s2, gamma_opt_cut_s2] = local_cut_series_from_max_to_min(psi_scan, gamma2_opt);
    [psi_sine_cut_s2, gamma_sine_cut_s2] = local_cut_series_from_max_to_min(psi_scan, gamma2_sine_pair);
end

% Requested: use the 6th-plot line data with x/y swapped after the cut.
if run_overlay_s1
    plot(gamma_sine_cut_s1, psi_sine_cut_s1, '-', 'LineWidth', 1.4, ...
        'Color', overlay_color_sin_cos, ...
        'DisplayName', '$s_1, z_\sin(\theta)$');
    plot(gamma_opt_cut_s1, psi_opt_cut_s1, '-', 'LineWidth', 1.8, ...
        'Color', overlay_color_w_opt, ...
        'DisplayName', '$s_1, z_{\mathrm{opt}}(\theta)$');
end
if run_overlay_s2
    plot(gamma_sine_cut_s2, psi_sine_cut_s2, '-', 'LineWidth', 1.4, ...
        'Color', overlay_color_sin_cos, ...
        'DisplayName', '$z_{\sin}(\theta), \mathrm{predicted}$');
    plot(gamma_opt_cut_s2, psi_opt_cut_s2, '-', 'LineWidth', 1.8, ...
        'Color', overlay_color_w_opt, ...
        'DisplayName', '$z_{\mathrm{opt}}(\theta), \mathrm{predicted}$');
end

% 9th-plot mean lines as-is.
plot(pm.xCosMean(:), -pm.yCosMean(:), 'o', 'MarkerSize', 8, ...
    'Color', overlay_color_sin_cos, 'MarkerFaceColor', overlay_color_sin_cos, ...
    'DisplayName', '$z_{\sin}(\theta), \mathrm{measured}$');
plot(pm.xW1Mean(:), -pm.yW1Mean(:), 's', 'MarkerSize', 8, ...
    'Color', overlay_color_w_opt, 'MarkerFaceColor', overlay_color_w_opt, ...
    'DisplayName', '$z_{\mathrm{opt}}(\theta), \mathrm{measured}$');

xlabel('$$\Delta\omega$$','Interpreter', 'latex');
ylabel('$$\psi$$','Interpreter', 'latex');
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
hold off;

if exist('tuneFigure', 'file') == 2
    tuneFigure;
end

lgd = legend(gca, 'show', 'Location', 'eastoutside');
local_force_legend_latex(lgd);
set(lgd, 'AutoUpdate', 'off');

local_expand_legend_box_for_pdf(lgd, 1.25, 1.10);
local_force_legend_latex(lgd);

if exist('saveFigure', 'file') == 2
    saveFigure;
end

fprintf('[RUN] Overlay figure created (6th lines swapped x/y + 9th mean lines).\n');
end

function [x_cut, y_cut] = local_cut_series_from_max_to_min(x, y)
x = x(:);
y = y(:);

n = min(numel(x), numel(y));
if n == 0
    x_cut = x;
    y_cut = y;
    return;
end
x = x(1:n);
y = y(1:n);

finite_mask = isfinite(x) & isfinite(y);
if nnz(finite_mask) < 2
    x_cut = x;
    y_cut = y;
    return;
end

idx_finite = find(finite_mask);
y_finite = y(idx_finite);

[~, idx_max_local] = max(y_finite);
[~, idx_min_local] = min(y_finite);

idx_max = idx_finite(idx_max_local);
idx_min = idx_finite(idx_min_local);

step = sign(idx_min - idx_max);
if step == 0
    x_cut = x(idx_max);
    y_cut = y(idx_max);
    return;
end

idx_range = idx_max:step:idx_min;
x_cut = x(idx_range);
y_cut = y(idx_range);
end

function local_expand_legend_box_for_pdf(lgd, width_scale, height_scale)
if nargin < 2 || isempty(width_scale)
    width_scale = 1.25;
end
if nargin < 3 || isempty(height_scale)
    height_scale = 1.10;
end

if isempty(lgd) || ~isvalid(lgd)
    return;
end

drawnow;

old_units = lgd.Units;
cleanup_obj = onCleanup(@() set(lgd, 'Units', old_units)); %#ok<NASGU>

lgd.Units = 'normalized';
pos = lgd.Position;

center_x = pos(1) + pos(3)/2;
center_y = pos(2) + pos(4)/2;

pos(3) = pos(3) * width_scale;
pos(4) = pos(4) * height_scale;

pos(1) = center_x - pos(3)/2;
pos(2) = center_y - pos(4)/2;

pos(1) = max(0, min(pos(1), 1 - pos(3)));
pos(2) = max(0, min(pos(2), 1 - pos(4)));

lgd.Position = pos;

drawnow;
end

function local_force_legend_latex(lgd)
if isempty(lgd) || ~isvalid(lgd)
    return;
end

set(lgd, 'Interpreter', 'latex');

legend_texts = findall(lgd, 'Type', 'Text');
if ~isempty(legend_texts)
    set(legend_texts, 'Interpreter', 'latex');
end

drawnow;
end
