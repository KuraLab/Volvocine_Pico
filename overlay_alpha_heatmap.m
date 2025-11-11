% Colormap view: vertical = kappa, horizontal = X, color = ratio (omega_agent / omega_base)
% Builds a kappa-by-X matrix by scanning OptimumAlpha/kappa* and computing ratios

%baseDir = fullfile(pwd, 'OptimumAlpha/Spring3');
baseDir = fullfile(pwd, 'SoftConnect');

% Detect if figure windows can be shown (VS Code batch often cannot)
canShow = true;
try
    canShow = feature('ShowFigureWindows') ~= 0;
catch
    canShow = true;
end

dirs = dir(fullfile(baseDir, 'Spring*'));
dirs = dirs([dirs.isdir]);
if isempty(dirs)
    error('No alpha* directories found under %s', baseDir);
end

% Sort directories by numeric kappa value (ascending)
kvalsExisting = nan(1, numel(dirs));
for ii = 1:numel(dirs)
    tok = regexp(dirs(ii).name, '([-+]?\d*\.?\d+)', 'match', 'once');
    if ~isempty(tok), kvalsExisting(ii) = str2double(tok); end
end
[kvalsExistingSorted, order] = sort(kvalsExisting);
dirs = dirs(order);

% Parameters
t_start = 20; t_end = 60;
% Horizontal axis (X) parameters — change these to adjust range/resolution easily
Xmin = 1.0;    % start of X range
Xmax = 4.0;    % end of X range
dX   = 0.1;    % X step
xTickStep = 0.5; % tick step along X axis
Xgrid = Xmin:dX:Xmax;  % grid along X
NX = numel(Xgrid);
% (obsolete) NK based on dirs not used beyond this point
% Vertical axis step for kappa (alpha) rows; adjust as needed (e.g., 0.5, 0.2, 2)
kStep = 1;  % CHANGE HERE to control vertical increment
% Display scale for vertical axis labels (e.g., 0.1 to show 1/10 values)
kDisplayScale = 1;  % CHANGE HERE to scale label numbers only
% Define vertical axis range based only on existing data (do not force include 0)
kMin = floor(min(kvalsExistingSorted));
kMax = ceil(max(kvalsExistingSorted));
kList = kMin:kStep:kMax;
NK = numel(kList);
M = NaN(NK, NX);
% Build matrix row by row
for i = 1:NK
    kTarget = kList(i);
    % Find matching existing dir for this kappa (exact match within tolerance)
    idx = find(abs(kvalsExistingSorted - kTarget) < 1e-9, 1, 'first');
    if isempty(idx)
        % No data yet for this kappa: leave row as NaN (blank)
        continue;
    end
    folder_path = fullfile(baseDir, dirs(idx).name);
    % Compute per-file ratio: omega(max agent in file) / omega(min agent in file)
    try
        T = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, NaN, false, false, 'per_file_max_vs_min');
    catch ME
        warning('Failed on %s: %s', folder_path, ME.message);
        continue;
    end
    if isempty(T)
        continue;
    end
    ratio_all = T.ratio_to_base;
    if isempty(ratio_all)
        continue;
    end
    if numel(ratio_all) == NX
        M(i, :) = ratio_all(:).';
    else
        % Map existing indices to X via linspace and interpolate onto Xgrid
        x_src = linspace(Xmin, Xmax, numel(ratio_all));
        try
            M(i, :) = interp1(x_src, ratio_all(:).', Xgrid, 'linear', 'extrap');
        catch
            % Fallback to nearest if interpolation fails
            M(i, :) = interp1(x_src, ratio_all(:).', Xgrid, 'nearest', 'extrap');
        end
    end
end

figVis = ternary(canShow, 'on', 'off');
f = figure('Color','w', 'Visible', figVis);
% Use row indices as Y and label them with kappa to avoid imagesc coordinate ambiguity
hImg = imagesc(Xgrid, 1:NK, M);
set(hImg, 'AlphaData', ~isnan(M));  % hide NaN as transparent (blank)
ax = gca;
set(ax, 'YDir', 'normal');
yticks(1:NK);
% Show vertical labels scaled by kDisplayScale (labels only; data mapping unchanged)
yticklabels(arrayfun(@(v) num2str(v * kDisplayScale), kList(:), 'UniformOutput', false));
% Compute symmetric color limits around 1 (baseline) for ratio values
vals = M(~isnan(M));
if isempty(vals)
    clim = [0.5 1.5];
else
    delta = max(abs(vals - 1));
    if ~isfinite(delta) || delta == 0
        delta = 0.25; % sensible default spread
    end
    deltaCap = 0.75; % cap extreme spreads to keep contrast
    delta = min(delta, deltaCap);
    clim = [1 - delta, 1 + delta];
end
% Use a colorblind-friendly diverging colormap centered at white
colormap(diverging_blue_white_red(256));
caxis(clim);
cb = colorbar;
cb.TickDirection = 'out';
% Neat colorbar ticks (include baseline 1.0 if within range)
[ticks, fmt] = compute_nice_ticks(clim, 5);
cb.Ticks = ticks;
cb.TickLabels = arrayfun(@(v) sprintf(fmt, v), ticks, 'UniformOutput', false);
% Apply label last to avoid any render quirks
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\bar{\dot{\phi_2}} / \bar{\dot{\phi_1}}$';
try
    cb.Label.FontSize = ax.FontSize;
catch
end
% Show missing data as light gray
set(ax, 'Color', [0.95 0.95 0.95]);
hold on;

xlabel('$$\omega_2/\pi$$'); ylabel('$$\tau/\pi$$');
% Set x-axis ticks using configured range and step
xlim([Xmin-dX/2, Xmax+dX/2]);
%xticks(Xmin-dX/2:xTickStep:Xmax);

set(gcf, 'Renderer', 'painters');
tuneFigure;
%saveFigure;
drawnow;
if canShow
    try, shg; catch, end
else
    % Auto-export when running in environments that cannot show figures
    outDir = fullfile(pwd, 'exports');
    if ~isfolder(outDir)
        try, mkdir(outDir); catch, end
    end
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    outPng = fullfile(outDir, sprintf('overlay_kappa_heatmap_%s.png', ts));
    try
        exportgraphics(f, outPng, 'Resolution', 200);
    catch
        try, saveas(f, outPng); catch, end
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end

% per-file agent selection is handled inside compute_avg_omega_id6_kappa5

function cmap = diverging_blue_white_red(n)
%DIVERGING_BLUE_WHITE_RED Simple blue-white-red diverging colormap.
%   cmap = diverging_blue_white_red(n) returns an n-by-3 colormap that
%   transitions from blue to white to red, suitable for data centered at 0
%   (or 1 after shifting). Designed to be reasonably colorblind-friendly.
    if nargin < 1 || isempty(n)
        n = 256;
    end
    n = max(3, round(n));
    % Create two gradients: blue->white and white->red
    n2 = ceil(n/2);
    up = [linspace(0,1,n2)', linspace(0,1,n2)', ones(n2,1)];   % blue -> white
    dn = [ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)'];   % white -> red (reverse green/blue)
    % Combine and ensure exact center white
    cmap = [up; dn];
    % Trim to exactly n rows
    if size(cmap,1) > n
        cmap = cmap(1:n, :);
    elseif size(cmap,1) < n
        cmap = interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,n), 'linear');
    end
    % Gamma correction for better perceptual balance
    gamma = 0.9;
    cmap = max(min(cmap.^gamma, 1), 0);
end

function [ticks, fmt] = compute_nice_ticks(clim, nTarget)
% Compute human-friendly tick positions within clim with roughly nTarget ticks.
% Ensures that 1.0 (baseline) is included if it lies in the range.
    if nargin < 2 || isempty(nTarget), nTarget = 5; end
    lo = clim(1); hi = clim(2);
    if ~isfinite(lo) || ~isfinite(hi) || lo >= hi
        ticks = linspace(0,1,nTarget);
        fmt = '%.2f';
        return;
    end
    span = hi - lo;
    rawStep = span / max(1, (nTarget-1));
    mag = 10.^floor(log10(rawStep));
    stepCand = [1, 2, 2.5, 5, 10] .* mag; % 1-2-2.5-5 series
    [~, idx] = min(abs(stepCand - rawStep));
    step = stepCand(idx);
    vmin = floor(lo/step)*step;
    vmax = ceil(hi/step)*step;
    ticks = vmin:step:vmax;
    % Clip within clim
    ticks = ticks(ticks >= lo - 1e-12 & ticks <= hi + 1e-12);
    % Ensure baseline 1 exists if inside range
    if 1 >= lo - 1e-12 && 1 <= hi + 1e-12
        if all(abs(ticks - 1) > 1e-9)
            ticks = sort([ticks 1]);
        end
    end
    % Choose format
    if step >= 1
        fmt = '%.0f';
    elseif step >= 0.1
        fmt = '%.1f';
    elseif step >= 0.01
        fmt = '%.2f';
    else
        fmt = '%.3f';
    end
end