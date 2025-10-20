% Colormap view: vertical = kappa, horizontal = X, color = ratio (omega_agent / omega_base)
% Builds a kappa-by-X matrix by scanning ArnoldPlot/kappa* and computing ratios

baseDir = fullfile(pwd, 'ArnoldPlot');

% Detect if figure windows can be shown (VS Code batch often cannot)
canShow = true;
try
    canShow = feature('ShowFigureWindows') ~= 0;
catch
    canShow = true;
end

dirs = dir(fullfile(baseDir, 'kappa*'));
dirs = dirs([dirs.isdir]);
if isempty(dirs)
    error('No kappa* directories found under %s', baseDir);
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
t_start = 20; t_end = 40; agent_id = 6;
Xgrid = 1.0:0.1:4.0;  % 31 points
NX = numel(Xgrid);
% (obsolete) NK based on dirs not used beyond this point
% Define vertical axis in integer kappa steps; rows without data remain blank (NaN)
kMin = floor(min(kvalsExistingSorted));
kMax = ceil(max(kvalsExistingSorted));
kList = kMin:1:kMax;
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
    try
        T = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, agent_id, false, false);
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
        x_src = linspace(1.0, 4.0, numel(ratio_all));
        try
            M(i, :) = interp1(x_src, ratio_all(:).', Xgrid, 'linear', 'extrap');
        catch
            % Fallback to nearest if interpolation fails
            M(i, :) = interp1(x_src, ratio_all(:).', Xgrid, 'nearest', 'extrap');
        end
    end
end

% Plot heatmap
figVis = ternary(canShow, 'on', 'off');
f = figure('Color','w', 'Visible', figVis);
% Use row indices as Y and label them with kappa to avoid imagesc coordinate ambiguity
hImg = imagesc(Xgrid, 1:NK, M);
set(hImg, 'AlphaData', ~isnan(M));  % hide NaN as transparent (blank)
ax = gca;
set(ax, 'YDir', 'normal');
yticks(1:NK);
yticklabels(arrayfun(@(v) num2str(v), kList(:), 'UniformOutput', false));
colormap(parula);
colorbar;
xlabel('X'); ylabel('kappa');
title(sprintf('Ratio heatmap across kappa folders (agent %d, [%g, %g] s)', agent_id, t_start, t_end));

set(gcf, 'Renderer', 'painters');
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
