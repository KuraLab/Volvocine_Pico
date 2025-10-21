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
t_start = 30; t_end = 60;
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
    % Determine agent_id dynamically: max agent id (excluding 99) across logs in this kappa folder
    [hasAgent, agent_id_kappa] = get_max_agent_id_excluding_99(folder_path);
    if ~hasAgent || isnan(agent_id_kappa)
        warning('No valid agent IDs (excluding 99) found in %s; skipping row.', folder_path);
        continue;
    end
    try
        T = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, agent_id_kappa, false, false);
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
xlabel('$$\omega_2/\pi$$'); ylabel('$$\kappa$$');
% Set x-axis ticks at 0.5 intervals over [1, 4]
xlim([1, 4]);
xticks(1:0.5:4);

set(gcf, 'Renderer', 'painters');
tuneFigure;
saveFigure;
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

function [hasAgent, max_id] = get_max_agent_id_excluding_99(folder_path)
    % Scan CSV logs under folder_path and return the maximum agent_id excluding 99.
    % hasAgent is true when at least one non-99 id is found; max_id is NaN otherwise.
    hasAgent = false; max_id = NaN;
    try
        files = dir(fullfile(folder_path, '*.csv'));
        if isempty(files), return; end
        ids_all = [];
        for iF = 1:numel(files)
            fpath = fullfile(files(iF).folder, files(iF).name);
            try
                T = readtable(fpath);
            catch
                continue;
            end
            if ismember('agent_id', T.Properties.VariableNames)
                ids = unique(T.agent_id);
                ids_all = [ids_all; ids(:)]; %#ok<AGROW>
            end
        end
        if isempty(ids_all), return; end
        ids_all = unique(ids_all);
        ids_all = ids_all(ids_all ~= 99);
        if isempty(ids_all), return; end
        max_id = max(ids_all);
        hasAgent = true;
    catch
        % leave as false/NaN
    end
end
