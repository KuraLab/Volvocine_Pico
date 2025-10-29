% Overlay ratio curves across multiple kappa directories under ArnoldPlot
% X-axis: 1.0:0.1:4.0 when N=31, else linspace fallback per series

baseDir = fullfile(pwd, 'ArnoldPlot');
dirs = dir(fullfile(baseDir, 'kappa*'));
dirs = dirs([dirs.isdir]);
if isempty(dirs)
    error('No kappa* directories found under %s', baseDir);
end

% Sort directories by numeric kappa value in the folder name (ascending)
kvals = nan(1, numel(dirs));
for ii = 1:numel(dirs)
    % Extract the first numeric token (supports integer/decimal)
    tok = regexp(dirs(ii).name, '([-+]?\d*\.?\d+)', 'match', 'once');
    if ~isempty(tok)
        kvals(ii) = str2double(tok);
    else
        kvals(ii) = NaN;
    end
end
[~, order] = sort(kvals);  % NaNs will go to the end
dirs = dirs(order);

t_start = 30; t_end = 60;  % agent_id will be set dynamically per kappa folder

figure('Color','w'); hold on;
colors = lines(max(1, 1+numel(dirs)));
labels = strings(0);

% Add kappa0 reference curve: y = x / 2.5 over [1.0:0.1:4.0]
x_ref = 1.0:0.1:4.0;
y_ref = x_ref ./ 2.5;
%h_ref = plot(x_ref, y_ref, '-o', 'LineWidth', 1.5,'MarkerSize', 4, 'MarkerFaceColor', 'auto');
%plotHandles(end+1) = h_ref; 
%labels(end+1) = "kappa0"; 

for i = 1:numel(dirs)
    folder_path = fullfile(baseDir, dirs(i).name);
    try
        % Per-file ratio: omega(max agent) / omega(min agent)
        T = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, NaN, false, false, 'per_file_max_vs_min');
    catch ME
        warning('Failed on %s: %s', dirs(i).name, ME.message);
        continue;
    end
    ratio_all = T.ratio_to_base;
    if isempty(ratio_all)
        continue;
    end
    if numel(ratio_all) == 31
        x = 1.0:0.1:4.0;
    else
        x = linspace(1.0, 4.0, numel(ratio_all));
    end
    plot(x, ratio_all, '-o', 'Color', colors(mod(i-1,size(colors,1))+1,:), ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
    labels(end+1) = string(sprintf('%s', dirs(i).name)); %#ok<SAGROW>
end



grid on; hold off;
xlabel('$$\omega_2/\pi$$');
ylabel('$$\bar{\dot{\phi_2}} / \bar{\dot{\phi_1}}$$');
ylim([0, 2]);
% Set x-axis ticks at 0.5 intervals over [1, 4]
xlim([1, 4]);
xticks(1:0.5:4);
legend(labels, 'Interpreter', 'none', 'Location', 'bestoutside');
tuneFigure;
saveFigure;