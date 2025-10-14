% Overlay ratio curves across multiple kappa directories under ArnoldPlot
% X-axis: 1.0:0.1:4.0 when N=31, else linspace fallback per series

baseDir = fullfile(pwd, 'ArnoldPlot');
dirs = dir(fullfile(baseDir, 'kappa*'));
dirs = dirs([dirs.isdir]);
if isempty(dirs)
    error('No kappa* directories found under %s', baseDir);
end

t_start = 15; t_end = 60; agent_id = 6;

figure('Color','w'); hold on;
colors = lines(max(1, numel(dirs)));
labels = strings(0);

for i = 1:numel(dirs)
    folder_path = fullfile(baseDir, dirs(i).name);
    try
        T = compute_avg_omega_id6_kappa5(folder_path, t_start, t_end, agent_id, false, false);
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
    labels(end+1) = string(dirs(i).name); %#ok<SAGROW>
end

grid on; hold off;
xlabel('X');
ylabel('Ratio: \omega_{agent} / \omega_{base}');
title(sprintf('Overlay ratios across kappa folders (agent %d, [%g, %g] s)', agent_id, t_start, t_end));
legend(labels, 'Interpreter', 'none', 'Location', 'bestoutside');
