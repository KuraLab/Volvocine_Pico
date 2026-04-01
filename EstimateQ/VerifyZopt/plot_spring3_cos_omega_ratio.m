% Plot omega ratio curves for both:
%   EstimateQ/VerifyZopt/Spring3/cos/**/merged_*.csv
%   EstimateQ/VerifyZopt/Spring3/w1/**/merged_*.csv
% Reuses compute_avg_omega_id6_kappa5 in per-file max-vs-min mode.

t_start = 35;
t_end = 90;

rootDirCos = fullfile(pwd, 'EstimateQ', 'VerifyZopt', 'Spring3', 'cos');
rootDirW1  = fullfile(pwd, 'EstimateQ', 'VerifyZopt', 'Spring3', 'w1');

[xCos, yCos] = local_collect_ratio_curve(rootDirCos, t_start, t_end);
[xW1, yW1]   = local_collect_ratio_curve(rootDirW1, t_start, t_end);

figure('Color', 'w');
plot(xCos, yCos, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'auto', 'DisplayName', 'Spring3/cos');
hold on;
plot(xW1, yW1, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'auto', 'DisplayName', 'Spring3/w1');
grid on;
xlim([min([xCos(:); xW1(:)]), max([xCos(:); xW1(:)])]);
xlabel('$$\Delta\omega$$');
ylabel('$$\bar{\dot{\phi_2}} / \bar{\dot{\phi_1}}$$');
legend('Location', 'best');
hold off;

if exist('tuneFigure', 'file') == 2
    tuneFigure;
end


function [x, y] = local_collect_ratio_curve(rootDir, t_start, t_end)
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

    r = T.ratio_to_base;
    r = r(~isnan(r));
    if isempty(r)
        continue;
    end

    % Use mean in case multiple CSV files exist in a folder.
    y(i) = mean(r);
end
end
