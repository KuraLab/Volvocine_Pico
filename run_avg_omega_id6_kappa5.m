% Runner script: computes and plots per-file ratio using dynamic agent IDs
% (target = max agent id in file, base = min agent id in file), between 15 and 60 seconds.

try
    % Set last argument to true to export PNG/CSV, false to skip saving
    % Dynamic like the heatmap: pass agent_id = NaN and agent_mode = 'per_file_max_vs_min'
    resultsTable = compute_avg_omega_id6_kappa5( ...
        fullfile(pwd,'SoftConnect','Spring1'), ...
        15, 60, ...
        NaN, ...            % agent_id (ignored in per_file_max_vs_min mode)
        false, ...          % export_outputs
        true, ...           % make_plot
        'per_file_max_vs_min');
    disp('--- Results (first few rows) ---');
    n = min(5, height(resultsTable));
    if n > 0
        disp(resultsTable(1:n, :));
    else
        disp('No rows to display.');
    end
catch ME
    disp(getReport(ME));
end

grid on; hold off;
xlabel('$$\omega_2/\pi$$');
ylabel('$$\bar{\dot{\phi_2}} / \bar{\dot{\phi_1}}$$');
ylim([0, 2]);
% Set x-axis ticks at 0.5 intervals over [1, 4]
xlim([1, 4]);
xticks(1:0.5:4);
tuneFigure;