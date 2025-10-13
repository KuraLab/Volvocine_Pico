% Runner script: computes and plots mean angular velocity for agent 6
% in ArnoldPlot/kappa5 between 10 and 50 seconds.

try
    % Set last argument to true to export PNG/CSV, false to skip saving
    resultsTable = compute_avg_omega_id6_kappa5(fullfile(pwd,'ArnoldPlot','kappa5'), 20, 60, 6, false);
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
