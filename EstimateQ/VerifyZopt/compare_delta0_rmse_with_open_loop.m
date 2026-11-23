function out = compare_delta0_rmse_with_open_loop(varargin)
% COMPARE_DELTA0_RMSE_WITH_OPEN_LOOP Match the manuscript's s2 RMSE table.
%   out = compare_delta0_rmse_with_open_loop('ShowFigures', false);
%
% Reproduces the six open-loop table values with the 240_0531 reference fit
% (M=N=10, log time 10--90 s), then evaluates the SAME model on closed-loop
% all numeric frequency folders in cos and w1 (k=1) by default. High-gain
% trials are opt-in through Conditions; HighGainDeltaOmegaOverPi=[] selects
% all numeric frequency folders, or a vector selects specific detunings.
% Default closed-loop evaluation: feedback time 0--115 s
% = original common-overlap log time 5--120 s, sampled every 0.01 s.
% The endpoint excludes the trial-end signal transition, which the user
% attributes to likely removal of the robots from the water.
% This is a comparison of normalized a2 (not exact unlogged controller flex).
% Each file/agent uses its full-log p10/p90 and clipping to [-0.5,0.5],
% exactly as in plot_spring3_a2_reference_fit_error. No closed-loop fitting.
% Table sign convention here: Delta omega = omega1 - omega2.
%
% A saved 240 model is accepted only if its s2 coefficients/mean match a
% fresh fit made by the original table code. Previous 255-model results are
% preserved in their original directory. Fixed-scale signals are not plotted.
% The comparison includes pre-feedback 1--4 s and 10 s bins as descriptive
% checks, not independent trials or a causal decomposition of error sources.
%
% Options: FeedbackDurationSec (115), Conditions (cos,w1),
% HighGainDeltaOmegaOverPi, OutputDir, ShowFigures (true),
% PlotTimeSeries (false), SaveTimeSeries (false), SaveZoomPanels (true),
% ShowReconstructedRate (true).
% High-gain trials at Delta omega=0.2*pi,0.3*pi,0.5*pi also export eight Module 1/2
% zoom overlays per condition into closed_loop_timeseries/overlaid_zoom,
% including with PlotTimeSeries=false or ReplotOnly=true.
% ReplotOnly=true redraws summary figures from OutputDir/rmse_comparison.mat
% without rerunning the log analysis. Detuning axes use rad/s, not Delta/pi.
% Outputs: combined CSV/LaTeX table, per-window/per-bin RMSE, matched-model
% phase/sensory time series and figures, MAT with provenance and coefficients.

    here = fileparts(mfilename('fullpath'));
    estimate_dir = fileparts(here);
    addpath(here, estimate_dir);
    p = inputParser;
    addParameter(p, 'FeedbackDurationSec', 115, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 90);
    addParameter(p, 'Conditions', {'cos','w1'}, @(x) iscellstr(x) && ~isempty(x) && all(ismember(x, {'cos','w1','k3w2','k5w2'})) && numel(unique(x)) == numel(x));
    addParameter(p, 'HighGainDeltaOmegaOverPi', [], ...
        @(x) isnumeric(x) && (isempty(x) || (isvector(x) && all(isfinite(x)) && numel(unique(x)) == numel(x))));
    addParameter(p, 'OutputDir', fullfile(here, 'k1_all_detunings_rmse_240'));
    addParameter(p, 'ShowFigures', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'PlotTimeSeries', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SaveTimeSeries', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SaveZoomPanels', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowReconstructedRate', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ReplotOnly', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    cfg = p.Results;
    if ~isfolder(cfg.OutputDir), mkdir(cfg.OutputDir); end
    if cfg.ReplotOnly
        saved = load(fullfile(cfg.OutputDir, 'rmse_comparison.mat'), 'out');
        out = saved.out;
        plot_cfg = out.config;
        plot_cfg.OutputDir = cfg.OutputDir;
        plot_cfg.ShowFigures = cfg.ShowFigures;
        plot_comparison(out.comparison, out.closed_loop_10s_bins, out.reference_fit.rmse, plot_cfg);
        if any(out.comparison.gain > 1)
            plot_detuned_windows(out.closed_loop_10s_bins, out.reference_fit.rmse, plot_cfg);
        end
        plot_closed_loop_detunings(out.comparison, out.reference_fit.rmse, plot_cfg);
        if cfg.SaveZoomPanels && any(out.comparison.gain > 1)
            evaluate_delta0_closed_loop('ReplotOnly', true, ...
                'OutputDir', fullfile(cfg.OutputDir, 'closed_loop_timeseries'), ...
                'Conditions', out.config.Conditions, 'PlotTimeSeries', false, ...
                'ShowFigures', cfg.ShowFigures, 'SaveZoomPanels', cfg.SaveZoomPanels, ...
                'ShowReconstructedRate', cfg.ShowReconstructedRate);
        end
        fprintf('Updated summary figures: %s\n', cfg.OutputDir);
        return;
    end

    % Some fitting helpers can create figures even with ShowFigures=false.
    previous_visibility = get(groot, 'DefaultFigureVisible');
    previous_figures = findall(groot, 'Type', 'figure');
    restore_figures = onCleanup(@() restore_fit_figures(previous_visibility, previous_figures));
    set(groot, 'DefaultFigureVisible', 'off');
    open_loop = plot_spring3_a2_reference_fit_error('ReferenceFolder', '240_0531', ...
        'M', 10, 'N', 10, 'ShowFigures', false);
    clear restore_figures;

    model_file = fullfile(estimate_dir, 'Spring3', '240_0531', 'gamma_exports', 'gamma_export_latest.mat');
    S = load(model_file, 'gamma_export');
    G = S.gamma_export;
    ix = find([G.agents.agent_id] == open_loop.z_agent_id, 1);
    saved_fit = G.agents(ix).a2_gamma;
    fit = open_loop.reference_fit;
    assert(isequal(G.phase_agent_ids(:), open_loop.phase_agent_ids(:)), 'Model agent order differs.');
    assert(saved_fit.M == fit.M && saved_fit.N == fit.N, 'Model Fourier orders differ.');
    assert(norm(saved_fit.coeff-fit.coeff) < 1e-8*max(1,norm(fit.coeff)) && ...
        abs(saved_fit.z_mean-fit.z_mean) < 1e-10, 'Saved model differs from manuscript reference fit.');

    folders = [220; 230; 240; 260; 270; 280];
    reported = [0.0471; 0.0367; 0.0230; 0.0556; 0.0793; 0.0994];
    reproduction = table(folders, reported, nan(6,1), 'VariableNames', {'folder_value','reported_rmse','reproduced_rmse'});
    rows = struct([]);
    for k = 1:numel(folders)
        folder_name = sprintf('%d_0531', folders(k));
        r = open_loop.summary_table(strcmp(open_loop.summary_table.folder_name, folder_name), :);
        assert(height(r) == 1, 'Expected one summary row for %s.', folder_name);
        reproduction.reproduced_rmse(k) = r.rmse;
        row = struct('operation', 'open_loop', 'condition', folder_name, ...
            'gain', 0, 'frequency_folder', folders(k), ...
            'delta_omega_over_pi', 2.5-folders(k)/100, ...
            'omega1_over_pi', 2.5, 'omega2_over_pi', folders(k)/100, ...
            'start_log_s', 10, 'end_log_s', 90, 'duration_s', 80, ...
            'n_samples', r.n_points, 'rmse', r.rmse, 'mae', r.mae, 'bias', r.bias);
        rows = [rows; row]; %#ok<AGROW>
    end
    assert(all(round(reproduction.reproduced_rmse,4) == reported), 'The supplied manuscript RMSE table was not reproduced.');

    onset = 5;
    trial_folders = cell(size(cfg.Conditions));
    high_gain_folders = 250-100*cfg.HighGainDeltaOmegaOverPi;
    assert(all(high_gain_folders > 0 & abs(high_gain_folders-round(high_gain_folders)) < 1e-8), ...
        'Detunings must correspond to positive integer folder names.');
    for k = 1:numel(cfg.Conditions)
        if condition_gain(cfg.Conditions{k}) > 1 && ~isempty(high_gain_folders)
            trial_folders{k} = round(high_gain_folders);
        else
            listing = dir(fullfile(here, 'Spring3', cfg.Conditions{k}));
            listing = listing([listing.isdir]);
            values = str2double({listing.name});
            trial_folders{k} = sort(values(isfinite(values)));
            assert(~isempty(trial_folders{k}), 'No numeric frequency folders for %s.', cfg.Conditions{k});
        end
    end
    closed_loop = evaluate_delta0_closed_loop('ModelFile', model_file, ...
        'AnalysisWindow', [onset onset+cfg.FeedbackDurationSec], ...
        'FeedbackStartSec', onset, 'ShowFigures', cfg.ShowFigures, ...
        'Conditions', cfg.Conditions, ...
        'FrequencyFolders', trial_folders, ...
        'PlotTimeSeries', cfg.PlotTimeSeries, 'SaveTimeSeries', cfg.SaveTimeSeries, ...
        'SaveZoomPanels', cfg.SaveZoomPanels, ...
        'ShowReconstructedRate', cfg.ShowReconstructedRate, ...
        'OutputDir', fullfile(cfg.OutputDir, 'closed_loop_timeseries'));
    windows = struct([]);
    bins = struct([]);
    for k = 1:numel(closed_loop.runs)
        run = closed_loop.runs(k);
        T = run.timeseries;
        fb_time = T.time_s-onset;
        error = T.s_error_file_2;
        main = fb_time >= 0 & fb_time <= cfg.FeedbackDurationSec & isfinite(error);
        assert(nnz(main) == round(cfg.FeedbackDurationSec/closed_loop.config.SampleDt)+1, ...
            'Not all requested closed-loop sensory samples are valid.');
        row = struct('operation', 'closed_loop', 'condition', run.condition, ...
            'gain', condition_gain(run.condition), 'frequency_folder', run.frequency_folder, ...
            'delta_omega_over_pi', run.delta_omega_over_pi, 'omega1_over_pi', 2.5, 'omega2_over_pi', run.frequency_folder/100, ...
            'start_log_s', onset, 'end_log_s', onset+cfg.FeedbackDurationSec, ...
            'duration_s', cfg.FeedbackDurationSec, 'n_samples', nnz(main), ...
            'rmse', sqrt(mean(error(main).^2)), 'mae', mean(abs(error(main))), 'bias', mean(error(main)));
        rows = [rows; row]; %#ok<AGROW>
        % 0--115 also records the alternative endpoint of original log 120 s.
        intervals = [1 4; onset onset+cfg.FeedbackDurationSec; onset 120; onset+10 onset+90; onset+60 onset+cfg.FeedbackDurationSec];
        names = {'before_feedback'; 'full_feedback'; 'through_original_log_120s'; 'feedback_10_to_90s'; 'feedback_60_to_end'};
        covered = intervals(:,2) <= onset+cfg.FeedbackDurationSec;
        intervals = intervals(covered,:);
        names = names(covered);
        for w = 1:size(intervals,1)
            mask = T.time_s >= intervals(w,1) & T.time_s <= intervals(w,2) & isfinite(error);
            wr = struct('condition', run.condition, 'window', names{w}, ...
                'frequency_folder', run.frequency_folder, 'delta_omega_over_pi', run.delta_omega_over_pi, ...
                'start_feedback_s', intervals(w,1)-onset, 'end_feedback_s', intervals(w,2)-onset, ...
                'n_samples', nnz(mask), 'rmse', sqrt(mean(error(mask).^2)), ...
                'mae', mean(abs(error(mask))), 'bias', mean(error(mask)));
            windows = [windows; wr]; %#ok<AGROW>
        end
        edges = unique([0:10:cfg.FeedbackDurationSec, cfg.FeedbackDurationSec]);
        for w = 1:numel(edges)-1
            mask = fb_time >= edges(w) & fb_time < edges(w+1) & isfinite(error);
            if w == numel(edges)-1, mask = mask | (fb_time == edges(end) & isfinite(error)); end
            br = struct('condition', run.condition, 'start_feedback_s', edges(w), ...
                'frequency_folder', run.frequency_folder, 'delta_omega_over_pi', run.delta_omega_over_pi, ...
                'end_feedback_s', edges(w+1), 'n_samples', nnz(mask), ...
                'rmse', sqrt(mean(error(mask).^2)));
            bins = [bins; br]; %#ok<AGROW>
        end
    end
    comparison = struct2table(rows);
    comparison.rmse_over_reference_training = comparison.rmse/fit.rmse;
    window_table = struct2table(windows);
    bin_table = struct2table(bins);
    out = struct('config', cfg, 'reference_model_file', model_file, ...
        'reference_fit', fit, 'table_reproduction', reproduction, ...
        'comparison', comparison, 'closed_loop_windows', window_table, ...
        'closed_loop_10s_bins', bin_table, 'closed_loop', closed_loop);
    out.notes = { ...
        'Delta omega = omega1 - omega2, matching the manuscript table (opposite the legacy plot labels).', ...
        'Same open-loop reference model and per-file normalization for every RMSE.', ...
        sprintf('Open-loop table windows are 80 s; primary closed-loop window is %.0f s.', cfg.FeedbackDurationSec), ...
        'Default cutoff is feedback time 115 s (log time 120 s), excluding the trial-end signal transition; user reports likely removal from water.', ...
        'Before-feedback check is only 3 s; phase-pair sampling also changes when feedback is active.', ...
        'RMSE differences alone do not identify a causal contribution of natural detuning versus feedback.'};
    writetable(comparison, fullfile(cfg.OutputDir, 'rmse_comparison.csv'));
    writetable(comparison(strcmp(comparison.operation,'closed_loop'),:), fullfile(cfg.OutputDir, 'closed_loop_rmse.csv'));
    writetable(reproduction, fullfile(cfg.OutputDir, 'open_loop_table_reproduction.csv'));
    writetable(window_table, fullfile(cfg.OutputDir, 'closed_loop_window_rmse.csv'));
    writetable(bin_table, fullfile(cfg.OutputDir, 'closed_loop_10s_rmse.csv'));
    save(fullfile(cfg.OutputDir, 'rmse_comparison.mat'), 'out');
    write_latex_table(comparison, fullfile(cfg.OutputDir, 'rmse_comparison.tex'));
    plot_comparison(comparison, bin_table, fit.rmse, cfg);
    if any(comparison.gain > 1)
        plot_detuned_windows(bin_table, fit.rmse, cfg);
    end
    plot_closed_loop_detunings(comparison, fit.rmse, cfg);
    fprintf('\nMatched-model s2 comparison (feedback time 0--%.0f s):\n', cfg.FeedbackDurationSec);
    disp(comparison(:, {'operation','condition','delta_omega_over_pi','n_samples','rmse','rmse_over_reference_training'}));
    fprintf('Additional windows:\n');
    disp(window_table);
    fprintf('Saved to: %s\n', cfg.OutputDir);
end

function plot_closed_loop_detunings(T, training_rmse, cfg)
    groups = {{'cos','w1'}, {'k3w2','k5w2'}};
    stems = {'k1_rmse_by_detuning', 'k3_k5_rmse_by_detuning'};
    titles = {'k=1', 'k=3,5'};
    for g = 1:numel(groups)
    conditions = groups{g};
    C = T(strcmp(T.operation,'closed_loop') & ismember(T.condition,conditions),:);
    if isempty(C), continue; end
    folders = sort(unique(C.frequency_folder));
    P = table(folders, (250-folders)/100, repmat(2.5,size(folders)), folders/100, ...
        'VariableNames', {'frequency_folder','delta_omega_over_pi','omega1_over_pi','omega2_over_pi'});
    visibility = 'off';
    if cfg.ShowFigures, visibility = 'on'; end
    f = figure('Color','w','Visible',visibility,'Position',[100 100 1000 650]);
    ax = axes(f);
    hold(ax,'on');
    for condition = conditions
        subset = C(strcmp(C.condition,condition{1}),:);
        if isempty(subset), continue; end
        values = nan(size(folders));
        for k = 1:numel(folders)
            mask = subset.frequency_folder==folders(k);
            assert(nnz(mask)==1, 'Expected one trial per condition and frequency folder.');
            values(k) = subset.rmse(mask);
        end
        P.(['rmse_' condition{1}]) = values;
        [x, order] = sort(P.delta_omega_over_pi*pi);
        plot(ax,x,values(order),'-o','DisplayName',condition_label(condition{1}));
    end
    yline(ax,training_rmse,'k--','Reference training RMSE','HandleVisibility','off');
    xlabel(ax,'$\Delta\omega\;(\mathrm{rad}\,\mathrm{s}^{-1})$','Interpreter','latex');
    ylabel(ax,'$\mathrm{RMSE}(s_2^{\mathrm{meas}}-\hat{s}_2)$','Interpreter','latex');
    title(ax,sprintf('%s, feedback time 0--%.0f s',titles{g},cfg.FeedbackDurationSec));
    ylim(ax,[0 Inf]);
    grid(ax,'on');
    legend(ax,'Location','best');
    tune_analysis_figure(f);
    writetable(P,fullfile(cfg.OutputDir,[stems{g} '.csv']));
    exportgraphics(f,fullfile(cfg.OutputDir,[stems{g} '.png']),'Resolution',160);
    exportgraphics(f,fullfile(cfg.OutputDir,[stems{g} '.pdf']),'ContentType','vector');
    if ~cfg.ShowFigures, close(f); end
    end
end

function restore_fit_figures(visibility, previous)
    set(groot, 'DefaultFigureVisible', visibility);
    current = findall(groot, 'Type', 'figure');
    created = current(~ismember(current, previous));
    close(created);
end

function gain = condition_gain(condition)
    switch condition
        case {'cos','w1'}, gain = 1;
        case 'k3w2', gain = 3;
        case 'k5w2', gain = 5;
        otherwise, error('Unknown gain for condition %s.', condition);
    end
end

function label = condition_label(condition)
    switch condition
        case 'cos', label = 'k=1 (cos)';
        case 'w1', label = 'k=1 (w1)';
        case 'k3w2', label = 'k=3 (w2)';
        case 'k5w2', label = 'k=5 (w2)';
        otherwise, label = condition;
    end
end

function write_latex_table(T, filename)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Cannot write %s.', filename);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%% Same 240_0531 s2 model; normalized a2; Delta omega = omega1 - omega2.\n');
    fprintf(fid, '%% Open-loop: log 10--90 s; closed-loop: see duration column.\n');
    fprintf(fid, '\\begin{tabular}{lccccc}\n\\toprule\n');
    fprintf(fid, 'Condition & $\\Delta\\omega$ & $\\omega_1$ & $\\omega_2$ & Duration (s) & RMSE \\\\\n\\midrule\n');
    for k = 1:height(T)
        condition = T.operation{k};
        if strcmp(condition, 'closed_loop'), condition = ['Closed-loop ' condition_label(T.condition{k})]; else, condition = 'Open-loop'; end
        fprintf(fid, '%s & $%.2f\\pi$ & $%.2f\\pi$ & $%.2f\\pi$ & %.0f & %.4f \\\\\n', ...
            condition, T.delta_omega_over_pi(k), T.omega1_over_pi(k), T.omega2_over_pi(k), T.duration_s(k), T.rmse(k));
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
end

function plot_comparison(T, bins, training_rmse, cfg)
    visibility = 'off';
    if cfg.ShowFigures, visibility = 'on'; end
    f = figure('Color', 'w', 'Visible', visibility, 'Position', [100 100 1500 550]);
    layout = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(layout);
    hold(ax, 'on');
    palette = [0 .35 .7; .85 .25 .1; .15 .6 .25; .55 .2 .65];
    closed_rows = find(strcmp(T.operation, 'closed_loop'));
    conditions = unique(T.condition(closed_rows), 'stable');
    open = T(strcmp(T.operation, 'open_loop'), :);
    [x, order] = sort(open.delta_omega_over_pi*pi);
    plot(ax, x, open.rmse(order), '-s', 'Color', [.4 .4 .4], 'DisplayName', 'Open loop');
    for k = 1:numel(conditions)
        mask = strcmp(T.operation, 'closed_loop') & strcmp(T.condition, conditions{k});
        subset = T(mask,:);
        [x, order] = sort(subset.delta_omega_over_pi*pi);
        plot(ax, x, subset.rmse(order), '-o', 'Color', palette(k,:), ...
            'DisplayName', condition_label(conditions{k}));
    end
    set(ax, 'TickLabelInterpreter', 'latex');
    ylabel(ax, '$\mathrm{RMSE}(s_2^{\mathrm{meas}}-\hat{s}_2)$', 'Interpreter', 'latex');
    xlabel(ax, '$\Delta\omega\;(\mathrm{rad}\,\mathrm{s}^{-1})$', 'Interpreter', 'latex');
    title(ax, 'Same model: identified at \Delta\omega=0.1\pi', 'Interpreter', 'tex');
    yline(ax, training_rmse, 'k--', 'Reference training RMSE', 'HandleVisibility', 'off');
    ylim(ax, [0 1.15*max(T.rmse)]);
    legend(ax, 'Location', 'best');
    grid(ax, 'on');
    ax = nexttile(layout);
    hold(ax, 'on');
    for k = 1:numel(conditions)
        subset = bins(strcmp(bins.condition, conditions{k}) & abs(bins.delta_omega_over_pi) < 1e-8, :);
        if isempty(subset), continue; end
        plot(ax, (subset.start_feedback_s+subset.end_feedback_s)/2, subset.rmse, '-o', ...
            'Color', palette(k,:), 'DisplayName', condition_label(conditions{k}));
    end
    yline(ax, training_rmse, 'k--', 'Reference training RMSE', 'HandleVisibility', 'off');
    xlabel(ax, 'Time since feedback onset $t$ (s)', 'Interpreter', 'latex');
    ylabel(ax, '$\mathrm{RMSE}(s_2^{\mathrm{meas}}-\hat{s}_2)$', 'Interpreter', 'latex');
    title(ax, 'Windowed RMSE at \Delta\omega=0', 'Interpreter', 'tex');
    xlim(ax, [0 cfg.FeedbackDurationSec]);
    ylim(ax, [0 Inf]);
    legend(ax, 'Location', 'best');
    grid(ax, 'on');
    tune_analysis_figure(f);
    exportgraphics(f, fullfile(cfg.OutputDir, 'rmse_comparison.png'), 'Resolution', 160);
    exportgraphics(f, fullfile(cfg.OutputDir, 'rmse_comparison.pdf'), 'ContentType', 'vector');
    if ~cfg.ShowFigures, close(f); end
end

function plot_detuned_windows(bins, training_rmse, cfg)
    visibility = 'off';
    if cfg.ShowFigures, visibility = 'on'; end
    f = figure('Color', 'w', 'Visible', visibility, 'Position', [100 100 1400 550]);
    layout = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    for condition = {'k3w2','k5w2'}
        ax = nexttile(layout);
        hold(ax, 'on');
        subset = bins(strcmp(bins.condition, condition{1}), :);
        deltas = sort(unique(subset.delta_omega_over_pi));
        colors = lines(numel(deltas));
        for k = 1:numel(deltas)
            rows = subset(abs(subset.delta_omega_over_pi-deltas(k)) < 1e-8,:);
            plot(ax, (rows.start_feedback_s+rows.end_feedback_s)/2, rows.rmse, '-o', ...
                'Color', colors(k,:), 'DisplayName', sprintf('$\\Delta\\omega=%+.1f\\pi$', deltas(k)));
        end
        yline(ax, training_rmse, 'k--', 'HandleVisibility', 'off');
        title(ax, condition_label(condition{1}));
        xlabel(ax, 'Time since feedback onset $t$ (s)', 'Interpreter', 'latex');
        ylabel(ax, '$\mathrm{RMSE}(s_2^{\mathrm{meas}}-\hat{s}_2)$', 'Interpreter', 'latex');
        xlim(ax, [0 cfg.FeedbackDurationSec]);
        ylim(ax, [0 Inf]);
        grid(ax, 'on');
        if ~isempty(deltas), legend(ax, 'Interpreter', 'latex', 'Location', 'best'); end
    end
    tune_analysis_figure(f);
    exportgraphics(f, fullfile(cfg.OutputDir, 'rmse_detuned_windows.png'), 'Resolution', 160);
    exportgraphics(f, fullfile(cfg.OutputDir, 'rmse_detuned_windows.pdf'), 'ContentType', 'vector');
    if ~cfg.ShowFigures, close(f); end
end
