function out = evaluate_delta0_closed_loop(varargin)
% EVALUATE_DELTA0_CLOSED_LOOP Phase rates and sensory-model errors in closed loop.
% Run from any directory after adding this function's directory to the path:
%   out = evaluate_delta0_closed_loop();
%   out = evaluate_delta0_closed_loop('ShowFigures', false);
%
% Uses every CSV in Spring3/{cos,w1}/250, without convergence selection.
% The saved OPEN-LOOP Spring3/255 model supplies both s1 and s2; there is no
% closed-loop refitting, time shifting, or error-minimizing rescaling.
% Summary window: 60--120 s; baseline: 1--4 s, before feedback starts at 5 s.
% Time zero is the start of common agent overlap (as in the existing plots).
% Figures instead use feedback onset (5 s by default) as time zero and show
% only post-onset data. CSV times and analysis windows retain log time.
%
% Important conventions:
% * April firmware encodes a0=floor(mod(phi,2*pi)*255/(2*pi)). Derivatives
%   therefore use a 255-code period. The shared loader's 256 unwrap is
%   undone before unwrapping with 255. Model INPUTS retain the 256 convention
%   used to identify the saved Fourier coefficients. These are separate.
% * Natural frequency of module 1 was sent as 7.85 rad/s (3.14*2.5 rounded
%   to two decimals); module 2 depends on its frequency folder. Estimate
%   each independently from the baseline as well.
% * a2 is compared with the model under two explicit normalizations:
%   (a) each closed-loop file's full-log p10/p90, clipped to [-.5,.5], as in
%       plot_spring3_a2_reference_fit_error (shape-normalized comparison);
%   (b) the OPEN-LOOP raw p10/p90, fixed across conditions, without clipping
%       (also report a clipped version). This retains amplitude/offset error.
% * A centered local linear phase slope suppresses phase-code quantization.
%   Results at 0.05/0.10/0.20 s are saved to quantify window sensitivity.
% * Samples across acquisition gaps > MaxGapSec are not interpolated.
% * The current saved model may have been exported after the April trials;
%   its path, export date, coefficients and source logs are retained in out.
%
% Name/value options: ModelFile, ReferenceDir, ClosedLoopRoot, OutputDir,
% AnalysisWindow, BaselineWindow, SampleDt, RateWindowSec, MaxGapSec,
% NominalOmegaRadS, PhaseCodePeriod, FeedbackStartSec, Conditions,
% ShowFigures, SaveOutputs. Conditions defaults to {'cos','w1'}; use
% {'cos','w1','k3w2','k5w2'} to include the high-gain trials.
% FrequencyFolders defaults to 250; a numeric vector applies to every
% condition, or a cell array supplies one numeric vector per condition.
% Delta omega = omega1 - omega2 = (2.5-folder/100)*pi. NominalOmegaRadS
% sets omega1; omega2 scales with folder/250 and is rounded to two decimals
% to match the historical parameter payload. At folder 250 both are 7.85.
% Outputs: MAT, per-file time-series CSV, summary/sensitivity CSV, PNG/PDF.
% ReplotOnly=true redraws the saved time series in OutputDir without
% recomputing the analysis or rewriting CSV/MAT results.
% Zoom figures compare Module 2 in two feedback-time windows, specified
% by ZoomFeedbackWindows (default [0 5; 60 65]).
% SaveZoomPanels=true also exports all six zoom panels as individual PNG/PDF.
% For k3w2/k5w2 at Delta omega=0.2*pi, 0.3*pi and 0.5*pi, it also saves six
% Module 1/2 overlays under OutputDir/overlaid_zoom/k{gain}_delta_omega_{d}pi,
% even with PlotTimeSeries=false or ReplotOnly=true. Phase-rate overlays
% compare the logged phase slopes with a reconstruction from logged a2 and W2.
% Set ShowReconstructedRate=false to show only the phase slopes from the log.

    here = fileparts(mfilename('fullpath'));
    estimate_dir = fileparts(here);
    repo_dir = fileparts(estimate_dir);
    addpath(repo_dir, estimate_dir);
    p = inputParser;
    addParameter(p, 'ModelFile', fullfile(estimate_dir, 'Spring3', '255', 'gamma_exports', 'gamma_export_latest.mat'));
    addParameter(p, 'ReferenceDir', '');
    addParameter(p, 'ClosedLoopRoot', fullfile(here, 'Spring3'));
    addParameter(p, 'OutputDir', fullfile(here, 'delta0_validation'));
    addParameter(p, 'AnalysisWindow', [60 120], @valid_window);
    addParameter(p, 'BaselineWindow', [1 4], @valid_window);
    addParameter(p, 'SampleDt', 0.01, @positive_scalar);
    addParameter(p, 'RateWindowSec', 0.10, @positive_scalar);
    addParameter(p, 'MaxGapSec', 0.05, @positive_scalar);
    addParameter(p, 'NominalOmegaRadS', 7.85, @positive_scalar);
    addParameter(p, 'PhaseCodePeriod', 255, @positive_scalar);
    addParameter(p, 'FeedbackStartSec', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(p, 'Conditions', {'cos','w1'}, @(x) iscellstr(x) && ~isempty(x) && numel(unique(x)) == numel(x));
    addParameter(p, 'FrequencyFolders', 250, @(x) isnumeric(x) || iscell(x));
    addParameter(p, 'ShowFigures', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'PlotTimeSeries', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SaveTimeSeries', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ReplotOnly', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ZoomFeedbackWindows', [0 5; 60 65], @(x) isnumeric(x) && isequal(size(x),[2 2]) && all(isfinite(x(:))) && all(x(:,1)>=0) && all(x(:,2)>x(:,1)));
    addParameter(p, 'SaveZoomPanels', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowReconstructedRate', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    cfg = p.Results;
    if cfg.ReplotOnly
        saved = load(fullfile(cfg.OutputDir, 'delta0_validation.mat'), 'out');
        out = saved.out;
        plot_cfg = out.config;
        defaults = fieldnames(cfg);
        for k = 1:numel(defaults)
            if ~isfield(plot_cfg, defaults{k}), plot_cfg.(defaults{k}) = cfg.(defaults{k}); end
        end
        plot_cfg.OutputDir = cfg.OutputDir;
        plot_cfg.ShowFigures = cfg.ShowFigures;
        plot_cfg.SaveOutputs = cfg.SaveOutputs;
        plot_cfg.ZoomFeedbackWindows = cfg.ZoomFeedbackWindows;
        plot_cfg.SaveZoomPanels = cfg.SaveZoomPanels;
        plot_cfg.ShowReconstructedRate = cfg.ShowReconstructedRate;
        n_redrawn = 0;
        for k = 1:numel(out.runs)
            run = out.runs(k);
            if ~ismember('Conditions', p.UsingDefaults) && ~ismember(run.condition,cfg.Conditions), continue; end
            [folder, stem] = fileparts(run.file);
            [~, folder_name] = fileparts(folder);
            folder_value = str2double(folder_name);
            assert(isfinite(folder_value), 'Cannot identify the frequency folder for %s.', run.file);
            delta_pi = (250-folder_value)/100;
            nominal = [plot_cfg.NominalOmegaRadS, round(plot_cfg.NominalOmegaRadS*folder_value/250,2)];
            label = [run.condition '_' stem];
            if folder_value ~= 250, label = sprintf('%s_%d_%s', run.condition, folder_value, stem); end
            if cfg.PlotTimeSeries
                plot_run(run.timeseries, plot_cfg, label, [0 plot_cfg.AnalysisWindow(2)], nominal, delta_pi);
                plot_run(run.timeseries, plot_cfg, [label '_zoom'], ...
                    [], nominal, delta_pi, true);
                n_redrawn = n_redrawn+1;
            end
        end
        selected_runs = out.runs;
        if ~ismember('Conditions', p.UsingDefaults)
            selected_runs = selected_runs(ismember({selected_runs.condition}, cfg.Conditions));
        end
        export_high_gain_zoom_panels(selected_runs, plot_cfg);
        out.dotphi_relative_deviation_stats = compute_dotphi_window_stats(out.runs, plot_cfg);
        if plot_cfg.SaveOutputs
            writetable(out.dotphi_relative_deviation_stats, ...
                fullfile(plot_cfg.OutputDir, 'dotphi_relative_deviation_stats.csv'));
        end
        fprintf('Updated full and zoom figures for %d saved runs: %s\n', n_redrawn, cfg.OutputDir);
        return;
    end
    assert(cfg.BaselineWindow(2) < cfg.AnalysisWindow(1), 'Baseline must precede analysis window.');
    S = load(cfg.ModelFile, 'gamma_export');
    G = S.gamma_export;
    assert(numel(G.phase_agent_ids) == 2, 'Model must have exactly two phase agents.');
    if isempty(cfg.ReferenceDir)
        cfg.ReferenceDir = fileparts(fileparts(cfg.ModelFile));
    end
    ref_files = dir(fullfile(cfg.ReferenceDir, '*.csv'));
    % Restrict calibration to the source files actually used for identification.
    used_names = {};
    for j = 1:numel(G.agents)
        for k = 1:numel(G.agents(j).used_files)
            [~, n, e] = fileparts(G.agents(j).used_files{k});
            used_names{end+1} = [n e]; %#ok<AGROW>
        end
    end
    ref_files = ref_files(ismember({ref_files.name}, used_names));
    assert(~isempty(ref_files), 'No original model source CSVs found in ReferenceDir.');
    models = cell(1, 2);
    ref_values = cell(1, 2);
    training_sse = zeros(1, 2);
    training_n = zeros(1, 2);
    for j = 1:2
        ix = find([G.agents.agent_id] == G.phase_agent_ids(j), 1);
        assert(~isempty(ix), 'Missing model agent.');
        models{j} = G.agents(ix).a2_gamma;
        assert(models{j}.available && ~isempty(models{j}.z_mean), 'Incomplete sensory model.');
    end
    for k = 1:numel(ref_files)
        a = load_corrected_agent_series_from_csv(fullfile(ref_files(k).folder, ref_files(k).name), ...
            G.phase_agent_ids, {'a0', 'a2'});
        ref_t0 = max([a(G.phase_agent_ids(1)).time(1), a(G.phase_agent_ids(2)).time(1)]);
        ref_end = min([a(G.phase_agent_ids(1)).time(end), a(G.phase_agent_ids(2)).time(end)]);
        ref_time = (ref_t0 + G.analysis_start_sec:cfg.SampleDt: ...
            min(ref_end, ref_t0 + G.analysis_start_sec + G.analysis_duration_sec)).';
        ref_phase = zeros(numel(ref_time), 2);
        for j = 1:2
            s = a(G.phase_agent_ids(j));
            ref_phase(:, j) = mod(sample_no_gaps(s.time, s.a0_corr, ref_time, cfg.MaxGapSec), 256)*(2*pi/256);
        end
        for j = 1:2
            s = a(G.phase_agent_ids(j));
            ref_values{j} = [ref_values{j}; s.a2(:)];
            bounds = prctile(s.a2, [10 90]);
            observed = clip_half((sample_no_gaps(s.time, s.a2, ref_time, cfg.MaxGapSec) - mean(bounds))/diff(bounds));
            predicted = predict_s(ref_phase(:, 1), ref_phase(:, 2), models{j});
            valid = isfinite(observed) & isfinite(predicted);
            training_sse(j) = training_sse(j) + sum((observed(valid)-predicted(valid)).^2);
            training_n(j) = training_n(j) + nnz(valid);
        end
    end
    reproduced_rmse = sqrt(training_sse./training_n);
    stored_rmse = cellfun(@(f) f.rmse, models);
    assert(all(abs(reproduced_rmse-stored_rmse) < 1e-4), ...
        'Source a2 RMSE does not reproduce saved fit: verify signal source and normalization.');
    ref_bounds = zeros(2, 2);
    for j = 1:2
        ref_bounds(j, :) = prctile(ref_values{j}, [10 90]);
    end
    assert(all(diff(ref_bounds, 1, 2) > 0), 'Open-loop sensor calibration has zero span.');
    if cfg.SaveOutputs && ~isfolder(cfg.OutputDir), mkdir(cfg.OutputDir); end
    out = struct('config', cfg, 'open_loop_model', G, 'reference_files', {ref_files}, ...
        'reference_p10_p90', ref_bounds, 'runs', struct([]));
    out.reference_fit_check = table(G.phase_agent_ids(:), stored_rmse(:), reproduced_rmse(:), ...
        'VariableNames', {'agent_id','saved_rmse','reproduced_rmse'});
    summary = struct([]);
    sensitivity = struct([]);
    conditions = cfg.Conditions(:).';
    folders_by_condition = cfg.FrequencyFolders;
    if isnumeric(folders_by_condition)
        folders_by_condition = repmat({folders_by_condition}, size(conditions));
    end
    assert(numel(folders_by_condition) == numel(conditions), 'Supply one folder vector per condition.');
    for c = 1:numel(conditions)
      folder_values = folders_by_condition{c};
      validateattributes(folder_values, {'numeric'}, {'vector','nonempty','positive','integer','finite'});
      assert(numel(unique(folder_values)) == numel(folder_values), 'Duplicate frequency folders.');
      for folder_value = folder_values(:).'
        delta_pi = (250-folder_value)/100;
        nominal_omega = [cfg.NominalOmegaRadS, round(cfg.NominalOmegaRadS*folder_value/250, 2)];
        files = dir(fullfile(cfg.ClosedLoopRoot, conditions{c}, sprintf('%d',folder_value), '*.csv'));
        assert(~isempty(files), 'No CSVs for %s/%d.', conditions{c}, folder_value);
        for k = 1:numel(files)
            fp = fullfile(files(k).folder, files(k).name);
            [a, ids] = load_corrected_agent_series_from_csv(fp, [], {'a0', 'a2'});
            ids = sort(ids(ids ~= 99));
            assert(numel(ids) == 2, 'Expected two non-99 agents in %s.', fp);
            assert(isequal(ids(:), G.phase_agent_ids(:)), ...
                'Closed/open-loop phase-agent IDs differ; explicitly establish correspondence first.');
            t0 = max([a(ids(1)).time(1), a(ids(2)).time(1)]);
            tend = min([a(ids(1)).time(end), a(ids(2)).time(end)]) - t0;
            assert(tend >= cfg.AnalysisWindow(2), 'Log does not cover the full analysis window.');
            % Include a derivative margin so the entire summary window has
            % centered slopes; the exported time series retains this margin.
            margin = actual_window(cfg.SampleDt, max(0.20, cfg.RateWindowSec))/2 + cfg.SampleDt;
            t = (0:cfg.SampleDt:min(tend, cfg.AnalysisWindow(2)+margin)).';
            phase = nan(numel(t), 2);
            model_phase = phase;
            raw = phase;
            file_bounds = zeros(2, 2);
            baseline = zeros(1, 2);
            for j = 1:2
                s = a(ids(j));
                native_t = s.time(:) - t0;
                code = mod(s.a0_corr(:), 256);
                physical_phase = unwrap(code * (2*pi/cfg.PhaseCodePeriod));
                phase(:, j) = sample_no_gaps(native_t, physical_phase, t, cfg.MaxGapSec);
                model_phase(:, j) = mod(sample_no_gaps(native_t, s.a0_corr(:), t, cfg.MaxGapSec), 256) * (2*pi/256);
                raw(:, j) = sample_no_gaps(native_t, s.a2(:), t, cfg.MaxGapSec);
                file_bounds(j, :) = prctile(s.a2, [10 90]);
                assert(diff(file_bounds(j, :)) > 0, 'Closed-loop sensor has zero span.');
                b = native_t >= cfg.BaselineWindow(1) & native_t <= cfg.BaselineWindow(2);
                assert(nnz(b) >= 20, 'Insufficient baseline samples.');
                q = polyfit(native_t(b), physical_phase(b), 1);
                baseline(j) = q(1);
            end
            rate = local_slope(phase, cfg.SampleDt, cfg.RateWindowSec);
            rate_raw = [gradient(phase(:, 1), cfg.SampleDt), gradient(phase(:, 2), cfg.SampleDt)];
            prediction = phase;
            s_file = phase;
            s_fixed = phase;
            for j = 1:2
                prediction(:, j) = predict_s(model_phase(:, 1), model_phase(:, 2), models{j});
                s_file(:, j) = clip_half((raw(:, j) - mean(file_bounds(j, :))) / diff(file_bounds(j, :)));
                s_fixed(:, j) = (raw(:, j) - mean(ref_bounds(j, :))) / diff(ref_bounds(j, :));
            end
            err_file = s_file - prediction;
            err_fixed = s_fixed - prediction;
            correction = rate - nominal_omega;
            correction_baseline = rate - baseline;
            T = table(t, 'VariableNames', {'time_s'});
            for j = 1:2
                names = {'phi_unwrapped_rad', 'phi_model_rad', 'dotphi_raw_rad_s', 'dotphi_rad_s', ...
                    'correction_nominal_rad_s', 'correction_baseline_rad_s', 'correction_over_nominal', ...
                    'a2_raw', 's_file_normalized', 's_fixed_scale', 's_predicted', 's_error_file', 's_error_fixed'};
                values = [phase(:, j), model_phase(:, j), rate_raw(:, j), rate(:, j), ...
                    correction(:, j), correction_baseline(:, j), correction(:, j)/nominal_omega(j), ...
                    raw(:, j), s_file(:, j), s_fixed(:, j), prediction(:, j), err_file(:, j), err_fixed(:, j)];
                for v = 1:numel(names)
                    T.(sprintf('%s_%d', names{v}, j)) = values(:, v);
                end
                main = t >= cfg.AnalysisWindow(1) & t <= cfg.AnalysisWindow(2);
                valid_s = main & isfinite(err_file(:, j)) & isfinite(err_fixed(:, j));
                valid_rate = main & isfinite(rate(:, j));
                assert(nnz(valid_s) > 0 && nnz(valid_rate) > 0, 'No valid evaluation samples.');
                r = struct('condition', conditions{c}, 'file', files(k).name, 'module', j, 'agent_id', ids(j), ...
                    'frequency_folder', folder_value, 'delta_omega_over_pi', delta_pi, ...
                    'nominal_omega_rad_s', nominal_omega(j), 'baseline_omega_rad_s', baseline(j), ...
                    'rate_window_actual_s', actual_window(cfg.SampleDt, cfg.RateWindowSec), ...
                    'n_rate', nnz(valid_rate), 'n_s', nnz(valid_s), ...
                    'mean_dotphi_rad_s', mean(rate(valid_rate, j)), ...
                    'mean_correction_rad_s', mean(correction(valid_rate, j)), ...
                    'rms_correction_rad_s', rms_value(correction(valid_rate, j)), ...
                    'rms_correction_over_omega', rms_value(correction(valid_rate, j))/nominal_omega(j), ...
                    'p95_abs_correction_over_omega', prctile(abs(correction(valid_rate, j)), 95)/nominal_omega(j), ...
                    'max_abs_correction_over_omega', max(abs(correction(valid_rate, j)))/nominal_omega(j), ...
                    'rms_correction_over_baseline', rms_value(correction_baseline(valid_rate, j))/baseline(j), ...
                    's_rmse_file', rms_value(err_file(valid_s, j)), ...
                    's_mae_file', mean(abs(err_file(valid_s, j))), 's_bias_file', mean(err_file(valid_s, j)), ...
                    's_rmse_fixed', rms_value(err_fixed(valid_s, j)), ...
                    's_rmse_fixed_clipped', rms_value(clip_half(s_fixed(valid_s, j)) - prediction(valid_s, j)), ...
                    's_open_loop_training_rmse', models{j}.rmse, ...
                    'closed_p10', file_bounds(j, 1), 'closed_p90', file_bounds(j, 2), ...
                    'reference_p10', ref_bounds(j, 1), 'reference_p90', ref_bounds(j, 2), ...
                    'sensor_span_ratio', diff(file_bounds(j, :))/diff(ref_bounds(j, :)), ...
                    'file_clip_fraction', mean(abs(s_file(valid_s, j)) >= 0.5), ...
                    'fixed_outside_half_fraction', mean(abs(s_fixed(valid_s, j)) > 0.5));
                summary = [summary; r]; %#ok<AGROW>
            end
            windows = unique([0.05, cfg.RateWindowSec, 0.20]);
            for w = windows
                rw = local_slope(phase, cfg.SampleDt, w);
                for j = 1:2
                    mask = main & isfinite(rw(:, j));
                    d = rw(mask, j) - nominal_omega(j);
                    sr = struct('condition', conditions{c}, 'file', files(k).name, 'module', j, ...
                        'frequency_folder', folder_value, 'delta_omega_over_pi', delta_pi, ...
                        'window_requested_s', w, 'window_actual_s', actual_window(cfg.SampleDt, w), ...
                        'n_samples', nnz(mask), 'mean_dotphi_rad_s', mean(rw(mask, j)), ...
                        'rms_correction_over_omega', rms_value(d)/nominal_omega(j), ...
                        'p95_abs_correction_over_omega', prctile(abs(d), 95)/nominal_omega(j));
                    sensitivity = [sensitivity; sr]; %#ok<AGROW>
                end
            end
            [~, stem] = fileparts(files(k).name);
            label = [conditions{c} '_' stem];
            if folder_value ~= 250, label = sprintf('%s_%d_%s', conditions{c}, folder_value, stem); end
            run = struct('condition', conditions{c}, 'file', fp, 'agent_ids', ids, ...
                'frequency_folder', folder_value, 'delta_omega_over_pi', delta_pi, 'nominal_omega_rad_s', nominal_omega, ...
                'common_start_abs', t0, 'baseline_omega_rad_s', baseline, 'file_p10_p90', file_bounds, ...
                'timeseries', T);
            out.runs = [out.runs; run];
            if cfg.SaveOutputs && cfg.SaveTimeSeries
                writetable(T, fullfile(cfg.OutputDir, [label '_timeseries.csv']));
            end
            if cfg.PlotTimeSeries && (cfg.ShowFigures || cfg.SaveOutputs)
                plot_run(T, cfg, label, [0 cfg.AnalysisWindow(2)], nominal_omega, delta_pi);
                plot_run(T, cfg, [label '_zoom'], [], nominal_omega, delta_pi, true);
            end
        end
      end
    end
    export_high_gain_zoom_panels(out.runs, cfg);
    out.summary = struct2table(summary);
    out.rate_window_sensitivity = struct2table(sensitivity);
    out.dotphi_relative_deviation_stats = compute_dotphi_window_stats(out.runs, cfg);
    out.notes = { ...
        'Rate correction dotphi-omega is observational; exact logged feedback is unavailable.', ...
        'Derivative code period=255 by default; Fourier model input code period=256.', ...
        'Per-file normalized sensor errors can hide amplitude and offset changes; see fixed-scale errors.', ...
        'Fixed scale pools full-log raw a2 from original model source files; no closed-loop fitting.', ...
        'Model is the chosen saved open-loop fit, not a verified copy of the model deployed in April.', ...
        'Each CSV is one trial; samples are serially correlated, not independent trials.'};
    if cfg.SaveOutputs
        writetable(out.summary, fullfile(cfg.OutputDir, 'summary.csv'));
        writetable(out.rate_window_sensitivity, fullfile(cfg.OutputDir, 'rate_window_sensitivity.csv'));
        writetable(out.dotphi_relative_deviation_stats, ...
            fullfile(cfg.OutputDir, 'dotphi_relative_deviation_stats.csv'));
        save(fullfile(cfg.OutputDir, 'delta0_validation.mat'), 'out');
    end
    fprintf('\nOpen-loop model: %s\nExported: %s\n', cfg.ModelFile, G.created_at);
    fprintf('Open-loop source RMSE reproduced (a2 normalization/phase conventions checked):\n');
    disp(out.reference_fit_check);
    fprintf('Summary window: %.1f--%.1f s; baseline: %.1f--%.1f s\n', cfg.AnalysisWindow, cfg.BaselineWindow);
    disp(out.summary(:, {'condition','delta_omega_over_pi','module','baseline_omega_rad_s','mean_dotphi_rad_s', ...
        'rms_correction_over_omega','s_rmse_file','s_rmse_fixed'}));
    fprintf('Results directory: %s\n', cfg.OutputDir);
end

function yes = valid_window(x)
    yes = isnumeric(x) && numel(x) == 2 && all(isfinite(x)) && x(1) >= 0 && x(2) > x(1);
end

function yes = positive_scalar(x)
    yes = isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
end

function y = sample_no_gaps(t, x, tq, max_gap)
    y = interp1(t, x, tq, 'linear', NaN);
    gaps = find(diff(t) > max_gap);
    for k = gaps(:).'
        y(tq > t(k) & tq < t(k+1)) = NaN;
    end
end

function duration = actual_window(dt, requested)
    duration = 2 * max(1, round(requested/(2*dt))) * dt;
end

function rate = local_slope(phase, dt, requested)
    h = max(1, round(requested/(2*dt)));
    offsets = (-h:h).' * dt;
    kernel = flipud(offsets/sum(offsets.^2));
    rate = nan(size(phase));
    for j = 1:size(phase, 2)
        x = phase(:, j);
        bad = ~isfinite(x);
        x(bad) = 0;
        d = conv(x, kernel, 'same');
        d(conv(double(bad), ones(2*h+1, 1), 'same') > 0) = NaN;
        d(1:h) = NaN;
        d(end-h+1:end) = NaN;
        rate(:, j) = d;
    end
end

function v = clip_half(x)
    v = min(max(x, -0.5), 0.5);
    v(~isfinite(x)) = NaN;
end

function v = rms_value(x)
    v = sqrt(mean(x.^2));
end

function y = predict_s(phi1, phi2, fit)
    assert(numel(fit.coeff) == numel(fit.basis_types), 'Invalid model dimensions.');
    y = zeros(size(phi1)) + fit.z_mean;
    for k = 1:numel(fit.coeff)
        m = fit.basis_phi1_order(k);
        n = fit.basis_phi2_order(k);
        switch fit.basis_types{k}
            case 'constant', b = ones(size(phi1));
            case 'phi1_cos', b = cos(m*phi1);
            case 'phi1_sin', b = sin(m*phi1);
            case 'phi2_cos', b = cos(n*phi2);
            case 'phi2_sin', b = sin(n*phi2);
            case 'mixed_cc', b = cos(m*phi1).*cos(n*phi2);
            case 'mixed_cs', b = cos(m*phi1).*sin(n*phi2);
            case 'mixed_sc', b = sin(m*phi1).*cos(n*phi2);
            case 'mixed_ss', b = sin(m*phi1).*sin(n*phi2);
            otherwise, error('Unsupported basis type: %s', fit.basis_types{k});
        end
        y = y + fit.coeff(k)*b;
    end
    y(~isfinite(phi1) | ~isfinite(phi2)) = NaN;
end

function plot_run(T, cfg, label, limits, nominal_omega, delta_pi, compare_windows)
    if nargin < 7, compare_windows = false; end
    plot_time = T.time_s - cfg.FeedbackStartSec;
    if compare_windows
        modules = [2 2];
        panel_limits = cfg.ZoomFeedbackWindows;
    else
        modules = [1 2];
        panel_limits = repmat(max(0,limits-cfg.FeedbackStartSec),2,1);
    end
    visibility = 'off';
    if cfg.ShowFigures, visibility = 'on'; end
    f = figure('Color', 'w', 'Visible', visibility, 'Position', [80 60 1400 900]);
    layout = tiledlayout(f, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, sprintf('%s | Delta omega=%+.1f pi | slope window %.3g s', strrep(label, '_', ' '), ...
        delta_pi, actual_window(cfg.SampleDt, cfg.RateWindowSec)), 'Interpreter', 'none');
    panel_axes = gobjects(3,2);
    for column = 1:2
        j = modules(column);
        bounds = panel_limits(column,:);
        assert(bounds(1)>=min(plot_time)-cfg.SampleDt/2 && bounds(2)<=max(plot_time)+cfg.SampleDt/2, ...
            'Plot window %.1f--%.1f s is outside the saved feedback-time range.',bounds);
        mask = plot_time>=bounds(1) & plot_time<=bounds(2);
        times = plot_time(mask);
        data = T(mask,:);
        ax = nexttile(layout, column);
        panel_axes(1,column) = ax;
        plot(ax, times, data.(sprintf('dotphi_rad_s_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        if j == 2
            yline(ax,nominal_omega(j),':','Color',[.1 .1 .1],'LineWidth',2.5, ...
                'Tag','NaturalOmega2','HandleVisibility','off');
        else
            yline(ax, nominal_omega(j), '--', 'Natural frequency', 'Color', [.3 .3 .3]);
        end
        ylabel(ax, sprintf('$\\dot{\\phi}_{%d}\\;(\\mathrm{rad/s})$', j), 'Interpreter', 'latex');
        ylim(ax, [0 3*pi]);
        if compare_windows
            title(ax,sprintf('Module 2: %.0f--%.0f s',bounds));
        else
            title(ax, sprintf('Module %d: phase rate', j));
        end
        ax = nexttile(layout, column+2);
        panel_axes(2,column) = ax;
        plot(ax, times, data.(sprintf('s_file_normalized_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        plot(ax, times, data.(sprintf('s_predicted_%d', j)), 'Color', [.85 .25 .1]);
        ylabel(ax, sprintf('$s_{%d}$', j), 'Interpreter', 'latex');
        ylim(ax, [-0.6 1.2]);
        yticks(ax, [-0.5 0 0.5 1]);
        legend(ax, 'Measured (file scale)', 'Open-loop prediction', 'Location', 'northeast', 'Interpreter', 'latex');
        ax = nexttile(layout, column+4);
        panel_axes(3,column) = ax;
        plot(ax, times, data.(sprintf('s_error_file_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        yline(ax, 0, 'k--', 'HandleVisibility', 'off');
        ylabel(ax, sprintf('$s_{%d}^{\\mathrm{meas}}-s_{%d}^{\\mathrm{pred}}$', j, j), 'Interpreter', 'latex');
        legend(ax, sprintf('$s_{%d}^{\\mathrm{meas}}-s_{%d}^{\\mathrm{pred}}$', j, j), 'Location', 'best', 'Interpreter', 'latex');
        for row = 1:3
            xlim(panel_axes(row,column),bounds);
        end
    end
    if compare_windows
        % Use identical vertical scales for the two time windows.
        for row = 1:3
            a = ylim(panel_axes(row,1));
            b = ylim(panel_axes(row,2));
            shared_limits = [min(a(1),b(1)),max(a(2),b(2))];
            ylim(panel_axes(row,1),shared_limits);
            ylim(panel_axes(row,2),shared_limits);
        end
    end
    axes_list = findall(f, 'Type', 'axes');
    for ax = axes_list(:).'
        grid(ax, 'on');
        xlabel(ax, '$t\;(\mathrm{s})$', 'Interpreter', 'latex');
    end
    tune_analysis_figure(f);
    for ax = axes_list(:).'
        ax.TickLabelInterpreter = 'latex';
    end
    for column = 1:2
        if modules(column)==2, style_module2_phase_rate(panel_axes(1,column),nominal_omega(2)); end
    end
    if cfg.SaveOutputs
        exportgraphics(f, fullfile(cfg.OutputDir, [label '.png']), 'Resolution', 160);
        exportgraphics(f, fullfile(cfg.OutputDir, [label '.pdf']), 'ContentType', 'vector');
        if compare_windows && cfg.SaveZoomPanels
            zoom_prefix = regexprep(label, '_merged_.*$', '');
            quantities = {'dotphi2','s2','s2_error'};
            for column = 1:2
                for row = 1:3
                    stem = sprintf('%s_%s_%g_%gs',zoom_prefix,quantities{row},panel_limits(column,:));
                    export_zoom_panel(panel_axes(row,column),cfg,stem);
                end
            end
            export_overlaid_zoom_panels(T,plot_time,cfg,zoom_prefix,panel_limits,nominal_omega);
        end
    end
    if ~cfg.ShowFigures, close(f); end
end

function export_zoom_panel(source_axes,cfg,stem)
    f = figure('Color','w','Visible','off','Position',[80 60 800 500]);
    ax = copyobj(source_axes,f);
    ax.Units = 'normalized';
    ax.Position = [0.16 0.18 0.80 0.68];
    title(ax,'');
    if ~isempty(source_axes.Legend) && isgraphics(source_axes.Legend)
        curves = flipud(findobj(ax,'Type','line'));
        legend(ax,curves,source_axes.Legend.String,'Location',source_axes.Legend.Location, ...
            'Interpreter',source_axes.Legend.Interpreter);
    end
    tune_analysis_figure(f);
    enlarge_panel_fonts(f, ax, 1.28);
    reference = findobj(ax,'Tag','NaturalOmega2');
    if ~isempty(reference), style_module2_phase_rate(ax,reference(1).Value); end
    export_panel_files(f,cfg,stem);
    close(f);
end

function style_module2_phase_rate(ax,omega2,omega_label,omega1)
    if nargin < 3, omega_label = '$\omega_2$'; end
    % Apply after tuneFigure, which otherwise resets grid opacity.
    ax.GridColor = [.8 .8 .8];
    ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.08;
    reference = findobj(ax,'Tag','NaturalOmega2');
    set(reference,'Color',[.1 .1 .1],'LineStyle',':','LineWidth',2.5,'Label','');
    set(findobj(ax,'Tag','NaturalOmega1'),'Color',[.1 .1 .1],'LineWidth',2.5);
    ticks = [0 pi 2*pi 3*pi];
    labels = {'$0$','$\pi$','$2\pi$','$3\pi$'};
    keep = ticks>=ax.YLim(1)-1e-9 & ticks<=ax.YLim(2)+1e-9 & abs(ticks-omega2)>0.3*pi;
    natural_ticks = omega2;
    natural_labels = {omega_label};
    if nargin >= 4 && abs(omega1-omega2) > 1e-6
        keep = keep & abs(ticks-omega1)>0.3*pi;
        natural_ticks = [omega2 omega1];
        natural_labels = {omega_label,'$\omega_1$'};
    end
    ticks = [ticks(keep),natural_ticks];
    labels = [labels(keep),natural_labels];
    [ticks,order] = sort(ticks);
    labels = labels(order);
    ax.YTick = ticks;
    ax.YTickLabel = labels;
    ax.TickLabelInterpreter = 'latex';
end

function export_high_gain_zoom_panels(runs, cfg)
    if ~cfg.SaveOutputs || ~cfg.SaveZoomPanels, return; end
    for k = 1:numel(runs)
        run = runs(k);
        switch run.condition
            case 'k3w2', gain = 3;
            case 'k5w2', gain = 5;
            otherwise, continue;
        end
        [folder, trial_stem] = fileparts(run.file);
        [~, folder_name] = fileparts(folder);
        folder_value = str2double(folder_name);
        delta_pi = (250-folder_value)/100;
        if ~any(abs(delta_pi-[0.2 0.3 0.5]) < 1e-8), continue; end
        delta_label = strrep(sprintf('%.1f',delta_pi),'.','p');
        prefix = sprintf('k%d_delta_omega_%spi',gain,delta_label);
        panel_cfg = cfg;
        panel_cfg.OutputDir = fullfile(cfg.OutputDir,'overlaid_zoom',prefix);
        % Preserve every trial if additional CSVs are added to a condition.
        same_condition = strcmp({runs.condition},run.condition);
        same_folder = cellfun(@(f) strcmp(fileparts(f),folder),{runs.file});
        if nnz(same_condition & same_folder) > 1
            panel_cfg.OutputDir = fullfile(panel_cfg.OutputDir,trial_stem);
        end
        if ~isfolder(panel_cfg.OutputDir), mkdir(panel_cfg.OutputDir); end
        nominal = [cfg.NominalOmegaRadS, round(cfg.NominalOmegaRadS*folder_value/250,2)];
        estimated_rate = [];
        if cfg.ShowReconstructedRate
            % Infer the controller settings from the experiment batches. May 27/28
            % used Spring3/255 W2 with no lag, while May 31/June 1 used 240_0531 and 1 s.
            reference = '240_0531';
            tau = 1;
            if startsWith(trial_stem, {'merged_20260527_', 'merged_20260528_'})
                reference = '255';
                tau = 0;
            end
            prc_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'Spring3', reference, 'gamma_exports', 'prc_snippet_ref_w2.txt');
            [prc_a, prc_b] = parse_prc_snippet_file(prc_file);
            estimated_rate = reconstruct_w2_rate(run.timeseries, cfg, nominal, gain, prc_a, prc_b, tau);
        end
        export_overlaid_zoom_panels(run.timeseries, ...
            run.timeseries.time_s-cfg.FeedbackStartSec, panel_cfg, prefix, ...
            cfg.ZoomFeedbackWindows, nominal, estimated_rate);
        fprintf('Saved six Module 1/2 zoom overlays (PNG/PDF): %s\n',panel_cfg.OutputDir);
    end
end

function rate = reconstruct_w2_rate(T, cfg, nominal_omega, gain, prc_a, prc_b, tau)
    % Each logged a2 is floor(rawADC/4), saturated at 255. Approximate the
    % firmware's 1000-sample ring buffer at 2 ms using 2 s of logged samples.
    window_samples = max(1, round(2/cfg.SampleDt));
    percentile_index = floor(window_samples*[0.1 0.9])+1;
    pre_feedback = T.time_s < cfg.FeedbackStartSec;
    slope_samples = round(actual_window(cfg.SampleDt,cfg.RateWindowSec)/cfg.SampleDt)+1;
    alpha = 1;
    if tau > 0, alpha = -expm1(-cfg.SampleDt/tau); end
    rate = nan(height(T),2);
    for j = 1:2
        a2 = T.(sprintf('a2_raw_%d',j));
        sensor = (4*a2+1.5)/4095;
        sensor(a2 >= 255) = NaN;
        % Use the center of each phase-code bin to evaluate the uploaded W2.
        phase = T.(sprintf('phi_unwrapped_rad_%d',j))+pi/cfg.PhaseCodePeriod;
        angles = mod(phase,2*pi)*(0:numel(prc_a)-1);
        z = cos(angles)*round(prc_a,6)+sin(angles)*round(prc_b,6);
        correction = nan(height(T),1);
        correction(pre_feedback) = 0;
        for k = window_samples:height(T)
            if pre_feedback(k) || ~isfinite(z(k)), continue; end
            samples = sort(sensor(k-window_samples+1:k),'MissingPlacement','last');
            if ~isfinite(samples(end)), continue; end
            bounds = samples(percentile_index);
            if bounds(1) == bounds(2), bounds = [0.10 0.16]; end
            flex = clip_half((sensor(k)-mean(bounds))/diff(bounds));
            correction(k) = gain*z(k)*flex;
        end
        % Compare at the same slope-window resolution as the logged phase.
        % Start a new filter segment when a gap makes the input unavailable.
        edges = diff([false; isfinite(correction); false]);
        starts = find(edges == 1);
        stops = find(edges == -1)-1;
        for k = 1:numel(starts)
            ix = starts(k):stops(k);
            if numel(ix) < slope_samples, continue; end
            filtered = filter(alpha,[1 alpha-1],correction(ix));
            phase_correction = cumsum(filtered)*cfg.SampleDt;
            rate(ix,j) = nominal_omega(j)+local_slope(phase_correction,cfg.SampleDt,cfg.RateWindowSec);
        end
    end
end

function export_overlaid_zoom_panels(T,plot_time,cfg,zoom_prefix,panel_limits,nominal_omega,estimated_rate)
    if nargin < 7, estimated_rate = []; end
    all_mask = false(size(plot_time));
    for column = 1:2
        bounds = panel_limits(column,:);
        assert(bounds(1)>=min(plot_time)-cfg.SampleDt/2 && bounds(2)<=max(plot_time)+cfg.SampleDt/2, ...
            'Plot window %.1f--%.1f s is outside the saved feedback-time range.',bounds);
        all_mask = all_mask | (plot_time>=bounds(1) & plot_time<=bounds(2));
    end
    err_all = [T.s_error_file_1(all_mask); T.s_error_file_2(all_mask)];
    err_min = floor(min(err_all)*10)/10;
    err_max = ceil(max(err_all)*10)/10;
    err_ylim = [min(-0.32, err_min), max(0.60, err_max + 0.30)];
    omega_label = '$\omega_{1,2}$';
    if abs(nominal_omega(1)-nominal_omega(2)) > 1e-6
        omega_label = '$\omega_2$';
    end
    quantities = {'dotphi12','s12','s12_error'};
    for column = 1:2
        bounds = panel_limits(column,:);
        mask = plot_time>=bounds(1) & plot_time<=bounds(2);
        times = plot_time(mask);
        data = T(mask,:);
        for row = 1:3
            f = figure('Color','w','Visible','off','Position',[80 60 800 500]);
            ax = axes(f,'Units','normalized','Position',[0.16 0.18 0.80 0.68]);
            hold(ax,'on');
            if row == 1
                p1 = plot(ax,times,data.dotphi_rad_s_1,'Color',[0 .35 .7]);
                p2 = plot(ax,times,data.dotphi_rad_s_2,'Color',[.85 .25 .1]);
                if abs(nominal_omega(1)-nominal_omega(2)) > 1e-6
                    yline(ax,nominal_omega(1),'--','Color',[.1 .1 .1],'LineWidth',2.5, ...
                        'Tag','NaturalOmega1','HandleVisibility','off');
                end
                yline(ax,nominal_omega(2),':','Color',[.1 .1 .1],'LineWidth',2.5, ...
                    'Tag','NaturalOmega2','HandleVisibility','off');
                ylabel(ax,'$\dot{\phi}_1,\dot{\phi}_2\;(\mathrm{rad/s})$','Interpreter','latex');
                ylim(ax,[0 3*pi]);
                if isempty(estimated_rate)
                    legend(ax,[p1 p2],{'$\dot{\phi}_1$','$\dot{\phi}_2$'}, ...
                        'Location','southeast','Interpreter','latex');
                else
                    p3 = plot(ax,times,estimated_rate(mask,1),'--','Color',[0 .60 .85]);
                    p4 = plot(ax,times,estimated_rate(mask,2),'--','Color',[.93 .60 0]);
                    legend(ax,[p1 p2 p3 p4],{'$\dot{\phi}_1$ (log)','$\dot{\phi}_2$ (log)', ...
                        '$\dot{\phi}_1$ from $s_1$','$\dot{\phi}_2$ from $s_2$'}, ...
                        'Location','southeast','NumColumns',2,'Interpreter','latex');
                end
            elseif row == 2
                p1 = plot(ax,times,data.s_file_normalized_1,'-','Color',[0 .35 .7]);
                p2 = plot(ax,times,data.s_predicted_1,'--','Color',[0 .60 .85]);
                p3 = plot(ax,times,data.s_file_normalized_2,'-','Color',[.85 .25 .1]);
                p4 = plot(ax,times,data.s_predicted_2,'--','Color',[.93 .60 0]);
                ylabel(ax,'$s_1,s_2$','Interpreter','latex');
                ylim(ax,[-0.6 1.3]);
                yticks(ax,[-0.5 0 0.5 1]);
                legend(ax,[p1 p2 p3 p4], ...
                    {'$s_1^{\mathrm{meas}}$','$s_1^{\mathrm{pred}}$','$s_2^{\mathrm{meas}}$','$s_2^{\mathrm{pred}}$'}, ...
                    'Location','northeast','NumColumns',2,'Interpreter','latex');
            else
                p1 = plot(ax,times,data.s_error_file_1,'Color',[0 .35 .7]);
                p2 = plot(ax,times,data.s_error_file_2,'Color',[.85 .25 .1]);
                yline(ax,0,'k--','HandleVisibility','off');
                ylabel(ax,'$s_j^{\mathrm{meas}}-s_j^{\mathrm{pred}}$','Interpreter','latex');
                ylim(ax,err_ylim);
                legend(ax,[p1 p2], ...
                    {'$s_1^{\mathrm{meas}}-s_1^{\mathrm{pred}}$','$s_2^{\mathrm{meas}}-s_2^{\mathrm{pred}}$'}, ...
                    'Location','northeast','Interpreter','latex');
            end
            xlim(ax,bounds);
            grid(ax,'on');
            xlabel(ax,'$t\;(\mathrm{s})$','Interpreter','latex');
            tune_analysis_figure(f);
            enlarge_panel_fonts(f, ax, 1.28);
            if row == 1
                style_module2_phase_rate(ax,nominal_omega(2),omega_label,nominal_omega(1));
            elseif row == 2
                ax.Legend.NumColumns = 2;
                ax.Legend.Location = 'northeast';
            end
            stem = sprintf('%s_%s_%g_%gs',zoom_prefix,quantities{row},bounds);
            export_panel_files(f,cfg,stem);
            close(f);
        end
    end
end

function enlarge_panel_fonts(fig, ax, scale)
    objs = findall(fig, '-property', 'FontSize');
    sizes = arrayfun(@(h) h.FontSize, objs);
    for k = 1:numel(objs)
        objs(k).FontSize = sizes(k) * scale;
    end
    ax.TickLabelInterpreter = 'latex';
    ax.Units = 'normalized';
    ax.Position = [0.21 0.25 0.75 0.68];
    drawnow;
end

function export_panel_files(f, cfg, stem)
    exportgraphics(f, fullfile(cfg.OutputDir, [stem '.png']), 'Resolution', 160);
    pdf_path = fullfile(cfg.OutputDir, [stem '.pdf']);
    try
        exportgraphics(f, pdf_path, 'ContentType', 'vector');
    catch err
        warning('Could not overwrite %s (%s).', pdf_path, err.message);
    end
end

function stats_table = compute_dotphi_window_stats(runs, cfg)
    rows = struct([]);
    fb_end = cfg.AnalysisWindow(2) - cfg.FeedbackStartSec;
    windows = [cfg.ZoomFeedbackWindows; 0 fb_end; -cfg.FeedbackStartSec fb_end];
    window_names = { ...
        sprintf('fb_%g_%gs', cfg.ZoomFeedbackWindows(1,:)), ...
        sprintf('fb_%g_%gs', cfg.ZoomFeedbackWindows(2,:)), ...
        sprintf('fb_0_%gs', fb_end), ...
        sprintf('log_0_%gs', cfg.AnalysisWindow(2))};
    for r = 1:numel(runs)
        run = runs(r);
        T = run.timeseries;
        fb_time = T.time_s - cfg.FeedbackStartSec;
        [folder, fname, fext] = fileparts(run.file);
        if isfield(run, 'nominal_omega_rad_s')
            nominal_omega = run.nominal_omega_rad_s;
        else
            [~, folder_name] = fileparts(folder);
            folder_value = str2double(folder_name);
            nominal_omega = [cfg.NominalOmegaRadS, round(cfg.NominalOmegaRadS*folder_value/250, 2)];
        end
        for w = 1:size(windows, 1)
            bounds = windows(w, :);
            mask = fb_time >= bounds(1) - 1e-9 & fb_time <= bounds(2) + 1e-9;
            for j = 1:2
                rel = T.(sprintf('correction_over_nominal_%d', j))(mask);
                rel = rel(isfinite(rel));
                abs_rel = abs(rel);
                row = struct( ...
                    'condition', run.condition, ...
                    'file', [fname fext], ...
                    'window', window_names{w}, ...
                    'feedback_start_s', bounds(1), ...
                    'feedback_end_s', bounds(2), ...
                    'log_start_s', bounds(1) + cfg.FeedbackStartSec, ...
                    'log_end_s', bounds(2) + cfg.FeedbackStartSec, ...
                    'module', j, ...
                    'nominal_omega_rad_s', nominal_omega(j), ...
                    'n_samples', numel(abs_rel), ...
                    'max_abs_correction_over_omega', max(abs_rel), ...
                    'var_abs_correction_over_omega', var(abs_rel), ...
                    'std_abs_correction_over_omega', std(abs_rel), ...
                    'mean_abs_correction_over_omega', mean(abs_rel), ...
                    'rms_correction_over_omega', rms_value(rel));
                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end
    stats_table = struct2table(rows);
end
