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
%   each independently from the baseline as well. dot(phi)-omega is an observed
%   rate correction, not a reconstruction of the unlogged k*z(phi)*flex.
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

    here = fileparts(mfilename('fullpath'));
    estimate_dir = fileparts(here);
    repo_dir = fileparts(estimate_dir);
    addpath(repo_dir);
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
    parse(p, varargin{:});
    cfg = p.Results;
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
                plot_run(T, cfg, [label '_zoom'], [cfg.AnalysisWindow(1), min(cfg.AnalysisWindow(1)+5, cfg.AnalysisWindow(2))], nominal_omega, delta_pi);
            end
        end
      end
    end
    out.summary = struct2table(summary);
    out.rate_window_sensitivity = struct2table(sensitivity);
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

function plot_run(T, cfg, label, limits, nominal_omega, delta_pi)
    plot_time = T.time_s - cfg.FeedbackStartSec;
    limits = max(0, limits - cfg.FeedbackStartSec);
    visibility = 'off';
    if cfg.ShowFigures, visibility = 'on'; end
    f = figure('Color', 'w', 'Visible', visibility, 'Position', [80 60 1400 900]);
    layout = tiledlayout(f, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, sprintf('%s | Delta omega=%+.1f pi | slope window %.3g s', strrep(label, '_', ' '), ...
        delta_pi, actual_window(cfg.SampleDt, cfg.RateWindowSec)), 'Interpreter', 'none');
    for j = 1:2
        ax = nexttile(layout, j);
        plot(ax, plot_time, T.(sprintf('dotphi_rad_s_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        yline(ax, nominal_omega(j), '--', 'Natural frequency', 'Color', [.3 .3 .3]);
        ylabel(ax, sprintf('$\\dot{\\phi}_{%d}\\;(\\mathrm{rad}\\,\\mathrm{s}^{-1})$', j), 'Interpreter', 'latex');
        ylim(ax, [0 Inf]);
        title(ax, sprintf('Module %d: phase rate', j));
        ax = nexttile(layout, j+2);
        plot(ax, plot_time, T.(sprintf('s_file_normalized_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        plot(ax, plot_time, T.(sprintf('s_predicted_%d', j)), 'Color', [.85 .25 .1]);
        ylabel(ax, sprintf('$s_{%d}$', j), 'Interpreter', 'latex');
        legend(ax, 'Measured (file scale)', 'Open-loop prediction', 'Location', 'best');
        ax = nexttile(layout, j+4);
        plot(ax, plot_time, T.(sprintf('s_error_file_%d', j)), 'Color', [0 .35 .7]);
        hold(ax, 'on');
        yline(ax, 0, 'k--');
        ylabel(ax, sprintf('$s_{%d}^{\\mathrm{meas}}-\\hat{s}_{%d}$', j, j), 'Interpreter', 'latex');
        legend(ax, 'File scale error', 'Location', 'best');
    end
    axes_list = findall(f, 'Type', 'axes');
    for ax = axes_list(:).'
        xlim(ax, limits);
        grid(ax, 'on');
        xlabel(ax, '$t\;(\mathrm{s})$', 'Interpreter', 'latex');
    end
    if cfg.SaveOutputs
        exportgraphics(f, fullfile(cfg.OutputDir, [label '.png']), 'Resolution', 160);
        exportgraphics(f, fullfile(cfg.OutputDir, [label '.pdf']), 'ContentType', 'vector');
    end
    if ~cfg.ShowFigures, close(f); end
end
