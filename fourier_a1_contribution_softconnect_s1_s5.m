% a2時系列のフーリエ級数寄与（Spring1 / Spring5）
% - 対象: ディレクトリ内のCSVをファイル名ソート順でn番目指定
%   SoftConnect\Spring1 (fileIndex1)
%   SoftConnect\Spring5 (fileIndex5)
% - 同一ファイル内の agent_id ごとに分割して個別解析
% - 前後10%をカット
% - FFTで連続スペクトラムを表示
% - 角周波数範囲: 0.3*2.5*pi 〜 1.2*2.5*pi (rad/s)
% - Spring1〜Spring5の最小/最大agent_idスペクトラムを重ね描き

clear; close all;

fileindex = 16;
fileIndex1 = fileindex;
fileIndex5 = fileindex;

valueCol = 'a2';

dir1 = fullfile('SoftConnect', 'Spring1');
dir5 = fullfile('SoftConnect', 'Spring5');

file1 = pick_nth_csv(dir1, fileIndex1);
file5 = pick_nth_csv(dir5, fileIndex5);

trimRatio = 0.10;
omegaBase = 2.5 * pi;
omegaRange = [0.0, 5.0] * omegaBase;
nFftFactor = 1;

fprintf('FFT omega range: %.6f .. %.6f rad/s\n', omegaRange(1), omegaRange(2));
fprintf('trim = %.0f%%%% each side\n\n', trimRatio * 100);

analyze_file_by_agent(file1, 'Spring1', trimRatio, omegaRange, valueCol, nFftFactor);
analyze_file_by_agent(file5, 'Spring5', trimRatio, omegaRange, valueCol, nFftFactor);

overlay_min_agent_spectra(fileindex, trimRatio, omegaRange, valueCol, nFftFactor);
overlay_max_agent_spectra(fileindex, trimRatio, omegaRange, valueCol, nFftFactor);

function csvPath = pick_nth_csv(parentDir, indexN)
    files = dir(fullfile(parentDir, '*.csv'));
    if isempty(files)
        error('CSVが見つかりません: %s', parentDir);
    end

    names = {files.name};
    names = sort(names);

    if indexN < 1 || indexN > numel(names)
        error('indexNが範囲外です: %s (1..%d)', parentDir, numel(names));
    end

    csvPath = fullfile(parentDir, names{indexN});
    fprintf('Selected CSV: %s (index %d)\n', csvPath, indexN);
end

function analyze_file_by_agent(csvPath, labelName, trimRatio, omegaRange, valueCol, nFftFactor)
    T = readtable(csvPath);

    if ~ismember('agent_id', T.Properties.VariableNames)
        error('agent_id列が見つかりません: %s', csvPath);
    end

    if ismember('time_local_sec', T.Properties.VariableNames)
        tAll = T.time_local_sec;
    elseif ismember('time_pc_sec_abs', T.Properties.VariableNames)
        tAll = T.time_pc_sec_abs;
    else
        error('time列が見つかりません: %s', csvPath);
    end

    if ~ismember(valueCol, T.Properties.VariableNames)
        error('%s列が見つかりません: %s', valueCol, csvPath);
    end
    xAll = T.(valueCol);
    agentAll = T.agent_id;

    tAll = double(tAll(:));
    xAll = double(xAll(:));
    agentAll = double(agentAll(:));

    validAll = isfinite(tAll) & isfinite(xAll) & isfinite(agentAll);
    tAll = tAll(validAll);
    xAll = xAll(validAll);
    agentAll = agentAll(validAll);

    agentList = unique(agentAll, 'sorted');

    fprintf('--- %s (%s) ---\n', labelName, csvPath);
    fprintf('detected agent_id: ');
    fprintf('%g ', agentList);
    fprintf('\n');

    for i = 1:numel(agentList)
        aid = agentList(i);
        idxAgent = (agentAll == aid);
        t = tAll(idxAgent);
        x = xAll(idxAgent);

        % 時刻順に並び替えて解析
        [t, order] = sort(t);
        x = x(order);

        [omegaVec, ampVec, tTrim, xTrim] = compute_fft_spectrum_from_vectors(t, x, trimRatio, csvPath, aid, nFftFactor);
        idxRange = (omegaVec >= omegaRange(1)) & (omegaVec <= omegaRange(2));

        figure('Color', 'w', 'Name', sprintf('%s FFT: %s agent %g', valueCol, labelName, aid));
        plot(omegaVec(idxRange), ampVec(idxRange), 'LineWidth', 1.2, 'Color', [0.2 0.45 0.8]);
        grid on;
        xlabel('Angular frequency \omega (rad/s)');
        ylabel('Amplitude');
        xlim(omegaRange);
        title(sprintf('%s agent_id=%g (FFT %.4f..%.4f rad/s)', ...
            labelName, aid, omegaRange(1), omegaRange(2)));
        tuneFigure;

        [tSeg, xSeg] = extract_middle_10sec(tTrim, xTrim);
        figure('Color', 'w', 'Name', sprintf('%s waveform: %s agent %g', valueCol, labelName, aid));
        plot(tSeg, xSeg, 'LineWidth', 1.2, 'Color', [0.15 0.15 0.15]);
        grid on;
        xlabel('Time (s, trimmed)');
        ylabel(valueCol);
        title(sprintf('%s agent_id=%g (middle ~10 s)', labelName, aid));
        tuneFigure;
        %saveFigure;

        %disp(table(nVals(:), contribNorm(:), 'VariableNames', {'n', sprintf('%s_agent_%g', labelName, aid)}));
    end
    fprintf('\n');
end

function [omegaVec, ampVec, tTrim, xTrim] = compute_fft_spectrum_from_vectors(t, x, trimRatio, csvPath, aid, nFftFactor)
    n = numel(x);
    if n < 20
        error('データ数が少なすぎます: %s (agent_id=%g)', csvPath, aid);
    end

    cut = floor(n * trimRatio);
    idx = (cut + 1):(n - cut);
    if isempty(idx) || numel(idx) < 10
        error('トリミング後のデータ数が不足: %s (agent_id=%g)', csvPath, aid);
    end
    tTrim = t(idx);
    xTrim = x(idx);

    % 10-90パーセンタイル幅で正規化
    xTrim = normalize_by_percentile_range(xTrim, 1, 99, csvPath, aid);

    % 基底関数の位相を安定化
    t = tTrim - tTrim(1);
    x = xTrim;

    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('時刻刻みが不正です: %s (agent_id=%g)', csvPath, aid);
    end

    n = numel(x);
    nfft = 2 ^ nextpow2(max(1, n * nFftFactor));
    X = fft(x, nfft);

    nHalf = floor(nfft / 2);
    Xh = X(1:(nHalf + 1));
    ampVec = abs(Xh) / n * 2;
    ampVec(1) = ampVec(1) / 2;
    if mod(nfft, 2) == 0
        ampVec(end) = ampVec(end) / 2;
    end

    fVec = (0:nHalf) / (nfft * dt);
    omegaVec = 2 * pi * fVec;
end

function xNorm = normalize_by_percentile_range(x, lowPct, highPct, csvPath, aid)
    if ~isvector(x)
        x = x(:);
    end
    x = double(x);
    x = x(isfinite(x));
    if numel(x) < 5
        warning('正規化に十分な点数がありません: %s (agent_id=%g)', csvPath, aid);
        xNorm = x;
        return;
    end

    pLow = prctile(x, lowPct);
    pHigh = prctile(x, highPct);
    denom = pHigh - pLow;
    center = 0.5 * (pLow + pHigh);
    if denom <= 0
        warning('正規化幅が0以下です: %s (agent_id=%g)', csvPath, aid);
        xNorm = x - center;
        return;
    end

    xNorm = (x - center) / denom;
end

function overlay_min_agent_spectra(fileIndex, trimRatio, omegaRange, valueCol, nFftFactor)
    springs = {'Spring1', 'Spring2', 'Spring3', 'Spring4', 'Spring5'};
    colors = lines(numel(springs));

    figure('Color', 'w', 'Name', sprintf('%s FFT overlay Spring1-5', valueCol));
    hold on;

    legendEntries = cell(1, numel(springs));
    for i = 1:numel(springs)
        springName = springs{i};
        dirPath = fullfile('SoftConnect', springName);
        csvPath = pick_nth_csv(dirPath, fileIndex);

        [omegaVec, ampVec, agentId] = min_agent_fft(csvPath, trimRatio, valueCol, nFftFactor);
        idxRange = (omegaVec >= omegaRange(1)) & (omegaVec <= omegaRange(2));

        plot(omegaVec(idxRange), ampVec(idxRange), 'LineWidth', 1.2, 'Color', colors(i, :));
        legendEntries{i} = sprintf('%s agent %g', springName, agentId);
    end

    grid on;
    xlabel('Angular frequency \omega (rad/s)');
    ylabel('Amplitude');
    xlim(omegaRange);
    title(sprintf('%s FFT overlay (min agent_id, %.4f..%.4f rad/s)', ...
        valueCol, omegaRange(1), omegaRange(2)));
    legend(legendEntries, 'Location', 'best');
    tuneFigure;
    hold off;
end

function overlay_max_agent_spectra(fileIndex, trimRatio, omegaRange, valueCol, nFftFactor)
    springs = {'Spring1', 'Spring2', 'Spring3', 'Spring4', 'Spring5'};
    colors = lines(numel(springs));

    figure('Color', 'w', 'Name', sprintf('%s FFT overlay Spring1-5 (max agent_id)', valueCol));
    hold on;

    legendEntries = cell(1, numel(springs));
    for i = 1:numel(springs)
        springName = springs{i};
        dirPath = fullfile('SoftConnect', springName);
        csvPath = pick_nth_csv(dirPath, fileIndex);

        [omegaVec, ampVec, agentId] = max_agent_fft(csvPath, trimRatio, valueCol, nFftFactor);
        idxRange = (omegaVec >= omegaRange(1)) & (omegaVec <= omegaRange(2));

        plot(omegaVec(idxRange), ampVec(idxRange), 'LineWidth', 1.2, 'Color', colors(i, :));
        legendEntries{i} = sprintf('%s agent %g', springName, agentId);
    end

    grid on;
    xlabel('Angular frequency \omega (rad/s)');
    ylabel('Amplitude');
    xlim(omegaRange);
    title(sprintf('%s FFT overlay (max agent_id, %.4f..%.4f rad/s)', ...
        valueCol, omegaRange(1), omegaRange(2)));
    legend(legendEntries, 'Location', 'best');
    tuneFigure;
    hold off;
end

function [omegaVec, ampVec, minAgentId] = min_agent_fft(csvPath, trimRatio, valueCol, nFftFactor)
    T = readtable(csvPath);

    if ~ismember('agent_id', T.Properties.VariableNames)
        error('agent_id列が見つかりません: %s', csvPath);
    end

    if ismember('time_local_sec', T.Properties.VariableNames)
        tAll = T.time_local_sec;
    elseif ismember('time_pc_sec_abs', T.Properties.VariableNames)
        tAll = T.time_pc_sec_abs;
    else
        error('time列が見つかりません: %s', csvPath);
    end

    if ~ismember(valueCol, T.Properties.VariableNames)
        error('%s列が見つかりません: %s', valueCol, csvPath);
    end
    xAll = T.(valueCol);
    agentAll = T.agent_id;

    tAll = double(tAll(:));
    xAll = double(xAll(:));
    agentAll = double(agentAll(:));

    validAll = isfinite(tAll) & isfinite(xAll) & isfinite(agentAll);
    tAll = tAll(validAll);
    xAll = xAll(validAll);
    agentAll = agentAll(validAll);

    minAgentId = min(agentAll);
    idxAgent = (agentAll == minAgentId);
    t = tAll(idxAgent);
    x = xAll(idxAgent);

    [t, order] = sort(t);
    x = x(order);

    [omegaVec, ampVec] = compute_fft_spectrum_from_vectors(t, x, trimRatio, csvPath, minAgentId, nFftFactor);
end

function [omegaVec, ampVec, maxAgentId] = max_agent_fft(csvPath, trimRatio, valueCol, nFftFactor)
    T = readtable(csvPath);

    if ~ismember('agent_id', T.Properties.VariableNames)
        error('agent_id列が見つかりません: %s', csvPath);
    end

    if ismember('time_local_sec', T.Properties.VariableNames)
        tAll = T.time_local_sec;
    elseif ismember('time_pc_sec_abs', T.Properties.VariableNames)
        tAll = T.time_pc_sec_abs;
    else
        error('time列が見つかりません: %s', csvPath);
    end

    if ~ismember(valueCol, T.Properties.VariableNames)
        error('%s列が見つかりません: %s', valueCol, csvPath);
    end
    xAll = T.(valueCol);
    agentAll = T.agent_id;

    tAll = double(tAll(:));
    xAll = double(xAll(:));
    agentAll = double(agentAll(:));

    validAll = isfinite(tAll) & isfinite(xAll) & isfinite(agentAll);
    tAll = tAll(validAll);
    xAll = xAll(validAll);
    agentAll = agentAll(validAll);

    maxAgentId = max(agentAll);
    idxAgent = (agentAll == maxAgentId);
    t = tAll(idxAgent);
    x = xAll(idxAgent);

    [t, order] = sort(t);
    x = x(order);

    [omegaVec, ampVec] = compute_fft_spectrum_from_vectors(t, x, trimRatio, csvPath, maxAgentId, nFftFactor);
end

function [tSeg, xSeg] = extract_middle_10sec(tTrim, xTrim)
    halfWindow = 5.0;

    tCenter = 0.5 * (tTrim(1) + tTrim(end));
    idxSeg = (tTrim >= (tCenter - halfWindow)) & (tTrim <= (tCenter + halfWindow));

    if sum(idxSeg) < 10
        warning('中央10秒区間の点数が少ないため、利用可能な全区間を表示します。');
        tSeg = tTrim;
        xSeg = xTrim;
        return;
    end

    tSeg = tTrim(idxSeg);
    xSeg = xTrim(idxSeg);
end
