function wavelet_session_power(dateStr, fmin, fmax)
% wavelet_session_power('2025-08-15', 1.0, 2.0)
% 依存: filtered_files_<date>.csv (list_merged_files.py が生成)
% 出力:
%   wavelet_<date>_1_2Hz_per_file.csv
%   wavelet_<date>_1_2Hz_agg.csv

if nargin < 1 || isempty(dateStr), dateStr = '2025-08-15'; end
if nargin < 2 || isempty(fmin),    fmin = 1.0; end
if nargin < 3 || isempty(fmax),    fmax = 1.5; end

requiredAgent = 99;
winStart = 25.0;  % 相対秒
winEnd   = 30.0;
minPoints = 5;
valueCols = {'a0','a1'};
resultsDir = fullfile('results', dateStr);
if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
keepListFile = fullfile(resultsDir, sprintf('filtered_files_%s.csv', dateStr));
legacyKeep = sprintf('filtered_files_%s.csv', dateStr); % 後方互換 (ルート)
baseDir = fullfile('merged_chunks_organized', dateStr);

% --- セッション境界 (Python と同じ)
sessionStarts = ["16:40","20:00","20:49","21:47","22:22"];
sessionIntervals = buildSessions(dateStr, sessionStarts);

% --- Keep リスト読み込み or 全ファイル
if isfile(keepListFile)
    Tkeep = readtable(keepListFile, 'VariableNamingRule','preserve');
    keepNames = string(Tkeep.filename);
    sessionMap = containers.Map(keepNames, string(Tkeep.session));
    fprintf('[INFO] keep list loaded: %d files\n', numel(keepNames));
elseif isfile(legacyKeep)
    Tkeep = readtable(legacyKeep, 'VariableNamingRule','preserve');
    keepNames = string(Tkeep.filename);
    sessionMap = containers.Map(keepNames, string(Tkeep.session));
    fprintf('[INFO] legacy keep list loaded: %d files\n', numel(keepNames));
else
    d = dir(fullfile(baseDir, 'merged_*.csv'));
    keepNames = string({d.name});
    sessionMap = containers.Map;
    for i=1:numel(keepNames)
        nm = keepNames(i);
        dt = parseFilenameDT(nm);
        sess = assignSession(dt, sessionIntervals);
        sessionMap(nm) = sess;
    end
    fprintf('[INFO] keep list not found -> using all %d files\n', numel(keepNames));
end

results = [];  % struct 配列
rCount = 0;

for i=1:numel(keepNames)
    fname = keepNames(i);
    fpath = fullfile(baseDir, fname);
    if ~isfile(fpath)
        fprintf('[SKIP missing file] %s\\n', fname);
        continue;
    end
    % 読み込み
    try
        T = readtable(fpath, 'VariableNamingRule','preserve');
    catch ME
        fprintf('[SKIP read error] %s (%s)\\n', fname, ME.message);
        continue;
    end
    neededCols = [{'time_local_sec','agent_id'}, valueCols];
    if ~all(ismember(neededCols, T.Properties.VariableNames))
        fprintf('[SKIP missing cols] %s\\n', fname);
        continue;
    end
    % agent フィルタ
    agentMask = T.agent_id == requiredAgent;
    if ~any(agentMask)
        % スキップ (Python 同様)
        continue;
    end
    Ta = T(agentMask, :);
    t_abs = Ta.time_local_sec;
    t0 = min(t_abs);
    t_rel = t_abs - t0;
    winMask = t_rel >= winStart & t_rel <= winEnd;
    if ~any(winMask)
        % ウィンドウ内データ無し -> スキップ
        continue;
    end
    Tw = Ta(winMask, :);
    if height(Tw) < minPoints
        continue;
    end
    % サンプリング間隔推定
    t = t_rel(winMask);
    dt_vec = diff(t);
    dt = median(dt_vec);
    if ~isfinite(dt) || dt <= 0
        continue;
    end
    fs = 1/dt;

    % CWT 計算 (MATLAB R2023 以降 cwt(x,fs) -> wt: freq x time)
    % Morlet がデフォ (analytic Morse になる場合は 'amor' 指定でも可)
    % 帯域抽出用にフル計算後フィルタ
    for vc = 1:numel(valueCols)
        col = valueCols{vc};
        if ~isnumeric(Tw.(col))
            y = double(table2array(Tw(:,col)));
        else
            y = double(Tw.(col));
        end
        if numel(y) < minPoints, continue; end
        try
            [wt,f] = cwt(y, fs);  % wt: numF x numT
        catch ME
            fprintf('[SKIP cwt error] %s %s (%s)\\n', fname, col, ME.message);
            continue;
        end
        bandIdx = f >= fmin & f <= fmax;
        if ~any(bandIdx), continue; end
        wtBand = wt(bandIdx, :);
        fBand = f(bandIdx);

        % --- 時間ベクトル昇順保証 ---
        [tSorted, tidx] = sort(t);
        wtBand = wtBand(:, tidx);   % 時間軸並び替え

        absC   = abs(wtBand);
        powerC = absC.^2;

        % 時間方向積分 (各周波数ごと)
        p_time = trapz(tSorted, powerC, 2);
        a_time = trapz(tSorted, absC,   2);

        % 周波数軸が降順なら反転
        if fBand(1) > fBand(end)
            fAsc  = flip(fBand);
            p_time = flip(p_time);
            a_time = flip(a_time);
        else
            fAsc = fBand;
        end

        % 周波数方向積分（昇順周波数で正の値）
        p_total = trapz(fAsc, p_time);
        a_total = trapz(fAsc, a_time);

        % 数値誤差で僅かに負の場合ゼロクリップ
        if p_total < 0 && p_total > -1e-9, p_total = 0; end
        if a_total < 0 && a_total > -1e-9, a_total = 0; end

        % ===== ここから平均化 (Option 3) へ変更 =====
        windowLen = tSorted(end) - tSorted(1);           % Δt
        deltaF = fAsc(end) - fAsc(1);                    % Δf
        if deltaF > 0 && windowLen > 0
            a_total_mean = a_total / (deltaF * windowLen);
        else
            a_total_mean = NaN; % 周波数1点や窓長0の異常ケース
        end
        % ===== ここまで =====

        if vc == 1
            rCount = rCount + 1;
            results(rCount).file = fname;
            if isKey(sessionMap, fname)
                sess = sessionMap(fname);
            else
                sess = assignSession(parseFilenameDT(fname), sessionIntervals);
            end
            results(rCount).session = sess;
        end
        results(rCount).(['power_' col]) = p_total;
        % 平均振幅に置換
        results(rCount).(['amp_' col])   = a_total_mean;
    end
end

if isempty(results)
    fprintf('[WARN] no results after processing\n');
    perFileTable = table();
    aggTable = table();
else
    perFileTable = struct2table(results);
    % 欠損列補完 (存在しない場合 NaN 列追加)
    allColsNeeded = [{'file','session'}, strcat('power_',valueCols), strcat('amp_',valueCols)];
    for c = 1:numel(allColsNeeded)
        if ~ismember(allColsNeeded{c}, perFileTable.Properties.VariableNames)
            perFileTable.(allColsNeeded{c}) = NaN(height(perFileTable),1);
        end
    end
    % 集計
    aggTable = aggregateSessions(perFileTable, valueCols);
end

outBase = fullfile(resultsDir, sprintf('wavelet_%s_%g_%gHz', dateStr, fmin, fmax));
perFileCSV = sprintf('%s_per_file.csv', outBase);
aggCSV = sprintf('%s_agg.csv', outBase);
if ~isempty(perFileTable)
    writetable(perFileTable, perFileCSV);
end
if ~isempty(aggTable)
    writetable(aggTable, aggCSV);
end
fprintf('[INFO] per-file saved: %s rows=%d\n', perFileCSV, height(perFileTable));
fprintf('[INFO] aggregated saved: %s rows=%d\n', aggCSV, height(aggTable));

% --- N,u グルーピング箱ひげ図 (S1除外, ampのみ) ---
% セッション→(N,u) マップ定義
sess2Nu = struct();
sess2Nu.S2 = [3 0];
sess2Nu.S3 = [4 0];
sess2Nu.S4 = [5 0];
sess2Nu.S5 = [6 0];
sess2Nu.S6 = [3 1];
% S1 は除外
if ~isempty(perFileTable)
    M = perFileTable;
    M = M(~strcmp(M.session,'S1'), :);
    if ~isempty(M)
        % N,u 列付与
        Ns = nan(height(M),1);
        Us = nan(height(M),1);
        for i=1:height(M)
            s = string(M.session(i));
            if isfield(sess2Nu, s)
                nu = sess2Nu.(s);
                Ns(i) = nu(1);
                Us(i) = nu(2);
            end
        end
        M.N = Ns;
        M.u = Us;
        % 有効行のみ
        M = M(~isnan(M.N) & ~isnan(M.u), :);

        if ~isempty(M)
            targetNs = [3 4 5 6];   % 並べる N
            uLevels  = [0 1];       % 左:0 右:1
            sep = 0.18;             % u=0/1 の左右オフセット
            widthBox = 0.30;
            featList = {
                'amp_a0','Amplitude a0','amp_Nu_boxplot_a0'
                'amp_a1','Amplitude a1','amp_Nu_boxplot_a1'
            };

            for fIdx = 1:size(featList,1)
                col = featList{fIdx,1};
                fnameBase = featList{fIdx,3};
                if ~ismember(col, M.Properties.VariableNames), continue; end
                % 値存在チェック
                if all(isnan(M.(col))), continue; end

                % データ収集
                allVals = [];
                groupIdx = [];
                positions = [];
                labelsNU = []; % 対応 (N,u) 記録
                gCounter = 0;
                for nIdx = 1:numel(targetNs)
                    nVal = targetNs(nIdx);
                    for ui = 1:numel(uLevels)
                        uVal = uLevels(ui);
                        rows = M.N==nVal & M.u==uVal;
                        if ~any(rows)
                            % 空 -> 単にスキップ (隙間保持のため位置だけ確保したい場合はダミーを入れないと詰まるが
                            % boxplot の positions 指定でスキップすると空白が残るので、ここは「追加しない」で OK
                            continue;
                        end
                        gCounter = gCounter + 1;
                        v = M.(col)(rows);
                        allVals = [allVals; v];
                        groupIdx = [groupIdx; repmat(gCounter, numel(v),1)];
                        pos = nIdx + (uVal==0)*(-sep) + (uVal==1)*(sep);
                        positions = [positions pos];
                        labelsNU = [labelsNU; nVal uVal];
                    end
                end

                if isempty(allVals)
                    continue;
                end

                % 両チャンネルで共通の Y 上限を計算
                rawA0 = M.amp_a0; rawA1 = M.amp_a1;
                globalYmax = max([rawA0; rawA1], [], 'omitnan');
                if isempty(globalYmax) || ~isfinite(globalYmax) || globalYmax <= 0
                    globalYmax = 1;
                end
                yMarginFactor = 0.15;                 % 上側余白割合
                yUpper = globalYmax * (1 + yMarginFactor);

                % 図を保存せず表示するため Visible をデフォルト (on) に
                fig = figure('Units','centimeters');

                fig.Position = [2 2 16 10.5];
                set(gca, 'FontName', 'Times New Roman');

                boxplot(allVals, groupIdx, 'Positions', positions, 'Widths', widthBox, ...
                    'Symbol','o','Whisker',1.5);

                ax = gca;
                ax.FontSize = 25;
                ax.LineWidth = 1.2;
                grid on; set(gcf,'Color','w');
                xlim([0.5 numel(targetNs)+0.5]);
                ylim([0 yUpper]);                         % 共通 Y 軸
                xticks(1:numel(targetNs));
                xticklabels(compose('%d', targetNs));
                xlabel('$N$','FontSize',25,'Interpreter', 'latex', 'FontName', 'Times');
                ylabel({'Mean Amplitude', ...
                    sprintf('(%.1f-%.1f Hz, %.0f-%.0f s)', fmin, fmax, winStart, winEnd)}, ...
                    'FontSize',25,'FontName', 'Times');


                % 線太さ調整
                set(findobj(ax,'Tag','Box'),'LineWidth',2.0);
                set(findobj(ax,'Tag','Median'),'LineWidth',2.0);
                set(findobj(ax,'Tag','Whisker'),'LineWidth',2.8);
                set(findobj(ax,'Tag','Outliers'),'MarkerSize',4);

                % 色付け (u=0 左 / u=1 右)
                boxes = findobj(ax,'Tag','Box');
                boxPos = arrayfun(@(h) mean(get(h,'XData')), boxes);
                [~,ord] = sort(boxPos);
                boxes = boxes(ord);
                for bi=1:numel(boxes)
                    xd = mean(get(boxes(bi),'XData'));
                    nIdxNear = round(xd);
                    if abs(xd - nIdxNear) < 0.5
                        if xd < nIdxNear
                            set(boxes(bi),'Color',[0.2 0.4 0.85]);  % u=0
                        else
                            set(boxes(bi),'Color',[0.9 0.35 0.15]); % u=1
                        end
                    end
                end
                % ---- 追加: 中央値線を赤にせず箱と同じ色へ ----
                meds = findobj(ax,'Tag','Median');
                medPos = arrayfun(@(h) mean(get(h,'XData')), meds);
                [~,mord] = sort(medPos);
                meds = meds(mord);
                for bi=1:min(numel(meds), numel(boxes))
                    set(meds(bi),'Color', get(boxes(bi),'Color'));
                end
                % -----------------------------------------------
                hold on;

                % 有意差ライン位置 (上部 12% 以内)
                yBase = yUpper * 0.88;
                yH    = (yUpper - yBase) * 0.35;
                starYOffset = (yUpper - yBase) * 0.55;

                for nIdx = 1:numel(targetNs)
                    nVal = targetNs(nIdx);
                    mask0 = M.N==nVal & M.u==0;
                    mask1 = M.N==nVal & M.u==1;
                    if ~any(mask0) || ~any(mask1), continue; end
                    v0 = M.(col)(mask0); v1 = M.(col)(mask1);
                    if numel(v0) > 1 && numel(v1) > 1
                        [~,p] = ttest2(v0, v1,'Vartype','unequal');
                        star = sigStar(p);
                        x0 = nIdx - sep; x1 = nIdx + sep;
                        plot([x0 x0 x1 x1],[yBase yBase+yH yBase+yH yBase],'k','LineWidth',1.5);
                        text(nIdx, yBase + yH*0.5 + starYOffset, star, ...
                            'HorizontalAlignment','center','FontWeight','bold','FontSize',16);
                    end
                end

                % 凡例
                % 既存 boxplot 線を凡例から除外
                set(findobj(ax,'Type','Line'),'HandleVisibility','off');
                % ダミー (NaN,NaN) で凡例用ハンドル作成
                hU0 = plot(nan,nan,'s','MarkerFaceColor',[0.2 0.4 0.85], ...
                    'MarkerEdgeColor','none','MarkerSize',10,'DisplayName','Baseline (u = 0)','HandleVisibility','on');
                hU1 = plot(nan,nan,'s','MarkerFaceColor',[0.9 0.35 0.15], ...
                    'MarkerEdgeColor','none','MarkerSize',10,'DisplayName','Proposed','HandleVisibility','on');
                legend([hU0 hU1],'Location','northoutside', 'Interpreter', 'latex');
                set(gca,'TickLabelInterpreter','latex');

                % 余白圧縮
                ti = ax.TightInset;
                ax.Position = [ti(1)+0.01 ti(2)+0.01 1 - ti(1) - ti(3) - 0.02 1 - ti(2) - ti(4) - 0.02];

                % 保存しない。必要なら手動で [File] -> [Save] してください。
                drawnow;  % 即時描画
                fprintf('[INFO] showed grouped boxplot %s (no file saved)\n', fnameBase);
                saveFigure;
            end
        end
    end
end
end  % ← 既存 end の直前に入れている場合は重複に注意

% ---------- 補助関数群 ----------

function intervals = buildSessions(dateStr, starts)
    % returns struct array with label,startTime,endTime (datetime or [])
    times = datetime.empty;
    for i=1:numel(starts)
        t = datetime(sprintf('%s %s:00', dateStr, starts(i)), 'InputFormat','yyyy-MM-dd HH:mm:ss');
        times(end+1) = t; %#ok<AGROW>
    end
    times = sort(times);
    intervals = struct('label', {}, 'start', {}, 'end', {});
    intervals(1).label = 'S1'; intervals(1).start = []; intervals(1).end = times(1);
    for k=1:numel(times)-1
        intervals(end+1).label = sprintf('S%d', k+1); %#ok<AGROW>
        intervals(end).start = times(k);
        intervals(end).end   = times(k+1);
    end
    intervals(end+1).label = sprintf('S%d', numel(times)+1);
    intervals(end).start = times(end);
    intervals(end).end = [];
end

function sess = assignSession(dt, intervals)
    if isempty(dt)
        sess = '';
        return;
    end
    sess = '';
    for k=1:numel(intervals)
        s = intervals(k).start;
        e = intervals(k).end;
        if isempty(s) && ~isempty(e)
            if dt < e, sess = intervals(k).label; return; end
        elseif ~isempty(s) && isempty(e)
            if dt >= s, sess = intervals(k).label; return; end
        else
            if dt >= s && dt < e, sess = intervals(k).label; return; end
        end
    end
end

function dt = parseFilenameDT(name)
    % name: merged_YYYYMMDD_HHMMSS.csv
    expr = 'merged_(\\d{4})(\\d{2})(\\d{2})_(\\d{2})(\\d{2})(\\d{2})\\.csv';
    tokens = regexp(name, expr, 'tokens', 'once');
    if isempty(tokens)
        dt = [];
        return;
    end
    y = str2double(tokens{1}); mo = str2double(tokens{2}); d = str2double(tokens{3});
    H = str2double(tokens{4}); M = str2double(tokens{5}); S = str2double(tokens{6});
    dt = datetime(y,mo,d,H,M,S);
end

function aggTable = aggregateSessions(perFileTable, valueCols)
    % features: power_a0, amp_a0, power_a1, amp_a1...
    sessions = unique(perFileTable.session);
    rows = [];
    for s = 1:numel(sessions)
        sess = sessions(s);
        mask = strcmp(perFileTable.session, sess);
        if ~any(mask), continue; end
        for vc = 1:numel(valueCols)
            base = valueCols{vc};
            feats = {['power_' base], ['amp_' base]};
            for f = 1:numel(feats)
                colName = feats{f};
                if ~ismember(colName, perFileTable.Properties.VariableNames), continue; end
                data = perFileTable.(colName)(mask);
                data = data(~isnan(data));
                if isempty(data), continue; end
                r.session = sess;
                r.feature = colName;
                r.count = numel(data);
                r.mean = mean(data);
                r.var = var(data, 1);   % Python の pandas var(ddof=1) なら var(data,0) に変更
                r.med = median(data);
                rows = [rows; r]; %#ok<AGROW>
            end
        end
    end
    if isempty(rows)
        aggTable = table();
    else
        aggTable = struct2table(rows);
        % pandas の var はデフォルト ddof=1 (母不偏) -> MATLAB var(data,0)。
        % 上で var(data,1) にしているので不偏に合わせるなら var(data,0) に変更。
        % 必要ならここを調整してください。
    end
end

function star = sigStar(p)
%if p < 1e-4
%    star = '****';
%elseif p < 1e-3
%    star = '***';
if p < 1e-2
    star = '**';
elseif p < 0.05
    star = '*';
else
    star = 'n.s.';
end
end