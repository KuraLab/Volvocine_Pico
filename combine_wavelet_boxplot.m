function combine_wavelet_boxplot(dates, fmin, fmax)
% 例: combine_wavelet_boxplot({'2025-08-15','2025-08-18'}, 1.0, 1.5)
if nargin<1 || isempty(dates), dates = {'2025-08-15','2025-08-18'}; end
if nargin<2, fmin = 1.0; end
if nargin<3, fmax = 1.5; end

allT = [];
for i=1:numel(dates)
    d = dates{i};
    csvPath = fullfile('results', d, sprintf('wavelet_%s_%g_%gHz_per_file.csv', d, fmin, fmax));
    if ~isfile(csvPath)
        fprintf('[WARN] missing per_file csv: %s\n', csvPath);
        continue;
    end
    T = readtable(csvPath, 'VariableNamingRule','preserve');
    T.date = repmat(string(d), height(T),1);
    allT = [allT; T]; %#ok<AGROW>
end
if isempty(allT)
    fprintf('[ERROR] no data loaded\n');
    return;
end

% (N,u) 割当
N = nan(height(allT),1);
U = nan(height(allT),1);
for i=1:height(allT)
    d = string(allT.date(i));
    s = string(allT.session(i));
    if d=="2025-08-15"
        switch s
            case "S2", N(i)=3; U(i)=0;
            case "S3", N(i)=4; U(i)=0;
            case "S4", N(i)=5; U(i)=0;
            case "S5", N(i)=6; U(i)=0;
            case "S6", N(i)=3; U(i)=1;
        end
    elseif d=="2025-08-18"
        switch s
            case "S2", N(i)=4; U(i)=1;
            case "S3", N(i)=5; U(i)=1;
            case "S4", N(i)=6; U(i)=1;
        end
    end
end
allT.N = N; allT.u = U;

% 有効データのみ
M = allT(~isnan(allT.N) & ~isnan(allT.u), :);
if isempty(M)
    fprintf('[ERROR] no rows after mapping\n');
    return;
end

% ==== 追加: 統計出力用準備 =========================================
outDir = fullfile('results','combined_stats');
if ~exist(outDir,'dir'), mkdir(outDir); end

targetNs = [3 4 5 6];
uLevels = [0 1];
sep = 0.18; widthBox = 0.30;

featList = {'amp_a0','Amplitude a0','box_full_amp_a0'; ...
            'amp_a1','Amplitude a1','box_full_amp_a1'};

% 追加: 集計テーブル用セル
groupStatRows = {};  % 各 (feature, N, u)
pairStatRows  = {};  % 各 (feature, N) の u=0 vs u=1

% ===== 各特徴について事前に group / pair 統計計算 =================
for fIdx_pre = 1:size(featList,1)
    col = featList{fIdx_pre,1};
    if ~ismember(col, M.Properties.VariableNames) || all(isnan(M.(col)))
        continue;
    end

    % --- group level (N,u)
    for nVal = targetNs
        for uVal = uLevels
            rows = M.N==nVal & M.u==uVal;
            v = M.(col)(rows);
            v = v(~isnan(v));
            if isempty(v), continue; end
            q1 = quantile(v,0.25);
            q3 = quantile(v,0.75);
            medv = median(v);
            iqrV = q3 - q1;
            mnv = mean(v);
            sdv = std(v,0);
            sev = sdv / sqrt(numel(v));
            groupStatRows(end+1,:) = {col, nVal, uVal, numel(v), medv, q1, q3, iqrV, mnv, sdv, sev}; %#ok<AGROW>
        end
    end

    % --- pair level (N 内で u=0 vs 1)
    for nVal = targetNs
        v0 = M.(col)(M.N==nVal & M.u==0); v0 = v0(~isnan(v0));
        v1 = M.(col)(M.N==nVal & M.u==1); v1 = v1(~isnan(v1));
        if numel(v0)>1 && numel(v1)>1
            % すべて Mann-Whitney U (ranksum)
            p = ranksum(v0, v1);
            testType   = 'ranksum';
            cliffs     = cliffs_delta(v1, v0);   % (u=1) - (u=0)
            effectType = 'cliffs_delta';
            effectVal  = cliffs;
            d = NaN; % cohen's d は計算しない (非パラメトリック統一)
            med0 = median(v0); med1 = median(v1);
            mean0= mean(v0);   mean1= mean(v1);
            q1_0 = quantile(v0,0.25); q3_0=quantile(v0,0.75);
            q1_1 = quantile(v1,0.25); q3_1=quantile(v1,0.75);
            redPctMedian = (med0 - med1)/med0 * 100;
            redPctMean   = (mean0 - mean1)/mean0 * 100;
            pairStatRows(end+1,:) = {col, nVal, numel(v0), numel(v1), ...
                med0, med1, q1_0, q3_0, q1_1, q3_1, ...
                mean0, mean1, p, testType, effectType, effectVal, d, cliffs, ...
                redPctMedian, redPctMean}; %#ok<AGROW>
        end
    end
end

% 追加: テーブル化 & CSV 書き出し (特徴ごとまとめて)
if ~isempty(groupStatRows)
    groupTbl = cell2table(groupStatRows, 'VariableNames', ...
        {'feature','N','u','n','median','q1','q3','iqr','mean','std','se'});
    outPathGroup = fullfile(outDir, sprintf('group_stats_%g_%gHz.csv', fmin, fmax));
    writetable(groupTbl, outPathGroup);
    fprintf('[INFO] wrote group stats: %s\n', outPathGroup);
else
    fprintf('[WARN] no group stats rows\n');
end

if ~isempty(pairStatRows)
    pairTbl = cell2table(pairStatRows, 'VariableNames', ...
        {'feature','N','n_u0','n_u1','median_u0','median_u1','q1_u0','q3_u0','q1_u1','q3_u1', ...
         'mean_u0','mean_u1','p_value','test','effect_type','effect_value','cohen_d','cliffs_delta', ...
         'reduction_pct_median','reduction_pct_mean'});
    outPathPair = fullfile(outDir, sprintf('pair_stats_%g_%gHz.csv', fmin, fmax));
    writetable(pairTbl, outPathPair);
    fprintf('[INFO] wrote pair stats: %s\n', outPathPair);
else
    fprintf('[WARN] no pair stats rows\n');
end
% ==================================================================

% ==== 有意差マーカー相対配置パラメータ（共通） ====
relGap = 0.05;       % 各ペア最大値の何倍上から括弧開始
barHeight = 0.04;    % 括弧高さ (pairMax 比)
starOffset = 0.05;   % 星位置 (括弧上 + pairMax*starOffset)
safetyMargin = 0.15; % 全体上限マージン (共通 y 上限に一度だけ適用)

% --- 第1パス: 全特徴で必要な最大高さを計算 (マージン未適用) ---
featInfo = struct('col',{},'pairMax',{},'rawTopNoMargin',{});
for fIdx = 1:size(featList,1)
    col = featList{fIdx,1};
    if ~ismember(col, M.Properties.VariableNames) || all(isnan(M.(col)))
        continue;
    end
    globalYmax = max(M.(col),[],'omitnan');
    if ~isfinite(globalYmax) || globalYmax<=0, globalYmax=1; end
    pairMax = zeros(numel(targetNs),1);
    for nIdx=1:numel(targetNs)
        v0 = M.(col)(M.N==targetNs(nIdx) & M.u==0);
        v1 = M.(col)(M.N==targetNs(nIdx) & M.u==1);
        if ~isempty(v0) && ~isempty(v1)
            pairMax(nIdx) = max([v0; v1], [], 'omitnan');
        else
            pairMax(nIdx) = 0;
        end
    end
    neededTop = 0;
    for nIdx=1:numel(targetNs)
        pm = pairMax(nIdx);
        if pm<=0, continue; end
        topThis = pm * (1 + relGap + barHeight + starOffset);
        if topThis > neededTop, neededTop = topThis; end
    end
    rawTopNoMargin = max(globalYmax, neededTop); % ← safetyMargin 未適用
    featInfo(end+1) = struct('col',col,'pairMax',pairMax,'rawTopNoMargin',rawTopNoMargin); %#ok<AGROW>
end
if isempty(featInfo)
    fprintf('[WARN] no valid features for plotting\n');
    return;
end
% 共通 y 上限 (全特徴で最大を取り safetyMargin を一度適用)
yUpperShared = max([featInfo.rawTopNoMargin]) * (1 + safetyMargin);

% --- 第2パス: 描画 (共通 yUpperShared 使用) ---
for fIdx = 1:numel(featInfo)
    col = featInfo(fIdx).col;
    pairMax = featInfo(fIdx).pairMax;

    % データ再構築 (ボックスプロット用)
    valsAll = []; grpAll = []; positions = []; gCounter = 0;
    for nIdx=1:numel(targetNs)
        nVal = targetNs(nIdx);
        for ui=1:numel(uLevels)
            uVal = uLevels(ui);
            rows = M.N==nVal & M.u==uVal;
            if ~any(rows), continue; end
            v = M.(col)(rows);
            gCounter = gCounter + 1;
            valsAll = [valsAll; v];
            grpAll  = [grpAll; repmat(gCounter, numel(v),1)];
            pos = nIdx + (uVal==0)*(-sep) + (uVal==1)*(sep);
            positions = [positions pos];
        end
    end
    if isempty(valsAll), continue; end

    yUpper = yUpperShared; % 共通上限

    figure('Units','centimeters','Position',[2 2 16 10.5],'Color','w');
    ax = gca; ax.FontName='Times'; ax.FontSize=25; ax.LineWidth=3;
    boxplot(valsAll, grpAll, 'Positions', positions, 'Widths', widthBox, ...
        'Symbol','o','Whisker',1.5);
    grid on; set(gcf,'Color','w');
    xlim([0.5 numel(targetNs)+0.5]);
    ylim([0 yUpper]);
    xticks(1:numel(targetNs));
    xticklabels(compose('%d', targetNs));
    xlabel('$N$','Interpreter','latex','FontSize',25,'FontName','Times');
    ylabel({'Mean Amplitude', sprintf('(%.1f-%.1f Hz)', fmin, fmax)}, ...
        'FontSize',25,'FontName','Times');

    set(findobj(ax,'Tag','Box'),'LineWidth',2.0);
    set(findobj(ax,'Tag','Median'),'LineWidth',2.0);
    set(findobj(ax,'Tag','Whisker'),'LineWidth',2.8);
    set(findobj(ax,'Tag','Outliers'),'MarkerSize',4);

    boxes = findobj(ax,'Tag','Box');
    boxPos = arrayfun(@(h) mean(get(h,'XData')), boxes);
    [~,ord]=sort(boxPos); boxes=boxes(ord);
    for bi=1:numel(boxes)
        xd = mean(get(boxes(bi),'XData'));
        nIdxNear = round(xd);
        if xd < nIdxNear
            set(boxes(bi),'Color',[0.2 0.4 0.85]); % u=0
        else
            set(boxes(bi),'Color',[0.9 0.35 0.15]); % u=1
        end
    end
    meds = findobj(ax,'Tag','Median');
    medPos = arrayfun(@(h) mean(get(h,'XData')), meds);
    [~,mord]=sort(medPos); meds=meds(mord);
    for bi=1:min(numel(meds), numel(boxes))
        set(meds(bi),'Color', get(boxes(bi),'Color'));
    end

    % 有意差 (u=0 vs u=1 同一N) を ranksum のみで描画
    for nIdx=1:numel(targetNs)
        nVal = targetNs(nIdx);
        v0 = M.(col)(M.N==nVal & M.u==0);
        v1 = M.(col)(M.N==nVal & M.u==1);
        if numel(v0)>1 && numel(v1)>1
            p = ranksum(v0,v1);
            if p<1e-2
                star='**';
            elseif p<0.05
                star='*';
            else
                star='n.s.';
            end
            pm = pairMax(nIdx);
            if pm<=0, continue; end
            yBasePair = pm * (1 + relGap);
            yTopPair  = yBasePair + pm * barHeight;
            starY     = yTopPair + pm * starOffset;
            x0 = nIdx - sep; x1 = nIdx + sep;
            hold on;
            plot([x0 x0 x1 x1],[yBasePair yTopPair yTopPair yBasePair],'k','LineWidth',1.3);
            text(nIdx, starY, star, 'HorizontalAlignment','center', ...
                'FontWeight','bold','FontSize',16);
        end
    end
    ax.FontSize=25;
    % 凡例
    % 既存 boxplot 線を凡例から除外
    set(findobj(ax,'Type','Line'),'HandleVisibility','off');
    % ダミー (NaN,NaN) で凡例用ハンドル作成
    hU0 = plot(nan,nan,'s','MarkerFaceColor',[0.2 0.4 0.85], ...
        'MarkerEdgeColor','none','MarkerSize',10,'DisplayName','Baseline ($$u = 0$$)','HandleVisibility','on');
    hU1 = plot(nan,nan,'s','MarkerFaceColor',[0.9 0.35 0.15], ...
        'MarkerEdgeColor','none','MarkerSize',10,'DisplayName','Proposed','HandleVisibility','on');
    legend([hU0 hU1],'Location','northeast', 'Interpreter', 'latex','FontSize',20);
    set(gca,'TickLabelInterpreter','latex');
    
    % 余白圧縮
    ti = ax.TightInset;
    ax.Position = [ti(1)+0.01 ti(2)+0.01 1 - ti(1) - ti(3) - 0.02 1 - ti(2) - ti(4) - 0.02];
    %saveFigure;
end
end

% ==== 追加: Cliff's delta 計算用ローカル関数 =======================
function d = cliffs_delta(b,a)
% d > 0 なら b が大きい傾向 (b = u=1, a = u=0)
a = a(:); b = b(:);
if isempty(a) || isempty(b)
    d = NaN; return;
end
% 効率化: ソートして累積 (O((n+m) log(n+m)))
aS = sort(a);
d_count = 0;
for bi = 1:numel(b)
    % b(bi) と a の大小関係
    d_count = d_count + sum(aS < b(bi)) - sum(aS > b(bi));
end
d = d_count / (numel(a)*numel(b));
end