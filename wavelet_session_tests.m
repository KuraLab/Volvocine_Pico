function wavelet_session_tests(dateStr, fmin, fmax)
% wavelet_session_tests('2025-08-15', 1.0, 2.0)
% 前提: wavelet_session_power を実行済みで results/<date>/ に per_file CSV があること
if nargin<1 || isempty(dateStr), dateStr = "2025-08-15"; end
if nargin<2, fmin = []; end
if nargin<3, fmax = []; end

resDir = fullfile("results", dateStr);
if ~isfolder(resDir)
    error("results ディレクトリが存在しません: %s (先に wavelet_session_power を実行してください)", resDir);
end

% 期待ファイル名（fmin/fmax が与えられた場合）
target = "";
if ~isempty(fmin) && ~isempty(fmax)
    % Python 側命名: wavelet_<date>_<fmin>_<fmax>Hz_per_file.csv （小数はそのまま出力している想定）
    % 小数点が含まれる場合は '_' に置換されているかもしれないので両方試す
    pat1 = sprintf("wavelet_%s_%g_%gHz_per_file.csv", dateStr, fmin, fmax);
    pat2 = strrep(pat1, ".", "_"); % 念のため
    c1 = dir(fullfile(resDir, pat1));
    c2 = dir(fullfile(resDir, pat2));
    cand = [c1; c2];
    if ~isempty(cand)
        target = fullfile(resDir, cand(1).name);
    end
end

% 汎用探索（ターゲット未決定時）
if target == ""
    listing = dir(fullfile(resDir, sprintf("wavelet_%s_*_Hz_per_file.csv", dateStr)));
    if isempty(listing)
        % さらに緩いパターン
        listing = dir(fullfile(resDir, sprintf("wavelet_%s*_per_file.csv", dateStr)));
    end
    if isempty(listing)
        % 代替: 利用可能ファイル列挙
        allFiles = dir(resDir);
        fprintf("利用可能ファイル一覧:\n");
        for k=1:numel(allFiles)
            if ~allFiles(k).isdir
                fprintf("  %s\n", allFiles(k).name);
            end
        end
        error("per-file CSV 見つからず (wavelet_session_power を実行したか確認)");
    end
    % 更新日時で最新
    [~, idx] = max([listing.datenum]);
    target = fullfile(resDir, listing(idx).name);
end

fprintf("[INFO] per-file CSV 使用: %s\n", target);
T = readtable(target, 'VariableNamingRule','preserve');

need = {'session','power_a0','power_a1','amp_a0','amp_a1'};
if ~all(ismember(need, T.Properties.VariableNames))
    error('必要列不足: %s', strjoin(setdiff(need, T.Properties.VariableNames), ','));
end

S2 = strcmp(T.session,'S2');
S6 = strcmp(T.session,'S6');

lines = {};
lines{end+1} = welchLine("power a0 S2 vs S6", T.power_a0(S2), T.power_a0(S6));
lines{end+1} = welchLine("amp a0 S2 vs S6",   T.amp_a0(S2),   T.amp_a0(S6));
lines{end+1} = welchLine("power a1 S2 vs S6", T.power_a1(S2), T.power_a1(S6));
lines{end+1} = welchLine("amp a1 S2 vs S6",   T.amp_a1(S2),   T.amp_a1(S6));
lines{end+1} = welchLine("power S2 a0 vs a1", T.power_a0(S2), T.power_a1(S2));
lines{end+1} = welchLine("amp   S2 a0 vs a1", T.amp_a0(S2),   T.amp_a1(S2));

% 追加: S3〜S6 各セッション内 a0 vs a1
for sess = ["S3","S4","S5","S6"]
    mask = strcmp(T.session, sess);
    if nnz(mask) >= 2
        lines{end+1} = welchLine("power " + sess + " a0 vs a1", T.power_a0(mask), T.power_a1(mask));
        lines{end+1} = welchLine("amp   " + sess + " a0 vs a1", T.amp_a0(mask),   T.amp_a1(mask));
    else
        lines{end+1} = sprintf("power %s a0 vs a1: サンプル不足", sess);
        lines{end+1} = sprintf("amp   %s a0 vs a1: サンプル不足", sess);
    end
end

% --- 既存 for i=1:numel(lines) 出力の直前あたりを以下に差し替え ---

% 生 p値抽出と効果量計算
rawP = [];
labels = [];
ds = [];
for i=1:numel(lines)
    tok = regexp(lines{i}, '^(.*): t=([-\d.]+) p=([0-9.eE-]+) n=\((\d+),(\d+)\)$', 'tokens','once');
    if ~isempty(tok)
        labels(end+1,1) = string(tok{1}); %#ok<AGROW>
        rawP(end+1,1) = str2double(tok{3}); %#ok<AGROW>
        % 効果量 d 再計算
        name = tok{1};
        parts = split(name);
        if contains(name,"S2 vs S6")
            sessA = "S2"; sessB = "S6";
            if contains(name,"power a0")
                xa = T.power_a0(strcmp(T.session,sessA));
                xb = T.power_a0(strcmp(T.session,sessB));
            elseif contains(name,"amp a0")
                xa = T.amp_a0(strcmp(T.session,sessA));
                xb = T.amp_a0(strcmp(T.session,sessB));
            elseif contains(name,"power a1")
                xa = T.power_a1(strcmp(T.session,sessA));
                xb = T.power_a1(strcmp(T.session,sessB));
            else
                xa = T.amp_a1(strcmp(T.session,sessA));
                xb = T.amp_a1(strcmp(T.session,sessB));
            end
        else
            % 同一セッション a0 vs a1
            sessTok = regexp(name,'(S[2-6])','tokens','once');
            if isempty(sessTok)
                xa = []; xb = [];
            else
                sessX = sessTok{1};
                xa = selectCol(T, sessX, "power", "a0", name);
                xb = selectCol(T, sessX, "power", "a1", name);
                if contains(name,"amp ")
                    xa = selectCol(T, sessX, "amp", "a0", name);
                    xb = selectCol(T, sessX, "amp", "a1", name);
                end
            end
        end
        ds(end+1,1) = cohend(xa, xb); %#ok<AGROW>
    end
end

% Holm 補正
if ~isempty(rawP)
    [pHolm, order] = holm_adjust(rawP);
    for k=1:numel(order)
        idx = order(k);
        % labels(idx) は string、lines は cell(char)。統一して検索:
        targetLabel = labels(idx) + ":";
        % 文字列配列に変換して検索
        lineStrs = string(lines);
        linesIdx = find(contains(lineStrs, targetLabel));
        if ~isempty(linesIdx)
            lines{linesIdx(1)} = lines{linesIdx(1)} + sprintf(" (Holm p=%.6g, d=%.3f)", pHolm(idx), ds(idx));
        end
    end
end

fprintf("\n[STATS]\n");
for i=1:numel(lines)
    disp(lines{i});
end

outFile = fullfile(resDir, sprintf("wavelet_stats_%s.txt", dateStr));
fid = fopen(outFile,'w');
for i=1:numel(lines)
    fprintf(fid,"%s\n", lines{i});
end
fclose(fid);
fprintf("[INFO] saved stats -> %s\n", outFile);
end

function line = welchLine(name, x, y)
x = x(~isnan(x)); y = y(~isnan(y));
if numel(x)<2 || numel(y)<2
    line = sprintf("%s: サンプル不足", name);
    return;
end
[~,p,~,st] = ttest2(x,y,'Vartype','unequal');
line = sprintf("%s: t=%.3f p=%s n=(%d,%d)", name, st.tstat, fmtP(p), numel(x), numel(y));
end

function s = fmtP(p)
if p < 1e-6
    s = sprintf("%.1e", p);
else
    s = sprintf("%.6f", p);
end
end

% 末尾に補助関数追加
function v = selectCol(T, sess, typ, ch, name)
mask = strcmp(T.session, sess);
col = typ + "_" + ch;
if ismember(col, T.Properties.VariableNames)
    v = T.(col)(mask);
else
    v = [];
end
end

function d = cohend(x,y)
x = x(~isnan(x)); y = y(~isnan(y));
nx = numel(x); ny = numel(y);
if nx<2 || ny<2, d = NaN; return; end
sx = var(x,1); sy = var(y,1);
sp = sqrt(((nx-1)*sx + (ny-1)*sy)/(nx+ny-2));
if sp==0, d = 0; else, d = (mean(x)-mean(y))/sp; end
end

function [adjP, order] = holm_adjust(p)
[ps, order] = sort(p);
m = numel(p);
adj = zeros(size(ps));
for i=1:m
    adj(i) = max(ps(i:end).* (m - (i-1)));
end
adj = min(adj,1);
% 元の順序に戻す
adjP = zeros(size(p));
adjP(order) = adj;
end