function plot_relative_phase_matlab_2modules(file_list, base_agent_id, n_seconds, plot_duration, allow_missing_agents, do_save_figure, apply_filter, filter_window_size, order_index)
% order_index : 出現 agent_id を昇順に並べた配列 agents に対するインデックス並び
%               例) agents = [3 7 9 12]; order_index = [1 3 2 4] -> 表示順 [3 9 7 12]
%               未指定なら元の昇順。インデックスは 1..numel(agents)

    % ファイルリストが空かどうかチェック
    if isempty(file_list)
        disp('[INFO] No files provided to plot.');
        return;
    end

    % 単一ファイルの場合はリストに変換
    if ischar(file_list)
        file_list = {file_list};
    end

    disp('[DEBUG] Plotting from files:');
    for i = 1:length(file_list)
        disp(['  - ', file_list{i}]);
    end

    % データを読み込む
    dfs = {};
    for i = 1:length(file_list)
        file = file_list{i};
        if ~isfile(file)
            disp(['[WARN] File not found: ', file]);
            continue;
        end

        try
            data = readtable(file);
            if all(ismember({'agent_id', 'chunk_id', 'time_pc_sec_abs', 'a0', 'a1', 'a2'}, data.Properties.VariableNames))
                dfs{end+1} = data(:, {'agent_id', 'chunk_id', 'time_pc_sec_abs', 'a0', 'a1', 'a2'});
            end
        catch ME
            disp(['[WARN] Failed to load ', file, ': ', ME.message]);
        end
    end

    if isempty(dfs)
        disp('[INFO] No valid data to plot.');
        return;
    end

    % データを結合
    df_all = vertcat(dfs{:});

    % --- Pythonと同じ定数 ---
    T_OVERFLOW = 2^32 / 1e6; % 約4294.967296秒
    T_TOL = 5.0;             % 許容誤差（秒）
    threshold_sec = T_OVERFLOW - T_TOL;
    jump_sec = T_OVERFLOW;

    % チャンク内ジャンプ補正
    df_all = correct_large_jump_matlab(df_all, threshold_sec, jump_sec);

    % チャンク開始時刻の未来ジャンプ補正
    df_all = correct_chunk_start_times_matlab(df_all, 4000.0, T_OVERFLOW);

    % agent_id==99のデータを分離
    df_99 = df_all(df_all.agent_id == 99, :);
    df_main = df_all(df_all.agent_id ~= 99, :);

    % 新しい時系列を定義 (100Hz)
    min_time = min(df_main.time_pc_sec_abs);
    max_time = max(df_main.time_pc_sec_abs);

    % 各エージェントの時間範囲を取得
    % エージェント ID を昇順取得し行ベクトル化
    agents = unique(df_main.agent_id,'sorted').';  % 行ベクトル強制

    agent_ranges = zeros(length(agents), 2);
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_main(df_main.agent_id == agent_id, :);
        agent_ranges(i, 1) = min(sub.time_pc_sec_abs);
        agent_ranges(i, 2) = max(sub.time_pc_sec_abs);
    end

    % ---- 追加: 表示順制御 (order_index) ----
    if nargin < 9 || isempty(order_index)
        order_index = 1:numel(agents);
    else
        order_index = order_index(:).';                 % 行ベクトル化
        valid_mask = order_index >= 1 & order_index <= numel(agents);
        if any(~valid_mask)
            warning('無効な order_index を無視: %s', mat2str(order_index(~valid_mask)));
            order_index = order_index(valid_mask);
        end
        % 重複除去
        [~, first_pos] = unique(order_index,'stable');
        if numel(first_pos) ~= numel(order_index)
            warning('order_index の重複を除去しました。');
            order_index = order_index(sort(first_pos));
        end
        if isempty(order_index)
            error('有効な order_index がありません。');
        end
    end
    display_agents = agents(order_index);              % 指定順
    display_agents = display_agents(:).';              % 念のため行ベクトル

    % 基準エージェント
    if nargin < 2 || isempty(base_agent_id)
        base_agent_id = min(agents);
    end
    if ~ismember(base_agent_id, agents)
        error('[ERROR] Base agent ID %d not found in data.', base_agent_id);
    end
    % 先頭へ（存在しない場合は追加）
    display_agents(display_agents == base_agent_id) = [];
    display_agents = [base_agent_id display_agents];

    % デバッグ表示（必要なら）
    % fprintf('[DEBUG] agents (asc): %s\n', mat2str(agents));
    % fprintf('[DEBUG] order_index : %s\n', mat2str(order_index));
    % fprintf('[DEBUG] display_agents(final): %s\n', mat2str(display_agents));
    % ---- 追加ここまで ----

    % --- ここから修正 ---
    if nargin < 5 || isempty(allow_missing_agents)
        allow_missing_agents = 0; % デフォルトは全員そろっている範囲
    end
    if nargin < 7 || isempty(apply_filter)
        apply_filter = false; % デフォルトはフィルタ無効
    end
    if nargin < 8 || isempty(filter_window_size)
        filter_window_size = 50; % デフォルトの窓サイズ
    end
    n = allow_missing_agents; % 欠損許容数

    % n_secondsを基準時刻として設定
    start_time_abs = min_time + n_seconds;
    new_time_series = (start_time_abs:0.01:max_time) - start_time_abs; % n_secondsを0とする時間軸

    % 各時刻で有効なエージェント数をカウント
    valid_counts = zeros(size(new_time_series));
    for t_idx = 1:length(new_time_series)
        t_abs = new_time_series(t_idx) + start_time_abs;
        valid_counts(t_idx) = sum(agent_ranges(:,1) <= t_abs & agent_ranges(:,2) >= t_abs);
    end

    % 「全体-n個」以上データがある時刻だけ抽出
    min_valid = length(agents) - n;
    valid_idx = find(valid_counts >= min_valid);

    % その範囲だけを使う
    new_time_series = new_time_series(valid_idx);
    % --- 以降、new_time_seriesを使って補間・プロット ---

    % 線形補間（全 agents 対象: 基準計算のため）
    interpolated_data = struct();
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_main(df_main.agent_id == agent_id, :);
        sub = sortrows(sub, 'time_pc_sec_abs');
        sub.a0 = correct_phase_discontinuity(sub.a0);
        [~, ia] = unique(sub.time_pc_sec_abs);
        sub = sub(ia, :);
        t_min = min(sub.time_pc_sec_abs) - start_time_abs;
        t_max = max(sub.time_pc_sec_abs) - start_time_abs;
        valid_mask = (new_time_series >= t_min) & (new_time_series <= t_max);
        interp_a0 = nan(size(new_time_series));
        interp_a0(valid_mask) = interp1(sub.time_pc_sec_abs - start_time_abs, sub.a0, new_time_series(valid_mask), 'linear', 'extrap');
        if apply_filter
            interp_a0 = movmean(interp_a0, filter_window_size, 'omitnan');
        end
        interpolated_data(agent_id).time = new_time_series;
        interpolated_data(agent_id).a0   = interp_a0;
    end

    % 基準エージェント位相
    base_agent_a0 = interpolated_data(base_agent_id).a0;

    % Agent99 a0/a1プロット用の最大時刻
    if ~isempty(df_99)
        max_time_99 = max(df_99.time_pc_sec_abs - start_time_abs);
    else
        max_time_99 = inf; % データがなければ無限大扱い
    end

    % 相対位相プロット用の最大時刻
    max_time_phase = max(new_time_series); % ここもデータの最大値

    % 両方の最大値の小さい方を採用
    common_xmax = min([max_time_phase, max_time_99, plot_duration - n_seconds]);

    % --- 相対位相プロット ---
    figure; hold on;
    colors = lines(length(display_agents));

    % 基準ライン
    h_base = plot(new_time_series, zeros(size(new_time_series)), ...
        'DisplayName', sprintf('Agent %d', base_agent_id), ...
        'Color', colors(1,:));
    legend_handles = h_base;
    legend_labels  = {sprintf('Agent %d', base_agent_id)};

    % 他エージェント
    for k = 2:length(display_agents)
        agent_id = display_agents(k);
        phase_diff = mod(interpolated_data(agent_id).a0 - base_agent_a0 + 128, 256) - 128;
        phase_diff = phase_diff * (2*pi/256);

        % ジャンプ除去
        phase_diff_with_nan = phase_diff;
        for j = 2:length(phase_diff)
            if abs(phase_diff(j) - phase_diff(j-1)) > pi
                phase_diff_with_nan(j) = NaN;
            end
        end

        h = plot(interpolated_data(agent_id).time, phase_diff_with_nan, ...
            'DisplayName', sprintf('Agent %d', agent_id), ...
            'Color', colors(k,:));
        legend_handles(end+1) = h; %#ok<AGROW>
        legend_labels{end+1}  = sprintf('Agent %d', agent_id); %#ok<AGROW>
    end

    ylim([-pi, pi]);
    yticks(-pi:pi:pi);
    yticklabels({'$-\pi$','0','$\pi$'});
    set(gca,'TickLabelInterpreter','latex');
    xlim([0 common_xmax]);
    xlabel('Time (s)', 'Interpreter','latex');
    ylabel('$$\phi_j - \phi_1$$','Interpreter','latex');
    grid on;
    tuneFigure;
    set(findall(gcf,'-property','FontSize'),'FontSize',28);
    ncol = min(3, numel(legend_labels));
    legend({'Module 1', 'Module 2', 'Module 3', 'Module 4', 'Module 5', 'Module 6'}, 'Location', 'northoutside', 'Interpreter', 'latex','NumColumns', ncol);
    if exist('do_save_figure','var') && do_save_figure
        saveFigure;
    end
    hold off;

    % --- Agent99 a0/a1プロットは不要のため削除 ---

    % --- Agent99 e1（a0）と e2（a1）のウェーブレット変換を別ウィンドウで表示 ---
    if ~isempty(df_99)
        % 元データの準備
        t99_all   = df_99.time_pc_sec_abs - start_time_abs;
        a0_99_all = correct_large_jump_99(df_99.a0); % e1
        a1_99_all = correct_large_jump_99(df_99.a1); % e2

        % スムージング（移動平均）
        windowsize = 1; % 必要に応じて調整
        a0_99_smooth = movmean(a0_99_all, windowsize);
        a1_99_smooth = movmean(a1_99_all, windowsize);

        % uint8 → 角度（-180〜180度）へ変換
        decode_angle = @(u) (double(u) * 360.0 / 255.0) - 180.0;
        a0_99_deg_all = decode_angle(a0_99_smooth);
        a1_99_deg_all = decode_angle(a1_99_smooth);

        % 表示区間 [0, common_xmax] に制限
        idx = t99_all >= 0 & t99_all <= common_xmax;
        t99      = t99_all(idx);
        a0_99_deg = a0_99_deg_all(idx);
        a1_99_deg = a1_99_deg_all(idx);

        % CWT パラメータ
        fs = 100;              % サンプリング周波数 [Hz]
        freq_range = [0.1 10]; % 表示する周波数範囲 [Hz]

        % ウェーブレット変換（e1 と e2）
        [wt_a0, f_a0] = cwt(a0_99_deg, fs, 'FrequencyLimits', freq_range); % e1
        [wt_a1, f_a1] = cwt(a1_99_deg, fs, 'FrequencyLimits', freq_range); % e2

        % カラースケール最大値
        cmax = 3;

        % プロット（e1: 別ウィンドウ）
        figure;
        surf(t99, f_a0, abs(wt_a0), 'EdgeColor', 'none');
        set(gca, 'YScale', 'log');
        axis tight; view(0, 90);
        xlim([0, common_xmax]);
        ylim([0.2, 10]); % 周波数範囲を設定
        xlabel('Time (s)');
        ylabel('Freq [Hz]');
        clim([0 cmax]);
        colorbar;
        tuneFigure
        set(findall(gcf,'-property','FontSize'),'FontSize',27);
        if exist('do_save_figure','var') && do_save_figure
            saveFigure;
        end

        % プロット（e2: 別ウィンドウ）
        figure;
        surf(t99, f_a1, abs(wt_a1), 'EdgeColor', 'none');
        set(gca, 'YScale', 'log');
        axis tight; view(0, 90);
        xlim([0, common_xmax]);
        ylim([0.2, 10]); % 周波数範囲を設定
        xlabel('Time (s)');
        ylabel('Freq [Hz]');
        clim([0 cmax]);
        colorbar;
        tuneFigure
        set(findall(gcf,'-property','FontSize'),'FontSize',27);
        if exist('do_save_figure','var') && do_save_figure
            saveFigure;
        end
    end


end

function corrected_phase = correct_phase_discontinuity(phase_data)
    % 位相データのジャンプを補正する関数
    corrected_phase = phase_data;
    for i = 2:length(corrected_phase)
        diff = corrected_phase(i) - corrected_phase(i - 1);
        if diff < -128
            corrected_phase(i:end) = corrected_phase(i:end) + 256;
        elseif diff > 128
            corrected_phase(i:end) = corrected_phase(i:end) - 256;
        end
    end
end

function df_all = correct_large_jump_matlab(df_all, threshold_sec, jump_sec)
    % グループ化（agent_id, chunk_id 単位）
    [G, ~] = findgroups(df_all.agent_id, df_all.chunk_id);
    fprintf('[INFO] Found %d unique chunks.\n', max(G));

    % 該当ブロックを修正
    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;
        end

        % 時系列を並び替え
        [~, sorted_idx_rel] = sort(df_all.time_pc_sec_abs(idx));
        idx = idx(sorted_idx_rel);
        
        time_series = df_all.time_pc_sec_abs(idx);
        time_diff = [0; diff(time_series)];

        jump_idx = find(time_diff > threshold_sec);
        for j = 1:length(jump_idx)
            fix_range = jump_idx(j):length(time_series);
            df_all.time_pc_sec_abs(idx(fix_range)) = df_all.time_pc_sec_abs(idx(fix_range)) - jump_sec;
            fprintf('[FIX] Corrected overflow at index %d, subtracted %.6f sec.\n', idx(jump_idx(j)), jump_sec);
        end
    end

end

function df_all = correct_chunk_start_times_matlab(df_all, threshold_sec, jump_sec)

    [G, chunk_keys] = findgroups(df_all.agent_id, df_all.chunk_id);
    chunk_start = splitapply(@(x) min(x), df_all.time_pc_sec_abs, G);
    median_start = median(chunk_start);

    has_chunk_id = size(chunk_keys, 2) >= 2;

    for i = 1:max(G)
        idx = find(G == i);
        if isempty(idx)
            continue;  % 空グループスキップ
        end
        start_time = df_all.time_pc_sec_abs(idx(1));
        if start_time - median_start > threshold_sec
            df_all.time_pc_sec_abs(idx) = df_all.time_pc_sec_abs(idx) - jump_sec;
            % agent_id, chunk_idの表示（chunk_idが無い場合も対応）
            if has_chunk_id
                aid = chunk_keys(i,1);
                cid = chunk_keys(i,2);
            else
                aid = chunk_keys(i);
                cid = -1;  % または NaN
            end
            fprintf('[FIX] Corrected chunk time for agent %d, chunk %d: %.3f → %.3f\n', ...
                aid, cid, start_time, start_time - jump_sec);
        end
    end
end

function corrected = correct_large_jump_99(data)
    % 200以上のジャンプがあれば、その方向に±255補正
    corrected = double(data);
    for i = 2:length(corrected)
        diff = corrected(i) - corrected(i-1);
        if diff > 200
            corrected(i:end) = corrected(i:end) - 255;
        elseif diff < -200
            corrected(i:end) = corrected(i:end) + 255;
        end
    end
end
