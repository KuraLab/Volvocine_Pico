function plot_relative_phase_matlab(file_list, base_agent_id)
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

    % 新しい時系列を定義 (100Hz)
    min_time = min(df_all.time_pc_sec_abs);
    max_time = max(df_all.time_pc_sec_abs);

    % 各エージェントの時間範囲を調整
    agents = unique(df_all.agent_id);
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_all(df_all.agent_id == agent_id, :);
        min_time = max(min_time, min(sub.time_pc_sec_abs));
        max_time = min(max_time, max(sub.time_pc_sec_abs));
    end

    if min_time >= max_time
        disp(['[INFO] No overlapping time range for agents. min_time=', num2str(min_time), ', max_time=', num2str(max_time)]);
        return;
    end

    new_time_series = (min_time:0.01:max_time) - min_time; % 最小値を基準にシフト

    % 線形補間で位相データを再定義
    interpolated_data = struct();
    for i = 1:length(agents)
        agent_id = agents(i);
        sub = df_all(df_all.agent_id == agent_id, :);
        sub = sortrows(sub, 'time_pc_sec_abs');
        sub.a0 = correct_phase_discontinuity(sub.a0);
        interpolated_data(agent_id).time = new_time_series;
        interpolated_data(agent_id).a0 = interp1(sub.time_pc_sec_abs - min_time, sub.a0, new_time_series, 'linear', 'extrap');
    end

    % 基準エージェントの選択
    if nargin < 2 || isempty(base_agent_id)
        base_agent_id = min(agents); % デフォルトで最小のエージェントIDを基準にする
    end
    if ~ismember(base_agent_id, agents)
        error('[ERROR] Base agent ID not found in the data.');
    end
    base_agent_a0 = interpolated_data(base_agent_id).a0;

    % プロット
    figure;
    hold on;
    colors = lines(length(agents));

    % 基準エージェント (id1 - id1) のプロット
    plot(new_time_series, zeros(size(new_time_series)), 'DisplayName', ['Agent ', num2str(base_agent_id), ' - Agent ', num2str(base_agent_id)], 'Color', colors(1, :));

    % 他のエージェントとの相対位相差をプロット
    for i = 1:length(agents)
        agent_id = agents(i);
        if agent_id == base_agent_id
            continue;
        end

        % 相対位相差を計算
        phase_diff = mod(interpolated_data(agent_id).a0 - base_agent_a0 + 128, 256) - 128;
        phase_diff = phase_diff * (2 * pi / 256); % 縦軸のデータを 2π/256 でスケール

        % NaNを挿入
        phase_diff_with_nan = phase_diff;
        for j = 2:length(phase_diff)
            if abs(phase_diff(j) - phase_diff(j - 1)) > pi % 128に相当するスケール
                phase_diff_with_nan(j) = NaN;
            end
        end

        % プロット
        plot(interpolated_data(agent_id).time, phase_diff_with_nan, 'DisplayName', ['Agent ', num2str(agent_id), ' - Agent ', num2str(base_agent_id)], 'Color', colors(i, :));
    end

    % 縦軸の目盛りをπ単位で設定し、範囲を -π から π に制限
    ylim([-pi, pi]);
    yticks(-pi:pi/2:pi);
    yticklabels({'-π', '-π/2', '0', 'π/2', 'π'});

    % プロットの設定
    xlabel('Time (s)');
    ylabel('Relative Phase (rad)');
    %legend('show', 'Location', 'best');
    grid on;
    tuneFigure;
    hold off;
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
