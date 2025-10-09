% filepath: main_plot_2modules.m

function main_plot_2modules()
    % 実験日付
    experiment_date = '2025-08-18';

    base_directory = 'merged_chunks_organized';

    % 最初の n 秒除外
    n_seconds_to_cut = 5;

    % プロット終了時刻
    plot_duration = 105;

    % 欠損許容
    allow_missing_agents = 1;
    do_save_figure = true;

    % 平滑化
    apply_filter = true;
    filter_window_size = 10;

    % 追加: 並び替え設定
    % agents 昇順 = エージェント番号が[3 7 9 12] のとき [1 3 2 4] -> [3 9 7 12]
    order_index = [1 2 3 4 5 6];   % 未指定なら [] でも可
    base_agent_id = [];        % 基準を自動(最小ID)にするなら []

    % 最新から n 番目
    n = 1;

    plot_nth_latest_file_by_date( ...
        experiment_date, n, base_directory, ...
        n_seconds_to_cut, plot_duration, ...
        allow_missing_agents, do_save_figure, ...
        apply_filter, filter_window_size, ...
        order_index, base_agent_id);
end

function plot_nth_latest_file_by_date(experiment_date, n, base_directory, n_seconds_to_cut, plot_duration, allow_missing_agents, do_save_figure, apply_filter, filter_window_size, order_index, base_agent_id)
    % 日付フォルダのパスを構築
    date_directory = fullfile(base_directory, experiment_date);
    
    % ディレクトリが存在するか確認
    if ~isfolder(date_directory)
        fprintf('[ERROR] Date directory not found: %s\n', date_directory);
        fprintf('[INFO] Available dates:\n');
        list_available_dates(base_directory);
        return;
    end

    % ディレクトリ内のCSVファイルを取得
    csv_files = dir(fullfile(date_directory, '*.csv'));
    if isempty(csv_files)
        fprintf('[INFO] No CSV files found in directory: %s\n', date_directory);
        return;
    end

    % ファイル名でソート（逆順：最新が先頭）
    file_names = string({csv_files.name}); % ファイル名を文字列配列に変換
    [~, idx] = sort(file_names, 'descend'); % ファイル名で降順にソート
    csv_files = csv_files(idx);

    % n番目のファイルを取得
    if n > length(csv_files) || n < 1
        fprintf('[ERROR] Invalid value for n: %d. There are only %d files for %s.\n', n, length(csv_files), experiment_date);
        fprintf('[INFO] Available files for %s:\n', experiment_date);
        for i = 1:min(5, length(csv_files)) % 最初の5つを表示
            fprintf('  %d: %s\n', i, csv_files(i).name);
        end
        if length(csv_files) > 5
            fprintf('  ... and %d more files\n', length(csv_files) - 5);
        end
        return;
    end

    nth_file = fullfile(date_directory, csv_files(n).name);
    fprintf('[INFO] Selected file (%s, %dth latest): %s\n', experiment_date, n, csv_files(n).name);

    % 呼び出しを order_index / base_agent_id 対応版へ
    if ~exist('order_index','var'); order_index = []; end
    if ~exist('base_agent_id','var'); base_agent_id = []; end
    plot_relative_phase_matlab_2modules( ...
        nth_file, base_agent_id, ...
        n_seconds_to_cut, plot_duration, ...
        allow_missing_agents, do_save_figure, ...
        apply_filter, filter_window_size, ...
        order_index);
end

function list_available_dates(base_directory)
    % 利用可能な日付フォルダを一覧表示
    if ~isfolder(base_directory)
        fprintf('[ERROR] Base directory not found: %s\n', base_directory);
        return;
    end
    
    % 日付パターン（YYYY-MM-DD）に一致するフォルダを取得
    folders = dir(base_directory);
    date_pattern = '^\d{4}-\d{2}-\d{2}$';
    
    date_folders = {};
    for i = 1:length(folders)
        if folders(i).isdir && ~ismember(folders(i).name, {'.', '..'})
            if ~isempty(regexp(folders(i).name, date_pattern, 'once'))
                date_folders{end+1} = folders(i).name;
            end
        end
    end
    
    if isempty(date_folders)
        fprintf('  No date folders found in %s\n', base_directory);
    else
        date_folders = sort(date_folders, 'descend'); % 最新日付を上に
        for i = 1:length(date_folders)
            % 各日付のファイル数も表示
            csv_count = length(dir(fullfile(base_directory, date_folders{i}, '*.csv')));
            fprintf('  %s (%d files)\n', date_folders{i}, csv_count);
        end
    end
end
