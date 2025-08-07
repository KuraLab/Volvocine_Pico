% filepath: main_plot.m

function main_plot()
    % 実験日付の設定（YYYY-MM-DD形式）
    experiment_date = '2025-08-07'; % ← ここで解析したい実験日を指定
    
    % ディレクトリの設定
    base_directory = 'merged_chunks_organized'; % 日付別に整理されたディレクトリ
    
    % 最初のn秒をカットする設定
    n_seconds_to_cut = 5; % プロット時に最初のn秒をカット

    % 何秒目までプロットするか
    plot_duration = 80; % 例: 183秒までプロット

    % 許容する欠損エージェント数
    allow_missing_agents = 1; % 欠損許容数
    do_save_figure = false;   % ← ここで保存有無を指定

    % 平均化フィルタの設定
    apply_filter = true;      % フィルタ適用の有無
    filter_window_size = 1;  % フィルタ窓サイズ

    % 最新からn番目のファイルをプロット
    n = 1; % ここでnを指定
    plot_nth_latest_file_by_date(experiment_date, n, base_directory, n_seconds_to_cut, plot_duration, allow_missing_agents, do_save_figure, apply_filter, filter_window_size);
end

function plot_nth_latest_file_by_date(experiment_date, n, base_directory, n_seconds_to_cut, plot_duration, allow_missing_agents, do_save_figure, apply_filter, filter_window_size)
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

    % プロット関数を呼び出し
    plot_relative_phase_matlab(nth_file, [], n_seconds_to_cut, plot_duration, allow_missing_agents, do_save_figure, apply_filter, filter_window_size);
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

% 便利関数: 利用可能な日付を一覧表示
function show_available_dates()
    base_directory = 'merged_chunks_organized';
    fprintf('=== Available experiment dates ===\n');
    list_available_dates(base_directory);
end

% 便利関数: 最新の日付のデータをプロット
function plot_latest_date(n)
    if nargin < 1
        n = 1; % デフォルトは最新ファイル
    end
    
    base_directory = 'merged_chunks_organized';
    
    % 最新の日付を取得
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
        fprintf('[ERROR] No date folders found\n');
        return;
    end
    
    latest_date = string(max(date_folders));
    fprintf('[INFO] Using latest date: %s\n', latest_date);
    
    % デフォルト設定で実行
    plot_nth_latest_file_by_date(latest_date, n, base_directory, 10, 183, 1, false, true, 10);
end