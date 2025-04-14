% filepath: main_plot.m

function main_plot()
    % ディレクトリの設定
    directory = 'merged_chunks'; % データが保存されているディレクトリ

    % ファイル選択モードを有効にするかどうか
    select_file_mode = false; % trueなら選択モード、falseなら最新ファイルをプロット

    if select_file_mode
        % 選択モードの場合
        selected_file = fullfile(directory, 'merged_chunks\merged_20250415_014757.csv'); % 選択したファイル名を指定
        if isfile(selected_file)
            fprintf('[INFO] Selected file: %s\n', selected_file);
            plot_relative_phase_matlab(selected_file);
        else
            fprintf('[ERROR] Selected file not found: %s\n', selected_file);
        end
    else
        % 最新からn番目のファイルをプロット
        n = 12; % ここでnを指定
        plot_nth_latest_file_in_merged_chunks(n, directory);
    end
end

function plot_nth_latest_file_in_merged_chunks(n, directory)
    % ディレクトリが存在するか確認
    if ~isfolder(directory)
        fprintf('[ERROR] Directory not found: %s\n', directory);
        return;
    end

    % ディレクトリ内のCSVファイルを取得
    csv_files = dir(fullfile(directory, '*.csv'));
    if isempty(csv_files)
        fprintf('[INFO] No CSV files found in directory: %s\n', directory);
        return;
    end

    % ファイル名でソート（逆順）
    file_names = string({csv_files.name}); % ファイル名を文字列配列に変換
    [~, idx] = sort(file_names, 'descend'); % ファイル名で降順にソート
    csv_files = csv_files(idx);

    % n番目のファイルを取得
    if n > length(csv_files) || n < 1
        fprintf('[ERROR] Invalid value for n: %d. There are only %d files.\n', n, length(csv_files));
        return;
    end

    nth_file = fullfile(directory, csv_files(n).name);
    fprintf('[INFO] %dth latest file found: %s\n', n, nth_file);

    % プロット関数を呼び出し
    plot_relative_phase_matlab(nth_file);
end