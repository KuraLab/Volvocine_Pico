import os
from Plotter import plot_chunks, plot_relativePhase

def plot_nth_latest_file_in_merged_chunks(n, directory="merged_chunks"):
    # ディレクトリが存在するか確認
    if not os.path.isdir(directory):
        print(f"[ERROR] Directory not found: {directory}")
        return

    # ディレクトリ内のCSVファイルを取得
    csv_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".csv")]
    if not csv_files:
        print(f"[INFO] No CSV files found in directory: {directory}")
        return

    # ファイル名でソート（逆順）
    csv_files.sort(key=lambda f: os.path.basename(f), reverse=True)

    # n番目のファイルを取得
    if n > len(csv_files) or n < 1:
        print(f"[ERROR] Invalid value for n: {n}. There are only {len(csv_files)} files.")
        return

    nth_file = csv_files[n - 1]
    print(f"[INFO] {n}th latest file found: {nth_file}")

    # プロット関数を呼び出し
    plot_relativePhase(nth_file)

if __name__ == "__main__":
    # ファイル選択モードを有効にするかどうか
    select_file_mode = False  # Trueなら選択モード、Falseなら最新ファイルをプロット

    if select_file_mode:
        # 選択モードの場合
        selected_file = "merged_chunks/merged_20250407_231718.csv"  # プロットしたいファイルを指定
        if os.path.isfile(selected_file):
            print(f"[INFO] Selected file: {selected_file}")
            plot_relativePhase(selected_file)
        else:
            print(f"[ERROR] Selected file not found: {selected_file}")
    else:
        # 最新からn番目のファイルをプロット
        n = 9  # ここでnを指定
        plot_nth_latest_file_in_merged_chunks(n)