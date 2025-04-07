import os
from Plotter import plot_chunks, plot_relativePhase

def plot_latest_file_in_merged_chunks(directory="merged_chunks"):
    # ディレクトリが存在するか確認
    if not os.path.isdir(directory):
        print(f"[ERROR] Directory not found: {directory}")
        return

    # ディレクトリ内のCSVファイルを取得
    csv_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".csv")]
    if not csv_files:
        print(f"[INFO] No CSV files found in directory: {directory}")
        return

    # 最新のファイルを取得
    latest_file = max(csv_files, key=os.path.getmtime)
    print(f"[INFO] Latest file found: {latest_file}")

    # プロット関数を呼び出し
    plot_relativePhase(latest_file)

if __name__ == "__main__":
    # ファイル選択モードを有効にするかどうか
    select_file_mode = False  # Trueなら選択モード、Falseなら最新ファイルをプロット

    if select_file_mode:
        # 選択モードの場合
        selected_file = "merged_chunks/example.csv"  # プロットしたいファイルを指定
        if os.path.isfile(selected_file):
            print(f"[INFO] Selected file: {selected_file}")
            plot_relativePhase(selected_file)
        else:
            print(f"[ERROR] Selected file not found: {selected_file}")
    else:
        # 最新ファイルをプロット
        plot_latest_file_in_merged_chunks()