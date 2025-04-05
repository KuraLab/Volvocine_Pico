import pandas as pd
from datetime import datetime
import os

SAVE_FOLDER = "."

def merge_and_save_chunks(chunk_files):
    # マージデータ専用フォルダの作成
    merged_folder = os.path.join(SAVE_FOLDER, "merged_chunks")
    if not os.path.exists(merged_folder):
        os.makedirs(merged_folder)

    # マージ処理
    merged_data = pd.DataFrame()
    for file in chunk_files:
        df = pd.read_csv(file)
        merged_data = pd.concat([merged_data, df], ignore_index=True)

    # 保存先ファイル名の生成
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    merged_file = os.path.join(merged_folder, f"merged_{timestamp}.csv")

    # マージデータを保存
    merged_data.to_csv(merged_file, index=False)
    print(f"[INFO] Merged data saved to {merged_file}")

    return merged_file
