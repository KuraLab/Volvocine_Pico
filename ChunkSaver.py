import pandas as pd
from datetime import datetime
import os

def merge_and_save_chunks(file_list, output_dir="."):
    if not file_list:
        print("[INFO] No chunk files to merge.")
        return None

    dfs = []
    for file in file_list:
        try:
            df = pd.read_csv(file)
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Failed to read {file}: {e}")

    if not dfs:
        print("[INFO] No valid dataframes loaded.")
        return None

    merged_df = pd.concat(dfs, ignore_index=True)

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outname = os.path.join(output_dir, f"merged_log_{timestamp}.csv")
    merged_df.to_csv(outname, index=False)
    print(f"[INFO] Merged log saved to {outname}")

    # 元ファイルを削除
    for file in file_list:
        try:
            os.remove(file)
            print(f"[INFO] Deleted chunk file: {file}")
        except Exception as e:
            print(f"[WARN] Failed to delete {file}: {e}")

    return outname
