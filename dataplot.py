import pandas as pd
import matplotlib.pyplot as plt

# ファイル名：例 "compressed_log_YYYYMMDD_HHMMSS.csv"
file_path = "compressed_log_20250403_114702.csv"
# ファイルの存在確認
import os
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File {file_path} does not exist.")
df = pd.read_csv(file_path)

# 1) 16ビットオーバーフロー補正（micros16）
time_us = []
wrap_offset_us = 0
prev_val = df["micros"].iloc[0]

for val in df["micros"]:
    # val: 0..65535 と想定 (16-bit)
    if val < prev_val:
        # 16-bit overflow => 65536カウント = 65536 * 16μs
        wrap_offset_us += 65536 * 16 * 16

    # "val * 16 + wrap_offset_us" => 実際のマイクロ秒
    time_us.append(val * 16 + wrap_offset_us)
    prev_val = val

# 2) 秒に変換
df["time_us"] = time_us
df["time_sec"] = (df["time_us"] - df["time_us"].iloc[0]) / 1e6

# 3) グラフ化
fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

# a0
axs[0].plot(df["time_sec"], df["a0"], label="a0", color="tab:blue")
axs[0].set_ylabel("a0")
axs[0].grid(True)
axs[0].legend()

# a1
axs[1].plot(df["time_sec"], df["a1"], label="a1", color="tab:green")
axs[1].set_ylabel("a1")
axs[1].grid(True)
axs[1].legend()

# a2
axs[2].plot(df["time_sec"], df["a2"], label="a2", color="tab:orange")
axs[2].set_ylabel("a2")
axs[2].set_xlabel("Time (sec)")
axs[2].grid(True)
axs[2].legend()

plt.tight_layout()
plt.show()
