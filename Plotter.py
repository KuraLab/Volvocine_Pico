import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import os
import numpy as np

def plot_chunks(file_list):
    # None チェックを追加
    if not file_list:
        print("[INFO] No files provided to plot.")
        return

    # 単一ファイルパスが来たらリストに変換
    if isinstance(file_list, str):
        file_list = [file_list]

    print("[DEBUG] Plotting from files:")
    for f in file_list:
        print(f"  - {f}")

    dfs = []
    for file in file_list:
        if not os.path.isfile(file):
            print(f"[WARN] File not found: {file}")
            continue

        try:
            df = pd.read_csv(file)
            if all(col in df.columns for col in ["agent_id", "chunk_id", "time_pc_sec_abs", "a0", "a1", "a2"]):
                dfs.append(df[["agent_id", "chunk_id", "time_pc_sec_abs", "a0", "a1", "a2"]])
        except Exception as e:
            print(f"[WARN] Failed to load {file}: {e}")

    if not dfs:
        print("[INFO] No valid data to plot.")
        return

    df_all = pd.concat(dfs, ignore_index=True)
    fig, axs = plt.subplots(3, 1, figsize=(9, 6), sharex=True)
    colors = {}
    color_cycle = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for (ag_id, _), sub in df_all.groupby(["agent_id", "chunk_id"]):
        if ag_id not in colors:
            colors[ag_id] = next(color_cycle)
        axs[0].plot(sub["time_pc_sec_abs"], sub["a0"], color=colors[ag_id])
        axs[1].plot(sub["time_pc_sec_abs"], sub["a1"], color=colors[ag_id])
        axs[2].plot(sub["time_pc_sec_abs"], sub["a2"], color=colors[ag_id])

    axs[0].set_ylabel("a0")
    axs[1].set_ylabel("a1")
    axs[2].set_ylabel("a2")
    axs[2].set_xlabel("PC time (sec)")
    for ax in axs:
        ax.grid(True)

    handles = [plt.Line2D([0], [0], color=color, lw=2, label=f"Agent {ag_id}")
               for ag_id, color in colors.items()]
    axs[0].legend(handles=handles, title="Agents")

    plt.tight_layout()
    plt.show()

def correct_phase_discontinuity(phase_data):
    """
    位相データのジャンプを補正する関数。
    急激な変化があった場合に 256 を加算または減算して連続性を保つ。
    """
    corrected_phase = phase_data.copy()
    for i in range(1, len(corrected_phase)):
        diff = corrected_phase[i] - corrected_phase[i - 1]
        if diff < -128:  # 急に128以上小さくなった場合
            corrected_phase[i:] += 256
        elif diff > 128:  # 急に128以上大きくなった場合
            corrected_phase[i:] -= 256
    return corrected_phase

def plot_relativePhase(file_list):
    # None チェックを追加
    if not file_list:
        print("[INFO] No files provided to plot.")
        return

    # 単一ファイルパスが来たらリストに変換
    if isinstance(file_list, str):
        file_list = [file_list]

    print("[DEBUG] Plotting from files:")
    for f in file_list:
        print(f"  - {f}")

    dfs = []
    for file in file_list:
        if not os.path.isfile(file):
            print(f"[WARN] File not found: {file}")
            continue

        try:
            df = pd.read_csv(file)
            if all(col in df.columns for col in ["agent_id", "chunk_id", "time_pc_sec_abs", "a0", "a1", "a2"]):
                dfs.append(df[["agent_id", "chunk_id", "time_pc_sec_abs", "a0", "a1", "a2"]])
        except Exception as e:
            print(f"[WARN] Failed to load {file}: {e}")

    if not dfs:
        print("[INFO] No valid data to plot.")
        return

    df_all = pd.concat(dfs, ignore_index=True)

    # 新しい時系列を定義 (100Hz)
    # 初期値を適切に設定
    min_time = df_all["time_pc_sec_abs"].min()  # 全体の最小値
    max_time = df_all["time_pc_sec_abs"].max()  # 全体の最大値

    for agent_id, sub in df_all.groupby("agent_id"):
        sub = sub.sort_values("time_pc_sec_abs")
        # 各エージェントのデータ範囲を考慮して更新
        min_time = max(min_time, sub["time_pc_sec_abs"].min())  # 各エージェントの最小値
        max_time = min(max_time, sub["time_pc_sec_abs"].max())  # 各エージェントの最大値

    # チェックを追加
    if min_time >= max_time:
        print(f"[INFO] No overlapping time range for agents. min_time={min_time}, max_time={max_time}")
        return

    new_time_series = np.arange(min_time, max_time, 0.01)

    # 線形補間で位相データを再定義
    interpolated_data = {}
    for agent_id, sub in df_all.groupby("agent_id"):
        sub = sub.sort_values("time_pc_sec_abs")
        
        # 位相データのジャンプを補正
        sub["a0"] = correct_phase_discontinuity(sub["a0"].values)
        
        # 線形補間
        interpolated_data[agent_id] = {
            "time": new_time_series,
            "a0": np.interp(new_time_series, sub["time_pc_sec_abs"], sub["a0"])
        }

    # 基準エージェントの選択
    base_agent_id = min(interpolated_data.keys())  # 最初のエージェントを基準とする
    base_agent_a0 = interpolated_data[base_agent_id]["a0"]

    fig, axs = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    # 元データをプロット (補正前のデータ)
    for agent_id, sub in df_all.groupby("agent_id"):
        sub = sub.sort_values("time_pc_sec_abs")
        axs[0].plot(sub["time_pc_sec_abs"], sub["a0"], label=f"Agent {agent_id} (raw)", color=next(colors))

    # 線形補間と補正後のデータを使用して相対位相差を計算
    for agent_id, data in interpolated_data.items():
        if agent_id == base_agent_id:
            continue  # 基準エージェントはスキップ

        # 相対位相差を計算 (単純な引き算)
        phase_diff = (data["a0"] - base_agent_a0 + 128) % 256 - 128  # 256でモッドを取る

        # 位相差のジャンプを検出してNaNを挿入
        phase_diff_with_nan = phase_diff.copy()
        for i in range(1, len(phase_diff)):
            if abs(phase_diff[i] - phase_diff[i - 1]) > 128:  # ジャンプを検出
                phase_diff_with_nan[i] = np.nan  # NaNを挿入

        print(f"[DEBUG] Agent {agent_id}: Phase diff with NaN (first 10 values) = {phase_diff_with_nan[:10]}")
        axs[1].plot(data["time"], phase_diff_with_nan, label=f"Agent {agent_id} - Agent {base_agent_id}", color=next(colors))

    # プロットの設定
    axs[0].set_ylabel("Phase (a0) (Raw)")
    axs[1].set_ylabel("Phase Diff")
    axs[1].set_xlabel("Time (s)")
    axs[0].legend(title="Agents (Raw Data)")
    axs[1].legend(title="Relative Phase")
    axs[0].grid(True)
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()
