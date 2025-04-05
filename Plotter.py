import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import os

def plot_chunks(file_list):
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
