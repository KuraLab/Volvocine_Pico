from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main():
    here = Path(__file__).resolve().parent
    csv_path = here / "data1.csv"

    if not csv_path.exists():
        raise FileNotFoundError(f"CSVが見つかりません: {csv_path}")

    df = pd.read_csv(csv_path).sort_values("t").reset_index(drop=True)

    df["Wh"] = df["Ws"] / 3600.0
    # Whをジュール(J)に変換（1 Wh = 3600 J）
    df["J"] = df["Ws"]

    def fit_slope(t, y):
        # y = a t + b
        a, b = np.polyfit(t, y, 1)
        return a, b


    ranges = {
        "water": (0.0, 20.0),
        "air":   (20.0, 40.0),
    }

    fits = {}

    for name, (t0, t1) in ranges.items():
        mask = (df["t"] >= t0) & (df["t"] < t1)
        d = df.loc[mask]

        a, b = fit_slope(d["t"].values, d["J"].values)
        fits[name] = (a, b)

        print(f"{name}: slope = {a:.6e} [W], intercept = {b:.6e}")

    # 消費電力を別途出力
    print("\n=== 消費電力 ===")
    for name, (a, b) in fits.items():
        print(f"{name}: {a:.6e} W ({a:.2f} W)")

    plots = [
        ("current.png", "A", "Current [A]"),
        ("voltage.png", "V", "Voltage [V]"),
        ("power.png",   "W", "Power [W]"),
        ("charge.png",  "As", "Charge [A·s]"),
        ("energy.png",  "Ws", "Energy [W·s]"),
        ("energy_J.png", "J", "Energy [J]"),
    ]

    for fname, col, ylabel in plots:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # energy_J.png の場合はラベルをつける
        if fname == "energy_J.png":
            ax.plot(df["t"], df[col], linewidth=4, label="Measured Energy")
        else:
            ax.plot(df["t"], df[col], linewidth=4)
        
        # "energy_J.png" の場合は近似直線を追加
        if fname == "energy_J.png":
            colors = {"water": "orange", "air": "lime"}
            for name, (t0, t1) in ranges.items():
                mask = (df["t"] >= t0) & (df["t"] < t1)
                d = df.loc[mask]
                a, b = fits[name]
                t_fit = np.array([t0, t1])
                y_fit = a * t_fit + b
                # 関数式を表示: y = ax + b
                fit_label = f"{name}: y = {a:.2e}·t + {b:.2e}"
                ax.plot(t_fit, y_fit, ":", label=fit_label, linewidth=4, color=colors[name])
            ax.legend(fontsize=24, loc='upper left')
            # データの最大値に余裕を持たせる（最大値の1.2倍）
            y_max = df[col].max()
            ax.set_ylim(top=y_max * 1.2)
        
        
        ax.set_xlabel("t [s]", fontsize=24)
        ax.set_ylabel(ylabel, fontsize=24)
        ax.tick_params(labelsize=24)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(here / fname, dpi=200)

    print("OK")


if __name__ == "__main__":
    main()
