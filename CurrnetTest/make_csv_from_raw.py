import re
from pathlib import Path

import numpy as np
import pandas as pd

# ========= ここにデータを貼り付ける =========
RAW = r"""
6.0,0.2626,1.195,0.186,0.90,4.988,23.9
6.1,0.2759,1.262,0.132,0.65,5.144,23.9
"""
# ==========================================


def parse_pasted_log(raw: str, n_cols: int = 7) -> np.ndarray:
    s = raw.strip()
    if not s:
        raise ValueError("RAW が空です。RAW の貼り付け欄にログを貼り付けてください。")

    recs = re.split(r"\s+", s)

    rows = []
    bad = []
    for i, rec in enumerate(recs):
        if not rec:
            continue
        if rec.startswith("#"):
            continue

        cols = rec.split(",")
        if len(cols) != n_cols:
            bad.append((i, rec, f"列数={len(cols)}"))
            continue
        try:
            rows.append([float(x) for x in cols])
        except ValueError:
            bad.append((i, rec, "float変換失敗"))

    if bad:
        msg = "\n".join([f"{b[0]}: {b[1]} ({b[2]})" for b in bad[:15]])
        raise ValueError(
            "壊れたレコードがありました（先頭15件だけ表示）:\n"
            f"{msg}\n"
            "対処:\n"
            "- 1行が '数値,数値,数値,数値,数値,数値,数値' の7列になっているか確認\n"
            "- 変な記号や全角スペースが混じっていないか確認\n"
        )

    if not rows:
        raise ValueError("有効なレコードが1件も読めませんでした。RAW の内容を確認してください。")

    return np.asarray(rows, dtype=float)


def main() -> None:
    data = parse_pasted_log(RAW, n_cols=7)

    cols = ["t", "As", "Ws", "A", "W", "V", "T"]
    df = pd.DataFrame(data, columns=cols).sort_values("t").reset_index(drop=True)

    out_path = Path("data.csv")
    df.to_csv(out_path, index=False)
    print(f"saved: {out_path.resolve()}")
    print(f"rows: {len(df)}")
    print(f"t range: {df['t'].iloc[0]} .. {df['t'].iloc[-1]}")


if __name__ == "__main__":
    main()
