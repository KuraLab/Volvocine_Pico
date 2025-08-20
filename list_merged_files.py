from pathlib import Path
import sys
import re
from datetime import datetime, time
import pandas as pd  # agent_id 集計用
from pathlib import Path

# 期待される agent_id ユニーク数 (セッションラベル -> 数)
#EXPECTED_AGENT_COUNTS = {"S1": 10, "S2": 4, "S3": 5, "S4": 6, "S5": 7, "S6": 4}
EXPECTED_AGENT_COUNTS = {"S1": 10, "S2": 5, "S3": 6, "S4": 7, "S5": 7, "S6": 4}
APPLY_AGENT_FILTER_DEFAULT = True  # 常にフィルタ適用
REQUIRED_AGENT_ID = 99  # この agent_id が含まれないファイルは除外

# 簡素化: オプション排除。日付引数 1 つ(省略可)のみ。
# 実行: python list_merged_files.py 2025-08-15
# 出力: セッション境界 / セッションごとの kept/dropped 統計

#DEFAULT_DATE = "2025-08-15"
DEFAULT_DATE = "2025-08-18"
#DEFAULT_SESSION_STARTS = ["16:40", "20:00", "20:49", "21:47", "22:22"]
DEFAULT_SESSION_STARTS = ["15:20", "16:00", "16:38", "21:47", "22:22"]
FILENAME_PATTERN = re.compile(r"merged_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})\.csv")

# -------------------------------------------------------------

def parse_filename_ts(name: str) -> datetime | None:
    m = FILENAME_PATTERN.match(name)
    if not m:
        return None
    y, mo, d, H, M, S = map(int, m.groups())
    return datetime(y, mo, d, H, M, S)


def list_files(target_date: str):
    day_dir = Path("merged_chunks_organized") / target_date
    if not day_dir.is_dir():
        print(f"[ERROR] Directory not found: {day_dir}")
        return []
    files = [p for p in day_dir.glob("merged_*.csv") if FILENAME_PATTERN.match(p.name)]
    files.sort()
    return files

# -------------------------------------------------------------

def build_session_intervals(target_date: str, starts: list[str]):
    if not starts:
        return []
    dt_list = []
    for s in starts:
        hh, mm = map(int, s.split(':'))
        dt_list.append(datetime.fromisoformat(f"{target_date} {hh:02d}:{mm:02d}:00"))
    dt_list.sort()
    intervals = []
    intervals.append({"label": "S1", "start": None, "end": dt_list[0]})
    for i in range(len(dt_list) - 1):
        intervals.append({"label": f"S{i+2}", "start": dt_list[i], "end": dt_list[i+1]})
    intervals.append({"label": f"S{len(dt_list)+1}", "start": dt_list[-1], "end": None})
    return intervals


def assign_session(file_dt: datetime, intervals):
    if not intervals:
        return None
    for iv in intervals:
        s = iv["start"]
        e = iv["end"]
        if s is None and e is not None:  # (-inf, e)
            if file_dt < e:
                return iv["label"]
        elif s is not None and e is None:  # [s, +inf)
            if file_dt >= s:
                return iv["label"]
        else:  # [s, e)
            if s <= file_dt < e:
                return iv["label"]
    return None

# -------------------------------------------------------------

def summarize_sessions(rows, intervals):
    from collections import defaultdict
    buckets = defaultdict(list)
    for r in rows:
        if r['session']:
            buckets[r['session']].append(r)
    print('\n[SESSION SUMMARY]')
    label_order = [iv['label'] for iv in intervals]
    for label in label_order:
        arr = buckets.get(label, [])
        times = [r['dt'] for r in arr if r['dt']]
        if times:
            start = min(times).strftime('%H:%M:%S')
            end = max(times).strftime('%H:%M:%S')
        else:
            start = end = '-'
        print(f"{label}: files={len(arr)} time_range={start}~{end}")

# -------------------------------------------------------------

def process(files, intervals, date_str: str):
    base_dir = Path("merged_chunks_organized") / date_str
    rows = []
    excluded = []
    for f in files:
        dt = parse_filename_ts(f.name)
        label = assign_session(dt, intervals) if (dt and intervals) else "-"
        agent_unique = None
        status_ok = True
        reason = None
        has_required = None
        if label and label.startswith('S') and label in EXPECTED_AGENT_COUNTS:
            try:
                df = pd.read_csv(base_dir / f.name, usecols=['agent_id'])
                col = df['agent_id'].dropna()
                agent_unique = col.nunique()
                has_required = any(x == REQUIRED_AGENT_ID or str(x) == str(REQUIRED_AGENT_ID) for x in col.unique())
            except Exception:
                agent_unique = 0
                has_required = False
            expected = EXPECTED_AGENT_COUNTS[label]
            if agent_unique != expected:
                status_ok = False
                reason = 'expected_mismatch'
            elif not has_required:
                status_ok = False
                reason = 'missing_agent99'
        if status_ok:
            rows.append({'filename': f.name, 'dt': dt, 'session': label, 'agent_unique': agent_unique})
        else:
            excluded.append({'filename': f.name, 'session': label, 'agent_unique': agent_unique, 'reason': reason})
    print(f"[INFO] kept {len(rows)} / {len(rows)+len(excluded)} files after agent filter")
    summarize_sessions(rows, intervals)
    # フィルタ結果
    from collections import defaultdict
    drop_counts = defaultdict(int)
    reason_counts = defaultdict(int)
    for ex in excluded:
        drop_counts[ex['session']] += 1
        reason_counts[ex['reason']] += 1
    print('\n[AGENT_ID FILTER RESULT]')
    for iv in intervals:
        lab = iv['label']
        exp = EXPECTED_AGENT_COUNTS.get(lab, '-')
        dropped = drop_counts.get(lab, 0)
        kept = sum(1 for r in rows if r['session']==lab)
        total = kept + dropped
        print(f"{lab}: expected={exp} total={total} kept={kept} dropped={dropped}")
    if reason_counts:
        print('[DROP REASONS] ' + ' '.join(f"{k}={v}" for k,v in reason_counts.items()))
    # --- 追加: kept ファイル一覧をCSV保存 (wavelet解析用) ---
    # 保存ディレクトリ results/<date>/
    try:
        base_out = Path('results') / date_str
        base_out.mkdir(parents=True, exist_ok=True)
        if rows:
            out_keep = base_out / f"filtered_files_{date_str}.csv"
            pd.DataFrame(rows)[['filename','session','agent_unique']].to_csv(out_keep, index=False)
            print(f"[INFO] kept file list saved -> {out_keep}")
        if excluded:
            out_drop = base_out / f"filtered_files_dropped_{date_str}.csv"
            pd.DataFrame(excluded).to_csv(out_drop, index=False)
            print(f"[INFO] dropped file list saved -> {out_drop}")
    except Exception as e:
        print(f"[WARN] failed to write filtered file lists: {e}")

# -------------------------------------------------------------

def main(argv=None):
    if argv is None:
        argv = sys.argv
    date_str = argv[1] if len(argv) > 1 else DEFAULT_DATE
    session_starts = DEFAULT_SESSION_STARTS
    files = list_files(date_str)
    if not files:
        return 1
    intervals = build_session_intervals(date_str, session_starts)
    disp = []
    for iv in intervals:
        s = iv['start'].strftime('%H:%M') if iv['start'] else '-inf'
        e = iv['end'].strftime('%H:%M') if iv['end'] else '+inf'
        disp.append(f"{iv['label']}:{s}->{e}")
    print('[INFO] intervals: ' + ' | '.join(disp))
    process(files, intervals, date_str)
    return 0

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
