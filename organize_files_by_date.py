#!/usr/bin/env python3
"""
merged_chunks内のファイルを日付ごとのフォルダに分けるスクリプト

ファイル名形式: merged_YYYYMMDD_HHMMSS.csv
出力: merged_chunks_organized/YYYY-MM-DD/ 以下にファイルを移動

使用例:
python organize_files_by_date.py
"""

import os
import shutil
import re
from pathlib import Path
from datetime import datetime

def extract_date_from_filename(filename):
    """
    ファイル名から日付を抽出する
    例: merged_20250406_041941.csv -> 2025-04-06
    """
    pattern = r'merged_(\d{4})(\d{2})(\d{2})_\d{6}\.csv'
    match = re.match(pattern, filename)
    if match:
        year, month, day = match.groups()
        return f"{year}-{month}-{day}"
    return None

def organize_files_by_date(source_dir, output_dir, move_files=False):
    """
    ファイルを日付ごとのフォルダに整理する
    
    Args:
        source_dir (str): 元のディレクトリ (merged_chunks)
        output_dir (str): 出力ディレクトリ (merged_chunks_organized)
        move_files (bool): Trueの場合は移動、Falseの場合はコピー
    """
    source_path = Path(source_dir)
    output_path = Path(output_dir)
    
    if not source_path.exists():
        print(f"エラー: {source_dir} が見つかりません")
        return
    
    # 出力ディレクトリを作成
    output_path.mkdir(exist_ok=True)
    
    # ファイルをスキャン
    csv_files = list(source_path.glob("merged_*.csv"))
    print(f"対象ファイル数: {len(csv_files)}")
    
    # 日付ごとにグループ化
    date_groups = {}
    unmatched_files = []
    
    for file_path in csv_files:
        filename = file_path.name
        date_str = extract_date_from_filename(filename)
        
        if date_str:
            if date_str not in date_groups:
                date_groups[date_str] = []
            date_groups[date_str].append(file_path)
        else:
            unmatched_files.append(file_path)
    
    print(f"日付が特定できたファイル: {sum(len(files) for files in date_groups.values())}")
    print(f"日付が特定できなかったファイル: {len(unmatched_files)}")
    print(f"対象日付数: {len(date_groups)}")
    
    # 各日付のフォルダを作成し、ファイルを移動/コピー
    for date_str, files in date_groups.items():
        date_dir = output_path / date_str
        date_dir.mkdir(exist_ok=True)
        
        print(f"\n{date_str}: {len(files)}ファイル")
        
        for file_path in files:
            dest_path = date_dir / file_path.name
            
            try:
                if move_files:
                    shutil.move(str(file_path), str(dest_path))
                    print(f"  移動: {file_path.name}")
                else:
                    shutil.copy2(str(file_path), str(dest_path))
                    print(f"  コピー: {file_path.name}")
            except Exception as e:
                print(f"  エラー ({file_path.name}): {e}")
    
    # マッチしなかったファイルの処理
    if unmatched_files:
        unmatched_dir = output_path / "unmatched"
        unmatched_dir.mkdir(exist_ok=True)
        print(f"\n未分類ファイル: {len(unmatched_files)}")
        
        for file_path in unmatched_files:
            dest_path = unmatched_dir / file_path.name
            try:
                if move_files:
                    shutil.move(str(file_path), str(dest_path))
                    print(f"  移動: {file_path.name}")
                else:
                    shutil.copy2(str(file_path), str(dest_path))
                    print(f"  コピー: {file_path.name}")
            except Exception as e:
                print(f"  エラー ({file_path.name}): {e}")

def show_date_summary(source_dir):
    """
    日付ごとのファイル数の要約を表示する
    """
    source_path = Path(source_dir)
    
    if not source_path.exists():
        print(f"エラー: {source_dir} が見つかりません")
        return
    
    csv_files = list(source_path.glob("merged_*.csv"))
    date_counts = {}
    
    for file_path in csv_files:
        date_str = extract_date_from_filename(file_path.name)
        if date_str:
            date_counts[date_str] = date_counts.get(date_str, 0) + 1
    
    print("日付ごとのファイル数:")
    print("-" * 30)
    for date_str in sorted(date_counts.keys()):
        print(f"{date_str}: {date_counts[date_str]}ファイル")
    
    print(f"\n合計: {sum(date_counts.values())}ファイル")
    print(f"日付数: {len(date_counts)}")

def main():
    """メイン関数"""
    # パス設定
    base_dir = Path(__file__).parent
    source_dir = base_dir / "merged_chunks"
    output_dir = base_dir / "merged_chunks_organized"
    
    print("=== ファイル整理スクリプト ===\n")
    
    # 1. 現在の状況を表示
    print("1. 現在の状況:")
    show_date_summary(source_dir)
    
    # 2. ユーザーに確認
    print(f"\n2. 整理設定:")
    print(f"   元フォルダ: {source_dir}")
    print(f"   出力フォルダ: {output_dir}")
    
    while True:
        choice = input("\nファイルを移動しますか？ (y=移動, c=コピー, n=キャンセル): ").lower().strip()
        if choice in ['y', 'yes']:
            move_files = True
            break
        elif choice in ['c', 'copy']:
            move_files = False
            break
        elif choice in ['n', 'no']:
            print("キャンセルしました。")
            return
        else:
            print("y, c, n のいずれかを入力してください。")
    
    # 3. ファイル整理実行
    print("\n3. ファイル整理を実行中...")
    organize_files_by_date(source_dir, output_dir, move_files)
    
    print("\n完了しました！")

if __name__ == "__main__":
    main()
