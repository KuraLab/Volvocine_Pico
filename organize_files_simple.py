#!/usr/bin/env python3
"""
merged_chunks内のファイルを日付ごとのフォルダに分ける簡易スクリプト

使用例:
python organize_files_simple.py        # コピーモード
python organize_files_simple.py move   # 移動モード
移動モードが一番いい
"""


import os
import shutil
import re
import sys
from pathlib import Path

def organize_files(move_mode=False):
    """ファイルを日付ごとに整理"""
    source_dir = Path("merged_chunks")
    output_dir = Path("merged_chunks_organized")
    
    if not source_dir.exists():
        print(f"エラー: {source_dir} が見つかりません")
        return
    
    # 出力ディレクトリ作成
    output_dir.mkdir(exist_ok=True)
    
    # ファイルパターン
    pattern = r'merged_(\d{4})(\d{2})(\d{2})_\d{6}\.csv'
    
    processed = 0
    errors = 0
    
    print(f"{'移動' if move_mode else 'コピー'}モードで実行中...")
    
    for file_path in source_dir.glob("merged_*.csv"):
        match = re.match(pattern, file_path.name)
        if match:
            year, month, day = match.groups()
            date_folder = f"{year}-{month}-{day}"
            
            # 日付フォルダを作成
            date_dir = output_dir / date_folder
            date_dir.mkdir(exist_ok=True)
            
            # ファイルを移動/コピー
            dest_path = date_dir / file_path.name
            
            try:
                if move_mode:
                    shutil.move(str(file_path), str(dest_path))
                else:
                    shutil.copy2(str(file_path), str(dest_path))
                
                print(f"{date_folder}: {file_path.name}")
                processed += 1
                
            except Exception as e:
                print(f"エラー ({file_path.name}): {e}")
                errors += 1
        else:
            print(f"スキップ (パターン不一致): {file_path.name}")
    
    print(f"\n完了: {processed}ファイル処理, {errors}エラー")

def show_summary():
    """日付別ファイル数の要約表示"""
    source_dir = Path("merged_chunks")
    pattern = r'merged_(\d{4})(\d{2})(\d{2})_\d{6}\.csv'
    
    date_counts = {}
    
    for file_path in source_dir.glob("merged_*.csv"):
        match = re.match(pattern, file_path.name)
        if match:
            year, month, day = match.groups()
            date_key = f"{year}-{month}-{day}"
            date_counts[date_key] = date_counts.get(date_key, 0) + 1
    
    print("日付別ファイル数:")
    for date_key in sorted(date_counts.keys()):
        print(f"  {date_key}: {date_counts[date_key]}ファイル")
    
    print(f"合計: {sum(date_counts.values())}ファイル")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].lower() in ['move', 'm']:
        mode = True
        print("移動モードで実行します")
    else:
        mode = False
        print("コピーモードで実行します")
    
    print("\n=== 現在の状況 ===")
    show_summary()
    
    print(f"\n=== ファイル整理開始 ===")
    organize_files(mode)
