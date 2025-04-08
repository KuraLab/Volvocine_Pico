import socket
import struct
import time
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from Plotter import plot_chunks
from keyinput import check_key
import os  # フォルダ作成用にosモジュールをインポート
from ServerResponse import handle_handshake, handle_parameter_request  # 新しいモジュールをインポート
from ChunkProcessor import build_dataframe_for_chunk, merge_and_save_chunks  # 新しいモジュールをインポート


# ---------------------------
# 設定
# ---------------------------
UDP_PORT = 5000
BUFFER_SIZE = 1024
SOCKET_TIMEOUT = 5.0
CHUNK_TIMEOUT = 10.0

STRUCT_FORMAT = "<6B"  # micros24 (3バイト), a0, a1, a2
RECORD_SIZE = struct.calcsize(STRUCT_FORMAT)  # 6 bytes

SAVE_FOLDER = "saved_chunks"  # 保存用フォルダ名

# 保存用フォルダを作成（存在しない場合のみ）
if not os.path.exists(SAVE_FOLDER):
    os.makedirs(SAVE_FOLDER)

agent_buffers = {}  # agent_id -> (chunk_data, send_micros_list, recv_time_list)
agent_lastrecv_time = {}
current_chunk_files = []

# サーバー側で管理するパラメータ
omega = 3.14 * 3  # 周波数
kappa = 1.5   # フィードバックゲイン
alpha = 0.2   # 位相遅れ定数

# ---------------------------
# チャンク処理
# ---------------------------
def is_valid_log_packet(data):
    # ダミー: agent_id==0 かつ payloadがない（最低1レコード＝5バイト未満）
    return len(data) >= 10 and data[0] != 0

# ---------------------------
# メイン受信ループ
# ---------------------------
def main():
    print(f"[INFO] Start listening UDP:{UDP_PORT}")
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind(("0.0.0.0", UDP_PORT))
    sock.settimeout(SOCKET_TIMEOUT)

    try:
        while True:
            try:
                data, addr = sock.recvfrom(BUFFER_SIZE)
                recv_time = time.time()

                # パラメータリクエストの処理
                if data.startswith(b"REQUEST_PARAMS"):  # パラメータリクエストの識別文字列
                    handle_parameter_request(sock, data, addr)
                    continue

                # ハンドシェイクメッセージの処理
                if data.startswith(b"HELLO"):  # バイト列で比較
                    handle_handshake(sock, data, addr)
                    continue

                if not is_valid_log_packet(data):
                    print(f"[INFO] Ignored dummy or malformed packet from {addr}, length={len(data)}")
                    continue
                elif len(data) < 5:
                    print(f"[WARN] Short packet from {addr}")
                    continue

                agent_id = data[0]
                send_micros = struct.unpack("<I", data[1:5])[0]
                raw = data[5:]

                if len(raw) % RECORD_SIZE != 0:
                    print(f"[WARN] Invalid record size from {addr}")
                    continue

                if agent_id in agent_lastrecv_time:
                    if (recv_time - agent_lastrecv_time[agent_id]) > CHUNK_TIMEOUT:
                        print(f"[INFO] Agent {agent_id} chunk timeout.")
                        chunk_data, send_list, recv_list = agent_buffers[agent_id]
                        if chunk_data:
                            _, saved_file = build_dataframe_for_chunk(agent_id, chunk_data, send_list, recv_list)
                            if saved_file:
                                current_chunk_files.append(saved_file)  # 保存されたファイルを追跡
                        agent_buffers[agent_id] = ([], [], [])

                if agent_id not in agent_buffers:
                    agent_buffers[agent_id] = ([], [], [])

                chunk_data, send_list, recv_list = agent_buffers[agent_id]

                offset_sec = recv_time - ((((send_micros >> 8) % 16777216) << 8) / 1e6)
                print(f"[DEBUG] Agent={agent_id}, send_micros={send_micros}, "
                      f"recv_time={recv_time:.6f}, offset_sec={offset_sec:.6f}")

                for i in range(len(raw) // RECORD_SIZE):
                    record = raw[i*RECORD_SIZE:(i+1)*RECORD_SIZE]
                    b0, b1, b2, a0, a1, a2 = struct.unpack(STRUCT_FORMAT, record)
                    micros24 = b0 | (b1 << 8) | (b2 << 16)
                    # 3バイトのタイムスタンプを32ビットに拡張
                    micros32 = micros24 & 0xFFFFFF  # 下位24ビットを取得
                    chunk_data.append((micros32, a0, a1, a2))

                send_list.append(send_micros)
                recv_list.append(recv_time)
                agent_buffers[agent_id] = (chunk_data, send_list, recv_list)
                # 最後のレコードの micros24 を取得してACK送信
                if len(raw) >= RECORD_SIZE:
                    last_record = raw[-RECORD_SIZE:]
                    b0, b1, b2, *_ = struct.unpack(STRUCT_FORMAT, last_record)
                    last_micros24 = b0 | (b1 << 8) | (b2 << 16)

                    ack = bytearray()
                    ack.append(agent_id)
                    ack += last_micros24.to_bytes(3, 'little')
                    sock.sendto(ack, addr)
                agent_lastrecv_time[agent_id] = recv_time

            except socket.timeout:
                pass

            key = check_key()
            if key in ('\r', '\n'):
                print("[INFO] Manual chunk flush.")
                for ag_id in list(agent_buffers.keys()):
                    data, send_list, recv_list = agent_buffers[ag_id]
                    if data:
                        _, saved_file = build_dataframe_for_chunk(ag_id, data, send_list, recv_list)
                        if saved_file:
                            current_chunk_files.append(saved_file)  # 保存されたファイルを追跡
                    agent_buffers[ag_id] = ([], [], [])
                merged_path = merge_and_save_chunks(current_chunk_files)
                print("[INFO] Merged and saved chunks.")
                plot_chunks(merged_path)
                current_chunk_files.clear()
                print("[DEBUG] current_chunk_files cleared.")

    except KeyboardInterrupt:
        print("[INFO] Interrupted by user.")

    finally:
        sock.close()
        print("[INFO] Socket closed.")

        for ag_id, (data, send_list, recv_list) in agent_buffers.items():
            if data:
                build_dataframe_for_chunk(ag_id, data, send_list, recv_list)

        plot_chunks(current_chunk_files)
        merge_and_save_chunks(current_chunk_files)
        current_chunk_files.clear()
        print("[DEBUG] current_chunk_files cleared.")
        print("[INFO] Exit complete.")

if __name__ == "__main__":
    main()