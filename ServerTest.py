import socket
import struct
import time
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from Plotter import plot_chunks
from ChunkSaver import merge_and_save_chunks
from keyinput import check_key
import os  # フォルダ作成用にosモジュールをインポート


# ---------------------------
# 設定
# ---------------------------
UDP_PORT = 5000
BUFFER_SIZE = 1024
SOCKET_TIMEOUT = 0.5
CHUNK_TIMEOUT = 1.0

STRUCT_FORMAT = "<HBBB"  # micros16, a0, a1, a2
RECORD_SIZE = struct.calcsize(STRUCT_FORMAT)  # 5 bytes

SAVE_FOLDER = "saved_chunks"  # 保存用フォルダ名

# 保存用フォルダを作成（存在しない場合のみ）
if not os.path.exists(SAVE_FOLDER):
    os.makedirs(SAVE_FOLDER)

agent_buffers = {}  # agent_id -> (chunk_data, send_micros_list, recv_time_list)
agent_lastrecv_time = {}
current_chunk_files = []

# ---------------------------
# チャンク処理
# ---------------------------
def build_dataframe_for_chunk(agent_id, chunk_data, chunk_send_micros, chunk_recv_times):
    global current_chunk_files

    if not chunk_data:
        return None

    # 平均オフセットを送信時刻の下位16bit再構成で算出
    wrapped_send_secs = [(((s >> 10) % 65536) << 10) / 1e6 for s in chunk_send_micros]
    offsets = [recv - send for send, recv in zip(wrapped_send_secs, chunk_recv_times)]
    offset = sum(offsets) / len(offsets)

    df = pd.DataFrame(chunk_data, columns=["micros16", "a0", "a1", "a2"])

    micros_list = df["micros16"].tolist()
    extended = [0] * len(micros_list)
    wrap_offset = 0
    prev = micros_list[0]
    extended[0] = prev
    for i in range(1, len(micros_list)):
        curr = micros_list[i]
        if curr < prev:
            wrap_offset += 65536
        extended[i] = curr + wrap_offset
        prev = curr

    df["micros32"] = extended
    df["micros32_raw"] = [val << 10 for val in extended]
    df["time_local_sec"] = [val / 1e6 for val in df["micros32_raw"]]
    df["time_pc_sec_abs"] = df["time_local_sec"] + offset

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    chunk_id = timestamp
    df["agent_id"] = agent_id
    df["chunk_id"] = chunk_id

    # 保存先を保存用フォルダに変更
    filename = os.path.join(SAVE_FOLDER, f"chunk_agent_{agent_id}_{timestamp}.csv")
    save_columns = [
        "time_pc_sec_abs", "micros32", "micros32_raw", "time_local_sec",
        "a0", "a1", "a2", "agent_id", "chunk_id"
    ]
    df.to_csv(filename, index=False, columns=save_columns)
    print(f"[INFO] Agent={agent_id}, chunk size={len(df)} -> Saved to {filename}")

    current_chunk_files.append(filename)
    print(f"[DEBUG] Added to current_chunk_files: {filename}")

    return df[["agent_id", "chunk_id", "time_pc_sec_abs", "a0", "a1", "a2"]]

def is_valid_log_packet(data):
    # ダミー: agent_id==0 かつ payloadがない（最低1レコード＝5バイト未満）
    return len(data) >= 10 and data[0] != 0

def handle_handshake(sock, data, addr):
    """
    クライアントからのハンドシェイクメッセージに応答する関数。
    """
    handshake_message = "HELLO"
    if data.decode('utf-8') == handshake_message:
        response = "READY"
        sock.sendto(response.encode('utf-8'), addr)
        print(f"[INFO] Handshake response sent to {addr}")

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

                # ハンドシェイクメッセージの処理
                if data.decode('utf-8') == "HELLO":
                    handle_handshake(sock, data, addr)
                    continue

                # デバッグログ: 受信データの内容を確認
                #print(f"[DEBUG] Received data from {addr}, length={len(data)}")

                if not is_valid_log_packet(data):
                    print(f"[INFO] Ignored dummy or malformed packet from {addr}, length={len(data)}")
                    continue
                elif len(data) < 5:
                    print(f"[WARN] Short packet from {addr}")
                    continue

                agent_id = data[0]
                send_micros = struct.unpack("<I", data[1:5])[0]
                raw = data[5:]

                # デバッグログ: パケット解析結果を確認
                #print(f"[DEBUG] Agent={agent_id}, send_micros={send_micros}, raw_length={len(raw)}")

                if len(raw) % RECORD_SIZE != 0:
                    print(f"[WARN] Invalid record size from {addr}")
                    continue

                if agent_id in agent_lastrecv_time:
                    if (recv_time - agent_lastrecv_time[agent_id]) > CHUNK_TIMEOUT:
                        print(f"[INFO] Agent {agent_id} chunk timeout.")
                        chunk_data, send_list, recv_list = agent_buffers[agent_id]
                        if chunk_data:
                            build_dataframe_for_chunk(agent_id, chunk_data, send_list, recv_list)
                        agent_buffers[agent_id] = ([], [], [])

                if agent_id not in agent_buffers:
                    agent_buffers[agent_id] = ([], [], [])

                chunk_data, send_list, recv_list = agent_buffers[agent_id]

                offset_sec = recv_time - ((((send_micros >> 10) % 65536) << 10) / 1e6)
                print(f"[DEBUG] Agent={agent_id}, send_micros={send_micros}, "
                      f"recv_time={recv_time:.6f}, offset_sec={offset_sec:.6f}")

                for i in range(len(raw) // RECORD_SIZE):
                    record = raw[i*RECORD_SIZE:(i+1)*RECORD_SIZE]
                    micros16, a0, a1, a2 = struct.unpack(STRUCT_FORMAT, record)
                    chunk_data.append((micros16, a0, a1, a2))

                send_list.append(send_micros)
                recv_list.append(recv_time)
                agent_buffers[agent_id] = (chunk_data, send_list, recv_list)
                agent_lastrecv_time[agent_id] = recv_time

            except socket.timeout:
                pass

            key = check_key()
            if key in ('\r', '\n'):
                print("[INFO] Manual chunk flush.")
                for ag_id in list(agent_buffers.keys()):
                    data, send_list, recv_list = agent_buffers[ag_id]
                    if data:
                        build_dataframe_for_chunk(ag_id, data, send_list, recv_list)
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