import socket
import threading
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# UDP設定
UDP_RECEIVE_PORT = 5000
BUFFER_SIZE = 1024
data_history = []
MAX_DATA_POINTS = 10000
PLOT_WINDOW = 10  # 秒

# データの各フィールド用
recv_times = []
digitals = []
analog0s = []
analog1s = []
analog2s = []
looptimes = []

# UDP受信スレッド
def udp_server_thread():
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind(("0.0.0.0", UDP_RECEIVE_PORT))
    print(f"[INFO] UDP server listening on port {UDP_RECEIVE_PORT}")

    while True:
        data, addr = sock.recvfrom(BUFFER_SIZE)
        msg = data.decode(errors='ignore').strip()
        try:
            fields = {kv.split('=')[0]: kv.split('=')[1] for kv in msg.split(',')}
            recv_time = time.time()
            digitals.append(int(fields["digital"]))
            analog0s.append(float(fields["analog0"]))
            analog1s.append(int(fields["analog1"]))
            analog2s.append(int(fields["analog2"]))
            looptimes.append(int(fields["looptime"]))
            recv_times.append(recv_time)

            # 古いデータを捨てる
            if len(recv_times) > MAX_DATA_POINTS:
                for lst in [recv_times, digitals, analog0s, analog1s, analog2s, looptimes]:
                    lst.pop(0)
        except Exception as e:
            print(f"[WARN] Parse error: {e}, msg={msg}")

# グラフ初期化
fig, axs = plt.subplots(3, 1, figsize=(10, 8))
lines = []

for ax, label in zip(axs, ["Analog0 (電圧)", "Analog1 / Analog2 (Raw)", "Loop time (μs)"]):
    ax.set_title(label)
    ax.grid(True)
    line1, = ax.plot([], [], label="Analog1", color="tab:blue")
    line2, = ax.plot([], [], label="Analog2", color="tab:orange")
    lines.append((line1, line2))
    ax.legend(loc="upper right")

# グラフ更新
def update(frame):
    if not recv_times:
        return []

    t0 = time.time()
    times = [t - t0 for t in recv_times]

    # 各データをセット
    lines[0][0].set_data(times, analog0s)  # Analog0
    lines[0][1].set_data(times, [0]*len(analog0s))  # Dummy second line

    lines[1][0].set_data(times, analog1s)
    lines[1][1].set_data(times, analog2s)

    lines[2][0].set_data(times, looptimes)
    lines[2][1].set_data(times, [0]*len(looptimes))  # Dummy

    for ax in axs:
        ax.set_xlim(-PLOT_WINDOW, 0)

    axs[0].set_ylim(0, 7)  # Analog0 (0~6.6V想定)
    axs[1].set_ylim(0, 4095)
    axs[2].set_ylim(0, 3000)

    return [l for pair in lines for l in pair]

# メイン
def main():
    threading.Thread(target=udp_server_thread, daemon=True).start()
    ani = animation.FuncAnimation(fig, update, interval=100, blit=False)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
