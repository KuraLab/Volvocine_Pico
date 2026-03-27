# サーバー側で管理するパラメータ
# エージェントIDごとのomega値を配列で管理
import math

omega_values = {
    4: 3.14 * 3.05   # エージェント1の周波数
}
default_omega = 3.14 * 2.95 # デフォルト周波数（未定義IDの場合）

kappa =0       # フィードバックゲイン
alpha = -3.14*0.6
servo_center = 90.0  # サーボ中心角度
servo_amplitude = 65.0 # サーボ振幅
stop_agent_id = 4      # 停止対象のエージェントID (0の場合はどのも停止しない等を意味づけることも可能)
stop_delay_seconds = 30000 # 停止までの秒数

# PRCのフーリエ係数（0..prc_harmonics を使用）
# z(psi) = Σ [ prc_a[n] * cos(n*psi) + prc_b[n] * sin(n*psi) ]
prc_harmonics = 10
prc_a = [0.0] * (prc_harmonics + 1)
prc_b = [0.0] * (prc_harmonics + 1)

# 既存の cos(psi-alpha) と等価な初期値（1次のみ）
prc_a[1] = math.cos(alpha)
prc_b[1] = math.sin(alpha)

# 例: 高調波を使う場合は以下を編集
# prc_a[2] = 0.10
# prc_b[3] = -0.05


def build_prc_payload():
    """PRC係数をレスポンス文字列へ展開"""
    fields = [f"prc_n:{prc_harmonics}"]
    for n in range(0, prc_harmonics + 1):
        fields.append(f"prc_a{n}:{prc_a[n]:.6f}")
        fields.append(f"prc_b{n}:{prc_b[n]:.6f}")
    return ",".join(fields)

def get_omega_for_agent(agent_id):
    """
    エージェントIDに応じたomega値を取得する関数
    """
    return omega_values.get(agent_id, default_omega)

def handle_handshake(sock, data, addr):
    """
    クライアントからのハンドシェイクメッセージに応答する関数。
    """
    handshake_message = "HELLO"
    try:
        if data.decode('utf-8') == handshake_message:
            response = "READY"
            sock.sendto(response.encode('utf-8'), addr)
            print(f"[INFO] Handshake response sent to {addr}")
    except UnicodeDecodeError:
        print(f"[WARN] Received non-UTF-8 data from {addr}, ignoring.")

def handle_parameter_request(sock, data, addr):
    """
    パラメータリクエストを処理し、デバッグ情報を表示
    """
    request_str = data.decode('utf-8')

    # リクエストデータを解析
    if request_str.startswith("REQUEST_PARAMS"):
        try:
            # デバッグ情報を解析
            parts = request_str.split(',')
            agent_id = int(parts[1].split(':')[1])  # id:<value>
            analog26 = int(parts[2].split(':')[1])  # analog26:<value>

            # analog26 を電圧値に変換
            voltage = (analog26 / 4095) * 3.3 * 2

            # エージェントIDに応じたomega値を取得
            omega = get_omega_for_agent(agent_id)

            # サーバー側のパラメータを送信
            response = (
                f"omega:{omega:.2f},kappa:{kappa:.2f},"
                f"center:{servo_center:.1f},amplitude:{servo_amplitude:.1f},"
                f"stop_id:{stop_agent_id},stop_delay:{stop_delay_seconds},"
                f"{build_prc_payload()}"
            )
            sock.sendto(response.encode('utf-8'), addr)
            print(f"[INFO] Sent parameters to Agent ID: {agent_id}, Omega: {omega:.2f}, Voltage: {voltage:.2f}: {response}")
            return agent_id

        except (IndexError, ValueError) as e:
            print(f"[ERROR] Failed to parse parameter request: {request_str}")
            print(f"[ERROR] {e}")
    else:
        print(f"[WARN] Invalid parameter request from {addr}: {request_str}")