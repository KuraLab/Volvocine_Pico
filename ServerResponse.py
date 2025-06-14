# サーバー側で管理するパラメータ
omega = 3.14 * 3  # 周波数
kappa = 10       # フィードバックゲイン
alpha = -3.14*0.5       # 位相遅れ定数
servo_center = 110.0  # サーボ中心角度
servo_amplitude = 60.0 # サーボ振幅
stop_agent_id = 2      # 停止対象のエージェントID (0の場合はどのも停止しない等を意味づけることも可能)
stop_delay_seconds = 100 # 停止までの秒数

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

            # サーバー側のパラメータを送信
            response = f"omega:{omega:.2f},kappa:{kappa:.2f},alpha:{alpha:.2f},center:{servo_center:.1f},amplitude:{servo_amplitude:.1f},stop_id:{stop_agent_id},stop_delay:{stop_delay_seconds}"
            sock.sendto(response.encode('utf-8'), addr)
            print(f"[INFO] Sent parameters to Agent ID: {agent_id}, Voltage: {voltage:.2f}: {response}")
            return agent_id

        except (IndexError, ValueError) as e:
            print(f"[ERROR] Failed to parse parameter request: {request_str}")
            print(f"[ERROR] {e}")
    else:
        print(f"[WARN] Invalid parameter request from {addr}: {request_str}")