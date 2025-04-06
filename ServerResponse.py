# サーバー側で管理するパラメータ
omega = 3.14 * 3  # 周波数
kappa = 1.5       # フィードバックゲイン
alpha = 0.2       # 位相遅れ定数

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

def handle_parameter_request(sock, addr):
    global omega, kappa, alpha
    response = f"omega:{omega:.2f},kappa:{kappa:.2f},alpha:{alpha:.2f}"
    sock.sendto(response.encode(), addr)
    print(f"[INFO] Sent parameters to {addr}: {response}")