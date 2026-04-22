import argparse
import queue
import re
import sys
import threading
import time
from collections import deque

import matplotlib
matplotlib.use("TkAgg", force=True)
matplotlib.rcParams["path.simplify"] = False
matplotlib.rcParams["agg.path.chunksize"] = 0
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
try:
    import serial  # type: ignore
    from serial.tools import list_ports  # type: ignore
except ModuleNotFoundError:
    serial = None
    list_ports = None


LINE_RE = re.compile(
    r"I=(?P<current>[-+]?\d*\.?\d+)\s*mA,\s*"
    r"Vbus=(?P<vbus>[-+]?\d*\.?\d+)\s*V,\s*"
    r"P=(?P<power_mw>[-+]?\d*\.?\d+)\s*mW"
    r"(?:\s*\((?P<power_w>[-+]?\d*\.?\d+)\s*W\))?"
)


def parse_line(line: str):
    # Fast path for compact CSV: current_mA,vbus_V,power_mW[,timestamp_us[,a1_raw,a2_raw]]
    if "," in line and "I=" not in line:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 3:
            try:
                current_ma = float(parts[0])
                vbus_v = float(parts[1])
                power_mw = float(parts[2])
                src_time_s = None
                a1_raw = None
                a2_raw = None
                if len(parts) >= 4:
                    src_time_s = float(parts[3]) * 1e-6
                if len(parts) >= 6:
                    a1_raw = float(parts[4])
                    a2_raw = float(parts[5])
                return current_ma, vbus_v, power_mw, src_time_s, a1_raw, a2_raw
            except ValueError:
                pass

    match = LINE_RE.search(line)
    if not match:
        return None

    current_ma = float(match.group("current"))
    vbus_v = float(match.group("vbus"))
    power_w_text = match.group("power_w")

    if power_w_text is not None:
        power_mw = float(power_w_text) * 1000.0
    else:
        power_mw = float(match.group("power_mw"))

    return current_ma, vbus_v, power_mw, None, None, None


def choose_port(port_arg: str | None, list_only: bool) -> str | None:
    if list_ports is None:
        print("[ERROR] pyserial is not installed in this Python environment.")
        print(f"[HINT] Current Python: {sys.executable}")
        print("[HINT] Install with: python -m pip install pyserial")
        return None

    ports = list(list_ports.comports())

    if list_only:
        if not ports:
            print("[INFO] No serial ports found")
            return None

        print("[INFO] Available serial ports:")
        for p in ports:
            print(f"  - {p.device}: {p.description}")
        return None

    if port_arg:
        return port_arg

    if len(ports) == 1:
        auto = ports[0].device
        print(f"[INFO] --port not provided. Auto-selected {auto}")
        return auto

    if not ports:
        print("[ERROR] No serial ports found. Connect device or specify --port COMx")
        return None

    def score_port(p):
        text = f"{p.device} {p.description} {getattr(p, 'hwid', '')}".lower()
        score = 0

        # Prefer USB CDC/UART adapters and de-prioritize legacy onboard serial ports.
        if "usb" in text or "ch340" in text or "cp210" in text or "ftdi" in text:
            score += 100
        if "arduino" in text or "serial device" in text:
            score += 50

        dev_lower = p.device.lower()
        if dev_lower.startswith("/dev/ttyacm") or dev_lower.startswith("/dev/ttyusb"):
            score += 120
        if dev_lower.startswith("/dev/ttys"):
            score -= 300
        if dev_lower.startswith("/dev/cu.usb"):
            score += 120

        if p.device.upper() == "COM1":
            score -= 200
        return score

    best = max(ports, key=score_port)

    print("[WARN] Multiple serial ports found. Auto-selecting the most likely USB port:")
    for p in ports:
        print(f"  - {p.device}: {p.description}")
    print(f"[INFO] Selected {best.device}")
    return best.device


def open_serial_port(
    port: str,
    baudrate: int,
    timeout: float,
    retries: int = 5,
    retry_delay: float = 0.4,
):
    if serial is None:
        print("[ERROR] pyserial is not installed in this Python environment.")
        print(f"[HINT] Current Python: {sys.executable}")
        print("[HINT] Install with: python -m pip install pyserial")
        return None

    attempts = max(1, retries + 1)
    for attempt in range(1, attempts + 1):
        try:
            return serial.Serial(port, baudrate, timeout=timeout)
        except serial.SerialException as exc:
            msg = str(exc)
            print(f"[ERROR] Failed to open {port}: {msg}")

            access_denied = (
                "PermissionError(13" in msg
                or "Access is denied" in msg
                or "access is denied" in msg
                or "\u30a2\u30af\u30bb\u30b9\u304c\u62d2\u5426" in msg
            )

            if access_denied:
                print("[HINT] The port is likely in use by another app.")
                print("[HINT] Close Arduino IDE/Serial Monitor/Thonny/TeraTerm and retry.")
                print("[HINT] If needed, unplug/replug USB and run with --port COMx.")

            busy_on_linux = "Device or resource busy" in msg or "Errno 16" in msg
            if busy_on_linux:
                print(f"[HINT] Linux detected port busy. Check holder with: lsof {port}")
                print(f"[HINT] You can also try: fuser -v {port}")
                print("[HINT] If ModemManager grabs the device, stop it temporarily and retry.")

            if busy_on_linux and attempt < attempts:
                print(
                    f"[INFO] Retrying open ({attempt}/{attempts - 1}) in {retry_delay:.1f}s..."
                )
                time.sleep(retry_delay)
                continue

            break

    return None


def serial_reader_loop(
    ser,
    stop_event: threading.Event,
    line_queue: queue.Queue[tuple[float, str]],
):
    while not stop_event.is_set():
        try:
            raw = ser.readline()
        except serial.SerialException:
            break

        if not raw:
            continue

        received_at = time.monotonic()
        text = raw.decode("utf-8", errors="ignore").strip()
        if text:
            line_queue.put((received_at, text))


def main():
    if serial is None:
        print("[ERROR] Missing dependency: pyserial")
        print(f"[HINT] Current Python: {sys.executable}")
        print("[HINT] Install with: python -m pip install pyserial")
        return

    parser = argparse.ArgumentParser(
        description="Plot INA226 current/voltage/power from serial logs in real time."
    )
    parser.add_argument("--port", default=None, help="Serial port (e.g. COM5)")
    parser.add_argument("--baudrate", type=int, default=921600, help="Baudrate")
    parser.add_argument(
        "--window",
        type=int,
        default=0,
        help="Max retained points (0 = unlimited; recommended for fixed x-span)",
    )
    parser.add_argument(
        "--interval-ms", type=int, default=10, help="Plot update interval in milliseconds"
    )
    parser.add_argument(
        "--x-span-sec",
        type=float,
        default=3.0,
        help="Fixed x-axis span in seconds",
    )
    parser.add_argument(
        "--avg-sec",
        type=float,
        default=0.0,
        help="Time window in seconds for moving-average overlay (0 to disable)",
    )
    parser.add_argument(
        "--autoscale-every",
        type=int,
        default=0,
        help="Recompute y-axis autoscale every N frames (0 to disable)",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print lines that fail to parse"
    )
    parser.add_argument(
        "--list-ports", action="store_true", help="List serial ports and exit"
    )
    args = parser.parse_args()

    port = choose_port(args.port, args.list_ports)
    if port is None:
        return

    ser = open_serial_port(port, args.baudrate, timeout=0.05)
    if ser is None:
        return

    ser.reset_input_buffer()

    print(f"[INFO] Opened {port} @ {args.baudrate} baud")
    print(f"[INFO] Matplotlib backend: {matplotlib.get_backend()}")

    point_limit = args.window if args.window > 0 else None
    t_data = deque(maxlen=point_limit)
    i_data = deque(maxlen=point_limit)
    v_data = deque(maxlen=point_limit)
    p_data = deque(maxlen=point_limit)
    dt_data = deque(maxlen=point_limit)
    a1_data = deque(maxlen=point_limit)
    a2_data = deque(maxlen=point_limit)
    i_avg_data = deque(maxlen=point_limit)
    v_avg_data = deque(maxlen=point_limit)
    p_avg_data = deque(maxlen=point_limit)
    dt_avg_data = deque(maxlen=point_limit)
    a1_avg_data = deque(maxlen=point_limit)
    a2_avg_data = deque(maxlen=point_limit)

    avg_window = deque()
    avg_sum_i = 0.0
    avg_sum_v = 0.0
    avg_sum_p = 0.0
    avg_sum_dt = 0.0
    avg_sum_a1 = 0.0
    avg_sum_a2 = 0.0
    line_queue: queue.Queue[tuple[float, str]] = queue.Queue()
    stop_event = threading.Event()
    reader_thread = threading.Thread(
        target=serial_reader_loop,
        args=(ser, stop_event, line_queue),
        daemon=True,
    )
    reader_thread.start()

    frame_count = 0
    start = time.monotonic()
    device_time_base_s = None
    prev_now = None

    fig, axes = plt.subplots(6, 1, figsize=(10, 12), sharex=True)
    (line_i,) = axes[0].plot([], [], lw=2, color="tab:blue")
    (line_v,) = axes[1].plot([], [], lw=2, color="tab:green")
    (line_p,) = axes[2].plot([], [], lw=2, color="tab:red")
    (line_dt,) = axes[3].plot([], [], lw=2, color="tab:purple")
    (line_a1,) = axes[4].plot([], [], lw=2, color="tab:brown")
    (line_a2,) = axes[5].plot([], [], lw=2, color="tab:gray")
    (line_i_avg,) = axes[0].plot([], [], lw=2, color="tab:orange", alpha=0.9)
    (line_v_avg,) = axes[1].plot([], [], lw=2, color="tab:orange", alpha=0.9)
    (line_p_avg,) = axes[2].plot([], [], lw=2, color="tab:orange", alpha=0.9)
    (line_dt_avg,) = axes[3].plot([], [], lw=2, color="tab:orange", alpha=0.9)
    (line_a1_avg,) = axes[4].plot([], [], lw=2, color="tab:orange", alpha=0.9)
    (line_a2_avg,) = axes[5].plot([], [], lw=2, color="tab:orange", alpha=0.9)

    axes[0].set_ylabel("Current [mA]")
    axes[1].set_ylabel("Vbus [V]")
    axes[2].set_ylabel("Power [mW]")
    axes[3].set_ylabel("dt [ms]")
    axes[4].set_ylabel("A1 raw")
    axes[5].set_ylabel("A2 raw")
    axes[5].set_xlabel("Time [s]")
    axes[1].set_ylim(0.0, 8.0)
    axes[4].set_ylim(0.0, 4095.0)
    axes[5].set_ylim(0.0, 4095.0)
    axes[5].set_xlim(0.0, max(1e-6, args.x_span_sec))

    for ax in axes:
        ax.grid(True, alpha=0.3)

    if args.avg_sec > 0:
        axes[0].legend([line_i, line_i_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")
        axes[1].legend([line_v, line_v_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")
        axes[2].legend([line_p, line_p_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")
        axes[3].legend([line_dt, line_dt_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")
        axes[4].legend([line_a1, line_a1_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")
        axes[5].legend([line_a2, line_a2_avg], ["Raw", f"Avg {args.avg_sec:g}s"], loc="upper left")

    fig.suptitle("INA226 Live Monitor")

    def update(_frame):
        nonlocal frame_count
        nonlocal avg_sum_i, avg_sum_v, avg_sum_p, avg_sum_dt, avg_sum_a1, avg_sum_a2
        nonlocal device_time_base_s
        nonlocal prev_now

        while True:
            try:
                received_at, text = line_queue.get_nowait()
            except queue.Empty:
                break

            if text.startswith("[PROF]"):
                print(text)
                continue

            parsed = parse_line(text)
            if parsed is None:
                if args.verbose:
                    print(f"[SKIP] {text}")
                continue

            current_ma, vbus_v, power_mw, source_time_s, a1_raw, a2_raw = parsed
            if source_time_s is not None:
                if device_time_base_s is None:
                    device_time_base_s = source_time_s
                now = source_time_s - device_time_base_s
                if now < -1.0:
                    # Handle wraparound/reboot by re-basing source time.
                    device_time_base_s = source_time_s
                    now = 0.0
            else:
                now = received_at - start

            # Guard against out-of-order timestamps that can create apparent gaps.
            if t_data and now <= t_data[-1]:
                now = t_data[-1] + 1e-6

            if prev_now is None:
                dt_ms = 0.0
            else:
                dt_ms = (now - prev_now) * 1000.0
            prev_now = now

            t_data.append(now)
            i_data.append(current_ma)
            v_data.append(vbus_v)
            p_data.append(power_mw)
            dt_data.append(dt_ms)
            a1_value = current_ma if a1_raw is None else a1_raw
            a2_value = vbus_v if a2_raw is None else a2_raw
            a1_data.append(a1_value)
            a2_data.append(a2_value)

            if args.avg_sec > 0:
                avg_window.append((now, current_ma, vbus_v, power_mw, dt_ms, a1_value, a2_value))
                avg_sum_i += current_ma
                avg_sum_v += vbus_v
                avg_sum_p += power_mw
                avg_sum_dt += dt_ms
                avg_sum_a1 += a1_value
                avg_sum_a2 += a2_value

                cutoff = now - args.avg_sec
                while avg_window and avg_window[0][0] < cutoff:
                    _, old_i, old_v, old_p, old_dt, old_a1, old_a2 = avg_window.popleft()
                    avg_sum_i -= old_i
                    avg_sum_v -= old_v
                    avg_sum_p -= old_p
                    avg_sum_dt -= old_dt
                    avg_sum_a1 -= old_a1
                    avg_sum_a2 -= old_a2

                count = len(avg_window)
                if count > 0:
                    i_avg_data.append(avg_sum_i / count)
                    v_avg_data.append(avg_sum_v / count)
                    p_avg_data.append(avg_sum_p / count)
                    dt_avg_data.append(avg_sum_dt / count)
                    a1_avg_data.append(avg_sum_a1 / count)
                    a2_avg_data.append(avg_sum_a2 / count)
                else:
                    i_avg_data.append(current_ma)
                    v_avg_data.append(vbus_v)
                    p_avg_data.append(power_mw)
                    dt_avg_data.append(dt_ms)
                    a1_avg_data.append(a1_value)
                    a2_avg_data.append(a2_value)
            else:
                i_avg_data.append(current_ma)
                v_avg_data.append(vbus_v)
                p_avg_data.append(power_mw)
                dt_avg_data.append(dt_ms)
                a1_avg_data.append(a1_value)
                a2_avg_data.append(a2_value)

        if not t_data:
            return (
                line_i,
                line_v,
                line_p,
                line_dt,
                line_a1,
                line_a2,
                line_i_avg,
                line_v_avg,
                line_p_avg,
                line_dt_avg,
                line_a1_avg,
                line_a2_avg,
            )

        x = list(t_data)
        line_i.set_data(x, list(i_data))
        line_v.set_data(x, list(v_data))
        line_p.set_data(x, list(p_data))
        line_dt.set_data(x, list(dt_data))
        line_a1.set_data(x, list(a1_data))
        line_a2.set_data(x, list(a2_data))
        if args.avg_sec > 0:
            line_i_avg.set_data(x, list(i_avg_data))
            line_v_avg.set_data(x, list(v_avg_data))
            line_p_avg.set_data(x, list(p_avg_data))
            line_dt_avg.set_data(x, list(dt_avg_data))
            line_a1_avg.set_data(x, list(a1_avg_data))
            line_a2_avg.set_data(x, list(a2_avg_data))
        else:
            line_i_avg.set_data([], [])
            line_v_avg.set_data([], [])
            line_p_avg.set_data([], [])
            line_dt_avg.set_data([], [])
            line_a1_avg.set_data([], [])
            line_a2_avg.set_data([], [])

        # Keep a fixed-width time window on the x-axis.
        x_max = x[-1]
        x_span = max(1e-6, args.x_span_sec)
        x_min = max(0.0, x_max - x_span)
        if x_max <= x_min:
            x_max = x_min + 1e-6

        # If point_limit is unlimited, keep at least the current visible span.
        if point_limit is None:
            while t_data and t_data[0] < x_min:
                t_data.popleft()
                i_data.popleft()
                v_data.popleft()
                p_data.popleft()
                dt_data.popleft()
                a1_data.popleft()
                a2_data.popleft()
                i_avg_data.popleft()
                v_avg_data.popleft()
                p_avg_data.popleft()
                dt_avg_data.popleft()
                a1_avg_data.popleft()
                a2_avg_data.popleft()

            x = list(t_data)
            line_i.set_data(x, list(i_data))
            line_v.set_data(x, list(v_data))
            line_p.set_data(x, list(p_data))
            line_dt.set_data(x, list(dt_data))
            line_a1.set_data(x, list(a1_data))
            line_a2.set_data(x, list(a2_data))
            if args.avg_sec > 0:
                line_i_avg.set_data(x, list(i_avg_data))
                line_v_avg.set_data(x, list(v_avg_data))
                line_p_avg.set_data(x, list(p_avg_data))
                line_dt_avg.set_data(x, list(dt_avg_data))
                line_a1_avg.set_data(x, list(a1_avg_data))
                line_a2_avg.set_data(x, list(a2_avg_data))

        axes[5].set_xlim(x_min, x_max)

        frame_count += 1
        if args.autoscale_every > 0 and frame_count % args.autoscale_every == 0:
            for ax in axes:
                ax.relim()
                ax.autoscale_view(scalex=False, scaley=True)
            axes[1].set_ylim(0.0, 8.0)
            axes[4].set_ylim(0.0, 4095.0)
            axes[5].set_ylim(0.0, 4095.0)

        return (
            line_i,
            line_v,
            line_p,
            line_dt,
            line_a1,
            line_a2,
            line_i_avg,
            line_v_avg,
            line_p_avg,
            line_dt_avg,
            line_a1_avg,
            line_a2_avg,
        )

    def cleanup_serial():
        stop_event.set()

        # Cancel any pending read first so the thread can exit quickly on Windows.
        try:
            ser.cancel_read()
        except Exception:
            pass

        if ser.is_open:
            ser.close()

        reader_thread.join(timeout=0.3)

    def handle_close(_event):
        cleanup_serial()
        if ser.is_open:
            return
        print("[INFO] Serial port closed")

    fig.canvas.mpl_connect("close_event", handle_close)
    plt.ioff()
    anim = FuncAnimation(
        fig,
        update,
        interval=args.interval_ms,
        blit=False,
        cache_frame_data=False,
    )
    plt.tight_layout()
    try:
        # VS Code debug launcher環境では block=True だけだと描画されない場合があるため、
        # 非ブロッキング表示 + 手動イベントループで確実にウィンドウ更新を回す。
        plt.show(block=False)
        while plt.fignum_exists(fig.number):
            plt.pause(0.05)
    except KeyboardInterrupt:
        plt.close(fig)
    finally:
        cleanup_serial()

    # Keep a reference until window closes to avoid animation GC warnings.
    _ = anim


if __name__ == "__main__":
    main()
