import argparse
import queue
import re
import threading
import time
from collections import deque

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import serial  # type: ignore
from serial.tools import list_ports  # type: ignore


LINE_RE = re.compile(
    r"I=(?P<current>[-+]?\d*\.?\d+)\s*mA,\s*"
    r"Vbus=(?P<vbus>[-+]?\d*\.?\d+)\s*V,\s*"
    r"P=(?P<power_mw>[-+]?\d*\.?\d+)\s*mW"
    r"(?:\s*\((?P<power_w>[-+]?\d*\.?\d+)\s*W\))?"
)


def parse_line(line: str):
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

    return current_ma, vbus_v, power_mw


def choose_port(port_arg: str | None, list_only: bool) -> str | None:
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

        # Prefer USB/serial adapter ports, avoid legacy motherboard COM ports.
        if "usb" in text or "ch340" in text or "cp210" in text or "ftdi" in text:
            score += 100
        if "arduino" in text or "serial device" in text:
            score += 50
        if p.device.upper() == "COM1":
            score -= 200
        return score

    best = max(ports, key=score_port)

    print("[WARN] Multiple serial ports found. Auto-selecting the most likely USB port:")
    for p in ports:
        print(f"  - {p.device}: {p.description}")
    print(f"[INFO] Selected {best.device}")
    return best.device


def open_serial_port(port: str, baudrate: int, timeout: float):
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

        return None


def serial_reader_loop(ser, stop_event: threading.Event, line_queue: queue.Queue[str]):
    while not stop_event.is_set():
        try:
            raw = ser.readline()
        except serial.SerialException:
            break

        if not raw:
            continue

        text = raw.decode("utf-8", errors="ignore").strip()
        if text:
            line_queue.put(text)


def main():
    parser = argparse.ArgumentParser(
        description="Plot INA226 current/voltage/power from serial logs in real time."
    )
    parser.add_argument("--port", default=None, help="Serial port (e.g. COM5)")
    parser.add_argument("--baudrate", type=int, default=115200, help="Baudrate")
    parser.add_argument(
        "--window", type=int, default=100, help="Number of points kept on the graph"
    )
    parser.add_argument(
        "--interval-ms", type=int, default=40, help="Plot update interval in milliseconds"
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

    t_data = deque(maxlen=args.window)
    i_data = deque(maxlen=args.window)
    v_data = deque(maxlen=args.window)
    p_data = deque(maxlen=args.window)
    line_queue: queue.Queue[str] = queue.Queue()
    stop_event = threading.Event()
    reader_thread = threading.Thread(
        target=serial_reader_loop,
        args=(ser, stop_event, line_queue),
        daemon=True,
    )
    reader_thread.start()

    frame_count = 0
    start = time.monotonic()

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    (line_i,) = axes[0].plot([], [], lw=2, color="tab:blue")
    (line_v,) = axes[1].plot([], [], lw=2, color="tab:green")
    (line_p,) = axes[2].plot([], [], lw=2, color="tab:red")

    axes[0].set_ylabel("Current [mA]")
    axes[1].set_ylabel("Vbus [V]")
    axes[2].set_ylabel("Power [mW]")
    axes[2].set_xlabel("Time [s]")

    for ax in axes:
        ax.grid(True, alpha=0.3)

    fig.suptitle("INA226 Live Monitor")

    def update(_frame):
        nonlocal frame_count

        lines_processed = 0
        while lines_processed < 1000:
            try:
                text = line_queue.get_nowait()
            except queue.Empty:
                break

            parsed = parse_line(text)
            if parsed is None:
                if args.verbose:
                    print(f"[SKIP] {text}")
                continue

            current_ma, vbus_v, power_mw = parsed
            now = time.monotonic() - start

            t_data.append(now)
            i_data.append(current_ma)
            v_data.append(vbus_v)
            p_data.append(power_mw)
            lines_processed += 1

        if not t_data:
            return line_i, line_v, line_p

        x = list(t_data)
        line_i.set_data(x, list(i_data))
        line_v.set_data(x, list(v_data))
        line_p.set_data(x, list(p_data))

        # Update x-limits each frame for scrolling view, but throttle y autoscale.
        x_min = x[0]
        x_max = x[-1] if x[-1] > x_min else x_min + 1e-6
        axes[2].set_xlim(x_min, x_max)

        frame_count += 1
        if frame_count % 5 == 0:
            for ax in axes:
                ax.relim()
                ax.autoscale_view(scalex=False, scaley=True)

        return line_i, line_v, line_p

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
    anim = FuncAnimation(
        fig,
        update,
        interval=args.interval_ms,
        blit=False,
        cache_frame_data=False,
    )
    plt.tight_layout()
    try:
        plt.show()
    except KeyboardInterrupt:
        plt.close(fig)
    finally:
        cleanup_serial()

    # Keep a reference until window closes to avoid animation GC warnings.
    _ = anim


if __name__ == "__main__":
    main()
