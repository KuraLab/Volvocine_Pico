# keyinput.py
import sys
import os

if os.name == 'nt':
    import msvcrt

    def check_key():
        if msvcrt.kbhit():
            return msvcrt.getwch()
        return None
else:
    import select
    import tty
    import termios

    def check_key():
        dr, _, _ = select.select([sys.stdin], [], [], 0)
        if dr:
            old_settings = termios.tcgetattr(sys.stdin)
            try:
                tty.setcbreak(sys.stdin.fileno())
                return sys.stdin.read(1)
            finally:
                termios.tcsetattr(sys.stdin, termios.TCSADRAIN, old_settings)
        return None
