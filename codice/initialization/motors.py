"""PI C-863 motors over serial (VISA): X = ASRL4 (COM4), Y = ASRL5 (COM5).

    motors = Motors()
    motors.move(92.3, 146.5)       # True if it got there
    motors.x.pos                   # current X in mm

Close PI MicroMove first: a serial port can be used by one program only.
"""
import time

import pyvisa

PORT_X, PORT_Y = "ASRL4::INSTR", "ASRL5::INSTR"


class Axis:
    def __init__(self, port, baud=19200, tolerance_mm=0.01):
        self.dev = pyvisa.ResourceManager().open_resource(port)
        self.dev.baud_rate, self.dev.timeout = baud, 2000
        self.dev.read_termination = self.dev.write_termination = "\n"
        self.tolerance = tolerance_mm
        self.query("ERR?")                          # reads and clears an old error, if any

    def query(self, command):
        return self.dev.query(command).strip()

    def value(self, command):
        """Number from replies like '1=+0012.3456'."""
        return float(self.query(command).split("=")[-1])

    @property
    def pos(self):
        return self.value("POS? 1")

    @property
    def limits(self):
        return self.value("TMN? 1"), self.value("TMX? 1")

    def move(self, target, timeout=60):
        """Absolute move in mm. Does not move if the target is outside the limits; True if it arrived."""
        low, high = self.limits
        if not low <= target <= high:
            return False
        self.dev.write(f"MOV 1 {target:.4f}")
        time.sleep(0.1)
        t0 = time.time()
        while not (self.query("ONT? 1").endswith("1") and abs(self.pos - target) <= self.tolerance):
            if time.time() - t0 > timeout:
                return False
            time.sleep(0.1)
        return int(self.query("ERR?")) == 0


class Motors:
    def __init__(self, port_x=PORT_X, port_y=PORT_Y):
        self.x, self.y = Axis(port_x), Axis(port_y)

    def move(self, x, y):
        return self.x.move(x) and self.y.move(y)
