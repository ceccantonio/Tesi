"""Tektronix MDO3104 oscilloscope through VISA.

    scope = Oscilloscope()
    scope.arm(); ...start the laser scan...; scope.wait(); data = scope.read()
"""
import time

import numpy as np
import pyvisa

VISA = "USB::0x0699::0x0408::C028813::INSTR"


class Oscilloscope:
    def __init__(self, visa=VISA, n_samples=None):
        """Open the scope and set it up to read the record (n_samples=None: all of it) as ASCII."""
        self.scope = pyvisa.ResourceManager().open_resource(visa)
        self.scope.timeout = 30000
        self.scope.write("HORizontal:MODE MANUAL")
        record = int(self.scope.query("HORizontal:RECORDLENGTH?"))
        self.scope.write("DATA:START 1")
        self.scope.write(f"DATA:STOP {min(record, n_samples or record)}")
        self.scope.write("DATA:WIDTH 1")
        self.scope.write("DATA:ENCdg ASCii")

    def _wait_until(self, condition, timeout):
        t0 = time.time()
        while not condition() and time.time() - t0 < timeout:
            time.sleep(0.1)

    def arm(self, timeout=120):
        """Single sequence; returns when the scope is READY (only then a trigger is accepted)."""
        self.scope.write("ACQ:STOPAFTER SEQUENCE")
        self.scope.write("ACQ:STATE RUN")
        self._wait_until(lambda: self.scope.query("TRIGger:STATE?").strip().upper() == "READY", timeout)

    def wait(self, timeout=600):
        """Return when the acquisition is over."""
        self._wait_until(lambda: int(self.scope.query("ACQ:STATE?")) == 0, timeout)

    def force_trigger(self):
        self.scope.write("TRIGger FORCe")

    def read(self, channels=("CH1", "CH2")):
        """Record of the channels in volts, array (N, n_channels)."""
        columns = []
        for ch in channels:
            self.scope.write(f"DATA:SOURCE {ch}")
            ymult, yoff, yzero = (float(self.scope.query(f"WFMPRE:{k}?")) for k in ("YMULT", "YOFF", "YZERO"))
            columns.append((np.array(self.scope.query("CURVE?").split(","), dtype=float) - yoff) * ymult + yzero)
            time.sleep(0.1)
        n = min(map(len, columns))
        return np.column_stack([c[:n] for c in columns])
