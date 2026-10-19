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

    def arm(self, timeout=120):
        """Single sequence; waits until the scope is READY (pre-trigger done: only then a trigger is accepted).
        Returns True if it got READY, False after the timeout."""
        self.scope.write("ACQ:STOPAFTER SEQUENCE")
        self.scope.write("ACQ:STATE RUN")
        print("Waiting for the scope to finish the pre-trigger (-> READY)...")
        last, t0 = None, time.time()
        while time.time() - t0 < timeout:
            state = self.scope.query("TRIGger:STATE?").strip().upper()
            if state != last:
                print(f"  trigger state: {state}")
                last = state
            if state == "READY":
                return True
            time.sleep(0.1)
        print(f"Scope not READY after {timeout} s (trigger state: {last})")
        return False

    def wait(self, timeout=600):
        """Wait until the acquisition is over; True if it finished in time."""
        print("Waiting for the acquisition to finish...")
        t0 = time.time()
        while int(self.scope.query("ACQ:STATE?")) != 0:
            if time.time() - t0 > timeout:
                return False
            time.sleep(0.1)
        return True

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
