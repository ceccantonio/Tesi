"""Newport TLB-6700 laser through UsbDllWrap.dll (pythonnet, Windows only).

    laser = Laser()
    laser.set_range(1520, 1570)
    laser.start()
"""
import sys
from pathlib import Path

DLL = Path(__file__).resolve().parent / "UsbDllWrap.dll"
WAVELENGTH_START, WAVELENGTH_STOP = 1520.0, 1570.0     # nm, as in trigger_scan.m


class Laser:
    def __init__(self, dll=DLL, addr=1, device_id=0x100A):
        """Load the DLL and open the device (once per session; no CloseDevices, as in the MATLAB code)."""
        sys.path.append(str(Path(dll).parent))
        import clr
        clr.AddReference("UsbDllWrap")
        import Newport
        from System.Text import StringBuilder
        self.usb, self.addr, self.buf = Newport.USBComm.USB(), addr, StringBuilder(64)
        self.usb.OpenDevices(device_id)

    def query(self, command):
        """Send a command and read the reply (also for writes: otherwise the 'OK' replies pile up)."""
        self.usb.Query(self.addr, command, self.buf)
        return self.buf.ToString()

    def set_range(self, start=WAVELENGTH_START, stop=WAVELENGTH_STOP):
        """Scan range in nm; also resets the scan counter."""
        self.query(f"SOUR:WAVE:START {start}")
        self.query(f"SOUR:WAVE:STOP {stop}")
        self.stop()
        self.query("SOUR:WAVE:DESSCANS 1")

    def start(self):
        self.query("OUTP:SCAN:START")

    def stop(self):
        self.query("OUTP:SCAN:STOP")
