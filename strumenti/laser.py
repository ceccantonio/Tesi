"""Laser Newport TLB-6700 via UsbDllWrap.dll (pythonnet, solo Windows). Uso: las = Laser(); las.imposta(1520, 1570); las.avvia()"""
import os
import sys

DLL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "laser", "UsbDllWrap.dll")


class Laser:
    def __init__(self, dll=DLL, addr=1, device_id=0x100A):
        """Carica la DLL e apre il dispositivo (una volta per sessione; niente CloseDevices, come in MATLAB)."""
        sys.path.append(os.path.dirname(os.path.abspath(dll)))
        import clr
        clr.AddReference("UsbDllWrap")
        import Newport
        from System.Text import StringBuilder
        self.usb, self.addr, self.buf = Newport.USBComm.USB(), addr, StringBuilder(64)
        self.usb.OpenDevices(device_id)

    def q(self, comando):
        """Invia un comando e legge la risposta (anche per i comandi di scrittura: altrimenti gli 'OK' si accumulano)."""
        self.usb.Query(self.addr, comando, self.buf)
        return self.buf.ToString()

    def imposta(self, start=1520.0, stop=1570.0, reset=True):
        """Intervallo della scansione in nm; reset=True azzera anche il contatore delle scansioni."""
        self.q(f"SOUR:WAVE:START {start}")
        self.q(f"SOUR:WAVE:STOP {stop}")
        if reset:
            self.ferma()
            self.q("SOUR:WAVE:DESSCANS 1")

    def avvia(self):
        self.q("OUTP:SCAN:START")

    def ferma(self):
        self.q("OUTP:SCAN:STOP")
