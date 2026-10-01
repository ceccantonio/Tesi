"""
Driver minimale per il laser sintonizzabile Newport/New Focus TLB-6700,
basato sulla DLL .NET proprietaria (UsbDllWrap.dll) installata insieme al
software Newport, non su VISA/seriale (il TLB-6700 non espone una porta
COM: usa un driver USB Jungo + wrapper .NET).

Richiede:
    pip install pythonnet

Il set di comandi SCPI (SOURce:WAVElength, OUTPut:SCAN:START, ecc.) è
preso da un esempio pubblico funzionante con lo stesso strumento:
https://github.com/marvync/Laser-automation-NewFocus
Verificalo comunque contro il manuale del tuo TLB-6700 se qualcosa non
si comporta come previsto (firmware diversi possono avere differenze).

Uso tipico:

    laser = TLB6700(device_key="6700 SN1234")  # trovala con trova_tlb6700.py
    laser.open()
    laser.set_scan_limits(1550.0, 1560.0)
    laser.set_scan_speeds(forward=0.1, backward=5.0)
    laser.set_output(True)
    laser.start_scan()
    ...
    laser.stop_scan()
    laser.set_output(False)
    laser.close()
"""

import glob
import os
import sys

PRODUCT_ID = 4106  # ProductID USB standard dei controller TLB-6700

CARTELLE_DA_CERCARE = [
    r"C:\Program Files\New Focus",
    r"C:\Program Files (x86)\New Focus",
    r"C:\Program Files\Newport",
    r"C:\Program Files (x86)\Newport",
]


def _trova_dll():
    for cartella in CARTELLE_DA_CERCARE:
        if not os.path.isdir(cartella):
            continue
        risultati = glob.glob(
            os.path.join(cartella, "**", "UsbDllWrap.dll"), recursive=True
        )
        if risultati:
            return risultati[0]
    return None


class TLB6700:
    def __init__(self, device_key, product_id=PRODUCT_ID):
        """
        device_key: stringa come "6700 SN1234", specifica del tuo laser.
                    Scoprila eseguendo trova_tlb6700.py.
        """
        self.device_key = device_key
        self.product_id = product_id
        self._tlb = None

    def open(self):
        dll_path = _trova_dll()
        if dll_path is None:
            raise RuntimeError(
                "UsbDllWrap.dll non trovata: verifica che il software Newport "
                "sia installato, o aggiungi il percorso a CARTELLE_DA_CERCARE."
            )
        sys.path.append(os.path.dirname(dll_path))

        import clr  # richiede pythonnet

        clr.AddReference("UsbDllWrap")
        import Newport
        from System.Text import StringBuilder

        self._tlb = Newport.USBComm.USB()
        self._tlb.OpenDevices(self.product_id, True)
        self._answer = StringBuilder(64)

    def close(self):
        if self._tlb is not None:
            self._tlb.CloseDevices()
            self._tlb = None

    def query(self, msg):
        self._answer.Clear()
        self._tlb.Query(self.device_key, msg, self._answer)
        return self._answer.ToString()

    def idn(self):
        return self.query("*IDN?")

    def reset(self):
        return self.query("*RST")

    def set_output(self, on):
        return self.query(f"OUTPut:STATe {1 if on else 0}")

    def get_output(self):
        return self.query("OUTPut:STATe?")

    def set_wavelength(self, nm):
        self.query(f"SOURce:WAVElength {nm}")
        self.query("OUTPut:TRACK 1")
        return self.query("SOURce:WAVElength?")

    def get_wavelength(self):
        return self.query("SOURce:WAVElength?")

    def set_scan_limits(self, start_nm, stop_nm):
        self.query(f"SOURce:WAVElength:START {start_nm}")
        self.query(f"SOURce:WAVElength:STOP {stop_nm}")
        return self.query("SOURce:WAVElength:START?"), self.query(
            "SOURce:WAVElength:STOP?"
        )

    def set_scan_speeds(self, forward, backward):
        """forward/backward in nm/s."""
        self.query(f"SOURce:WAVE:SLEW:FORWard {forward}")
        self.query(f"SOURce:WAVE:SLEW:RETurn {backward}")
        return self.query("SOURce:WAVE:SLEW:FORWard?"), self.query(
            "SOURce:WAVE:SLEW:RETurn?"
        )

    def start_scan(self, num_scans=1):
        self.query(f"SOUR:WAVE:DESSCANS {num_scans}")
        return self.query("OUTPut:SCAN:START")

    def stop_scan(self):
        return self.query("OUTPut:SCAN:STOP")

    def operation_complete(self):
        return self.query("*OPC?")
