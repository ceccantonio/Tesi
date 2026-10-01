"""
Driver per il laser sintonizzabile Newport/New Focus TLB-6700, basato su
un driver pubblico testato (Rhys Povey, 2021) che usa la stessa DLL .NET
proprietaria (UsbDllWrap.dll) installata insieme al software Newport
(il TLB-6700 non usa VISA/porte seriali: driver USB Jungo + wrapper .NET).

Richiede:
    pip install pythonnet

Differenze rispetto alla prima versione (basata solo sul notebook
marvync/Laser-automation-NewFocus):
- la DeviceKey viene trovata automaticamente con GetDeviceKeys(), non va
  più scritta a mano (deve però esserci un solo laser collegato)
- il laser viene messo esplicitamente in modalità "remote control"
  (SYSTem:MCONtrol REMote) alla connessione: senza questo, i comandi di
  scan mandati da qui potrebbero essere accettati ma non avere alcun
  effetto reale se il laser resta in modalità locale/manuale
- i comandi di scrittura (write) controllano che la risposta sia "OK" e
  sollevano un errore se il laser rifiuta il comando, invece di
  ignorarlo in silenzio

Uso tipico:

    laser = TLB6700()
    laser.open()
    laser.setup_sweep(start_nm=1550.0, stop_nm=1560.0, speed=0.1)
    laser.on()
    laser.start_sweep()
    ...
    laser.stop_sweep()
    laser.off()
    laser.close()
"""

import glob
import os
import sys
import time

PRODUCT_ID = 4106  # ProductID USB dei controller TLB-6700 (decimale di 0x100A)

CARTELLE_DA_CERCARE = [
    r"C:\Program Files\New Focus",
    r"C:\Program Files (x86)\New Focus",
    r"C:\Program Files\Newport",
    r"C:\Program Files (x86)\Newport",
]

# Range di sintonizzazione tipico dei TLB-6700 in banda C (verifica contro
# l'etichetta/manuale del tuo esemplare se diverso).
WAVELENGTH_MIN_NM = 1520
WAVELENGTH_MAX_NM = 1570


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
    def __init__(self, product_id=PRODUCT_ID):
        self.product_id = product_id
        self._tlb = None
        self._answer = None
        self.device_key = None

    def open(self):
        """Trova la DLL, apre il dispositivo e scopre la sua DeviceKey."""
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

        ndevices, keystrings = self._tlb.GetDeviceKeys("")
        if ndevices != 1:
            raise RuntimeError(
                f"Trovati {ndevices} dispositivi laser (ne serve esattamente 1). "
                f"Chiavi viste: {list(keystrings) if keystrings else keystrings}"
            )
        self.device_key = keystrings[0]
        print(f"{self.device_key} connesso.")

        self.reset()
        self.write("syst:mcon rem")  # modalità controllo remoto

    def close(self):
        if self._tlb is not None:
            try:
                self.write("syst:mcon loc")  # torna in modalità locale
            finally:
                self._tlb.CloseDevices()
                self._tlb = None

    # --- I/O di base -----------------------------------------------------
    def query(self, text):
        self._answer.Clear()
        self._tlb.Query(self.device_key, text, self._answer)
        return self._answer.ToString()

    def write(self, text):
        """Manda un comando (non una query) e controlla che risponda 'OK'."""
        if text.strip().endswith("?"):
            raise ValueError(f"write() non va usato per query: {text!r}")
        risposta = self.query(text)
        if risposta != "OK":
            raise RuntimeError(f"Comando {text!r} fallito, risposta: {risposta!r}")

    def idn(self):
        return self.query("*IDN?")

    def reset(self):
        return self.query("*RST")

    def operation_complete(self):
        return self.query("*OPC?")

    # --- Output on/off -----------------------------------------------------
    def set_output(self, on):
        if self.query("outp:stat?") != ("1" if on else "0"):
            self.write(f"outp:stat {1 if on else 0}")

    def get_output(self):
        return self.query("outp:stat?")

    def on(self, wait=True):
        self.set_output(True)
        if wait:
            time.sleep(float(self.query("ondelay?")) / 1000)

    def off(self):
        self.set_output(False)

    # --- Wavelength / tracking --------------------------------------------
    def set_track(self, on):
        """Necessario perché il laser segua davvero i comandi di wavelength/scan."""
        self.write(f"outp:trac {1 if on else 0}")

    def set_wavelength(self, nm):
        self.set_track(True)
        if float(self.query("sour:wave?")) != nm:
            self.write(f"sour:wave {nm}")
        return self.get_wavelength()

    def get_wavelength(self):
        return self.query("sens:wave?")

    # --- Scan / sweep --------------------------------------------------------
    def set_scan_limits(self, start_nm, stop_nm):
        for nm in (start_nm, stop_nm):
            if not (WAVELENGTH_MIN_NM <= nm <= WAVELENGTH_MAX_NM):
                raise ValueError(
                    f"Wavelength {nm} nm fuori dal range tipico "
                    f"[{WAVELENGTH_MIN_NM}, {WAVELENGTH_MAX_NM}] nm."
                )
        self.write(f"sour:wave:start {start_nm}")
        self.write(f"sour:wave:stop {stop_nm}")

    def set_scan_speeds(self, forward, backward):
        """forward/backward in nm/s."""
        self.write(f"sour:wave:slew:forw {forward}")
        self.write(f"sour:wave:slew:ret {backward}")

    def setup_sweep(self, start_nm, stop_nm, speed, return_speed=None, scans=1):
        """Configura uno sweep completo: limiti, velocità, numero di scan."""
        self.stop_sweep()
        self.write("sour:wave:scancfg 0")  # laser acceso anche nel ritorno
        self.write(f"sour:wave:desscans {scans}")
        self.set_scan_limits(start_nm, stop_nm)
        self.set_scan_speeds(speed, return_speed if return_speed is not None else speed)

    def start_sweep(self):
        self.write("outp:scan:start")

    def stop_sweep(self):
        self.write("outp:scan:stop")

    # Alias per compatibilità con gli script che già usano questi nomi
    def start_scan(self, num_scans=1):
        self.write(f"sour:wave:desscans {num_scans}")
        self.start_sweep()

    def stop_scan(self):
        self.stop_sweep()
