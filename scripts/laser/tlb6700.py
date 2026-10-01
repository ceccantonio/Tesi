"""
Driver per il laser sintonizzabile Newport/New Focus TLB-6700.

Basato sul codice MATLAB usato nel tuo laboratorio (NP_USB_connect_velocity.m,
TLB6700_init.m, trigger_scan.m - cartella laser_fix/), non su esempi generici
trovati online: usa la stessa DLL .NET proprietaria (UsbDllWrap.dll) del
software Newport, con lo stesso schema di comunicazione provato funzionante
su questo laser specifico.

Punto importante (diverso dai tentativi precedenti): il dispositivo si
indirizza con un semplice **indirizzo USB numerico** (USBADDR, 1 nel codice
MATLAB del laboratorio), non con una DeviceKey testuale tipo "6700 SN1234".
Niente comandi SYSTem:MCONtrol: nel codice del laboratorio non ci sono, e
quindi non li mandiamo nemmeno qui.

Richiede:
    pip install pythonnet

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

PRODUCT_ID = 0x100A  # = 4106 decimale, come in deviceID = hex2dec('100A')
USB_ADDR = 1  # indirizzo USB del laser sul bus (0-31), come USBADDR nel .m

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

# Interpretazione dei codici di errore del firmware USB Newport, presa da
# NP_USB_reperror.m
_MESSAGGI_ERRORE = {
    0: "operazione eseguita correttamente",
    1: "operazione eseguita correttamente",
    -2: "USBADDRESSNOTFOUND: device ID non trovato tra i dispositivi aperti sul bus",
    -3: "USBINVALIDADDRESS: device ID fuori dal range valido 0-31",
}


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


def _descrivi_errore(codice):
    return _MESSAGGI_ERRORE.get(codice, f"codice di errore sconosciuto ({codice})")


class TLB6700:
    def __init__(self, usb_addr=USB_ADDR, product_id=PRODUCT_ID):
        self.usb_addr = usb_addr
        self.product_id = product_id
        self._tlb = None
        self._answer = None

    def open(self):
        """Trova la DLL e apre il dispositivo (indirizzo USB fisso, niente DeviceKey)."""
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
        self._answer = StringBuilder(64)

        # Overload a UN argomento (come nel codice del laboratorio): restituisce
        # un codice di errore intero, non un bool.
        esito = self._tlb.OpenDevices(self.product_id)
        print(f"OpenDevices: {_descrivi_errore(esito)} (codice {esito})")
        if esito not in (0, 1):
            raise RuntimeError(f"OpenDevices fallito: {_descrivi_errore(esito)}")

    def close(self):
        if self._tlb is not None:
            self._tlb.CloseDevices()
            self._tlb = None

    # --- I/O di base -----------------------------------------------------
    def query(self, text):
        self._answer.Clear()
        esito = self._tlb.Query(self.usb_addr, text, self._answer)
        if esito not in (0, 1):
            raise RuntimeError(f"Query {text!r} fallita: {_descrivi_errore(esito)}")
        return self._answer.ToString()

    def write(self, text):
        """Manda un comando (non una query)."""
        if text.strip().endswith("?"):
            raise ValueError(f"write() non va usato per query: {text!r}")
        esito = self._tlb.Write(self.usb_addr, text)
        if esito not in (0, 1):
            raise RuntimeError(f"Comando {text!r} fallito: {_descrivi_errore(esito)}")

    def idn(self):
        return self.query("*IDN?")

    def reset(self):
        return self.query("*RST")

    def operation_complete(self):
        return self.query("*OPC?")

    # --- Output on/off -----------------------------------------------------
    def set_output(self, on):
        self.write(f"OUTP:STAT {1 if on else 0}")

    def get_output(self):
        return self.query("OUTP:STAT?")

    def on(self):
        self.set_output(True)
        print("Laser acceso.")

    def off(self):
        self.set_output(False)
        print("Laser spento.")

    # --- Wavelength / tracking --------------------------------------------
    def set_track(self, on):
        """Necessario perché il laser segua davvero i comandi di wavelength/scan."""
        self.write(f"OUTP:TRAC {1 if on else 0}")

    def set_wavelength(self, nm):
        self.set_track(True)
        self.write(f"SOUR:WAVE {nm}")
        return self.get_wavelength()

    def get_wavelength(self):
        return self.query("SENS:WAVE?")

    # --- Scan / sweep --------------------------------------------------------
    def set_scan_limits(self, start_nm, stop_nm):
        for nm in (start_nm, stop_nm):
            if not (WAVELENGTH_MIN_NM <= nm <= WAVELENGTH_MAX_NM):
                raise ValueError(
                    f"Wavelength {nm} nm fuori dal range tipico "
                    f"[{WAVELENGTH_MIN_NM}, {WAVELENGTH_MAX_NM}] nm."
                )
        self.write(f"SOUR:WAVE:START {start_nm}")
        self.write(f"SOUR:WAVE:STOP {stop_nm}")

    def set_scan_speeds(self, forward, backward):
        """forward/backward in nm/s."""
        self.write(f"SOUR:WAVE:SLEW:FORW {forward}")
        self.write(f"SOUR:WAVE:SLEW:RET {backward}")

    def setup_sweep(self, start_nm, stop_nm, speed, return_speed=None, scans=1):
        """Configura uno sweep completo: limiti, velocità, numero di scan."""
        self.stop_sweep()
        self.write(f"SOUR:WAVE:DESSCANS {scans}")
        self.set_scan_limits(start_nm, stop_nm)
        self.set_scan_speeds(speed, return_speed if return_speed is not None else speed)

    def start_sweep(self):
        self.write("OUTP:SCAN:START")

    def stop_sweep(self):
        self.write("OUTP:SCAN:STOP")

    # Alias per compatibilità con gli script che già usano questi nomi
    def start_scan(self, num_scans=1):
        self.write(f"SOUR:WAVE:DESSCANS {num_scans}")
        self.start_sweep()

    def stop_scan(self):
        self.stop_sweep()
