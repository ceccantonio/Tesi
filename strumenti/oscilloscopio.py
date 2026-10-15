"""Oscilloscopio Tektronix MDO3104 via VISA. Uso: osc = Oscilloscopio(); osc.arma(); ...; dati = osc.leggi()"""
import time

import numpy as np
import pyvisa

VISA = "USB::0x0699::0x0408::C028813::INSTR"


class Oscilloscopio:
    def __init__(self, visa=VISA, timeout_ms=30000):
        """Apre lo scope e lo prepara a leggere tutto il record in ASCII."""
        self.osc = pyvisa.ResourceManager().open_resource(visa)
        self.osc.timeout = timeout_ms
        self.osc.write("HORizontal:MODE MANUAL")
        self.osc.write("DATA:START 1")
        self.osc.write(f"DATA:STOP {int(self.osc.query('HORizontal:RECORDLENGTH?'))}")
        self.osc.write("DATA:WIDTH 1")
        self.osc.write("DATA:ENCdg ASCii")

    def _fino_a(self, cond, timeout):
        t0 = time.time()
        while not cond():
            if time.time() - t0 > timeout:
                return False
            time.sleep(0.1)
        return True

    def arma(self, timeout=120):
        """Single sequence; ritorna True quando e' READY (solo allora il trigger viene accettato)."""
        self.osc.write("ACQ:STOPAFTER SEQUENCE")
        self.osc.write("ACQ:STATE RUN")
        return self._fino_a(lambda: self.osc.query("TRIGger:STATE?").strip().upper() == "READY", timeout)

    def attendi(self, timeout=600):
        """Aspetta la fine dell'acquisizione; True se finita nel tempo."""
        return self._fino_a(lambda: int(self.osc.query("ACQ:STATE?")) == 0, timeout)

    def forza_trigger(self):
        self.osc.write("TRIGger FORCe")

    def leggi(self, canali=("CH1", "CH2")):
        """Record completo dei canali in volt, array (N, n_canali)."""
        col = []
        for ch in canali:
            self.osc.write(f"DATA:SOURCE {ch}")
            ymult, yoff, yzero = (float(self.osc.query(f"WFMPRE:{k}?")) for k in ("YMULT", "YOFF", "YZERO"))
            col.append((np.array(self.osc.query("CURVE?").split(","), dtype=float) - yoff) * ymult + yzero)
            time.sleep(0.1)
        n = min(map(len, col))
        return np.column_stack([c[:n] for c in col])
