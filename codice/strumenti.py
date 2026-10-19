"""Macchine del laboratorio: laser Newport TLB-6700, oscilloscopio Tektronix MDO3104, motori PI C-863.

    laser, osc, mot = init()                       # collega tutto (init(laser=False, ...) per escludere qualcosa)
    cartella = scansiona_punti(punti, laser, osc, mot)
"""
import datetime
import sys
import time
from pathlib import Path

import numpy as np
import pyvisa

RADICE = Path(__file__).resolve().parent.parent
DATI = RADICE / "Chiara_Bru" / "runs_oscilloscope"          # qui finiscono AAAA-MM-GG/run_NNN
DLL = RADICE / "scripts" / "laser" / "UsbDllWrap.dll"

OSC_VISA = "USB::0x0699::0x0408::C028813::INSTR"
MOTORE_X, MOTORE_Y = "ASRL4::INSTR", "ASRL5::INSTR"         # X = COM4, Y = COM5
LIMITI_MM = (0.0, 150.0)                                    # corsa dei motori
LAMBDA_START, LAMBDA_STOP = 1520.0, 1570.0                  # nm, come trigger_scan.m
ATTESA_RITORNO_S = 5.0                                      # dopo una scansione, perche' il laser torni a 1520 nm


class Laser:
    def __init__(self, dll=DLL, addr=1, device_id=0x100A):
        """Carica la DLL e apre il dispositivo (una volta per sessione; niente CloseDevices, come in MATLAB)."""
        sys.path.append(str(Path(dll).parent))
        import clr
        clr.AddReference("UsbDllWrap")
        import Newport
        from System.Text import StringBuilder
        self.usb, self.addr, self.buf = Newport.USBComm.USB(), addr, StringBuilder(64)
        self.usb.OpenDevices(device_id)

    def q(self, comando):
        """Comando + lettura della risposta (anche per le scritture: altrimenti gli 'OK' si accumulano)."""
        self.usb.Query(self.addr, comando, self.buf)
        return self.buf.ToString()

    def imposta(self, start=LAMBDA_START, stop=LAMBDA_STOP):
        self.q(f"SOUR:WAVE:START {start}")
        self.q(f"SOUR:WAVE:STOP {stop}")
        self.ferma()
        self.q("SOUR:WAVE:DESSCANS 1")                      # azzera il contatore delle scansioni

    def avvia(self):
        self.q("OUTP:SCAN:START")

    def ferma(self):
        self.q("OUTP:SCAN:STOP")


class Oscilloscopio:
    def __init__(self, visa=OSC_VISA, n_campioni=None):
        """Apre lo scope e lo prepara a leggere il record (n_campioni=None: tutto) in ASCII."""
        self.osc = pyvisa.ResourceManager().open_resource(visa)
        self.osc.timeout = 30000
        self.osc.write("HORizontal:MODE MANUAL")
        record = int(self.osc.query("HORizontal:RECORDLENGTH?"))
        self.osc.write("DATA:START 1")
        self.osc.write(f"DATA:STOP {min(record, n_campioni or record)}")
        self.osc.write("DATA:WIDTH 1")
        self.osc.write("DATA:ENCdg ASCii")

    def _attendi(self, condizione, timeout):
        t0 = time.time()
        while not condizione() and time.time() - t0 < timeout:
            time.sleep(0.1)

    def arma(self, timeout=120):
        """Single sequence; torna quando lo scope e' READY (solo allora il trigger viene accettato)."""
        self.osc.write("ACQ:STOPAFTER SEQUENCE")
        self.osc.write("ACQ:STATE RUN")
        self._attendi(lambda: self.osc.query("TRIGger:STATE?").strip().upper() == "READY", timeout)

    def attendi(self, timeout=600):
        """Torna a acquisizione finita."""
        self._attendi(lambda: int(self.osc.query("ACQ:STATE?")) == 0, timeout)

    def forza_trigger(self):
        self.osc.write("TRIGger FORCe")

    def leggi(self, canali=("CH1", "CH2")):
        """Record dei canali in volt, array (N, n_canali)."""
        col = []
        for ch in canali:
            self.osc.write(f"DATA:SOURCE {ch}")
            ymult, yoff, yzero = (float(self.osc.query(f"WFMPRE:{k}?")) for k in ("YMULT", "YOFF", "YZERO"))
            col.append((np.array(self.osc.query("CURVE?").split(","), dtype=float) - yoff) * ymult + yzero)
            time.sleep(0.1)
        n = min(map(len, col))
        return np.column_stack([c[:n] for c in col])


class Asse:
    def __init__(self, porta, baud=19200, tolleranza_mm=0.01):
        self.i = pyvisa.ResourceManager().open_resource(porta)
        self.i.baud_rate, self.i.timeout = baud, 2000
        self.i.read_termination = self.i.write_termination = "\n"
        self.tol = tolleranza_mm
        self.q("ERR?")                                      # legge e azzera un eventuale errore vecchio

    def q(self, comando):
        return self.i.query(comando).strip()

    def v(self, comando):
        """Numero da risposte tipo '1=+0012.3456'."""
        return float(self.q(comando).split("=")[-1])

    @property
    def pos(self):
        return self.v("POS? 1")

    @property
    def limiti(self):
        return self.v("TMN? 1"), self.v("TMX? 1")

    def muovi(self, target, timeout=60):
        """Movimento assoluto in mm. Non muove se il target e' fuori dai limiti; True se arrivato."""
        tmin, tmax = self.limiti
        if not tmin <= target <= tmax:
            return False
        self.i.write(f"MOV 1 {target:.4f}")
        time.sleep(0.1)
        t0 = time.time()
        while not (self.q("ONT? 1").endswith("1") and abs(self.pos - target) <= self.tol):
            if time.time() - t0 > timeout:
                return False
            time.sleep(0.1)
        return int(self.q("ERR?")) == 0


class Motori:
    def __init__(self, x=MOTORE_X, y=MOTORE_Y):
        self.x, self.y = Asse(x), Asse(y)

    def muovi(self, x, y):
        return self.x.muovi(x) and self.y.muovi(y)


def init(laser=True, osc=True, motori=True, n_campioni=None):
    """Collega le macchine e le restituisce come (laser, osc, motori); quelle escluse sono None."""
    return (Laser() if laser else None,
            Oscilloscopio(n_campioni=n_campioni) if osc else None,
            Motori() if motori else None)


def scansione(laser, osc):
    """Una scansione completa 1520-1570 nm sincronizzata con lo scope. Ritorna (N, 2): CH1, CH2."""
    laser.imposta()
    osc.arma()
    laser.avvia()
    osc.attendi()
    time.sleep(1.0)
    dati = osc.leggi()
    time.sleep(ATTESA_RITORNO_S)
    return dati


def acquisizione_ferma(osc):
    """Acquisizione senza scansione (laser fermo, trigger forzato). Ritorna (N, 2): CH1, CH2."""
    osc.arma()
    osc.forza_trigger()
    osc.attendi()
    time.sleep(1.0)
    return osc.leggi()


def nuova_cartella(prefisso="run"):
    """DATI/AAAA-MM-GG/<prefisso>_NNN, creata."""
    base = DATI / datetime.date.today().isoformat()
    base.mkdir(parents=True, exist_ok=True)
    cartella = base / f"{prefisso}_{len(list(base.glob(prefisso + '_*'))) + 1:03d}"
    cartella.mkdir()
    return cartella


def scansiona_punti(punti, laser, osc, mot, assestamento=1.0):
    """Per ogni (x, y): muove i motori, fa una scansione e salva un CSV. Registro: griglia.csv. Ritorna la cartella."""
    cartella = nuova_cartella("run")
    registro = cartella / "griglia.csv"
    registro.write_text("indice,x_mm,y_mm,file\n")
    try:
        for i, (x, y) in enumerate(punti):
            print(f"[{i + 1}/{len(punti)}] X={x:.4f}  Y={y:.4f}")
            if not mot.muovi(x, y):
                print("   motori non arrivati o fuori limiti: punto saltato")
                continue
            time.sleep(assestamento)
            nome = f"punto_{i:03d}_x{x:.3f}_y{y:.3f}.csv"
            np.savetxt(cartella / nome, scansione(laser, osc), delimiter=",", header="CH1,CH2", comments="")
            with registro.open("a") as f:
                f.write(f"{i},{x:.4f},{y:.4f},{nome}\n")
    except KeyboardInterrupt:
        laser.ferma()
        print("Interrotto: scansione del laser fermata")
    return cartella
