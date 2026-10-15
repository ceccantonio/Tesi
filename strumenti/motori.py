"""Motori PI C-863 via seriale (VISA). Uso: mot = Motori(); mot.muovi(x_mm, y_mm); mot.x.pos"""
import time

import pyvisa


class Asse:
    def __init__(self, porta, baud=19200, tolleranza_mm=0.01):
        self.i = pyvisa.ResourceManager().open_resource(porta)
        self.i.baud_rate, self.i.timeout = baud, 2000
        self.i.read_termination = self.i.write_termination = "\n"
        self.tol = tolleranza_mm
        self.q("ERR?")                      # legge e azzera un eventuale errore vecchio

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
        """Movimento assoluto in mm. Non muove se il target e' fuori dai limiti; ritorna True se arrivato."""
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
    def __init__(self, x="ASRL4::INSTR", y="ASRL5::INSTR"):
        self.x, self.y = Asse(x), Asse(y)

    def muovi(self, x, y):
        return self.x.muovi(x) and self.y.muovi(y)
