"""
Test isolato del laser TLB-6700, senza oscilloscopio: avvia uno scan e
stampa la wavelength corrente ogni secondo, per verificare "a occhio" se
il laser sta davvero scansionando prima di rimetterlo nel flusso completo
con l'oscilloscopio.

La DeviceKey viene trovata automaticamente da TLB6700.open() (deve
esserci un solo laser TLB-6700 collegato).
"""

import os
import sys
import time

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from tlb6700 import TLB6700

# ---------------------------------------------------------------- CONFIGURAZIONE
WAVELENGTH_START = 1550.0  # nm (range piccolo, solo per test rapido isolato)
WAVELENGTH_STOP = 1560.0   # nm
SCAN_SPEED_FORWARD = 0.5   # nm/s (più rapido del default, solo per il test)
SCAN_SPEED_BACKWARD = 5.0  # nm/s

DURATA_TEST_S = 20
INTERVALLO_LETTURA_S = 1.0
# ---------------------------------------------------------------------------

laser = TLB6700()
laser.open()
try:
    print("IDN:", laser.idn())

    laser.setup_sweep(
        start_nm=WAVELENGTH_START,
        stop_nm=WAVELENGTH_STOP,
        speed=SCAN_SPEED_FORWARD,
        return_speed=SCAN_SPEED_BACKWARD,
    )
    print(f"Scan configurato: {WAVELENGTH_START}-{WAVELENGTH_STOP} nm, "
          f"forward={SCAN_SPEED_FORWARD} nm/s, backward={SCAN_SPEED_BACKWARD} nm/s")

    laser.on()
    print("Laser acceso (output ON).")

    print("\nAvvio scan...")
    laser.start_sweep()

    print(f"\nMonitoro la wavelength per {DURATA_TEST_S} s (ogni {INTERVALLO_LETTURA_S} s):")
    t0 = time.time()
    while time.time() - t0 < DURATA_TEST_S:
        wl = laser.get_wavelength()
        print(f"  t={time.time() - t0:5.1f}s  wavelength = {wl} nm")
        time.sleep(INTERVALLO_LETTURA_S)

    print("\nFermo lo scan...")
    laser.stop_sweep()
    laser.off()
finally:
    # Chiude sempre la connessione, anche se qualcosa sopra fallisce:
    # altrimenti il dispositivo resta "agganciato" per il prossimo tentativo.
    laser.close()
    print("Connessione al laser chiusa.")
