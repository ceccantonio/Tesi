"""
Test isolato del laser TLB-6700, senza oscilloscopio: avvia uno scan e
stampa la wavelength corrente ogni secondo, per verificare "a occhio" se
il laser sta davvero scansionando prima di rimetterlo nel flusso completo
con l'oscilloscopio.

Riempi DEVICE_KEY con il valore trovato tramite trova_tlb6700.py.
"""

import time

from tlb6700 import TLB6700

# ---------------------------------------------------------------- CONFIGURAZIONE
DEVICE_KEY = "6700 SN0000"  # <-- sostituisci con la tua DeviceKey

WAVELENGTH_START = 1550.0  # nm
WAVELENGTH_STOP = 1560.0   # nm
SCAN_SPEED_FORWARD = 0.5   # nm/s (più rapido del default, solo per il test)
SCAN_SPEED_BACKWARD = 5.0  # nm/s

DURATA_TEST_S = 20
INTERVALLO_LETTURA_S = 1.0
# ---------------------------------------------------------------------------

laser = TLB6700(DEVICE_KEY)
laser.open()
print("IDN:", laser.idn())

limiti = laser.set_scan_limits(WAVELENGTH_START, WAVELENGTH_STOP)
print("Limiti scan impostati (start, stop):", limiti)

velocita = laser.set_scan_speeds(SCAN_SPEED_FORWARD, SCAN_SPEED_BACKWARD)
print("Velocità scan impostate (forward, backward):", velocita)

print("Output ON:", laser.set_output(True))

print("\nAvvio scan...")
print("Risposta OUTPut:SCAN:START ->", laser.start_scan())

print(f"\nMonitoro la wavelength per {DURATA_TEST_S} s (ogni {INTERVALLO_LETTURA_S} s):")
t0 = time.time()
while time.time() - t0 < DURATA_TEST_S:
    wl = laser.get_wavelength()
    print(f"  t={time.time() - t0:5.1f}s  wavelength = {wl}")
    time.sleep(INTERVALLO_LETTURA_S)

print("\nFermo lo scan...")
print("Risposta OUTPut:SCAN:STOP ->", laser.stop_scan())
laser.set_output(False)
laser.close()
print("Fatto.")
