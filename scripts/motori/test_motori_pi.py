"""
Test dei due controller PI C-863.11 (Mercury) via porta seriale.

Per default NON muove niente: apre ogni porta e legge identificativo, errori,
stato servo, riferimento, limiti, velocità e posizione, così puoi controllare
che la comunicazione funzioni e capire quale porta è quale motore.

Per provare un movimento piccolo imposta MOVE_MM (relativo, in mm) a un valore
diverso da 0: lo script controlla i limiti, chiede conferma scrivendo "si", e
aspetta che il motore arrivi a destinazione.

Prima di lanciarlo CHIUDI PI MicroMove: la porta seriale può essere usata da
un solo programma alla volta.

Comandi PI GCS usati: *IDN?, ERR?, SVO?, FRF?, TMN?, TMX?, VEL?, POS?, ONT?, MVR.
"""

import time

import pyvisa

# ---------------------------------------------------------------- CONFIGURAZIONE
PORTE = ["ASRL4::INSTR", "ASRL5::INSTR"]
BAUD_RATE = 19200
TERMINATORE = "\n"
TIMEOUT_MS = 2000

AXIS = "1"  # ID dell'asse sul controller C-863 (quasi sempre 1)

MOVE_MM = -3.0  # spostamento relativo di prova in mm; 0 = nessun movimento
PORTA_DA_MUOVERE = "ASRL5::INSTR"  # usata solo se MOVE_MM != 0
TIMEOUT_MOVIMENTO_S = 30
# ---------------------------------------------------------------------------


def apri(rm, porta):
    inst = rm.open_resource(porta)
    inst.baud_rate = BAUD_RATE
    inst.timeout = TIMEOUT_MS
    inst.read_termination = TERMINATORE
    inst.write_termination = TERMINATORE
    return inst


def chiedi(inst, comando):
    try:
        return inst.query(comando).strip()
    except Exception as e:
        return f"ERRORE: {e}"


def leggi_numero(risposta):
    """Estrae il valore da risposte tipo '1=+0012.3456'."""
    try:
        return float(risposta.split("=")[-1])
    except ValueError:
        return None


def stampa_stato(porta, inst):
    print(f"\n=== {porta} ===")
    for etichetta, comando in [
        ("IDN", "*IDN?"),
        ("Errore", "ERR?"),
        ("Servo (1=on)", f"SVO? {AXIS}"),
        ("Referenziato (1=si)", f"FRF? {AXIS}"),
        ("Limite min", f"TMN? {AXIS}"),
        ("Limite max", f"TMX? {AXIS}"),
        ("Velocita", f"VEL? {AXIS}"),
        ("Posizione", f"POS? {AXIS}"),
        ("On target (1=si)", f"ONT? {AXIS}"),
    ]:
        print(f"  {etichetta:<20} {comando:<10} -> {chiedi(inst, comando)}")


def muovi_relativo(inst, mm):
    pos = leggi_numero(chiedi(inst, f"POS? {AXIS}"))
    tmin = leggi_numero(chiedi(inst, f"TMN? {AXIS}"))
    tmax = leggi_numero(chiedi(inst, f"TMX? {AXIS}"))
    if None in (pos, tmin, tmax):
        print("Non riesco a leggere posizione/limiti: nessun movimento.")
        return

    target = pos + mm
    print(f"\nPosizione {pos} mm -> target {target} mm (limiti {tmin} .. {tmax} mm)")
    if not (tmin <= target <= tmax):
        print("Il target è fuori dai limiti: nessun movimento.")
        return

    if input("Muovo il motore? scrivi 'si' per confermare: ").strip().lower() != "si":
        print("Annullato.")
        return

    inst.write(f"MVR {AXIS} {mm}")
    t0 = time.time()
    while time.time() - t0 < TIMEOUT_MOVIMENTO_S:
        if chiedi(inst, f"ONT? {AXIS}").endswith("1"):
            break
        time.sleep(0.2)
    else:
        print("Timeout: il motore non ha segnalato 'on target'.")

    print(f"Posizione finale: {chiedi(inst, f'POS? {AXIS}')} | errore: {chiedi(inst, 'ERR?')}")


def main():
    rm = pyvisa.ResourceManager()
    aperte = {}
    try:
        for porta in PORTE:
            aperte[porta] = apri(rm, porta)
            stampa_stato(porta, aperte[porta])

        if MOVE_MM != 0:
            muovi_relativo(aperte[PORTA_DA_MUOVERE], MOVE_MM)
        else:
            print("\nMOVE_MM = 0: nessun movimento eseguito (solo lettura).")
    finally:
        for inst in aperte.values():
            inst.close()


if __name__ == "__main__":
    main()
