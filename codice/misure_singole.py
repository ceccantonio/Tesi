"""Misure di singoli punti, da confrontare dopo con analisi.py (python analisi.py confronto 1 2 3).

    python misure_singole.py

Per ogni misura chiede un nome, delle note (scrivi cosa hai cambiato nel setup) e, se USA_MOTORI, la posizione.
Invio a vuoto sul nome = esci. Ogni misura finisce in misure_singole/ con un numero progressivo e una riga in misure.csv.
"""
import csv
import datetime

import numpy as np

from strumenti import DATI, acquisizione_ferma, init, scansione

CON_SCANSIONE = False   # True: scansione del laser (spettro); False: laser fermo, trigger forzato (livello)
USA_MOTORI = False      # True: chiede la posizione (X Y) e muove i motori prima di misurare
N_CAMPIONI = 100000     # senza scansione: campioni letti per misura (None = tutto il record, ~240 MB)

CARTELLA = DATI / "misure_singole"
CAMPI = ["id", "ora", "nome", "x_mm", "y_mm", "con_scansione", "file", "note"]


def main():
    laser, osc, mot = init(laser=CON_SCANSIONE, motori=USA_MOTORI, n_campioni=None if CON_SCANSIONE else N_CAMPIONI)
    CARTELLA.mkdir(exist_ok=True)
    registro = CARTELLA / "misure.csv"
    if not registro.exists():
        registro.write_text(",".join(CAMPI) + "\n")
    while True:
        nome = input("\nNome della misura (invio per uscire): ").strip()
        if not nome:
            break
        note = input("Note (cosa hai cambiato): ").strip()
        x = y = ""
        if USA_MOTORI:
            pos = input("Posizione 'X Y' in mm (invio = non muovere): ").split()
            if len(pos) == 2:
                x, y = pos
                if not mot.muovi(float(x), float(y)):
                    print("Motori non arrivati o fuori limiti: misura annullata")
                    continue
        dati = scansione(laser, osc) if CON_SCANSIONE else acquisizione_ferma(osc)
        with registro.open(newline="") as f:
            n = max([int(r["id"]) for r in csv.DictReader(f)] + [0]) + 1
        file = f"m{n:03d}_{''.join(c if c.isalnum() or c in '-_' else '_' for c in nome)}.csv"
        np.savetxt(CARTELLA / file, dati, delimiter=",", header="CH1,CH2", comments="")
        with registro.open("a", newline="") as f:
            csv.writer(f).writerow([n, datetime.datetime.now().isoformat(timespec="seconds"), nome, x, y, CON_SCANSIONE, file, note])
        print(f"Misura {n} salvata ({file}), mediana di CH2 = {np.median(dati[:, 1]):.5g} V")


if __name__ == "__main__":
    main()
