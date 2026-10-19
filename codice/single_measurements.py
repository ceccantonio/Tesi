"""Single-point measurements, to be compared afterwards with analysis.py (python analysis.py compare 1 2 3).

    python single_measurements.py

For each measurement it asks for a name, some notes (write what you changed in the setup) and, if USE_MOTORS, the position.
An empty name quits. Each measurement goes to misure_singole/ with a progressive number and a row in misure.csv
(same folder and columns as misure_singole.ipynb).
"""
import csv
import datetime

import numpy as np

from acquisition import DATA, scan, still_acquisition
from initialization import init

SCAN = False         # True: laser scan (spectrum); False: laser still, forced trigger (a single level)
USE_MOTORS = False   # True: asks for the position (X Y) and moves the motors before measuring
N_SAMPLES = 100000   # without scan: samples read per measurement (None = the whole record, ~240 MB)

FOLDER = DATA / "misure_singole"
FIELDS = ["id", "ora", "nome", "x_mm", "y_mm", "con_scansione", "file", "note"]   # kept as in the notebook


def main():
    laser, scope, motors = init(laser=SCAN, motors=USE_MOTORS, n_samples=None if SCAN else N_SAMPLES)
    FOLDER.mkdir(exist_ok=True)
    registry = FOLDER / "misure.csv"
    if not registry.exists():
        registry.write_text(",".join(FIELDS) + "\n")
    while True:
        name = input("\nMeasurement name (empty to quit): ").strip()
        if not name:
            break
        notes = input("Notes (what you changed): ").strip()
        x = y = ""
        if USE_MOTORS:
            pos = input("Position 'X Y' in mm (empty = do not move): ").split()
            if len(pos) == 2:
                x, y = pos
                if not motors.move(float(x), float(y)):
                    print("Motors did not arrive or target out of limits: measurement cancelled")
                    continue
        data = scan(laser, scope) if SCAN else still_acquisition(scope)
        if data is None:
            print("Scope not READY: measurement cancelled")
            continue
        with registry.open(newline="") as f:
            n = max([int(r["id"]) for r in csv.DictReader(f)] + [0]) + 1
        file = f"m{n:03d}_{''.join(c if c.isalnum() or c in '-_' else '_' for c in name)}.csv"
        np.savetxt(FOLDER / file, data, delimiter=",", header="CH1,CH2", comments="")
        with registry.open("a", newline="") as f:
            csv.writer(f).writerow([n, datetime.datetime.now().isoformat(timespec="seconds"), name, x, y, SCAN, file, notes])
        print(f"Measurement {n} saved ({file}), CH2 median = {np.median(data[:, 1]):.5g} V")


if __name__ == "__main__":
    main()
