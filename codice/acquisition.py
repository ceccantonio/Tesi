"""Acquisition sequences shared by the scan scripts. They take instruments already initialized with initialization.init()."""
import datetime
import time
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "Chiara_Bru" / "runs_oscilloscope"   # YYYY-MM-DD/run_NNN folders go here
MOTOR_LIMITS_MM = (0.0, 150.0)                                                          # motor travel
RETURN_WAIT_S = 5.0                                                                     # after a scan, so that the laser goes back to 1520 nm


def scan(laser, scope):
    """One full scan of the laser synchronized with the scope. Returns (N, 2): CH1, CH2."""
    laser.set_range()
    scope.arm()
    laser.start()
    scope.wait()
    time.sleep(1.0)
    data = scope.read()
    time.sleep(RETURN_WAIT_S)
    return data


def still_acquisition(scope):
    """Acquisition without a scan (laser still, forced trigger). Returns (N, 2): CH1, CH2."""
    scope.arm()
    scope.force_trigger()
    scope.wait()
    time.sleep(1.0)
    return scope.read()


def new_folder(prefix="run"):
    """DATA/YYYY-MM-DD/<prefix>_NNN, created."""
    day = DATA / datetime.date.today().isoformat()
    day.mkdir(parents=True, exist_ok=True)
    folder = day / f"{prefix}_{len(list(day.glob(prefix + '_*'))) + 1:03d}"
    folder.mkdir()
    return folder


def scan_points(points, laser, scope, motors, settle=1.0):
    """For each (x, y): move the motors, do one scan, save one CSV. Registry: griglia.csv (indice,x_mm,y_mm,file,
    the same format as the earlier notebooks). Returns the folder."""
    folder = new_folder("run")
    registry = folder / "griglia.csv"
    registry.write_text("indice,x_mm,y_mm,file\n")
    try:
        for i, (x, y) in enumerate(points):
            print(f"[{i + 1}/{len(points)}] X={x:.4f}  Y={y:.4f}")
            if not motors.move(x, y):
                print("   motors did not arrive or target out of limits: point skipped")
                continue
            time.sleep(settle)
            name = f"point_{i:03d}_x{x:.3f}_y{y:.3f}.csv"
            np.savetxt(folder / name, scan(laser, scope), delimiter=",", header="CH1,CH2", comments="")
            with registry.open("a") as f:
                f.write(f"{i},{x:.4f},{y:.4f},{name}\n")
    except KeyboardInterrupt:
        laser.stop()
        print("Interrupted: laser scan stopped")
    return folder
