"""Raster scan of the whole sample: N_X x N_Y membranes, one laser scan on each, serpentine path.

    python raster_scan.py

The sample is tilted: the membrane centers lie on a rotated lattice, with the angle taken from FIRST_Y and FIRST_ROW_END_Y.
Close PI MicroMove before starting (a serial port can be used by one program only).
"""
import numpy as np

from acquisition import MOTOR_LIMITS_MM, scan_points
from initialization import init

FIRST_X = None                # <-- X of the center of the first membrane (top left), in mm
FIRST_Y = 133.978             # Y of the center of the first membrane
FIRST_ROW_END_Y = 134.078     # Y of the center of the last membrane of the first row (gives the tilt)
DESCRIPTION = ""            # what this run is (sample, setup...): saved in the run folder as description.txt, shown by analysis.py
N_X, N_Y = 11, 11             # membranes along X and along Y
PITCH = 0.30                  # mm between neighbouring membranes
RIGHT_IS_X_INCREASING = True  # "right" = increasing X
DOWN_IS_Y_DECREASING = True   # "down" = decreasing Y


def path():
    """Membrane centers along a serpentine: first row to the right, then down one row and back, and so on. Returns (points, tilt in degrees)."""
    sx = 1.0 if RIGHT_IS_X_INCREASING else -1.0
    sy = -1.0 if DOWN_IS_Y_DECREASING else 1.0
    uy = (FIRST_ROW_END_Y - FIRST_Y) / ((N_X - 1) * PITCH)       # sine of the tilt
    ux = sx * np.sqrt(1 - uy ** 2)                                  # unit vector along a row
    vx, vy = uy, -ux                                                # unit vector towards the next row
    if np.sign(vy) != sy:
        vx, vy = -uy, ux
    points = []
    for j in range(N_Y):
        for i in (range(N_X) if j % 2 == 0 else range(N_X - 1, -1, -1)):
            points.append((FIRST_X + PITCH * (i * ux + j * vx), FIRST_Y + PITCH * (i * uy + j * vy)))
    return points, np.degrees(np.arcsin(uy))


def main():
    if FIRST_X is None:
        print("Set FIRST_X (X of the first membrane) at the top of this file")
        return
    points, tilt = path()
    x, y = np.array(points).T
    print(f"{N_X} x {N_Y} membranes, pitch {PITCH} mm, tilt {tilt:.3f} degrees")
    print(f"Start at X={x[0]:.4f} Y={y[0]:.4f}, end at X={x[-1]:.4f} Y={y[-1]:.4f}")
    print(f"X from {x.min():.3f} to {x.max():.3f}, Y from {y.min():.3f} to {y.max():.3f} (limits {MOTOR_LIMITS_MM[0]}..{MOTOR_LIMITS_MM[1]} mm)")
    print(f"Disk space: about {len(points) * 0.24:.1f} GB")
    if min(x.min(), y.min()) < MOTOR_LIMITS_MM[0] or max(x.max(), y.max()) > MOTOR_LIMITS_MM[1]:
        print("Points outside the motor limits: check FIRST_X, FIRST_Y and the pitch")
        return
    if input(f"Move the motors and run {len(points)} scans? type 'yes': ").strip().lower() != "yes":
        return
    laser, scope, motors = init()
    print("Data in", scan_points(points, laser, scope, motors, description=DESCRIPTION))


if __name__ == "__main__":
    main()
