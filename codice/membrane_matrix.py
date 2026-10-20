"""Matrix of points on ONE membrane: N_X x N_Y points from edge to edge, one laser scan on each, serpentine path.

    python membrane_matrix.py

Close PI MicroMove before starting (a serial port can be used by one program only).
"""
import numpy as np

from acquisition import MOTOR_LIMITS_MM, scan_points
from initialization import init

CENTER_X, CENTER_Y = 91.814, 146.508   # mm, center of the membrane
SIDE_UM = 500.0                        # the points go from -side/2 to +side/2 (edges included)
DESCRIPTION = ""            # what this run is (sample, setup...): saved in the run folder as description.txt, shown by analysis.py
N_X, N_Y = 8, 8
RIGHT_IS_X_INCREASING = True           # "right" = increasing X
DOWN_IS_Y_DECREASING = True            # "down" = decreasing Y; the first row is the top one


def path():
    """Serpentine points: top row to the right, then down one row and back, and so on."""
    xs = CENTER_X + np.linspace(-SIDE_UM / 2, SIDE_UM / 2, N_X) / 1000
    ys = CENTER_Y + np.linspace(-SIDE_UM / 2, SIDE_UM / 2, N_Y) / 1000
    rows = ys[::-1] if DOWN_IS_Y_DECREASING else ys
    rightwards = xs if RIGHT_IS_X_INCREASING else xs[::-1]
    return [(float(x), float(y)) for j, y in enumerate(rows) for x in (rightwards if j % 2 == 0 else rightwards[::-1])]


def main():
    points = path()
    x, y = np.array(points).T
    print(f"Matrix {N_X} x {N_Y}, X step {SIDE_UM / (N_X - 1):.2f} um, Y step {SIDE_UM / (N_Y - 1):.2f} um")
    print(f"X from {x.min():.4f} to {x.max():.4f}, Y from {y.min():.4f} to {y.max():.4f} (limits {MOTOR_LIMITS_MM[0]}..{MOTOR_LIMITS_MM[1]} mm)")
    print(f"Disk space: about {len(points) * 0.24:.1f} GB")
    if min(x.min(), y.min()) < MOTOR_LIMITS_MM[0] or max(x.max(), y.max()) > MOTOR_LIMITS_MM[1]:
        print("Points outside the motor limits: check the center")
        return
    if input(f"Move the motors and run {len(points)} scans? type 'yes': ").strip().lower() != "yes":
        return
    laser, scope, motors = init()
    print("Data in", scan_points(points, laser, scope, motors, description=DESCRIPTION))


if __name__ == "__main__":
    main()
