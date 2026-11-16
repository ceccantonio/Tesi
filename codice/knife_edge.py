"""Knife-edge measurement of the beam waist: the X motor moves the cutting edge across the beam in steps of STEP_UM; at each position
the scope reads ONE channel (laser still, no scan) and the signal versus position is fitted with an error function.

    python knife_edge.py                 # measure (asks z and a description, then for confirmation); CSV and plot in runs_oscilloscope/knife_edge/
    python knife_edge.py file.csv        # only the fit of an old measurement, no instruments

For a gaussian beam, intensity exp(-2 r^2 / w^2), the power that passes the edge is P(x) = B + A/2 (1 + erf(sqrt(2) (x - x0) / w)),
where w is the beam radius (1/e^2) at that z. A can be negative (signal that decreases as the edge advances).
Before starting, put the edge where the signal is zero (beam completely blocked) or full (beam completely free): the scan starts from the current X and goes
in ONE direction (DIRECTION), for SPAN_UM, until the beam is completely free (or blocked); then the edge goes back to the starting X.
Every measurement is saved with its z (in um) and a description (first lines of the CSV, starting with #): beam_profile.py uses them to plot w(z) and fit the divergence.
Close PI MicroMove first (a serial port can be used by one program only).
"""
import datetime
import sys
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

from acquisition import DATA

CHANNEL = "CH3"      # the only channel read
SPAN_UM = 600.0      # length of the scan from the current position: it must be larger than the beam (about 4 w)
DIRECTION = 1        # +1 = X increasing, -1 = X decreasing
STEP_UM = 1.0
N_SAMPLES = 100000   # samples read at each point: their mean is the signal
SETTLE_S = 0.2       # pause after each move
FOLDER = DATA / "knife_edge"


def edge(x, B, A, x0, w):
    return B + A / 2 * (1 + erf(np.sqrt(2) * (x - x0) / abs(w)))


def fit_edge(x_um, y):
    """Error function fitted to (position in um, signal): parameters (B, A, x0, w), the waist w (um) and its error."""
    A0 = y[-1] - y[0]
    x0 = x_um[np.argmin(abs(y - (y[0] + A0 / 2)))]
    p, cov = curve_fit(edge, x_um, y, p0=[y[0], A0, x0, abs(x_um[-1] - x_um[0]) / 4])
    return p, abs(p[3]), np.sqrt(cov[3, 3])


def read_edge(path):
    """(x_mm, mean_V, std_V, z_um or None, description) from a CSV of this script (lines starting with # are the z, the description and the header)."""
    z, notes, rows = None, [], []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            text = line.lstrip("# ").strip()
            if text.startswith("z_um"):
                z = float(text.split("=")[1])
            elif text.startswith("z_mm"):                                # files written in mm by an older version
                z = float(text.split("=")[1]) * 1000
            elif text and not text.startswith("x_mm"):
                notes.append(text)
        else:
            try:
                rows.append([float(v) for v in line.split(",")])
            except ValueError:
                pass                                                     # a header line without #
    x_mm, mean, std = np.array(rows).T
    return x_mm, mean, std, z, " ".join(notes)


def report(x_um, y, y_err, name, title=""):
    """Fit of the error function: prints the waist and saves <name>.png."""
    p, w, w_err = fit_edge(x_um, y)
    print(f"waist w = {w:.2f} +- {w_err:.2f} um (1/e^2 radius), edge centre x0 = {p[2]:.2f} um, "
          f"10-90% distance {1.2816 * w:.2f} um, FWHM of the intensity {1.1774 * w:.2f} um")
    plt.figure(figsize=(7, 4.5))
    plt.errorbar(x_um, y, y_err, fmt="o", ms=3, color="k", label="data")
    xs = np.linspace(x_um.min(), x_um.max(), 500)
    plt.plot(xs, edge(xs, *p), "r", label=f"error function: w = {w:.1f} $\\pm$ {w_err:.1f} um")
    plt.title(title, fontsize=9)
    plt.xlabel("edge position (um)")
    plt.ylabel(f"{CHANNEL} (V)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(name.with_suffix(".png"), dpi=150)
    print("Saved", name.with_suffix(".png"))
    plt.show()


def acquire(scope):
    """Mean and standard deviation of CHANNEL in one acquisition, or None if the scope did not finish it."""
    if not scope.acquire():
        return None
    data = scope.read((CHANNEL,))[:, 0]
    return data.mean(), data.std()


def ask_z():
    while True:
        try:
            return float(input("z of this measurement (um): ").replace(",", "."))
        except ValueError:
            print("  write a number")


def measure():
    from initialization import init                                    # only here, so that the fit of an old file needs no instruments
    z = ask_z()
    description = input("Description (setup, what you changed; Enter for none): ").strip()
    _, scope, motors = init(laser=False, n_samples=N_SAMPLES)
    axis = motors.x
    axis.tolerance = 0.002                                             # mm; the real position is saved anyway
    start = axis.pos
    targets = start + DIRECTION * np.arange(0, SPAN_UM + STEP_UM / 2, STEP_UM) / 1000
    print(f"{len(targets)} points from X = {targets[0]:.4f} to {targets[-1]:.4f} mm (step {STEP_UM} um), channel {CHANNEL}, z = {z:g} um")
    if input("Move the X motor and measure? type 'yes': ").strip().lower() != "yes":
        return
    rows = []
    try:
        for i, target in enumerate(targets):
            if not axis.move(float(target)):
                print(f"[{i + 1}/{len(targets)}] X = {target:.4f}: motor did not arrive or target out of limits, skipped")
                continue
            time.sleep(SETTLE_S)
            result = acquire(scope)
            if result is None:
                print(f"[{i + 1}/{len(targets)}] the scope did not finish the acquisition, skipped")
                continue
            rows.append((axis.pos, *result))
            print(f"[{i + 1}/{len(targets)}] X = {rows[-1][0]:.4f} mm   {CHANNEL} = {result[0]:.5g} V")
    except KeyboardInterrupt:
        print("Interrupted")
    axis.move(float(start))                                            # the edge goes back where it was
    if len(rows) < 10:
        print("Too few points to fit")
        return
    FOLDER.mkdir(parents=True, exist_ok=True)
    name = FOLDER / f"edge_{datetime.datetime.now():%Y-%m-%d_%H%M%S}"
    np.savetxt(name.with_suffix(".csv"), rows, delimiter=",", header=f"z_um = {z:g}\n{description}\nx_mm,mean_V,std_V", comments="# ")
    print("Saved", name.with_suffix(".csv"))
    x_mm, mean, std = np.array(rows).T
    report((x_mm - x_mm[0]) * 1000, mean, std / np.sqrt(N_SAMPLES), name, f"z = {z:g} um  {description}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        x_mm, mean, std, z, description = read_edge(sys.argv[1])
        report((x_mm - x_mm[0]) * 1000, mean, std / np.sqrt(N_SAMPLES), Path(sys.argv[1]),
               (f"z = {z:g} um  " if z is not None else "") + description)
    else:
        measure()
