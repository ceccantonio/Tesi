"""Beam radius w versus z and fit of the caustic, from the knife-edge measurements (one CSV per z, taken with knife_edge.py).

    python beam_profile.py                        # all the edge_*.csv of runs_oscilloscope/knife_edge/
    python beam_profile.py folder                 # all the edge_*.csv of that folder
    python beam_profile.py a.csv b.csv c.csv      # only these files

Every CSV is fitted with the error function (w at its z, written in the file). Then the points (z, w) are fitted with the caustic of a gaussian beam
    w(z) = w0 sqrt(1 + ((z - z0) / zR)^2)
with w0 = waist, z0 = its position, zR = Rayleigh length. Results:
    divergence (far-field half angle)  theta = w0 / zR
    beam quality                       M^2 = pi w0^2 / (lambda zR) = pi w0 theta / lambda      (M^2 = 1: ideal gaussian beam)
For a good fit take at least 6 z, some of them farther from the focus than zR (w larger than about 1.4 w0).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from knife_edge import FOLDER, fit_edge, read_edge

WAVELENGTH_NM = 1550.0   # for M^2 (the laser scans 1520-1570 nm)


def caustic(z, w0, z0, zR):
    return w0 * np.sqrt(1 + ((z - z0) / zR) ** 2)


def with_error(f, p, cov):
    """f(*p) and its error from the covariance of the parameters (first order, correlations included)."""
    step = 1e-6 * np.maximum(abs(p), 1e-12)
    grad = np.array([(f(*(p + d)) - f(*(p - d))) / (2 * d[i]) for i, d in enumerate(np.diag(step))])
    return f(*p), float(np.sqrt(grad @ cov @ grad))


def main(args):
    files = []
    for a in (args or [FOLDER]):
        a = Path(a)
        files += sorted(a.glob("edge_*.csv")) if a.is_dir() else [a]
    z, w, w_err, names = [], [], [], []
    for f in files:
        x_mm, mean, std, z_mm, description = read_edge(f)
        if z_mm is None:
            print(f"{f.name}: no z in the file, skipped")
            continue
        _, width, error = fit_edge((x_mm - x_mm[0]) * 1000, mean)
        z.append(z_mm)
        w.append(width)
        w_err.append(error)
        names.append(description)
    if len(z) < 4:
        sys.exit(f"{len(z)} measurements with z: at least 4 are needed to fit the caustic")
    order = np.argsort(z)
    z, w, w_err = np.array(z)[order], np.array(w)[order], np.array(w_err)[order]
    names = [names[i] for i in order]
    print("   z (mm)    w (um)   description")
    for zi, wi, ei, n in zip(z, w, w_err, names):
        print(f"{zi:9.3f}  {wi:7.2f} +- {ei:.2f}   {n}")

    p, cov = curve_fit(caustic, z, w, p0=[w.min(), z[np.argmin(w)], (z.max() - z.min()) / 2], sigma=np.maximum(w_err, 1e-3 * w), absolute_sigma=False)
    p[0], p[2] = abs(p[0]), abs(p[2])
    w0, z0, zR = p
    err = np.sqrt(np.diag(cov))
    theta, theta_err = with_error(lambda a, b, c: a / c, p, cov)                  # um / mm = mrad
    M2, M2_err = with_error(lambda a, b, c: np.pi * a ** 2 / (WAVELENGTH_NM * c), p, cov)       # w0 in um, lambda in nm, zR in mm
    print(f"\nwaist w0 = {w0:.2f} +- {err[0]:.2f} um at z0 = {z0:.3f} +- {err[1]:.3f} mm, Rayleigh length zR = {zR:.2f} +- {err[2]:.2f} mm")
    print(f"divergence (half angle) = {theta:.2f} +- {theta_err:.2f} mrad (full angle {2 * theta:.2f} mrad)")
    print(f"beam quality M^2 = {M2:.2f} +- {M2_err:.2f}  (lambda = {WAVELENGTH_NM:g} nm)")

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(z, w, w_err, fmt="o", color="k", label="knife edge")
    zs = np.linspace(z.min(), z.max(), 400)
    plt.plot(zs, caustic(zs, *p), "r", label=f"fit: $w_0$ = {w0:.1f} um, $z_R$ = {zR:.1f} mm\n$\\theta$ = {theta:.2f} mrad, $M^2$ = {M2:.2f}")
    plt.xlabel("z (mm)")
    plt.ylabel("beam radius w (um)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = files[0].parent / "beam_profile.png"                                   # next to the CSV files
    plt.savefig(out, dpi=150)
    print("Saved", out)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
