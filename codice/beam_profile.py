"""Beam radius w versus z and fit of the caustic, from the knife-edge measurements (one CSV per z, taken with knife_edge.py).

    python beam_profile.py                        # all the edge_*.csv of runs_oscilloscope/knife_edge/
    python beam_profile.py folder                 # all the edge_*.csv of that folder
    python beam_profile.py a.csv b.csv c.csv      # only these files

Every CSV is fitted with the error function (w at its z, written in the file). Then the points (z, w) are fitted with the caustic of a gaussian beam
    w(z) = w0 sqrt(1 + (M^2 lambda (z - z0) / (pi w0^2))^2)
with w0 = waist, z0 = its position, M^2 = beam quality (1 = ideal gaussian beam); lambda is fixed (WAVELENGTH_NM). Fit parameters: w0, z0, M^2. Derived:
    Rayleigh length                    zR = pi w0^2 / (M^2 lambda)
    divergence (far-field half angle)  theta = M^2 lambda / (pi w0) = w0 / zR
For a good fit take at least 6 z, some of them farther from the focus than zR (w larger than about 1.4 w0).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from knife_edge import FOLDER, fit_edge, read_edge

WAVELENGTH_NM = 1550.0   # for M^2 (the laser scans 1520-1570 nm)


def caustic(z, w0, z0, M2):
    """w in um, z in mm, lambda in nm: the units cancel in M^2 lambda (z - z0) / (pi w0^2)."""
    return w0 * np.sqrt(1 + (M2 * WAVELENGTH_NM * (z - z0) / (np.pi * w0 ** 2)) ** 2)


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

    zR0 = (z.max() - z.min()) / 2                                                  # first guess: a Rayleigh length of half the range
    p, cov = curve_fit(caustic, z, w, p0=[w.min(), z[np.argmin(w)], np.pi * w.min() ** 2 / (WAVELENGTH_NM * zR0)],
                       sigma=np.maximum(w_err, 1e-3 * w), absolute_sigma=False)
    p[0], p[2] = abs(p[0]), abs(p[2])
    w0, z0, M2 = p
    err = np.sqrt(np.diag(cov))
    zR, zR_err = with_error(lambda a, b, c: np.pi * a ** 2 / (c * WAVELENGTH_NM), p, cov)         # mm
    theta, theta_err = with_error(lambda a, b, c: c * WAVELENGTH_NM / (np.pi * a), p, cov)        # lambda in nm / w0 in um = mrad
    print(f"\nwaist w0 = {w0:.2f} +- {err[0]:.2f} um at z0 = {z0:.3f} +- {err[1]:.3f} mm")
    print(f"beam quality M^2 = {M2:.2f} +- {err[2]:.2f}  (lambda = {WAVELENGTH_NM:g} nm)")
    print(f"divergence (half angle) = {theta:.2f} +- {theta_err:.2f} mrad (full angle {2 * theta:.2f} mrad), Rayleigh length zR = {zR:.2f} +- {zR_err:.2f} mm")

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
