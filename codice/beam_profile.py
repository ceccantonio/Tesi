"""Beam radius w versus z and fit of the caustic, from the knife-edge measurements (one CSV per z, taken with knife_edge.py).

    python beam_profile.py                        # all the edge_*.csv of runs_oscilloscope/knife_edge/
    python beam_profile.py folder                 # all the edge_*.csv of that folder
    python beam_profile.py a.csv b.csv c.csv      # only these files
    python beam_profile.py --z 7600 8000          # only the measurements with z between 7600 and 8000 um (Z_RANGE in the code)
    python beam_profile.py --no-m2                # fit without M^2 as a free parameter: M^2 = M2_FIXED (1 = ideal gaussian beam), only w0 and z0 are fitted

Every CSV is fitted with the error function (w at its z, written in the file). Then the points (z, w) are fitted with the caustic of a gaussian beam
    w(z) = w0 sqrt(1 + (M^2 lambda (z - z0) / (pi w0^2))^2)
with w0 = waist, z0 = its position, M^2 = beam quality (1 = ideal gaussian beam); lambda is fixed (WAVELENGTH_NM). Fit parameters: w0, z0, M^2. Derived:
    Rayleigh length                    zR = pi w0^2 / (M^2 lambda)
    divergence (far-field half angle)  theta = M^2 lambda / (pi w0) = w0 / zR
z in um. For a good fit take at least 6 z, some of them farther from the focus than zR (w larger than about 1.4 w0): the z range must be at least about 2 zR.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from knife_edge import FOLDER, fit_edge, read_edge

WAVELENGTH_NM = 1550.0   # for M^2 (the laser scans 1520-1570 nm)
FIT_M2 = True            # True: M^2 is fitted (w0, z0, M^2); False: M^2 is fixed to M2_FIXED and only w0 and z0 are fitted (also with the option --no-m2)
M2_FIXED = 1.0           # M^2 used when it is not fitted
Z_RANGE = None           # (z_min, z_max) in um: only the measurements with z in this range are used; None = all (also with the option --z MIN MAX)


def caustic(z, w0, z0, M2):
    """w and z in um, lambda in nm (1e-3 converts it to um)."""
    return w0 * np.sqrt(1 + (M2 * WAVELENGTH_NM * 1e-3 * (z - z0) / (np.pi * w0 ** 2)) ** 2)


def with_error(f, p, cov):
    """f(*p) and its error from the covariance of the parameters (first order, correlations included)."""
    step = 1e-6 * np.maximum(abs(p), 1e-12)
    grad = np.array([(f(*(p + d)) - f(*(p - d))) / (2 * d[i]) for i, d in enumerate(np.diag(step))])
    return f(*p), float(np.sqrt(grad @ cov @ grad))


def main(args):
    z_range = Z_RANGE
    if "--z" in args:                                                            # --z MIN MAX
        i = args.index("--z")
        z_range = (float(args[i + 1]), float(args[i + 2]))
        del args[i:i + 3]
    fit_m2 = FIT_M2 and "--no-m2" not in args                                     # option: fit without M^2
    args = [a for a in args if not a.startswith("--")]                           # the options are not files
    files = []
    for a in (args or [FOLDER]):
        a = Path(a)
        files += sorted(a.glob("edge_*.csv")) if a.is_dir() else [a]
    z, w, w_err, names = [], [], [], []
    for f in files:
        x_mm, mean, std, z_um, description = read_edge(f)
        if z_um is None:
            print(f"{f.name}: no z in the file, skipped")
            continue
        _, width, error = fit_edge((x_mm - x_mm[0]) * 1000, mean)
        z.append(z_um)
        w.append(width)
        w_err.append(error)
        names.append(description)
    if len(z) < 4:
        sys.exit(f"{len(z)} measurements with z: at least 4 are needed to fit the caustic")
    order = np.argsort(z)
    z, w, w_err = np.array(z)[order], np.array(w)[order], np.array(w_err)[order]
    names = [names[i] for i in order]
    if z_range is not None:                                                        # only the z in the range
        keep = (z >= z_range[0]) & (z <= z_range[1])
        print(f"z between {z_range[0]:g} and {z_range[1]:g} um: {keep.sum()} of {len(z)} measurements")
        z, w, w_err = z[keep], w[keep], w_err[keep]
        names = [n for n, k in zip(names, keep) if k]
        if len(z) < 4:
            sys.exit("fewer than 4 measurements in that z range")
    print("   z (um)    w (um)   description")
    for zi, wi, ei, n in zip(z, w, w_err, names):
        print(f"{zi:9.1f}  {wi:7.2f} +- {ei:.2f}   {n}")

    sigma = np.maximum(w_err, 1e-3 * w)                                            # error of every point, at least 0.1% of w
    if fit_m2:
        zR0 = (z.max() - z.min()) / 2                                              # first guess: a Rayleigh length of half the range
        p, cov = curve_fit(caustic, z, w, p0=[w.min(), z[np.argmin(w)], np.pi * w.min() ** 2 / (WAVELENGTH_NM * 1e-3 * zR0)],
                           sigma=sigma, absolute_sigma=False)
        p[0], p[2] = abs(p[0]), abs(p[2])
        w0, z0, M2 = p
        err = np.sqrt(np.diag(cov))
        model = lambda zz: caustic(zz, *p)
        zR, zR_err = with_error(lambda a, b, c: np.pi * a ** 2 / (c * WAVELENGTH_NM * 1e-3), p, cov)         # um
        theta, theta_err = with_error(lambda a, b, c: c * WAVELENGTH_NM / (np.pi * a), p, cov)        # lambda in nm / w0 in um = mrad
    else:                                                                          # M^2 fixed: only w0 and z0 are fitted
        M2 = M2_FIXED
        fixed = lambda zz, w0, z0: caustic(zz, w0, z0, M2)
        p, cov = curve_fit(fixed, z, w, p0=[w.min(), z[np.argmin(w)]], sigma=sigma, absolute_sigma=False)
        p[0] = abs(p[0])
        w0, z0 = p
        err = np.sqrt(np.diag(cov))
        model = lambda zz: fixed(zz, *p)
        zR, zR_err = with_error(lambda a, b: np.pi * a ** 2 / (M2 * WAVELENGTH_NM * 1e-3), p, cov)           # um
        theta, theta_err = with_error(lambda a, b: M2 * WAVELENGTH_NM / (np.pi * a), p, cov)          # mrad
    chi2 = np.sum(((w - model(z)) / sigma) ** 2) / (len(z) - len(p))               # chi^2 per degree of freedom: about 1 or less = the curve describes the points
    print(f"\nwaist w0 = {w0:.2f} +- {err[0]:.2f} um at z0 = {z0:.1f} +- {err[1]:.2f} um")
    if fit_m2:
        print(f"beam quality M^2 = {M2:.2f} +- {err[2]:.2f}  (lambda = {WAVELENGTH_NM:g} nm)")
    else:
        print(f"beam quality M^2 = {M2:g} FIXED, not fitted  (lambda = {WAVELENGTH_NM:g} nm)")
    print(f"divergence (half angle) = {theta:.2f} +- {theta_err:.2f} mrad (full angle {2 * theta:.2f} mrad), Rayleigh length zR = {zR:.1f} +- {zR_err:.1f} um")
    print(f"chi^2 / dof = {chi2:.1f}")

    if z.max() - z.min() < 2 * zR:
        print(f"WARNING: the z range ({z.max() - z.min():.0f} um) is smaller than 2 zR ({2 * zR:.0f} um): w hardly changes, so M^2 and the divergence are poorly determined "
              "(see their errors). Measure at farther z (a longer translation, even coarse) to see w grow.")

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(z, w, w_err, fmt="o", color="k", label="knife edge")
    zs = np.linspace(z.min(), z.max(), 400)
    plt.plot(zs, model(zs), "r", label=f"fit: $w_0$ = {w0:.1f} um, $z_R$ = {zR:.1f} um\n$\\theta$ = {theta:.2f} mrad, $M^2$ = {M2:.2f}" + ("" if fit_m2 else " (fixed)") + f", chi2/dof = {chi2:.1f}")
    plt.xlabel("z (um)")
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
