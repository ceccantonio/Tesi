"""Inverse design: finds the periodicity a and the hole radius of each membrane by fitting the simulated spectra to the measured ones.

    python inverse_design.py            # fits what is written in the settings below
    python inverse_design.py 12 40      # only the membranes 12 and 40 (and their plots)

The simulation (GHz_aniso_metasurface/theory/PPML_apr2026/sweep_a_r.m) gives the transmissions Txx, Tyy, Txy, Tyx on a grid of
(a, radius) in the measured window; the thickness is fixed. Each membrane has two runs, co-pol and cross-pol. A measured spectrum
(CH2, in volts) is modelled as  S = alpha * T(a, radius, wavelength + shift) + beta,  with one (alpha, beta) for each polarization
and the same (a, radius, shift) for both. The shape of the spectra decides: scale and offset are free, the shift absorbs an error
of the wavelength axis. The curves are compared directly, so several resonances in the window are fine.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline
from scipy.io import loadmat
from scipy.optimize import least_squares, nnls

import analysis

SWEEP = Path(__file__).resolve().parent.parent / "GHz_aniso_metasurface" / "theory" / "PPML_apr2026" / "sweep_a_r.mat"
RUN_CO = ""                  # run folder with the co-pol spectra, e.g. "2026-10-20/run_003"
RUN_CROSS = ""               # run folder with the cross-pol spectra of the same membranes (same indices in griglia.csv)
USE_REFERENCE = False        # True: each run is divided by its transmission reference (reference/copol.csv for the co-pol run, reference/crosspol.csv for the cross-pol one)
CO, CROSS = "Txx", "Txy"     # simulated components of the two runs ("Tyy", "Tyx" for the other pair)
MEMBRANES = None             # indices of the membranes to fit (e.g. [12, 40]); None = all the points of RUN_CO
MAX_SHIFT_NM = 2.0           # the fit may shift the measured wavelengths by up to +- this (error of the wavelength axis)
N_CANDIDATES = 5             # the best grid points (far from each other) from which the fit starts: the misfit has several minima


def load_sweep(path):
    m = loadmat(path, squeeze_me=True)
    return {"a": np.atleast_1d(m["a_values"]).astype(float), "r": np.atleast_1d(m["r_values"]).astype(float),
            "lam": np.atleast_1d(m["lambdav"]).astype(float), "T": {k: m[k] for k in ("Txx", "Tyy", "Txy", "Tyx")}}


class Model:
    """Simulated spectrum of one component as a function of (a, radius), by spline interpolation of the sweep."""
    def __init__(self, sweep, component):
        self.a, self.r, self.lam, self.grid = sweep["a"], sweep["r"], sweep["lam"], sweep["T"][component]
        self.splines = [RectBivariateSpline(sweep["a"], sweep["r"], self.grid[:, :, k]) for k in range(len(self.lam))]

    def __call__(self, a, r, lam, shift=0.0):
        return np.interp(lam + shift, self.lam, [s(a, r)[0, 0] for s in self.splines])


def measured(folder, row, lam):
    """Measured spectrum of a point of a run on the wavelengths lam (nm): block median of CH2 in the CH1 ramp, averaged in bins."""
    l, s = analysis.spectrum(*analysis.load_point(folder / row.file))
    step = np.median(np.diff(lam))
    edges = np.r_[lam - step / 2, lam[-1] + step / 2]
    count = np.histogram(l, edges)[0]
    return np.where(count > 0, np.histogram(l, edges, weights=s)[0] / np.maximum(count, 1), np.nan)


def scale_offset(T, s):
    """alpha >= 0 and beta of the best  s = alpha * T + beta,  and the residual normalized by the amplitude of s."""
    columns = np.column_stack([T, np.ones_like(T), -np.ones_like(T)])           # beta = second - third coefficient (any sign)
    c = nnls(columns, s)[0]
    return c[0], c[1] - c[2], (s - columns @ c) / np.std(s)


def fit_membrane(models, lam, data):
    """(a, radius, shift) that best reproduce the measured spectra (one per model) + the scale and offset of each. Returns a dict."""
    a_grid, r_grid = models[0].a, models[0].r

    def residuals(p):
        return np.concatenate([scale_offset(m(p[0], p[1], lam, p[2]), s)[2] for m, s in zip(models, data)])

    cost = np.full((len(a_grid), len(r_grid)), np.inf)                           # 1) misfit on the whole sweep grid (global search)
    shift_of = np.zeros(cost.shape)
    for ia in range(len(a_grid)):
        for ir in range(len(r_grid)):
            for shift in (-MAX_SHIFT_NM, 0.0, MAX_SHIFT_NM):
                c = sum(np.sum(scale_offset(np.interp(lam + shift, m.lam, m.grid[ia, ir]), s)[2] ** 2) for m, s in zip(models, data))
                if c < cost[ia, ir]:
                    cost[ia, ir], shift_of[ia, ir] = c, shift
    starts = []                                                                  # the best minima, at least 3 grid steps apart
    for flat in np.argsort(cost, axis=None):
        ia, ir = np.unravel_index(flat, cost.shape)
        if all(abs(ia - ja) >= 3 or abs(ir - jr) >= 3 for ja, jr, _ in starts):
            starts.append((ia, ir, shift_of[ia, ir]))
        if len(starts) == N_CANDIDATES:
            break
    best = None                                                                  # 2) refinement from each of them
    for ia, ir, shift in starts:
        fit = least_squares(residuals, [a_grid[ia], r_grid[ir], shift], loss="soft_l1", f_scale=0.3,
                            bounds=([a_grid[0], r_grid[0], -MAX_SHIFT_NM], [a_grid[-1], r_grid[-1], MAX_SHIFT_NM]))
        if best is None or fit.cost < best.cost:
            best = fit
    a, r, shift = best.x
    result = {"a_nm": a, "r_nm": r, "shift_nm": shift, "cost": best.cost}
    for name, m, s in zip(("co", "cross"), models, data):
        result["alpha_" + name], result["beta_" + name] = scale_offset(m(a, r, lam, shift), s)[:2]
    return result


def plot_fit(result, models, lam, data, label):
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for ax, key, name, m, s in zip(axes, ("co", "cross"), ("co-pol", "cross-pol"), models, data):
        ax.plot(lam, s, "k.", ms=3, label="measurement")
        ax.plot(lam, result["alpha_" + key] * m(result["a_nm"], result["r_nm"], lam, result["shift_nm"]) + result["beta_" + key],
                color="tab:green", lw=2, label="simulation (scale and offset fitted)")
        ax.set_ylabel(f"{name}: CH2 (V)")
        ax.legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("wavelength (nm)")
    axes[0].set_title(f"{label}: a = {result['a_nm']:.1f} nm, radius = {result['r_nm']:.1f} nm, wavelength shift = {result['shift_nm']:+.2f} nm")
    plt.tight_layout()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("membranes", nargs="*", type=int, help="indices of the membranes (empty = MEMBRANES)")
    args = parser.parse_args()
    wanted = args.membranes or MEMBRANES
    sweep = load_sweep(SWEEP)
    inside = (sweep["lam"] >= sweep["lam"][0] + MAX_SHIFT_NM) & (sweep["lam"] <= sweep["lam"][-1] - MAX_SHIFT_NM)
    lam = sweep["lam"][inside]                                                     # the wavelengths where a shift stays inside the sweep
    models = [Model(sweep, CO), Model(sweep, CROSS)]
    folders = [analysis.run_path(RUN_CO), analysis.run_path(RUN_CROSS)]
    rows_co = analysis.read_registry(folders[0])
    rows_cross = analysis.read_registry(folders[1]).set_index("indice")
    results = []
    for row in rows_co.itertuples():
        if (wanted and row.indice not in wanted) or row.indice not in rows_cross.index:
            continue
        analysis.POLARIZATION = "co" if USE_REFERENCE else ""                      # each run is divided by the reference of its polarization
        co = measured(folders[0], row, lam)
        analysis.POLARIZATION = "cross" if USE_REFERENCE else ""
        data = [co, measured(folders[1], rows_cross.loc[row.indice], lam)]
        ok = ~np.isnan(data[0]) & ~np.isnan(data[1])
        if ok.sum() < 10:
            print(f"membrane {row.indice}: not enough measured wavelengths in the window, skipped")
            continue
        result = fit_membrane(models, lam[ok], [d[ok] for d in data])
        result.update(index=row.indice, x_mm=row.x_mm, y_mm=row.y_mm)
        print(f"membrane {row.indice}: a = {result['a_nm']:.1f} nm, radius = {result['r_nm']:.1f} nm, shift = {result['shift_nm']:+.2f} nm, "
              f"misfit = {result['cost']:.3g}")
        results.append(result)
        if wanted:
            plot_fit(result, models, lam[ok], [d[ok] for d in data], f"Membrane {row.indice}")
    if not results:
        print("no membrane fitted: check RUN_CO, RUN_CROSS, MEMBRANES and the sweep file")
        return
    table = pd.DataFrame(results)
    table.to_csv(folders[0] / "inverse_design.csv", index=False)
    print("Saved", folders[0] / "inverse_design.csv")
    if not wanted and len(table) > 3:
        analysis.heatmap(table, "a_nm", "fitted periodicity a (nm)", folders[0] / "heatmap_a.png", fmt="{:.1f}")
        analysis.heatmap(table, "r_nm", "fitted hole radius (nm)", folders[0] / "heatmap_radius.png", fmt="{:.1f}")
    plt.show()


if __name__ == "__main__":
    main()
