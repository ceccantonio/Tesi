"""Data analysis: modulation amplitude, quality factor, heatmaps and comparison of single measurements.

    python analysis.py run 2026-10-06/run_007       # modulation and Q of every point -> analysis.csv + two heatmaps
    python analysis.py point 2026-10-06/run_007 3 12 40   # detail plots of some points of a run (no heatmaps)
    python analysis.py view misure_singole/m004_x.csv     # plot of single acquisitions (CSV files), no run folder needed
    python analysis.py fano 2026-10-06/run_007/point_011_x91.818_y146.526.csv   # Fano fit of the resonances of one CSV
    python analysis.py list                         # list of the single measurements
    python analysis.py compare 1 2 3 --ref 1        # compare single measurements and divide them by measurement 1
    python analysis.py                              # no arguments: does what is written in "WHAT TO LOOK AT" below

Method (as in movemin.py): moving average of CH2 (np.convolve) inside the CH1 ramp.
 - exact start and end of the ramp: fit of low plateau + linear rise + high plateau (fit_ramp), which also uses points
   before the start and after the end; the x axis is the wavelength: LAMBDA_0 at the start and LAMBDA_1 at the end of the
   ramp, by linear proportion (assumption to verify: the scan goes from 1520 to 1570 nm over the whole ramp);
 - modulation = max - min of the moving average inside the ramp;
 - Q of each peak = height (prominence) / width at half height, in volts per sample;
   each peak also gets its wavelength and its width in nm.

Fano fit (around each located peak): y = y0 + A (q + eps)^2 / (1 + eps^2), eps = 2 (x - x0) / Gamma, with x in nm: x0 = resonance wavelength,
Gamma = width, q = asymmetry, A >= 0 (so that q is unique); Q = x0 / Gamma. The resonances are located on the moving average; the fit is
always done in a window around each located peak (FANO_WINDOW widths on each side, or FANO_RANGE chosen by hand), on the raw
CH2 data or, with FANO_ON_AVERAGE, on the moving average. A quick first fit on block-median data finds the starting values, then the
fit is refined on all the samples of the window (FANO_ALL_SAMPLES).

Run folders contain griglia.csv (columns indice, x_mm, y_mm, file) and one CSV (columns CH1, CH2) per point:
that format is the one written by the earlier notebooks and by acquisition.py.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.optimize import least_squares
from scipy.signal import fftconvolve, find_peaks, peak_widths

DATA = Path(__file__).resolve().parent.parent / "Chiara_Bru" / "runs_oscilloscope"

# ---------------------------------------------------------------- WHAT TO LOOK AT
# Used when you run `python analysis.py` with no arguments (or call main() from your own code / a notebook).
RUN = ""                    # run folder (e.g. "2026-10-06/run_007"); "" = skip
HEATMAPS = True             # True: draw the heatmaps of RUN (modulation and Q); False: skip them (e.g. to look only at some points)
POINT = None                # point(s) of RUN to plot in detail: an index (12) or a list ([3, 12, 40]); None = skip
VIEW = None                 # CSV file(s) of single acquisitions to plot, no run folder needed: a path or a list of paths; None = skip
COMPARE = None              # ids of the single measurements to compare (e.g. [1, 2, 3], [] = all); None = skip
REFERENCE = None            # id of the reference measurement for the ratios (None = the first one)
RECOMPUTE = False           # True: ignore analysis.csv and read all the CSV files again
FANO = ""                   # CSV file to fit with a Fano line shape (e.g. "2026-10-06/run_007/point_011_x91.818_y146.526.csv"); "" = skip

# ---------------------------------------------------------------- parameters of the method
WINDOW = 10000              # moving-average length, in samples
FAST_MOVING_AVERAGE = False # True: same moving average computed with an FFT (identical result, much faster on long records); False: np.convolve as in movemin.py
V_START = 0.9               # the CH1 ramp starts when CH1 exceeds this value...
END_TOLERANCE = 0.3         # ...and ends when CH1 is within this of its maximum (starting guess for the ramp fit)
RAMP_MARGIN = 0.25          # the ramp fit also uses points before and after: this fraction of the duration on each side
LAMBDA_0, LAMBDA_1 = 1520.0, 1570.0   # nm at the start and at the end of the CH1 ramp, i.e. of the laser scan
PROMINENCE_THRESHOLD = 0.1  # a peak counts if taller than this fraction of the excursion (modulation)...
MIN_SNR = 10                # ...and than this multiple of the residual noise of the moving average...
RUN_THRESHOLD = 0.1         # ...and than this fraction of the tallest peak of the whole run (no resonance = n.a.)
FANO_PEAKS = 1              # Fano fit: how many resonances to fit (the tallest ones); None = all the located ones
FANO_WINDOW = 3.0           # the fit is always done around the located peak, on the data within +- this many widths (FWHM) of it
FANO_RANGE = None           # (min_nm, max_nm) to choose the fit window by hand; None = automatic
FANO_ON_AVERAGE = False     # False: fit the raw CH2 data; True: fit the moving average of CH2
FANO_SLOPE = False          # True: adds a sloped baseline s * (x - x0) to the model
FANO_ROBUST = True          # True: soft_l1 loss (resists spurious spikes); False: least squares
FANO_POINTS = 5000          # the first, quick fit (it finds the starting values) uses the window data reduced to about this many points (block median)
FANO_ALL_SAMPLES = True     # True: then the fit is refined on ALL the samples of the window; False: stop at the reduced data (as in the old notebook)
FANO_DIP = False            # True: the resonances are minima of CH2 (to locate them); False: maxima
Q_UNIT = 1e-6               # Q is shown in units of this (V per sample)
DARK_LEVEL = 0.0            # CH2 level without light (single measurements)


# ---------------------------------------------------------------- reading and basic method

def run_path(folder):
    p = Path(folder)
    return p if p.is_absolute() or p.exists() else DATA / p


def read_csv(path):
    df = pd.read_csv(path, dtype=np.float32)
    return df["CH1"].to_numpy(float), df["CH2"].to_numpy(float)


def moving_average(y, n=WINDOW):
    """Same as movemin.py: odd window, np.convolve in 'same' mode (edges padded with zeros).
    With FAST_MOVING_AVERAGE the same convolution is done with an FFT: np.convolve costs N x window operations (minutes on a
    5M-sample record), the FFT a couple of seconds."""
    n = int(n) | 1
    kernel = np.ones(n) / n
    return fftconvolve(y, kernel, mode="same") if FAST_MOVING_AVERAGE else np.convolve(y, kernel, mode="same")


def ramp(ch1):
    """(start, end) of the CH1 ramp in samples, from the threshold, or None."""
    above = np.flatnonzero(ch1 > V_START)
    if above.size == 0:
        return None
    i0 = above[0]
    i1 = i0 + np.flatnonzero(ch1[i0:] >= ch1.max() - END_TOLERANCE)[0]
    return (i0, i1) if i1 > i0 else None


def fit_ramp(ch1):
    """Exact start and end of the CH1 ramp, in samples (real numbers), or None.
    Fit of low plateau + linear rise + high plateau, also using points before the start and after the end
    (seen from the V_START threshold the rise has already begun)."""
    r = ramp(ch1)
    if r is None:
        return None
    i0, i1 = r
    length = i1 - i0
    first = max(0, i0 - int(RAMP_MARGIN * length) - int(0.15 * length))
    last = min(len(ch1) - 1, i1 + int(RAMP_MARGIN * length))
    falling = np.flatnonzero(ch1[i1:last + 1] < ch1.max() - 1.0)       # the high plateau ends when CH1 starts to fall
    if falling.size:
        last = i1 + falling[0]
    step = max(1, (last - first + 1) // 20000)
    t = np.arange(first, last + 1, step, dtype=float)
    v = ch1[first:last + 1:step]

    def model(p):
        return p[0] + (p[1] - p[0]) * np.clip((t - p[2]) / p[3], 0.0, 1.0)

    p0 = [np.median(v[:max(3, len(v) // 10)]), ch1.max(), i0 - 0.1 * length, 1.1 * length]
    fit = least_squares(lambda p: model(p) - v, p0, loss="soft_l1", f_scale=0.3,
                        bounds=([-np.inf, -np.inf, first, 10.0], [np.inf, np.inf, last, 3.0 * length]))
    return float(fit.x[2]), float(fit.x[2] + fit.x[3])


def wavelength(samples, f):
    """Wavelength (nm) of the samples by linear proportion: LAMBDA_0 at the ramp start f[0], LAMBDA_1 at the ramp end f[1]."""
    return LAMBDA_0 + (np.asarray(samples, dtype=float) - f[0]) / (f[1] - f[0]) * (LAMBDA_1 - LAMBDA_0)


def valid_zone(f, length):
    """Samples of the ramp (from fit_ramp) without WINDOW/2 at the record edges, where the moving average uses zeros."""
    if f is None:
        return None
    half = (int(WINDOW) | 1) // 2
    a, b = max(int(np.ceil(f[0])), half), min(int(np.floor(f[1])), length - 1 - half)
    return (a, b) if b > a else None


def analyze_point(ch1, ch2, smooth=None):
    """(modulation in V, peaks, ramp (start, end) or None). Each peak is a dict: index, height, width, Q, left, right,
    level, lambda_nm, width_nm. smooth = moving_average(ch2) if you already have it (it is the slow step)."""
    f = fit_ramp(ch1)
    zone = valid_zone(f, len(ch1))
    if zone is None:
        return np.nan, [], f
    a, b = zone
    curve = (moving_average(ch2) if smooth is None else smooth)[a:b + 1]
    modulation = curve.max() - curve.min()
    d = np.diff(ch2[a:b + 1])
    sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.sqrt(WINDOW)   # residual noise of the moving average
    found, prop = find_peaks(curve, prominence=max(PROMINENCE_THRESHOLD * modulation, MIN_SNR * sigma), distance=WINDOW)
    if len(found) == 0:
        return modulation, [], f
    widths, levels, left, right = peak_widths(curve, found, rel_height=0.5,
                                              prominence_data=(prop["prominences"], prop["left_bases"], prop["right_bases"]))
    nm_per_sample = (LAMBDA_1 - LAMBDA_0) / (f[1] - f[0])
    return modulation, [{"index": a + int(p), "height": float(h), "width": float(w), "Q": float(h / w),
                         "left": a + float(lo), "right": a + float(hi), "level": float(lv),
                         "lambda_nm": float(wavelength(a + p, f)), "width_nm": float(w * nm_per_sample)}
                        for p, h, w, lv, lo, hi in zip(found, prop["prominences"], widths, levels, left, right)], f


# ---------------------------------------------------------------- analysis of a run

def analyze_run(folder, recompute=False):
    """Modulation and Q of every point of the run (folder with griglia.csv). The table is saved in analysis.csv."""
    folder = run_path(folder)
    if (folder / "analysis.csv").exists() and not recompute:
        return pd.read_csv(folder / "analysis.csv")
    points = []
    for r in pd.read_csv(folder / "griglia.csv").itertuples():
        if (folder / r.file).exists():
            modulation, peaks, f = analyze_point(*read_csv(folder / r.file))
            print(f"point {r.indice}: modulation {modulation:.4f} V, {len(peaks)} peaks")
            points.append((r, modulation, peaks, f))
    threshold = RUN_THRESHOLD * max([p["height"] for _, _, peaks, _ in points for p in peaks], default=0)
    rows = []
    for r, modulation, peaks, f in points:
        best = max([p for p in peaks if p["height"] >= threshold], key=lambda p: p["height"], default=None)   # the tallest peak
        rows.append({"index": r.indice, "x_mm": r.x_mm, "y_mm": r.y_mm, "modulation_V": modulation,
                     "Q": best["Q"] if best else np.nan, "height_V": best["height"] if best else np.nan,
                     "width_samples": best["width"] if best else np.nan,
                     "lambda_nm": best["lambda_nm"] if best else np.nan, "width_nm": best["width_nm"] if best else np.nan,
                     "ramp_start": f[0] if f else np.nan, "ramp_end": f[1] if f else np.nan})
    table = pd.DataFrame(rows)
    table.to_csv(folder / "analysis.csv", index=False)
    return table


def lattice(x, y):
    """(row, column) of every point of a serpentine path; None if it cannot be reconstructed."""
    jumps = np.abs(np.diff(y))
    if len(x) < 4 or jumps.max() == 0:
        return None
    rows = np.concatenate([[0], np.cumsum(jumps > 0.5 * jumps.max())])         # the change of row is a big jump in Y
    lengths = np.bincount(rows)
    if len(set(lengths)) != 1:
        return None
    columns = np.concatenate([np.arange(n) if r % 2 == 0 else np.arange(n)[::-1] for r, n in enumerate(lengths)])
    return rows, columns


def heatmap(table, column, label, file=None, scale=1.0, fmt="{:.2f}"):
    """Heatmap with the value written on each point: on a lattice (rows/columns of the serpentine) or at the real X, Y positions."""
    x, y = table["x_mm"].to_numpy(), table["y_mm"].to_numpy()
    z = table[column].to_numpy() / scale
    valid = np.isfinite(z)
    if not valid.any():
        print(f"{column}: no valid value, heatmap not drawn")
        return
    norm, cmap = Normalize(z[valid].min(), z[valid].max()), plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(9, 7))
    rc = lattice(x, y)
    if rc is not None:
        rows, columns = rc
        grid = np.full((rows.max() + 1, columns.max() + 1), np.nan)
        grid[rows, columns] = z
        ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, norm=norm, origin="upper")
        cx, cy = columns, rows
        ax.set_xlabel("column (along the first row)")
        ax.set_ylabel("row (from the top)")
    else:
        d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
        pitch = d[d > 0].min() if len(x) > 1 else 0.3
        for xi, yi, zi in zip(x, y, z):
            ax.add_patch(Rectangle((xi - pitch / 2, yi - pitch / 2), pitch, pitch, edgecolor="white",
                                   facecolor=cmap(norm(zi)) if np.isfinite(zi) else "lightgray"))
        ax.set_xlim(x.min() - pitch, x.max() + pitch)
        ax.set_ylim(y.min() - pitch, y.max() + pitch)
        ax.set_aspect("equal")
        cx, cy = x, y
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
    for xi, yi, zi in zip(cx, cy, z):
        ax.text(xi, yi, fmt.format(zi) if np.isfinite(zi) else "n.a.", ha="center", va="center", fontsize=9,
                color="gray" if not np.isfinite(zi) else "black" if norm(zi) > 0.55 else "white")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=label)
    ax.set_title(label)
    plt.tight_layout()
    if file:
        fig.savefig(file, dpi=200)
        print("Saved", file)
    plt.show()


def run_heatmaps(folder, recompute=False):
    """Heatmaps of the modulation and of the quality factor of a run."""
    folder = run_path(folder)
    table = analyze_run(folder, recompute)
    heatmap(table, "modulation_V", "modulation amplitude: max - min of the CH2 moving average in the ramp (V)",
            folder / "heatmap_modulation.png", fmt="{:.3f}")
    heatmap(table, "Q", f"quality factor Q = height / width  [{Q_UNIT:g} V per sample]",
            folder / "heatmap_quality_factor.png", scale=Q_UNIT)


def point_plot(folder, index, show=True):
    """Detail plot of one point of a run (see csv_plot). show=False: draw it without opening it yet."""
    folder = run_path(folder)
    r = pd.read_csv(folder / "griglia.csv").query("indice == @index").iloc[0]
    csv_plot(folder / r.file, f"Point {index} (X={r.x_mm}, Y={r.y_mm})", show)


def csv_plot(path, label=None, show=True):
    """Plot of ONE acquisition (any CSV with CH1, CH2: a point of a run, a single measurement...): CH1 ramp with its fit,
    raw CH2, moving average, measurement zone and peaks, in wavelength.
    show=False: draw the figure without opening it yet (to open several together with one plt.show())."""
    path = run_path(path)
    ch1, ch2 = read_csv(path)
    smooth = moving_average(ch2)
    modulation, peaks, f = analyze_point(ch1, ch2, smooth)
    if f is None:
        print("no ramp in CH1: nothing to plot")
        return
    step = max(1, len(ch2) // 20000)
    lam = wavelength(np.arange(len(ch2)), f)[::step]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    a1.plot(lam, ch1[::step], color="tab:red", lw=1)
    a1.axvline(LAMBDA_0, color="k", ls="--", lw=0.8, label=f"ramp start (fit) = {LAMBDA_0:g} nm")
    a1.axvline(LAMBDA_1, color="k", ls=":", lw=0.8, label=f"ramp end (fit) = {LAMBDA_1:g} nm")
    a1.set_ylabel("CH1 (V)")
    a1.legend(loc="best", fontsize=8)
    a2.plot(lam, ch2[::step], color="gray", lw=0.6, alpha=0.6, label="raw CH2")
    a2.plot(lam, smooth[::step], color="red", lw=1.5, label=f"moving average (N={WINDOW})")
    zone = valid_zone(f, len(ch1))
    if zone:
        a2.axvspan(*wavelength(zone, f), color="orange", alpha=0.12, label="measurement zone")
    for n, p in enumerate(peaks, start=1):
        a2.hlines(p["level"], *wavelength([p["left"], p["right"]], f), color="blue", lw=2)
        a2.annotate(f"{n}: {p['lambda_nm']:.2f} nm, Q={p['Q'] / Q_UNIT:.3g}", (p["lambda_nm"], smooth[p["index"]]),
                    textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9)
    a2.set_xlim(LAMBDA_0 - 0.3 * (LAMBDA_1 - LAMBDA_0), LAMBDA_1 + 0.3 * (LAMBDA_1 - LAMBDA_0))
    a2.set_xlabel("wavelength (nm)")
    a2.set_ylabel("CH2 (V)")
    a2.legend(loc="best", fontsize=8)
    a1.set_title(f"{label or path.name}: modulation {modulation:.4f} V, {len(peaks)} peaks")
    plt.tight_layout()
    if show:
        plt.show()


# ---------------------------------------------------------------- Fano fit

def fano(x, y0, A, x0, gamma, q, s=0.0):
    """Fano line shape: y0 + s*(x - x0) + A * (q + eps)^2 / (1 + eps^2), eps = 2 (x - x0) / gamma."""
    eps = 2.0 * (x - x0) / gamma
    return y0 + s * (x - x0) + A * (q + eps) ** 2 / (1.0 + eps ** 2)


def block_median(x, y, n_points):
    """Reduce (x, y) to about n_points with a block median: removes the spurious spikes and makes the fit fast."""
    n = max(1, len(x) // int(n_points))
    k = len(x) // n
    return x[:k * n].reshape(k, n).mean(axis=1), np.median(y[:k * n].reshape(k, n), axis=1)


def fano_fit_peak(x, y, x0_guess, gamma_guess, x_all=None, y_all=None):
    """Fano fit of one resonance. Stage 1 on (x, y), the data reduced for speed: starts from several values of q and keeps the best
    one (the Fano fit has several local minima). Stage 2, if (x_all, y_all) are given: the fit is refined on ALL the raw samples,
    starting from the stage-1 result. A >= 0 is imposed: (A, q) and (-A q^2, -1/q) describe the same curve, so without it q
    would not be unique. Returns a dict with the parameters, their uncertainties, Q = x0 / Gamma and the model."""
    scale = float(gamma_guess)
    span = float(np.max(y) - np.min(y))
    edge = max(3, len(y) // 5)
    y0_guess = float(np.median(np.concatenate([y[:edge], y[-edge:]])))
    slope = bool(FANO_SLOPE)
    lower = [-np.inf, 0.0, -2.0, 0.05, -50.0] + ([-np.inf] if slope else [])
    upper = [np.inf, 20 * span, 2.0, 20.0, 50.0] + ([np.inf] if slope else [])
    f_scale = max(1.4826 * float(np.median(np.abs(y - np.median(y)))) * 0.1, span * 0.02, 1e-12)
    floor = float(np.percentile(y, 5))                                         # the minimum of the structure
    A_guess = max(y0_guess - floor, 0.1 * span)

    def model(p, xx):                                                            # xx normalized: the resonance is near 0, width ~1
        return fano((xx - x0_guess) / scale, p[0], p[1], p[2], p[3], p[4], p[5] if slope else 0.0)

    def jacobian(p, xx):                                                         # analytic: much faster than finite differences
        eps = 2.0 * ((xx - x0_guess) / scale - p[2]) / p[3]
        den = 1.0 + eps ** 2
        d_eps = 2.0 * (p[4] + eps) * (1.0 - p[4] * eps) / den ** 2             # d[(q + eps)^2 / (1 + eps^2)] / d eps
        columns = [np.ones_like(eps), (p[4] + eps) ** 2 / den, p[1] * d_eps * (-2.0 / p[3]) - (p[5] if slope else 0.0),
                   p[1] * d_eps * (-eps / p[3]), 2.0 * p[1] * (p[4] + eps) / den]
        if slope:
            columns.append((xx - x0_guess) / scale - p[2])
        return np.column_stack(columns)

    def solve(xx, yy, starts):
        residuals = lambda p: model(p, xx) - yy
        fits = [least_squares(residuals, p0, jac=lambda p: jacobian(p, xx), bounds=(lower, upper),
                              loss="soft_l1" if FANO_ROBUST else "linear", f_scale=f_scale) for p0 in starts]
        return min(fits, key=lambda fit: fit.cost), residuals

    best, residuals = solve(x, y, [[floor, A_guess, 0.0, 1.0, q0] + ([0.0] if slope else []) for q0 in (-3.0, -1.0, -0.3, 0.3, 1.0, 3.0)])
    if x_all is not None:
        x, y = x_all, y_all
        best, residuals = solve(x, y, [best.x])
    p = best.x
    sigma2 = float(np.sum(residuals(p) ** 2)) / max(1, len(y) - len(p))
    err = np.sqrt(np.clip(np.diag(np.linalg.pinv(best.jac.T @ best.jac) * sigma2), 0, None))
    x0, gamma, s = x0_guess + scale * p[2], scale * p[3], (p[5] / scale) if slope else 0.0
    result = {"y0": p[0], "A": p[1], "x0": x0, "gamma": gamma, "q": p[4], "s": s,
              "err_y0": err[0], "err_A": err[1], "err_x0": scale * err[2], "err_gamma": scale * err[3], "err_q": err[4],
              "rms": float(np.sqrt(np.mean(residuals(p) ** 2))), "Q": abs(x0) / gamma,
              "model": lambda xx: fano(xx, p[0], p[1], x0, gamma, p[4], s)}
    result["err_Q"] = result["Q"] * np.hypot(result["err_x0"] / x0 if x0 else 0.0, result["err_gamma"] / gamma)
    return result


def fano_fit(path):
    """Fano fit of the resonances of one CSV (CH1, CH2) inside the CH1 ramp, with the x axis in wavelength (nm).
    The resonances are located on the moving average; each one is fitted in a window around the located peak, on the raw data (or on the moving average if FANO_ON_AVERAGE).
    Prints the parameters, plots data + fit + residuals, and returns the list of results (one dict per resonance)."""
    path = run_path(path)
    ch1, ch2 = read_csv(path)
    f = fit_ramp(ch1)
    zone = valid_zone(f, len(ch1))
    if zone is None:
        print("no ramp in CH1: nothing to fit")
        return []
    a, b = zone
    lam = wavelength(np.arange(len(ch2)), f)
    smooth = moving_average(ch2)
    _, peaks, _ = analyze_point(ch1, -ch2 if FANO_DIP else ch2, -smooth if FANO_DIP else smooth)   # only to locate the resonances
    if not peaks:
        print("no resonance located on the moving average: lower PROMINENCE_THRESHOLD or MIN_SNR, or check FANO_DIP")
        return []
    peaks.sort(key=lambda p: -p["height"])                                   # tallest first
    chosen = peaks[:FANO_PEAKS]
    chosen.sort(key=lambda p: p["index"])                                    # in order of position
    print(f"{path.name}: {len(peaks)} resonances located, fitting {len(chosen)}")
    fits = []
    for n, p in enumerate(chosen, start=1):
        if FANO_RANGE is not None:
            inside = np.flatnonzero((lam >= min(FANO_RANGE)) & (lam <= max(FANO_RANGE)) & (np.arange(len(lam)) >= a) & (np.arange(len(lam)) <= b))
            if inside.size == 0:
                print(f"resonance {n}: FANO_RANGE is outside the ramp")
                continue
            ia, ib = inside[0], inside[-1]
        else:
            half = FANO_WINDOW * p["width"]
            ia, ib = max(a, int(p["index"] - half)), min(b, int(p["index"] + half))
        xw, yw = lam[ia:ib + 1], (smooth if FANO_ON_AVERAGE else ch2)[ia:ib + 1]
        xr, yr = block_median(xw, yw, FANO_POINTS)
        fit = fano_fit_peak(xr, yr, lam[p["index"]], p["width_nm"], *((xw, yw) if FANO_ALL_SAMPLES else (None, None)))
        fit.update(n=n, xw=xw, yw=yw, xr=xw if FANO_ALL_SAMPLES else xr, yr=yw if FANO_ALL_SAMPLES else yr)
        fits.append(fit)

    print(f"{'n':>2} {'x0 (nm)':>12} {'Gamma (nm)':>12} {'q':>9} {'A':>11} {'y0':>11} {'Q = x0/Gamma':>14} {'rms':>10}")
    for r in fits:
        print(f"{r['n']:>2} {r['x0']:>12.5f} {r['gamma']:>12.5g} {r['q']:>9.3f} {r['A']:>11.4g} {r['y0']:>11.4g} {r['Q']:>14.5g} {r['rms']:>10.3g}")
    for r in fits:
        print(f"resonance {r['n']}: x0 = {r['x0']:.5f} +- {r['err_x0']:.2g} nm | Gamma = {r['gamma']:.5g} +- {r['err_gamma']:.2g} nm | "
              f"q = {r['q']:.3f} +- {r['err_q']:.2g} | A = {r['A']:.4g} +- {r['err_A']:.2g} | y0 = {r['y0']:.4g} +- {r['err_y0']:.2g} | "
              f"Q = {r['Q']:.5g} +- {r['err_Q']:.2g}")
    if not fits:
        return fits

    step = max(1, len(ch2) // 20000)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(lam[::step], ch2[::step], color="gray", lw=0.6, alpha=0.7, label="raw CH2")
    ax.plot(lam[::step], smooth[::step], color="red", lw=1.0, alpha=0.7, label=f"moving average (N={WINDOW})")
    ax.axvspan(lam[a], lam[b], color="orange", alpha=0.1)
    for r in fits:
        ax.plot(r["xw"], r["model"](r["xw"]), color="blue", lw=2)
        ax.annotate(str(r["n"]), (r["x0"], r["model"](np.array([r["x0"]]))[0]), textcoords="offset points", xytext=(0, 10), ha="center")
    ax.set_xlim(LAMBDA_0 - 0.3 * (LAMBDA_1 - LAMBDA_0), LAMBDA_1 + 0.3 * (LAMBDA_1 - LAMBDA_0))
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("CH2 (V)")
    ax.legend(loc="best")
    plt.tight_layout()

    fig, axes = plt.subplots(2, len(fits), figsize=(5.2 * len(fits), 5.5), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1]}, squeeze=False)
    for k, r in enumerate(fits):
        a1, a2 = axes[0, k], axes[1, k]
        a1.plot(r["xw"][::max(1, len(r["xw"]) // 20000)], r["yw"][::max(1, len(r["xw"]) // 20000)], color="lightgray", lw=0.5, label="moving average" if FANO_ON_AVERAGE else "raw data")
        a1.plot(r["xr"][::max(1, len(r["xr"]) // 5000)], r["yr"][::max(1, len(r["xr"]) // 5000)], "k.", ms=2, label="fitted data")
        xx = np.linspace(r["xw"][0], r["xw"][-1], 2000)
        a1.plot(xx, r["model"](xx), color="blue", lw=2, label="Fano fit")
        a1.set_title(f"resonance {r['n']}: x0={r['x0']:.5g} nm, Gamma={r['gamma']:.4g} nm, q={r['q']:.2f}", fontsize=10)
        a1.set_ylabel("CH2 (V)")
        a1.legend(loc="best", fontsize=8)
        sub = slice(None, None, max(1, len(r["xr"]) // 5000))
        a2.plot(r["xr"][sub], (r["yr"] - r["model"](r["xr"]))[sub], "k.", ms=2)
        a2.axhline(0, color="blue", lw=1)
        a2.set_ylabel("residuals")
        a2.set_xlabel("wavelength (nm)")
    plt.tight_layout()
    plt.show()
    return fits


# ---------------------------------------------------------------- single measurements

def measurements():
    """Registry of the single measurements: misure_singole/misure.csv (columns id, ora, nome, x_mm, y_mm, con_scansione, file, note)."""
    return pd.read_csv(DATA / "misure_singole" / "misure.csv", keep_default_na=False)


def spectrum(ch1, ch2, n_points=5000):
    """CH2 inside the CH1 ramp as (wavelength in nm, block median of CH2)."""
    f = fit_ramp(ch1)
    i0, i1 = int(np.ceil(f[0])), int(np.floor(f[1]))
    n = max(1, (i1 - i0 + 1) // n_points)
    k = (i1 - i0 + 1) // n
    return wavelength(i0 + (np.arange(k) + 0.5) * n - 0.5, f), np.median(ch2[i0:i0 + k * n].reshape(k, n), axis=1) - DARK_LEVEL


def compare(ids=None, ref=None):
    """Compare single measurements (same kind as the reference) and divide them by measurement ref (default: the first)."""
    m = measurements()
    m = m[m["id"].isin(ids)] if ids else m
    ref = ref or int(m["id"].iloc[0])
    scanned = str(m[m["id"] == ref]["con_scansione"].iloc[0]) == "True"
    m = m[m["con_scansione"].astype(str) == str(scanned)]
    data = {r.id: read_csv(DATA / "misure_singole" / r.file) for r in m.itertuples()}
    labels = {r.id: f"{r.id} {r.nome}" for r in m.itertuples()}
    if scanned:
        curves = {i: spectrum(*d) for i, d in data.items()}
        lam = np.linspace(max(c[0][0] for c in curves.values()), min(c[0][-1] for c in curves.values()), 2000)
        s = {i: np.interp(lam, *c) for i, c in curves.items()}
        for i in s:
            plt.plot(lam, s[i], lw=2.2 if i == ref else 1.2, label=labels[i])
        plt.xlabel("wavelength (nm)")
        plt.ylabel("CH2 (V)")
        plt.legend()
        plt.tight_layout()
        values = {i: s[i].mean() for i in s}
        ratios = {i: (s[i] / s[ref]).mean() for i in s}
        errors = {i: np.nan for i in s}
    else:
        level, err = {}, {}
        for i, (_, ch2) in data.items():
            blocks = np.median(ch2[:len(ch2) // 1000 * 1000].reshape(-1, 1000), axis=1)
            level[i], err[i] = np.median(ch2) - DARK_LEVEL, blocks.std(ddof=1) / np.sqrt(len(blocks))
        values = level
        ratios = {i: level[i] / level[ref] for i in level}
        errors = {i: abs(ratios[i]) * np.hypot(err[i] / level[i], err[ref] / level[ref]) for i in level}
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
        x = np.arange(len(level))
        a1.bar(x, [level[i] for i in level], yerr=[err[i] for i in level], capsize=3)
        a2.bar(x, [ratios[i] for i in level], yerr=[errors[i] for i in level], capsize=3)
        for a, t in ((a1, "CH2 level (V)"), (a2, f"ratio to measurement {ref}")):
            a.set_xticks(x)
            a.set_xticklabels([labels[i] for i in level], rotation=30, ha="right")
            a.set_ylabel(t)
        a2.axhline(1, color="gray", lw=0.8, ls=":")
        plt.tight_layout()
    table = pd.DataFrame({"id": list(values), "name": [labels[i] for i in values], "value_V": list(values.values()),
                          f"ratio_to_{ref}": list(ratios.values()), "ratio_error": list(errors.values()),
                          "notes": [m[m["id"] == i]["note"].iloc[0] for i in values]})
    print(table.to_string(index=False))
    plt.show()
    return table


# ---------------------------------------------------------------- from the terminal

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command")
    a = sub.add_parser("run", help="modulation, Q and heatmaps of a run")
    a.add_argument("folder")
    a.add_argument("--recompute", action="store_true", help="ignore analysis.csv and read all the CSV files again")
    a = sub.add_parser("point", help="check plot of one or more points of the run (no heatmaps)")
    a.add_argument("folder")
    a.add_argument("indices", nargs="+", type=int, help="point indices, e.g. 3 12 40")
    a = sub.add_parser("view", help="plot one or more single acquisitions (CSV files), no heatmaps")
    a.add_argument("csv", nargs="+")
    a = sub.add_parser("fano", help="Fano fit of the resonances of one CSV")
    a.add_argument("csv")
    sub.add_parser("list", help="list of the single measurements")
    a = sub.add_parser("compare", help="compare single measurements")
    a.add_argument("ids", nargs="*", type=int, help="measurement ids (empty = all)")
    a.add_argument("--ref", type=int, help="id of the reference measurement for the ratios")
    args = p.parse_args()
    if args.command is None:                                         # no arguments: use WHAT TO LOOK AT
        if not (RUN and (HEATMAPS or POINT is not None)) and VIEW is None and COMPARE is None and not FANO:
            print("Nothing to do: set RUN (with HEATMAPS and/or POINT), VIEW, COMPARE and/or FANO in WHAT TO LOOK AT at the top of analysis.py, "
                  "or use a command (python analysis.py -h)")
        if RUN and HEATMAPS:
            run_heatmaps(RUN, RECOMPUTE)
        if RUN and POINT is not None:
            for i in ([POINT] if isinstance(POINT, int) else POINT):
                point_plot(RUN, i, show=False)
            plt.show()                                                   # all the points together
        if VIEW is not None:
            for path in ([VIEW] if isinstance(VIEW, (str, Path)) else VIEW):
                csv_plot(path, show=False)
            plt.show()
        if COMPARE is not None:
            compare(COMPARE, REFERENCE)
        if FANO:
            fano_fit(FANO)
    elif args.command == "view":
        for path in args.csv:
            csv_plot(path, show=False)
        plt.show()
    elif args.command == "fano":
        fano_fit(args.csv)
    elif args.command == "run":
        run_heatmaps(args.folder, args.recompute)
    elif args.command == "point":
        for i in args.indices:
            point_plot(args.folder, i, show=False)
        plt.show()
    elif args.command == "list":
        print(measurements().to_string(index=False))
    else:
        compare(args.ids, args.ref)


if __name__ == "__main__":
    main()
