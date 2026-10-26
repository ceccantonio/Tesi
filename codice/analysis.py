"""Data analysis: modulation amplitude, quality factor, heatmaps and comparison of single measurements.

    python analysis.py run 2026-10-06/run_007       # modulation and Q of every point -> analysis.csv + two heatmaps
    python analysis.py mosaic 2026-10-06/run_007    # ONE png with the spectra of all the membranes, each in the position of its membrane
    python analysis.py point 2026-10-06/run_007 3 12 40   # detail plots of some points of a run (no heatmaps)
    python analysis.py view misure_singole/m004_x.csv     # plot of single acquisitions (CSV files), no run folder needed
    python analysis.py fano 2026-10-06/run_007/point_011_x91.818_y146.526.csv   # Fano fit of the resonances of one CSV
    python analysis.py list                         # list of the single measurements
    python analysis.py compare 1 2 3                # compare single measurements, in volts (add --ratio --ref 1 to divide them by measurement 1)
    python analysis.py compare 2026-10-06/run_007:12 reference/copol.csv   # a point of a run and the reference, from different folders
    python analysis.py                              # no arguments: does what is written in "WHAT TO LOOK AT" below

Method (as in movemin.py): moving average of CH2 (np.convolve) inside the CH1 ramp.
 - exact start and end of the ramp: fit of low plateau + linear rise + high plateau (fit_ramp), which also uses points
   before the start and after the end; the x axis is the wavelength: LAMBDA_0 at the start and LAMBDA_1 at the end of the
   ramp, by linear proportion (assumption to verify: the scan goes from 1520 to 1570 nm over the whole ramp);
 - modulation = max - min of the moving average inside the ramp;
 - Q of each peak = height (prominence) / width at half height, in volts per sample;
   each peak also gets its wavelength and its width in nm.

Fano fit: y = y0 + A (q + eps)^2 / (1 + eps^2), eps = 2 (x - x0) / Gamma, with x in nm (scipy curve_fit, A >= 0 so that q is
unique, Q = x0 / Gamma). The resonances are located on the moving average; each one is fitted in a window around it, on the raw
CH2 data or on the moving average (FANO_ON_AVERAGE).

Run folders contain griglia.csv (columns indice, x_mm, y_mm, file) and one CSV (columns CH1, CH2) per point:
that format is the one written by the earlier notebooks and by acquisition.py.
"""
import argparse
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit, least_squares
from scipy.signal import fftconvolve, find_peaks, peak_widths

DATA = Path(__file__).resolve().parent.parent / "Chiara_Bru" / "runs_oscilloscope"

# ---------------------------------------------------------------- WHAT TO LOOK AT
# Used when you run `python analysis.py` with no arguments (or call main() from your own code / a notebook).
RUN = ""                    # run folder (e.g. "2026-10-06/run_007"); "" = skip
HEATMAPS = True             # True: draw the heatmaps of RUN (modulation and Q); False: skip them (e.g. to look only at some points)
POINT = None                # point(s) of RUN to plot in detail: an index (12) or a list ([3, 12, 40]); None = skip
VIEW = None                 # CSV file(s) of single acquisitions to plot, no run folder needed: a path or a list of paths; None = skip
COMPARE = None              # what to compare (in volts), from ANY folder, mixed as you like: ids of the single measurements (1), CSV files ("2026-10-06/run_007/point_012.csv", "reference/copol.csv"), points of a run ("2026-10-06/run_007:12"); [] = all the single measurements; None = skip. Or a dict of lists for one panel per group, one above the other: {"co": [...], "cross": [...]}
COMPARE_RATIO = False       # False: only the spectra and their mean values; True: also divide each one by REFERENCE
REFERENCE = None            # used only with COMPARE_RATIO = True: the one the others are divided by, written in the same way (None = the first of each group), or a dict like COMPARE, e.g. the point without the membrane
RECOMPUTE = False           # True: ignore analysis.csv and read all the CSV files again
SIMULATION = ""             # simulation file (CSV or text: wavelength or frequency in one column, the spectrum in another): drawn above the measurement in the point plots
SIM_COLUMNS = (0, 1)        # (column of the wavelength/frequency, column of the simulated spectrum) in the SIMULATION file, counting from 0
MOSAIC = False              # True: also draw ONE png with the spectra (1520-1570 nm) of all the membranes of RUN, each one in the position of its membrane
POLARIZATION = ""           # "co" or "cross": every measured spectrum is divided by the transmission reference of that polarization (reference/copol.csv or reference/crosspol.csv, measured once), so you get T like in the simulation; "" = no division
FANO = ""                   # CSV file to fit with a Fano line shape (e.g. "2026-10-06/run_007/point_011_x91.818_y146.526.csv"); "" = skip

# ---------------------------------------------------------------- parameters of the method
WINDOW = 10000              # moving-average length, in samples
FAST_MOVING_AVERAGE = False # True: same moving average computed with an FFT (identical result, much faster on long records); False: np.convolve as in movemin.py
V_START = 0.9               # the CH1 ramp starts when CH1 exceeds this value...
END_TOLERANCE = 0.3         # ...and ends when CH1 is within this of its maximum (starting guess for the ramp fit)
RAMP_MARGIN = 0.25          # the ramp fit also uses points before and after: this fraction of the duration on each side
LAMBDA_0, LAMBDA_1 = 1520.0, 1570.0   # nm at the start and at the end of the CH1 ramp, i.e. of the laser scan
REFERENCE_DIR = DATA / "reference"    # the transmission references (scan without the sample), measured once: copol.csv and crosspol.csv
PROMINENCE_THRESHOLD = 0.1  # a peak counts if taller than this fraction of the excursion (modulation)...
MIN_SNR = 10                # ...and than this multiple of the residual noise of the moving average...
RUN_THRESHOLD = 0.1         # ...and than this fraction of the tallest peak of the whole run (no resonance = n.a.)
FANO_PEAKS = 1              # Fano fit: how many resonances to fit (the tallest ones); None = all the located ones
FANO_WINDOW = 3.0           # each fit uses the data within +- this many widths (FWHM) around the resonance
FANO_ON_AVERAGE = False     # False: fit the raw CH2 data; True: fit the moving average of CH2
Q_UNIT = 1e-6               # Q is shown in units of this (V per sample)
DARK_LEVEL = 0.0            # CH2 level without light (single measurements)


# ---------------------------------------------------------------- reading and basic method

def run_path(folder):
    p = Path(folder)
    return p if p.is_absolute() or p.exists() else DATA / p


def read_table(path):
    """Small CSV (a registry) as a table. Ignored: the lines starting with # (e.g. a description written at the top) and
    whatever follows a # in the header line (e.g. "indice,x_mm,y_mm,file  #my description")."""
    with open(path, encoding="utf-8-sig") as f:
        lines = [line for line in f if not line.startswith("#")]
    if lines:
        lines[0] = lines[0].split("#")[0].rstrip() + "\n"
    return pd.read_csv(io.StringIO("".join(lines)))


def description(folder):
    """Text of description.txt in the run folder (what the run is; written by the acquisition scripts, never parsed), or ""."""
    path = Path(folder) / "description.txt"
    return path.read_text(encoding="utf-8").strip() if path.exists() else ""


def read_registry(folder):
    """griglia.csv of a run as a table with the columns indice, x_mm, y_mm, file. Tolerant to a BOM, spaces, capital letters and to
    another name of the file column (it is then the column whose values end with .csv)."""
    table = read_table(folder / "griglia.csv")
    table.columns = [c.strip().lower() for c in table.columns]
    if "file" not in table.columns:
        names = [c for c in table.columns if table[c].astype(str).str.strip().str.endswith(".csv").all()]
        if names:
            table = table.rename(columns={names[0]: "file"})
        else:
            print(f"griglia.csv has the columns {list(table.columns)}: no column with the CSV file names")
    return table


def read_csv(path):
    df = pd.read_csv(path, dtype=np.float32, comment="#")                  # what follows a # (a description) is ignored
    df.columns = df.columns.str.strip()                                    # "CH2  #description" would leave spaces in the name
    return df["CH1"].to_numpy(float), df["CH2"].to_numpy(float)


def suffix():
    """"" or "_T_<reference name>": added to the names of the results, so that the ones computed with a transmission reference do not mix with the others."""
    return f"_T_{POLARIZATION}pol" if POLARIZATION else ""


def unit():
    return "T" if POLARIZATION else "V"


def signal_name():
    return "transmission T" if POLARIZATION else "CH2 (V)"


_reference = {}


def reference_spectrum():
    """(wavelength, spectrum) of the transmission reference of POLARIZATION (REFERENCE_DIR/copol.csv or crosspol.csv), or (None, level) if
    it has no scan; None if POLARIZATION is not set. The big CSV is read once: its spectrum is saved next to it (copol.npz, crosspol.npz)."""
    if not POLARIZATION:
        return None
    path = REFERENCE_DIR / f"{POLARIZATION}pol.csv"
    if path not in _reference:
        cache = path.with_suffix(".npz")
        if cache.exists() and cache.stat().st_mtime >= path.stat().st_mtime:
            saved = np.load(cache)
            _reference[path] = (saved["lam"] if len(saved["lam"]) else None, saved["s"])
        else:
            ch1, ch2 = read_csv(path)
            _reference[path] = spectrum(ch1, ch2) if fit_ramp(ch1) is not None else (None, np.median(ch2) - DARK_LEVEL)
            np.savez(cache, lam=np.array([] if _reference[path][0] is None else _reference[path][0]), s=_reference[path][1])
    return _reference[path]


def load_point(path):
    """CH1 and CH2 of one acquisition. With POLARIZATION, CH2 is divided by the reference measurement at the same wavelength: T."""
    ch1, ch2 = read_csv(path)
    ref = reference_spectrum()
    f = fit_ramp(ch1) if ref is not None else None
    if ref is None or (ref[0] is not None and f is None):
        return ch1, ch2
    if ref[0] is None:
        den = ref[1]
    else:
        den = np.maximum(np.interp(wavelength(np.arange(len(ch2)), f), *ref), 0.02 * ref[1].max())   # where the reference is dark T is not defined
    return ch1, (ch2 - DARK_LEVEL) / den


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
    cache = folder / f"analysis{suffix()}.csv"
    if cache.exists() and not recompute:
        try:
            saved = pd.read_csv(cache)
        except pd.errors.EmptyDataError:
            saved = pd.DataFrame()
        if "x_mm" in saved.columns and len(saved):                  # an empty or broken analysis.csv is ignored
            return saved
    points, missing, spectra = [], [], {}
    registry = read_registry(folder)
    print(f"{folder}: {len(registry)} points in griglia.csv")
    if description(folder):
        print("description:", description(folder))
    for r in registry.itertuples():
        if (folder / r.file).exists():
            ch1, ch2 = load_point(folder / r.file)
            modulation, peaks, f = analyze_point(ch1, ch2)
            print(f"point {r.indice}: modulation {modulation:.4f} {unit()}, {len(peaks)} peaks")
            points.append((r, modulation, peaks, f))
            spectra[r.indice] = point_spectrum(ch1, ch2)                       # for the mosaic, in the same pass (the files are big)
        else:
            missing.append(r.file)
    if missing:
        present = sorted(q.name for q in folder.glob("*.csv") if q.name not in ("griglia.csv", "analysis.csv"))[:3]
        print(f"{len(missing)} files listed in griglia.csv are not in the folder.")
        print(f"  first one looked for: {folder / missing[0]}")
        print(f"  CSV files actually in the folder: {present if present else 'none'}")
    if not points:
        print("no data file found: nothing to analyze (check the run folder and griglia.csv)")
        return pd.DataFrame(columns=["index", "x_mm", "y_mm", "modulation_V", "Q"])
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
    table.to_csv(cache, index=False)
    np.savez(folder / f"spectra{suffix()}.npz", index=list(spectra), spectra=np.array(list(spectra.values())))
    return table


def lattice(x, y):
    """(row, column) of every point of a serpentine path; None if it cannot be reconstructed. Works for any tilt of the sample:
    inside a row every step is +u or -u (u = the first step), a step that is neither of them is the turn to the next row."""
    steps = np.diff(np.column_stack([x, y]), axis=0)
    if len(x) < 4:
        return None
    u = steps[0]
    along = (np.linalg.norm(steps - u, axis=1) < 0.3 * np.linalg.norm(u)) | (np.linalg.norm(steps + u, axis=1) < 0.3 * np.linalg.norm(u))
    rows = np.concatenate([[0], np.cumsum(~along)])
    lengths = np.bincount(rows)
    if len(set(lengths)) != 1:
        return None
    columns = np.concatenate([np.arange(n) if r % 2 == 0 else np.arange(n)[::-1] for r, n in enumerate(lengths)])
    return rows, columns


_checked = set()


def layout(x, y):
    """Where each point goes in a grid drawn like the real sample: (row, column, n_rows, n_columns), or None if the points are not a
    regular serpentine. It also checks that the grid agrees with the real X, Y positions (position = origin + column * u + row * v,
    fitted on all the points) and prints the check once. The grid is oriented so that columns grow towards +X and rows towards -Y (down)."""
    rc = lattice(x, y)
    key = (len(x), round(float(np.sum(x) + np.sum(y)), 4))
    if rc is None:
        if key not in _checked:
            print("lattice check: the points are NOT a regular serpentine (missing or out-of-order points): drawn at the real X, Y positions")
        _checked.add(key)
        return None
    rows, columns = rc
    design = np.column_stack([np.ones(len(x)), columns, rows])
    coef = np.linalg.lstsq(design, np.column_stack([x, y]), rcond=None)[0]
    u, v = coef[1], coef[2]
    error = np.hypot(*(np.column_stack([x, y]) - design @ coef).T).max()
    flip_columns, flip_rows = u[0] < 0, v[1] > 0
    if key not in _checked:
        print(f"lattice check: {columns.max() + 1} columns x {rows.max() + 1} rows, pitch {np.linalg.norm(u):.3f} x {np.linalg.norm(v):.3f} mm, "
              f"tilt {np.degrees(np.arctan2(u[1], u[0])):+.2f} deg, largest distance of a point from the ideal lattice {error * 1000:.1f} um")
        print(f"  the column index grows towards {'-X' if flip_columns else '+X'} and the row index towards {'+Y' if flip_rows else '-Y'} in the real positions: "
              f"grid drawn with columns towards +X (right) and rows towards -Y (down)" + (" after flipping" if flip_columns or flip_rows else ""))
        if error > 0.2 * min(np.linalg.norm(u), np.linalg.norm(v)):
            print("  WARNING: some points are far from the ideal lattice: the picture may not match the real positions")
        else:
            print("  OK: the picture matches the real positions")
    _checked.add(key)
    n_rows, n_columns = rows.max() + 1, columns.max() + 1
    return (n_rows - 1 - rows if flip_rows else rows), (n_columns - 1 - columns if flip_columns else columns), n_rows, n_columns


def heatmap(table, column, label, file=None, scale=1.0, fmt="{:.2f}", subtitle=""):
    """Heatmap with the value written on each point: on a lattice (rows/columns of the serpentine) or at the real X, Y positions."""
    if column not in table.columns or table[column].isna().all():
        print(f"{column}: no valid value, heatmap not drawn")
        return
    x, y = table["x_mm"].to_numpy(), table["y_mm"].to_numpy()
    z = table[column].to_numpy() / scale
    valid = np.isfinite(z)
    if not valid.any():
        print(f"{column}: no valid value, heatmap not drawn")
        return
    norm, cmap = Normalize(z[valid].min(), z[valid].max()), plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(9, 7))
    lay = layout(x, y)
    if lay is not None:
        rows, columns, n_rows, n_columns = lay
        grid = np.full((n_rows, n_columns), np.nan)
        grid[rows, columns] = z
        ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, norm=norm, origin="upper")
        cx, cy = columns, rows
        ax.set_xlabel("column (towards +X)")
        ax.set_ylabel("row (towards -Y, down)")
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
    ax.set_title(label + ("\n" + subtitle if subtitle else ""))
    plt.tight_layout()
    if file:
        fig.savefig(file, dpi=200)
        print("Saved", file)
    plt.show()


def run_heatmaps(folder, recompute=False):
    """Heatmaps of the modulation and of the quality factor of a run."""
    folder = run_path(folder)
    table = analyze_run(folder, recompute)
    heatmap(table, "modulation_V", f"modulation amplitude: max - min of the moving average in the ramp ({unit()})",
            folder / f"heatmap_modulation{suffix()}.png", fmt="{:.3f}", subtitle=description(folder))
    heatmap(table, "Q", f"quality factor Q = height / width  [{Q_UNIT:g} {unit()} per sample]",
            folder / f"heatmap_quality_factor{suffix()}.png", scale=Q_UNIT, subtitle=description(folder))


LAM_GRID = np.linspace(LAMBDA_0, LAMBDA_1, 501)


def point_spectrum(ch1, ch2):
    """CH2 (or T) of one acquisition on LAM_GRID, from the ramp fit (NaN if there is no ramp)."""
    if fit_ramp(ch1) is None:
        return np.full(len(LAM_GRID), np.nan)
    return np.interp(LAM_GRID, *spectrum(ch1, ch2))


def run_spectra(folder, recompute=False):
    """(indices, spectra on LAM_GRID) of the points of a run, cached in spectra.npz (analyze_run writes it in the same pass)."""
    folder = run_path(folder)
    cache = folder / f"spectra{suffix()}.npz"
    if cache.exists() and not recompute:
        saved = np.load(cache)
        return saved["index"], saved["spectra"]
    index, spectra = [], []
    for r in read_registry(folder).itertuples():
        if (folder / r.file).exists():
            index.append(r.indice)
            spectra.append(point_spectrum(*load_point(folder / r.file)))
            print(f"spectrum of point {r.indice}")
    np.savez(cache, index=index, spectra=np.array(spectra))
    return np.array(index), np.array(spectra)


def mosaic(folder, recompute=False):
    """ONE png with the spectrum (LAMBDA_0 to LAMBDA_1) of every membrane, each one in the position of its membrane: the top left of the
    png is the top left of the sample (same orientation as the heatmaps). Same vertical scale everywhere."""
    folder = run_path(folder)
    index, spectra = run_spectra(folder, recompute)
    registry = read_registry(folder).set_index("indice").loc[index]
    lay = layout(registry["x_mm"].to_numpy(), registry["y_mm"].to_numpy())
    if lay is None or len(index) == 0:
        print("mosaic: the points are not a regular lattice (or there is no data): the spectra cannot be placed")
        return
    rows, columns, n_rows, n_columns = lay
    low, high = np.nanpercentile(spectra, [1, 99])
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(2.0 * n_columns + 0.8, 1.7 * n_rows + 1.3), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for k, i in enumerate(index):
        ax = axes[rows[k], columns[k]]
        ax.axis("on")
        ax.plot(LAM_GRID, spectra[k], lw=0.7)
        ax.set_xlim(LAMBDA_0, LAMBDA_1)
        ax.set_ylim(low - 0.05 * (high - low), high + 0.05 * (high - low))
        ax.set_xticks([LAMBDA_0, (LAMBDA_0 + LAMBDA_1) / 2, LAMBDA_1])             # axes with numbers on every panel
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.tick_params(labelsize=7, length=3, width=0.8, pad=2)
        labels = ax.get_xticklabels()                                              # the end labels stay inside their own panel
        labels[0].set_ha("left")
        labels[-1].set_ha("right")
        ax.grid(alpha=0.3, lw=0.3)
        ax.text(0.04, 0.93, str(i), transform=ax.transAxes, fontsize=8, va="top", fontweight="bold")
    fig.supxlabel("wavelength (nm)", fontsize=11)
    fig.supylabel(signal_name(), fontsize=11)
    fig.suptitle(f"{folder.name}: {signal_name()} from {LAMBDA_0:g} to {LAMBDA_1:g} nm, each panel in the position of its membrane, same vertical scale"
                 + (f"\n{description(folder)}" if description(folder) else ""), fontsize=9)
    plt.tight_layout(rect=(0.01, 0.01, 1, 0.96), h_pad=1.5, w_pad=1.5)
    fig.savefig(folder / f"mosaic_spectra{suffix()}.png", dpi=150)
    print("Saved", folder / f"mosaic_spectra{suffix()}.png")
    plt.show()


def point_plot(folder, index, show=True):
    """Detail plot of one point of a run (see csv_plot). show=False: draw it without opening it yet."""
    folder = run_path(folder)
    r = read_registry(folder).query("indice == @index").iloc[0]
    csv_plot(folder / r.file, f"Point {index} (X={r.x_mm}, Y={r.y_mm})", show)


def read_simulation(path):
    """Simulation file (CSV or text, e.g. a COMSOL export; the lines starting with # or % are comments) as (wavelength in nm, values),
    sorted by wavelength. SIM_COLUMNS = (column of x, column of the spectrum). The unit of x is guessed from its size:
    below 5 it is um, below 1000 it is THz (converted with c / f), otherwise nm."""
    path = Path(path) if Path(path).exists() else DATA / path
    with open(path, encoding="utf-8-sig") as f:
        text = "\n".join(line.strip() for line in f if line.strip() and not line.lstrip().startswith(("#", "%")))
    table = pd.read_csv(io.StringIO(text), sep=r"[,;\s]+", header=None, engine="python").apply(pd.to_numeric, errors="coerce")
    x, y = table.iloc[:, SIM_COLUMNS[0]], table.iloc[:, SIM_COLUMNS[1]]       # a text header row becomes NaN and is dropped
    ok = (x.notna() & y.notna()).to_numpy()
    x, y = x.to_numpy()[ok], y.to_numpy()[ok]
    x = x * 1000 if np.median(x) < 5 else 299792.458 / x if np.median(x) < 1000 else x
    order = np.argsort(x)
    return x[order], y[order]


def draw_simulation(ax, xmin=LAMBDA_0, xmax=LAMBDA_1):
    """Draws the SIMULATION spectrum between xmin and xmax (nm) in the panel ax."""
    x, y = read_simulation(SIMULATION)
    inside = (x >= xmin) & (x <= xmax)
    ax.plot(x[inside], y[inside], color="tab:green", lw=1.5, label=f"simulation ({Path(SIMULATION).name})")
    ax.set_ylabel("simulation")
    ax.legend(loc="best", fontsize=8)


def csv_plot(path, label=None, show=True):
    """Plot of ONE acquisition (any CSV with CH1, CH2: a point of a run, a single measurement...): raw CH2, moving average,
    measurement zone and peaks, in wavelength, only between LAMBDA_0 and LAMBDA_1. If SIMULATION is set, its spectrum is drawn
    in a panel above the measurement, with the same wavelength axis.
    show=False: draw the figure without opening it yet (to open several together with one plt.show())."""
    path = run_path(path)
    ch1, ch2 = load_point(path)
    smooth = moving_average(ch2)
    modulation, peaks, f = analyze_point(ch1, ch2, smooth)
    if f is None:
        print("no ramp in CH1: nothing to plot")
        return
    lam = wavelength(np.arange(len(ch2)), f)
    inside = np.flatnonzero((lam >= LAMBDA_0) & (lam <= LAMBDA_1))[::max(1, len(ch2) // 20000)]       # only the points of the scan
    fig, axes = plt.subplots(2 if SIMULATION else 1, 1, figsize=(10, 7 if SIMULATION else 4.5), sharex=True, squeeze=False)
    axes = axes[:, 0]
    a = axes[-1]
    a.plot(lam[inside], ch2[inside], color="gray", lw=0.6, alpha=0.6, label="transmission T" if POLARIZATION else "raw CH2")
    a.plot(lam[inside], smooth[inside], color="red", lw=1.5, label=f"moving average (N={WINDOW})")
    zone = valid_zone(f, len(ch1))
    if zone:
        a.axvspan(*wavelength(zone, f), color="orange", alpha=0.12, label="measurement zone")
    for n, p in enumerate(peaks, start=1):
        a.hlines(p["level"], *wavelength([p["left"], p["right"]], f), color="blue", lw=2)
        a.annotate(f"{n}: {p['lambda_nm']:.2f} nm, Q={p['Q'] / Q_UNIT:.3g}", (p["lambda_nm"], smooth[p["index"]]),
                   textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9)
    a.set_xlim(LAMBDA_0, LAMBDA_1)
    a.set_xlabel("wavelength (nm)")
    a.set_ylabel(signal_name())
    a.legend(loc="best", fontsize=8)
    if SIMULATION:
        draw_simulation(axes[0])
    axes[0].set_title(f"{label or path.name}: modulation {modulation:.4f} {unit()}, {len(peaks)} peaks")
    plt.tight_layout()
    if show:
        plt.show()


# ---------------------------------------------------------------- Fano fit

def fano(x, y0, A, x0, gamma, q):
    """Fano line shape: y0 + A (q + eps)^2 / (1 + eps^2), eps = 2 (x - x0) / gamma."""
    eps = 2.0 * (x - x0) / gamma
    return y0 + A * (q + eps) ** 2 / (1.0 + eps ** 2)


def fit_fano(x, y, x0, gamma):
    """curve_fit of the Fano shape on (x, y), starting from a resonance at x0 with width gamma. A >= 0 is imposed ((A, q) and
    (-A q^2, -1/q) are the same curve, so without it q would not be unique). It starts from several values of q and keeps the
    best one, because the Fano shape has several local minima. Returns (parameters, errors) or None if it never converges."""
    lo, hi = np.percentile(y, [5, 95])
    edge = np.median(np.r_[y[:len(y) // 5], y[-len(y) // 5:]])
    bounds = ([-np.inf, 0, x0 - 2 * gamma, gamma / 20, -50], [np.inf, 20 * (hi - lo), x0 + 2 * gamma, 20 * gamma, 50])
    best = None
    for q0 in (-3, -1, -0.3, 0.3, 1, 3):
        try:
            p, cov = curve_fit(fano, x, y, p0=[lo, max(edge - lo, 0.1 * (hi - lo)), x0, gamma, q0], bounds=bounds,
                               loss="soft_l1", f_scale=0.02 * (hi - lo), x_scale="jac", max_nfev=3000)   # soft_l1: resists spikes
        except RuntimeError:
            continue
        error = np.abs(fano(x, *p) - y).sum()                       # the best start = the smallest total residual
        if best is None or error < best[0]:
            best = error, p, np.sqrt(np.diag(cov))
    return best[1:] if best else None


def fano_fit(path):
    """Fano fit of the resonances of one CSV (CH1, CH2), in wavelength (nm). The FANO_PEAKS tallest resonances found on the
    moving average are fitted one by one, each in a window of +- FANO_WINDOW widths around it, on the raw CH2 data or, with
    FANO_ON_AVERAGE, on the moving average. Prints the parameters and plots data, fit and residuals."""
    path = run_path(path)
    ch1, ch2 = load_point(path)
    smooth = moving_average(ch2)
    _, peaks, f = analyze_point(ch1, ch2, smooth)
    if not peaks:
        print("no resonance found on the moving average: lower PROMINENCE_THRESHOLD or MIN_SNR")
        return []
    zone_start, zone_end = valid_zone(f, len(ch1))
    data = smooth if FANO_ON_AVERAGE else ch2
    results = []
    for n, peak in enumerate(sorted(sorted(peaks, key=lambda p: -p["height"])[:FANO_PEAKS], key=lambda p: p["index"]), start=1):
        half = int(FANO_WINDOW * peak["width"])
        a, b = max(zone_start, peak["index"] - half), min(zone_end, peak["index"] + half)
        x, y = wavelength(np.arange(a, b + 1), f), data[a:b + 1]
        fit = fit_fano(x, y, peak["lambda_nm"], peak["width_nm"])
        if fit is None:
            print(f"resonance {n}: the fit did not converge")
            continue
        p, err = fit
        print(f"resonance {n}: x0 = {p[2]:.5f} +- {err[2]:.2g} nm | Gamma = {p[3]:.5g} +- {err[3]:.2g} nm | q = {p[4]:.3f} +- {err[4]:.2g} | "
              f"A = {p[1]:.4g} | y0 = {p[0]:.4g} | Q = x0/Gamma = {p[2] / p[3]:.5g}")
        results.append({"x0": p[2], "gamma": p[3], "q": p[4], "A": p[1], "y0": p[0], "Q": p[2] / p[3], "params": p, "errors": err, "x": x, "y": y})
    if results:
        lam_all = wavelength(np.arange(len(ch2)), f)
        inside = np.flatnonzero((lam_all >= LAMBDA_0) & (lam_all <= LAMBDA_1))[::max(1, len(ch2) // 20000)]   # only the points of the scan
        fig, axes = plt.subplots(2 if SIMULATION else 1, 1, figsize=(10, 6.5 if SIMULATION else 4), sharex=True, squeeze=False)
        axes = axes[:, 0]                                                         # 1) the whole scan, with the fit of every resonance
        ax = axes[-1]
        if SIMULATION:
            draw_simulation(axes[0])                                              # the simulation above
        ax.plot(lam_all[inside], ch2[inside], color="gray", lw=0.6, alpha=0.7, label="raw CH2")
        ax.plot(lam_all[inside], smooth[inside], color="red", lw=1.0, alpha=0.7, label=f"moving average (N={WINDOW})")
        ax.axvspan(*wavelength((zone_start, zone_end), f), color="orange", alpha=0.1, label="measurement zone")
        for k, r in enumerate(results):
            ax.plot(r["x"], fano(r["x"], *r["params"]), color="blue", lw=2, label="Fano fit" if k == 0 else None)
            ax.annotate(str(k + 1), (r["x0"], fano(r["x0"], *r["params"])), textcoords="offset points", xytext=(0, 10), ha="center")
        ax.set_xlim(LAMBDA_0, LAMBDA_1)
        ax.set_xlabel("wavelength (nm)")
        ax.set_ylabel(signal_name())
        ax.legend(loc="best")
        plt.tight_layout()
        plt.show()
        for k, r in enumerate(results):                                           # 2) one figure per resonance: data, fit and residuals
            n = max(1, len(r["x"]) // 5000)                                       # block median (as in the old notebook), only to draw fewer points
            m = len(r["x"]) // n
            xb, yb = r["x"][:m * n].reshape(m, n).mean(axis=1), np.median(r["y"][:m * n].reshape(m, n), axis=1)
            fig, axes = plt.subplots(3 if SIMULATION else 2, 1, figsize=(7, 7.5 if SIMULATION else 5.5), sharex=True,
                                     gridspec_kw={"height_ratios": [2, 3, 1] if SIMULATION else [3, 1]})
            a1, a2 = axes[-2], axes[-1]
            if SIMULATION:
                draw_simulation(axes[0], r["x"][0], r["x"][-1])                   # the simulation above, same wavelengths as the window
            a1.plot(r["x"], r["y"], color="lightgray", lw=0.5, label="moving average" if FANO_ON_AVERAGE else "raw data")
            a1.plot(xb, yb, "k.", ms=2, label="block median")
            xx = np.linspace(r["x"][0], r["x"][-1], 2000)
            a1.plot(xx, fano(xx, *r["params"]), color="blue", lw=2, label="Fano fit")
            axes[0].set_title(f"resonance {k + 1}: x0={r['x0']:.5g} nm, Gamma={r['gamma']:.4g} nm, q={r['q']:.2f}", fontsize=10)
            a1.set_ylabel(signal_name())
            a1.legend(loc="best", fontsize=8)
            a2.plot(xb, yb - fano(xb, *r["params"]), "k.", ms=2)
            a2.axhline(0, color="blue", lw=1)
            a2.set_ylabel("residuals")
            a2.set_xlabel("wavelength (nm)")
            plt.tight_layout()
        plt.show()
    return results


# ---------------------------------------------------------------- single measurements

def measurements():
    """Registry of the single measurements: misure_singole/misure.csv (columns id, ora, nome, x_mm, y_mm, con_scansione, file, note)."""
    table = read_table(DATA / "misure_singole" / "misure.csv")
    return table.fillna("")


def spectrum(ch1, ch2, n_points=5000):
    """CH2 inside the CH1 ramp as (wavelength in nm, block median of CH2)."""
    f = fit_ramp(ch1)
    i0, i1 = int(np.ceil(f[0])), int(np.floor(f[1]))
    n = max(1, (i1 - i0 + 1) // n_points)
    k = (i1 - i0 + 1) // n
    return wavelength(i0 + (np.arange(k) + 0.5) * n - 0.5, f), np.median(ch2[i0:i0 + k * n].reshape(k, n), axis=1) - DARK_LEVEL


def load_item(item):
    """What to compare can be: the id of a single measurement (misure_singole/misure.csv), a CSV file in ANY folder, or a point of a run,
    written "run folder:point index" (e.g. "2026-10-06/run_007:12"). Returns (label, ch1, ch2, notes)."""
    if isinstance(item, (int, np.integer)):
        r = measurements().query("id == @item").iloc[0]
        ch1, ch2 = load_point(DATA / "misure_singole" / r.file)
        return f"{r.id} {r.nome}", ch1, ch2, r.note
    folder_text, _, index = str(item).rpartition(":")
    if folder_text and index.isdigit():                                       # "run folder:point index"
        folder, n = run_path(folder_text), int(index)
        row = read_registry(folder).query("indice == @n").iloc[0]
        ch1, ch2 = load_point(folder / row.file)
        return f"{folder.parent.name}/{folder.name} point {n}", ch1, ch2, description(folder)
    path = run_path(str(item))
    ch1, ch2 = load_point(path)
    return "/".join(path.parts[-3:-1] + (path.stem,)), ch1, ch2, ""             # e.g. "2026-10-06/run_007/point_012": two folders tell them apart


def compare_group(items, ref=None, ratio=False, ax=None, title=""):
    """Compares the spectra of one list of items (see load_item). By default only the spectra are drawn and their mean value printed, in the
    units of the signal (volts, or T with POLARIZATION). With ratio=True each one is also divided by the reference ref (default: the first item).
    The items of another kind than the first one (laser scan / no scan) are skipped. ax: panel where the spectra are drawn."""
    items = list(items) if items else list(measurements()["id"])
    ref = items[0] if ref is None else ref
    items = items if ref in items else [ref] + items
    loaded = [load_item(item) for item in items]
    r = items.index(ref)
    scanned = fit_ramp(loaded[r][1]) is not None
    keep = [k for k in range(len(items)) if (fit_ramp(loaded[k][1]) is not None) == scanned]
    for k in set(range(len(items))) - set(keep):
        print(f"skipped {loaded[k][0]}: it is not the same kind of acquisition as the first one (scan / no scan)")
    labels = {k: loaded[k][0] + (f" ({k})" if [loaded[j][0] for j in keep].count(loaded[k][0]) > 1 else "") for k in keep}
    if scanned:
        curves = {k: spectrum(loaded[k][1], loaded[k][2]) for k in keep}
        lam = np.linspace(max(c[0][0] for c in curves.values()), min(c[0][-1] for c in curves.values()), 2000)
        s = {k: np.interp(lam, *c) for k, c in curves.items()}
        if ax is None:
            plt.figure(figsize=(10, 4.5))
            ax = plt.gca()
        for k in keep:
            ax.plot(lam, s[k], lw=2.2 if ratio and k == r else 1.2, label=labels[k])
        ax.set_xlabel("wavelength (nm)")
        ax.set_ylabel(signal_name())
        ax.legend(loc="best", fontsize=8)
        if title:
            ax.set_title(title)
        values = {k: s[k].mean() for k in keep}
        errors = {k: np.nan for k in keep}
        ratios = {k: (s[k] / s[r]).mean() for k in keep}
    else:
        level, err = {}, {}
        for k in keep:
            ch2 = loaded[k][2]
            blocks = np.median(ch2[:len(ch2) // 1000 * 1000].reshape(-1, 1000), axis=1)
            level[k], err[k] = np.median(ch2) - DARK_LEVEL, blocks.std(ddof=1) / np.sqrt(len(blocks))
        values = level
        ratios = {k: level[k] / level[r] for k in keep}
        errors = {k: abs(ratios[k]) * np.hypot(err[k] / level[k], err[r] / level[r]) for k in keep}
        fig, axes = plt.subplots(1, 2 if ratio else 1, figsize=(12 if ratio else 7, 4.5), squeeze=False)
        x = np.arange(len(keep))
        axes[0, 0].bar(x, [level[k] for k in keep], yerr=[err[k] for k in keep], capsize=3)
        panels = [(axes[0, 0], f"{signal_name()} level")]
        if ratio:
            axes[0, 1].bar(x, [ratios[k] for k in keep], yerr=[errors[k] for k in keep], capsize=3)
            axes[0, 1].axhline(1, color="gray", lw=0.8, ls=":")
            panels.append((axes[0, 1], f"ratio to {labels[r]}"))
        for a, t in panels:
            a.set_xticks(x)
            a.set_xticklabels([labels[k] for k in keep], rotation=30, ha="right")
            a.set_ylabel(t)
        plt.tight_layout()
    table = pd.DataFrame({"name": [labels[k] for k in keep], f"mean_{unit()}": [values[k] for k in keep]})
    if scanned is False:
        table["error"] = [err[k] for k in keep]
    if ratio:
        table[f"ratio_to_{labels[r]}"] = [ratios[k] for k in keep]
        table["ratio_error"] = [errors[k] for k in keep]
    table["notes"] = [loaded[k][3] for k in keep]
    print((f"--- {title} ---\n" if title else "") + table.to_string(index=False))
    return table


def compare(items=None, ref=None, ratio=False):
    """Compare spectra that can be anywhere: ids of the single measurements, CSV files in any folder (also the transmission references:
    "reference/copol.csv"), points of a run ("folder:index"), mixed as you like (see load_item). By default the signals are compared as they
    are (volts) with no ratio; ratio=True also divides each one by the reference ref.
    items can also be a dict of lists, e.g. {"co": [...], "cross": [...]}: one panel for each group, one above the other (ref can be a dict too)."""
    groups = items if isinstance(items, dict) else {"": items}
    fig, axes = plt.subplots(len(groups), 1, figsize=(10, 4.5 * len(groups)), sharex=True, squeeze=False)
    tables = {}
    for ax, (name, group) in zip(axes[:, 0], groups.items()):
        tables[name] = compare_group(group, ref[name] if isinstance(ref, dict) else ref, ratio, ax, name)
    if not any(ax.lines for ax in axes[:, 0]):
        plt.close(fig)                                                              # the groups were not scans: they drew their own figures
    else:
        plt.figure(fig.number)
        plt.tight_layout()
    plt.show()
    return tables if isinstance(items, dict) else tables[""]


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
    a = sub.add_parser("mosaic", help="one png with the spectra of all the membranes of a run, each in its position")
    a.add_argument("folder")
    a.add_argument("--recompute", action="store_true", help="read all the CSV files again")
    a = sub.add_parser("view", help="plot one or more single acquisitions (CSV files), no heatmaps")
    a.add_argument("csv", nargs="+")
    a = sub.add_parser("fano", help="Fano fit of the resonances of one CSV")
    a.add_argument("csv")
    sub.add_parser("list", help="list of the single measurements")
    a = sub.add_parser("compare", help="compare spectra from any folder: single-measurement ids, CSV files, points of a run (folder:index)")
    a.add_argument("items", nargs="*", help="ids, CSV files or run_folder:index (empty = all the single measurements)")
    a.add_argument("--ratio", action="store_true", help="also divide each item by the reference (default: only the spectra, in volts)")
    a.add_argument("--ref", help="with --ratio: the item the others are divided by (default: the first)")
    args = p.parse_args()
    if args.command is None:                                         # no arguments: use WHAT TO LOOK AT
        if not (RUN and (HEATMAPS or MOSAIC or POINT is not None)) and VIEW is None and COMPARE is None and not FANO:
            print("Nothing to do: set RUN (with HEATMAPS, MOSAIC and/or POINT), VIEW, COMPARE and/or FANO in WHAT TO LOOK AT at the top of analysis.py, "
                  "or use a command (python analysis.py -h)")
        if RUN and HEATMAPS:
            run_heatmaps(RUN, RECOMPUTE)
        if RUN and MOSAIC:
            mosaic(RUN, RECOMPUTE)
        if RUN and POINT is not None:
            for i in ([POINT] if isinstance(POINT, int) else POINT):
                point_plot(RUN, i, show=False)
            plt.show()                                                   # all the points together
        if VIEW is not None:
            for path in ([VIEW] if isinstance(VIEW, (str, Path)) else VIEW):
                csv_plot(path, show=False)
            plt.show()
        if COMPARE is not None:
            compare(COMPARE, REFERENCE, COMPARE_RATIO)
        if FANO:
            fano_fit(FANO)
    elif args.command == "mosaic":
        mosaic(args.folder, args.recompute)
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
        number = lambda text: int(text) if text.isdigit() else text           # "3" is the id of a single measurement, anything else a file or run:index
        compare([number(i) for i in args.items], number(args.ref) if args.ref else None, args.ratio)


if __name__ == "__main__":
    main()
