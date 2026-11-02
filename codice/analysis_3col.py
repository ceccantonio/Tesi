"""Analysis of files with 3 columns (CH1 = laser ramp, transmission detector, reflection detector): same analysis as analysis.py
(which is not modified and is used as it is), but you choose how many columns to read and what to analyze: transmission and/or reflection.

    python analysis_3col.py run 2026-10-06/run_007                  # every command of analysis.py works (run, point, view, nm, mosaic, fano, compare...)
    python analysis_3col.py --channels R run 2026-10-06/run_007     # only the reflection
    python analysis_3col.py --channels TR --columns 3 mosaic 2026-10-06/run_007
    python analysis_3col.py                                         # no arguments: does what is written in "WHAT TO LOOK AT" below

The method parameters (WINDOW, V_START, LAMBDA_0/1, FANO_*, ...) stay in analysis.py. The results of the reflection have R in the name
(analysis_R.csv, heatmap_modulation_R.png, mosaic_spectra_R.png...), the ones of the transmission are the same as analysis.py.
With both channels ("TR") the plots of the points (point, view) and the mosaic show T, R and their sum T + R on the same axes,
and the means are printed; with POLARIZATION the absorption (and losses) is estimated as A = 1 - (T + R). The run also gets the heatmaps of
the mean T + R and of A (sum_TR*.csv). The reflection divided by the reference (POLARIZATION) is R = (reflection - DARK_R) / (transmission of the reference without the sample) * R_SCALE:
the light reflected is compared with the light sent, up to the relative gain of the two detectors (R_SCALE, to be calibrated, e.g. with a mirror).
"""
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text
from pathlib import Path

import analysis as base

# ---------------------------------------------------------------- WHAT TO LOOK AT
COLUMNS = 3                 # 3: read the 3 columns of the files; 2: read only CH1 and the transmission (the third column is ignored; also for old 2-column files)
CHANNELS = "TR"             # what to analyze: "T" = transmission, "R" = reflection, "TR" = both (one after the other)
T_COLUMN = "CH2"            # name of the transmission column in the files (if the name is not there, the 2nd column is used)
R_COLUMN = "CH3"            # name of the reflection column in the files (if the name is not there, the 3rd column is used)
DARK_T = 0.0                # level (V) of the transmission detector without light
DARK_R = 0.0                # level (V) of the reflection detector without light
R_SCALE = 1.0               # reflection with POLARIZATION: R is multiplied by this (relative gain of the reflection detector with respect to the transmission one)
SIMULATION_R = ""           # simulation file of the reflection, drawn above the reflection plots ("" = none); the one of the transmission is SIMULATION
SIM_COLUMNS_R = (0, 1)      # (column of the wavelength/frequency, column of the simulated reflection) in SIMULATION_R

# same as the ones in analysis.py: what to look at (see there for the explanation of each one)
RUN = ""
HEATMAPS = True
POINT = None
VIEW = None
COMPARE = None
COMPARE_RATIO = False
REFERENCE = None
RECOMPUTE = False
SAVE_NM = True
SAVE_NM_POINTS = None
SIMULATION = ""
SIM_COLUMNS = (0, 1)
MOSAIC = False
MOSAIC_YSCALE = "common"
MOSAIC_BLOCKS = (1, 1)
POLARIZATION = ""
REFERENCE_SET = "reference"
FANO = ""

LOOK_AT = ("RUN", "HEATMAPS", "POINT", "VIEW", "COMPARE", "COMPARE_RATIO", "REFERENCE", "RECOMPUTE", "SAVE_NM", "SAVE_NM_POINTS",
           "SIMULATION", "SIM_COLUMNS", "MOSAIC", "MOSAIC_YSCALE", "MOSAIC_BLOCKS", "POLARIZATION", "REFERENCE_SET", "FANO")

_channel = "T"              # the channel that analysis.py is working on now
_file = {}                  # the last file read, with all its columns: the second channel does not read it again
_saved = set()              # spectra_nm files already written in this run of the program


# ---------------------------------------------------------------- reading
def read_columns(path):
    """{"CH1", "T", "R", "n"} of one CSV: the arrays (float32) and the number of columns of the file. "R" is None if the file or COLUMNS do not have it."""
    key = (str(path), COLUMNS)
    if _file.get("key") != key:
        df = pd.read_csv(path, dtype=np.float32, comment="#", usecols=range(2) if COLUMNS == 2 else None)   # what follows a # is ignored
        df.columns = df.columns.str.strip()

        def column(name, position):
            if name in df.columns:
                return df[name].to_numpy()
            return df.iloc[:, position].to_numpy() if position < df.shape[1] else None

        _file.clear()
        _file.update(key=key, data={"CH1": column("CH1", 0), "T": column(T_COLUMN, 1), "R": column(R_COLUMN, 2) if COLUMNS >= 3 else None,
                                    "n": df.shape[1]})
    return _file["data"]


def read_csv(path):
    """(CH1, signal of the current channel) of one acquisition, as float arrays."""
    data = read_columns(path)
    if data[_channel] is None:
        raise ValueError(f"{Path(path).name}: no reflection column to read (COLUMNS = {COLUMNS}, the file has {data['n']} columns)")
    return data["CH1"].astype(float), data[_channel].astype(float)


def save_nm(folder, name, ch1, volts):
    """spectra_nm/name with the wavelength (nm) and all the detector columns (V), only from LAMBDA_0 to LAMBDA_1: the same file whatever the channel."""
    key = (str(folder), name)
    if key in _saved:
        return
    f = base.fit_ramp(ch1)
    if f is None:
        print(f"{name}: no ramp in CH1, not saved in spectra_nm")
        return
    data = read_columns(Path(folder) / name)
    i0, i1 = int(np.ceil(f[0])), int(np.floor(f[1]))
    columns = {"wavelength_nm": base.wavelength(np.arange(i0, i1 + 1), f)}
    for channel, label in (("T", T_COLUMN), ("R", R_COLUMN)):
        if data[channel] is not None:
            columns[f"{label}_V"] = data[channel][i0:i1 + 1].astype(float)
    table = pd.DataFrame(columns)
    if base.SAVE_NM_POINTS:
        n = max(1, len(table) // int(base.SAVE_NM_POINTS))
        k = len(table) // n
        blocks = table.iloc[:k * n].to_numpy().reshape(k, n, -1)
        table = pd.DataFrame(np.column_stack([blocks[:, :, 0].mean(axis=1), np.median(blocks[:, :, 1:], axis=1)]), columns=table.columns)
    out = Path(folder) / "spectra_nm"
    out.mkdir(exist_ok=True)
    table.to_csv(out / name, index=False, float_format="%.5f")
    for extra in ("griglia.csv", "description.txt"):
        if (Path(folder) / extra).exists() and not (out / extra).exists():
            base.shutil.copy(Path(folder) / extra, out / extra)
    _saved.add(key)


# ---------------------------------------------------------------- what depends on the channel
_reference_spectrum, _normalize, _spectrum, _run_heatmaps = base.reference_spectrum, base.normalize, base.spectrum, base.run_heatmaps
_csv_plot, _mosaic = base.csv_plot, base.mosaic
_in_reference = False


def reference_spectrum():
    """The reference is always the transmission measured without the sample, whatever the channel being analyzed."""
    global _in_reference
    current = _channel
    set_channel("T")
    _in_reference = True
    try:
        return _reference_spectrum()
    finally:
        _in_reference = False
        set_channel(current)


def spectrum(ch1, ch2, n_points=5000):
    """analysis.spectrum subtracts the dark level: with POLARIZATION the data already had it subtracted (normalize), so it is added back."""
    lam, y = _spectrum(ch1, ch2, n_points)
    return lam, (y + base.DARK_LEVEL if base.POLARIZATION and not _in_reference else y)


def normalize(ch1, ch2):
    out = _normalize(ch1, ch2)
    return out * R_SCALE if _channel == "R" and base.POLARIZATION else out


def suffix():
    return ("_R" if _channel == "R" else "") + (f"_T_{base.POLARIZATION}pol" if base.POLARIZATION else "")


def unit():
    return _channel if base.POLARIZATION else "V"


def signal_name():
    if _channel == "R":
        return "reflection R" if base.POLARIZATION else f"{R_COLUMN} (V)"
    return "transmission T" if base.POLARIZATION else f"{T_COLUMN} (V)"


def set_channel(channel):
    """Points analysis.py to the transmission ("T") or to the reflection ("R"): column, dark level and simulation."""
    global _channel
    _channel = channel
    base.DARK_LEVEL = DARK_R if channel == "R" else DARK_T
    if channel == "R":
        base.SIMULATION, base.SIM_COLUMNS = SIMULATION_R, SIM_COLUMNS_R
    else:
        base.SIMULATION, base.SIM_COLUMNS = SIMULATION, SIM_COLUMNS


# ---------------------------------------------------------------- both channels: T and R together, and their sum
DUAL = False                # True when T and R are both analyzed
_done = set()               # plots already drawn (analysis.py asks for them once per channel)


def names():
    """(label of T, label of R, label of the sum, y axis)"""
    if base.POLARIZATION:
        return "T", "R", "T + R", "T, R"
    return f"{T_COLUMN} (V)", f"{R_COLUMN} (V)", f"{T_COLUMN} + {R_COLUMN} (V)", "signal (V)"


def summary(t, r):
    """Text with the mean values; with POLARIZATION also the absorption (and losses) A = 1 - (T + R)."""
    text = f"mean T = {t:.4f}, R = {r:.4f}, T + R = {t + r:.4f}"
    if base.POLARIZATION:
        return text + f"  ->  A = 1 - (T + R) = {1 - t - r:.4f}"
    return text + "  (volts: set POLARIZATION to compare them with the light sent and estimate the absorption)"


def pair_spectra(path):
    """(wavelength, T, R) of one CSV, each one normalized if POLARIZATION; None if there is no ramp."""
    current, out = _channel, {}
    try:
        for c in "TR":
            set_channel(c)
            ch1, y = base.load_point(path)
            if base.fit_ramp(ch1) is None:
                return None
            out[c] = base.spectrum(ch1, y)
    finally:
        set_channel(current)
    return out["T"][0], out["T"][1], np.interp(out["T"][0], *out["R"])


def overlay_csv_plot(path, label=None, show=True):
    """Plot of ONE acquisition with T, R and T + R on the same axes (replaces csv_plot when both channels are analyzed)."""
    path = base.run_path(path)
    if ("plot", str(path)) in _done:
        return
    _done.add(("plot", str(path)))
    pair = pair_spectra(path)
    if pair is None:
        print("no ramp in CH1: nothing to plot")
        return
    lam, t, r = pair
    sims = []
    for channel, file, columns in (("T", SIMULATION, SIM_COLUMNS), ("R", SIMULATION_R, SIM_COLUMNS_R)):
        if file:
            base.SIM_COLUMNS = columns
            x, y = base.read_simulation(file)
            sims.append((channel, x, y))
    fig, axes = plt.subplots(2 if sims else 1, 1, figsize=(10, 7 if sims else 4.8), sharex=True, squeeze=False)
    axes = axes[:, 0]
    name_t, name_r, name_sum, ylabel = names()
    a = axes[-1]
    a.plot(lam, t, color="tab:blue", lw=1.2, label=name_t)
    a.plot(lam, r, color="tab:red", lw=1.2, label=name_r)
    a.plot(lam, t + r, color="black", lw=1.5, label=name_sum)
    if base.POLARIZATION:
        a.axhline(1, color="gray", lw=0.8, ls=":")
    for channel, x, y in sims:
        inside = (x >= base.LAMBDA_0) & (x <= base.LAMBDA_1)
        axes[0].plot(x[inside], y[inside], color="tab:blue" if channel == "T" else "tab:red", lw=1.3, label=f"simulation {channel}")
    if sims:
        axes[0].set_ylabel("simulation")
        axes[0].legend(loc="best", fontsize=8)
    a.set_xlim(base.LAMBDA_0, base.LAMBDA_1)
    a.set_xlabel("wavelength (nm)")
    a.set_ylabel(ylabel)
    a.legend(loc="best", fontsize=8)
    text = summary(t.mean(), r.mean())
    axes[0].set_title(f"{label or path.name}\n{text}", fontsize=9)
    print(f"{label or path.name}: {text}")
    plt.tight_layout()
    if show:
        plt.show()


def overlay_mosaic(folder, recompute=False):
    """ONE png with T, R and T + R of every membrane, each panel in the position of its membrane (replaces mosaic when both channels are analyzed)."""
    folder = base.run_path(folder)
    if ("mosaic", str(folder)) in _done:
        return
    _done.add(("mosaic", str(folder)))
    current, data = _channel, {}
    try:
        for c in "TR":
            set_channel(c)
            data[c] = dict(zip(*base.run_spectra(folder, recompute)))
    finally:
        set_channel(current)
    index = [i for i in data["T"] if i in data["R"]]
    registry = base.read_registry(folder).set_index("indice").loc[index]
    lay = base.layout(registry["x_mm"].to_numpy(), registry["y_mm"].to_numpy())
    if lay is None or not index:
        print("mosaic: the points are not a regular lattice (or there is no data): the spectra cannot be placed")
        return
    rows, columns, n_rows, n_columns = lay
    curves = [(data["T"][i], data["R"][i], data["T"][i] + data["R"][i]) for i in index]
    low, high = np.nanmin([c for k in curves for c in k]), np.nanmax([c for k in curves for c in k])
    name_t, name_r, name_sum, ylabel = names()

    def draw(ax, k):
        ax.axis("on")
        t, r, total = curves[k]
        ax.plot(base.LAM_GRID, t, color="tab:blue", lw=0.7)
        ax.plot(base.LAM_GRID, r, color="tab:red", lw=0.7)
        ax.plot(base.LAM_GRID, total, color="black", lw=0.8)
        ax.set_xlim(base.LAMBDA_0, base.LAMBDA_1)
        lo, hi = (min(np.nanmin(c) for c in curves[k]), max(np.nanmax(c) for c in curves[k])) if base.MOSAIC_YSCALE == "each" else (low, high)
        ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
        ax.set_xticks([base.LAMBDA_0, (base.LAMBDA_0 + base.LAMBDA_1) / 2, base.LAMBDA_1])
        ax.yaxis.set_major_locator(base.MaxNLocator(3))
        ax.tick_params(labelsize=7, length=3, width=0.8, pad=2)
        labels = ax.get_xticklabels()
        labels[0].set_ha("left")
        labels[-1].set_ha("right")
        ax.grid(alpha=0.3, lw=0.3)
        ax.text(0.04, 0.93, str(index[k]), transform=ax.transAxes, fontsize=8, va="top", fontweight="bold",
                bbox={"facecolor": "white", "alpha": 0.85, "pad": 1.5, "lw": 0})

    row_blocks = np.array_split(np.arange(n_rows), max(1, min(base.MOSAIC_BLOCKS[0], n_rows)))
    column_blocks = np.array_split(np.arange(n_columns), max(1, min(base.MOSAIC_BLOCKS[1], n_columns)))
    total_blocks = len(row_blocks) * len(column_blocks)
    width, height = (2.0, 1.7) if total_blocks == 1 else (2.8, 2.3)
    n = 0
    for rb in row_blocks:
        for cb in column_blocks:
            n += 1
            fig, axes = plt.subplots(len(rb), len(cb), figsize=(width * len(cb) + 0.8, height * len(rb) + 1.3), squeeze=False)
            for ax in axes.flat:
                ax.axis("off")
            for k in range(len(index)):
                if rows[k] in rb and columns[k] in cb:
                    draw(axes[rows[k] - rb[0], columns[k] - cb[0]], k)
            fig.supxlabel("wavelength (nm)", fontsize=11)
            fig.supylabel(ylabel, fontsize=11)
            where = "" if total_blocks == 1 else f"\nblock {n} of {total_blocks}: rows {rb[0] + 1}-{rb[-1] + 1}, columns {cb[0] + 1}-{cb[-1] + 1} of {n_rows} x {n_columns}"
            fig.suptitle(f"{folder.name}: blue = {name_t}, red = {name_r}, black = {name_sum}, {base.LAMBDA_0:g} to {base.LAMBDA_1:g} nm" + where
                         + (f"\n{base.description(folder)}" if base.description(folder) else ""), fontsize=9)
            plt.tight_layout(rect=(0.01, 0.01, 1, 0.93), h_pad=1.5, w_pad=1.5)
            name = "mosaic_spectra_TR" + (f"_T_{base.POLARIZATION}pol" if base.POLARIZATION else "") + ("" if total_blocks == 1 else f"_part{n}of{total_blocks}") + ".png"
            fig.savefig(folder / name, dpi=150)
            print("Saved", folder / name)
    plt.show()


def sum_maps(folder, recompute=False):
    """Mean T, R, T + R (and A = 1 - T - R with POLARIZATION) of every point of a run: sum_TR*.csv and the heatmaps of T + R and of A."""
    folder = base.run_path(folder)
    current, data = _channel, {}
    try:
        for c in "TR":
            set_channel(c)
            data[c] = dict(zip(*base.run_spectra(folder, recompute)))
    finally:
        set_channel(current)
    registry = base.read_registry(folder).set_index("indice")
    rows = []
    for i in data["T"]:
        if i in data["R"]:
            t, r = np.nanmean(data["T"][i]), np.nanmean(data["R"][i])
            rows.append({"index": i, "x_mm": registry.loc[i, "x_mm"], "y_mm": registry.loc[i, "y_mm"], "T": t, "R": r, "T_plus_R": t + r,
                         "absorption": 1 - t - r if base.POLARIZATION else np.nan})
    table = pd.DataFrame(rows)
    pol = f"_T_{base.POLARIZATION}pol" if base.POLARIZATION else ""
    table.to_csv(folder / f"sum_TR{pol}.csv", index=False)
    print("Saved", folder / f"sum_TR{pol}.csv")
    volts = "" if base.POLARIZATION else " (V)"
    base.heatmap(table, "T_plus_R", f"mean T + R between {base.LAMBDA_0:g} and {base.LAMBDA_1:g} nm{volts}", folder / f"heatmap_sum_TR{pol}.png",
                 fmt="{:.3f}", subtitle=base.description(folder))
    if base.POLARIZATION:
        base.heatmap(table, "absorption", "absorption + losses A = 1 - (T + R), mean value", folder / f"heatmap_absorption{pol}.png",
                     fmt="{:.3f}", subtitle=base.description(folder))


def run_heatmaps(folder, recompute=False):
    """The heatmaps of the channel being analyzed; after the last channel, also the ones of the sum."""
    _run_heatmaps(folder, recompute)
    if DUAL and _channel == "R":
        sum_maps(folder)


class _Pyplot:
    """analysis.py writes "raw CH2" and "transmission T" in the legends: with the reflection they are renamed just before showing the figures."""

    def __getattr__(self, name):
        return getattr(plt, name)

    def show(self, *args, **kwargs):
        if _channel == "R":
            for number in plt.get_fignums():
                for text in plt.figure(number).findobj(Text):
                    text.set_text(text.get_text().replace("raw CH2", f"raw {R_COLUMN}").replace("transmission T", "reflection R"))
        return plt.show(*args, **kwargs)


# ---------------------------------------------------------------- from the terminal
def main():
    global COLUMNS, CHANNELS, DUAL
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--columns", type=int, choices=(2, 3), help="2: only CH1 and transmission; 3: also the reflection column")
    p.add_argument("--channels", help="T, R or TR")
    own, rest = p.parse_known_args()
    sys.argv = [sys.argv[0]] + rest
    COLUMNS = own.columns or COLUMNS
    CHANNELS = (own.channels or CHANNELS).upper()
    channels = [c for c in "TR" if c in CHANNELS]
    if COLUMNS == 2 and "R" in channels:
        print("COLUMNS = 2: the reflection column is not read, only the transmission is analyzed (use COLUMNS = 3 for the reflection)")
        channels.remove("R")
    if not channels:
        print("Nothing to analyze: CHANNELS must contain T and/or R")
        return
    for name in LOOK_AT:                                                  # what is written here goes to analysis.py
        setattr(base, name, globals()[name])
    base.REFERENCE_DIR = base.DATA / base.REFERENCE_SET
    base.read_csv, base.save_nm, base.reference_spectrum, base.normalize, base.spectrum = read_csv, save_nm, reference_spectrum, normalize, spectrum
    DUAL = len(channels) == 2
    # with both channels: T and R together in the plots of the points, in the mosaic and with their sum
    base.csv_plot, base.mosaic, base.run_heatmaps = (overlay_csv_plot, overlay_mosaic, run_heatmaps) if DUAL else (_csv_plot, _mosaic, _run_heatmaps)
    base.suffix, base.unit, base.signal_name, base.plt = suffix, unit, signal_name, _Pyplot()
    for channel in channels:
        print(f"=== {'transmission' if channel == 'T' else 'reflection'} ({T_COLUMN if channel == 'T' else R_COLUMN}) ===")
        set_channel(channel)
        base.main()


if __name__ == "__main__":
    main()
