"""Mosaic to compare the spectra of two runs, membrane by membrane (legible in a presentation).

    python compare_runs.py 2026-10-27/run_009 2026-10-27/run_010            # one pair of panels per membrane: run A and run B side by side
    python compare_runs.py runA runB --pol co --mode overlay --blocks 4 4   # both runs in the same panel, 16 pngs
    python compare_runs.py runA runB --pol co --ref-a reference_2026-10-27 --ref-b reference_2026-11-03   # runs of different samples: each with its own reference
    python compare_runs.py                                                  # no arguments: does what is written in "WHAT TO COMPARE" below

The two runs are matched by position in the grid (row, column of the membrane), not by index. Every panel has the same axes (wavelength from
LAMBDA_0 to LAMBDA_1 of analysis.py, y from YLIM with the ticks YTICKS), the two runs have different colors. The spectra are the ones of
analysis.py (run_spectra, cached in spectra*.npz of each run), so the run is read again only the first time.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import analysis as base

# ---------------------------------------------------------------- WHAT TO COMPARE
RUN_A = ""                  # first run (folder, as in analysis.py: "2026-10-27/run_009"); "" = nothing to do
RUN_B = ""                  # second run
LABEL_A = ""                # name of run A in the legend; "" = the folder name
LABEL_B = ""
REFERENCE_A = ""            # transmission reference of run A (folder with copol.csv / crosspol.csv, inside runs_oscilloscope: "reference_2026-10-27"); "" = REFERENCE_SET of analysis.py
REFERENCE_B = ""            # the same for run B: different samples / set-ups have each its own reference (used with POLARIZATION)
COLORS = ("#0072B2", "#D55E00")   # color of run A and of run B (blue and orange: distinguishable also by color-blind people and in grey)
POLARIZATION = ""           # as in analysis.py: "co" or "cross" = spectra divided by the transmission reference (T); "" = volts
MODE = "side"               # "side": every membrane has two panels next to each other (A on the left, B on the right); "overlay": the two runs in the same panel
BLOCKS = (3, 3)             # split in several pngs with bigger panels: (blocks along the rows, blocks along the columns). (1, 1) = one png with everything
YLIM = (0, 1)               # vertical axis of every panel
YTICKS = (0, 0.5, 1)        # ticks (with numbers) on the vertical axis of every panel
FONT = 13                   # size of the numbers on the axes (points); titles and labels are a bit bigger
LINEWIDTH = 1.6
RECOMPUTE = False           # True: read all the points of both runs again
OUT = ""                    # folder of the pngs; "" = runs_oscilloscope/confronti/<A>_vs_<B>


def load(run, recompute, reference=""):
    """({(row, column): (point index, spectrum on LAM_GRID)}, n_rows, n_columns) of a run, with the grid drawn like the real sample.
    The spectra are divided by the reference of that run (only if they still have to be computed: the ones in spectra*.npz are not computed again)."""
    base.REFERENCE_DIR = base.run_path(reference or base.REFERENCE_SET)
    folder = base.run_path(run)
    index, spectra = base.run_spectra(folder, recompute)
    registry = base.read_registry(folder).set_index("indice").loc[index]
    lay = base.layout(registry["x_mm"].to_numpy(), registry["y_mm"].to_numpy())
    if lay is None or len(index) == 0:
        raise SystemExit(f"{run}: the points are not a regular lattice (or there is no data): they cannot be placed")
    rows, columns, n_rows, n_columns = lay
    return {(int(r), int(c)): (int(i), s) for r, c, i, s in zip(rows, columns, index, spectra)}, n_rows, n_columns


def style(ax):
    ax.set_xlim(base.LAMBDA_0, base.LAMBDA_1)
    ax.set_ylim(*YLIM)
    ax.set_xticks([base.LAMBDA_0, (base.LAMBDA_0 + base.LAMBDA_1) / 2, base.LAMBDA_1])
    ax.set_yticks(YTICKS)
    ax.set_yticklabels([f"{t:g}" for t in YTICKS])
    ax.tick_params(labelsize=FONT, length=4, width=1.0, pad=2)
    labels = ax.get_xticklabels()                                              # the end labels stay inside their own panel
    labels[0].set_ha("left")
    labels[-1].set_ha("right")
    ax.grid(alpha=0.35, lw=0.5)


def compare_mosaic(run_a, run_b, recompute=False):
    a, rows_a, columns_a = load(run_a, recompute, REFERENCE_A)
    b, rows_b, columns_b = load(run_b, recompute, REFERENCE_B)
    if (rows_a, columns_a) != (rows_b, columns_b):
        print(f"WARNING: the grids are different ({rows_a} x {columns_a} and {rows_b} x {columns_b}): matched by position, missing panels stay empty")
    n_rows, n_columns = max(rows_a, rows_b), max(columns_a, columns_b)
    names = (LABEL_A or Path(run_a).name, LABEL_B or Path(run_b).name)
    if names[0] == names[1]:                                                   # runs with the same name in different days
        names = tuple(Path(r).parent.name + "/" + Path(r).name for r in (run_a, run_b))
    safe = lambda text: str(text).replace("\\", "_").replace("/", "_").replace(":", "")
    out = Path(OUT) if OUT else base.DATA / "confronti" / f"{safe(run_a)}_vs_{safe(run_b)}"
    out.mkdir(parents=True, exist_ok=True)

    row_blocks = np.array_split(np.arange(n_rows), max(1, min(BLOCKS[0], n_rows)))
    column_blocks = np.array_split(np.arange(n_columns), max(1, min(BLOCKS[1], n_columns)))
    total = len(row_blocks) * len(column_blocks)
    side = MODE == "side"
    panel_w, panel_h = 1.8, 1.7                                                # inches of a panel (the numbers are in points: bigger panels = relatively smaller numbers)
    cell_w = 2 * panel_w + 1.4 if side else panel_w + 0.9
    cell_h = 2.5
    margin = {"left": 1.0, "right": 0.3, "bottom": 1.0, "top": 1.5}            # inches around the panels: axis labels, title and legend
    n = 0
    for rb in row_blocks:
        for cb in column_blocks:
            n += 1
            width = margin["left"] + margin["right"] + cell_w * len(cb)
            height = margin["bottom"] + margin["top"] + cell_h * len(rb)
            fig = plt.figure(figsize=(width, height))
            outer = fig.add_gridspec(len(rb), len(cb), left=margin["left"] / width, right=1 - margin["right"] / width,
                                     bottom=margin["bottom"] / height, top=1 - margin["top"] / height, wspace=0.2, hspace=0.45)
            for r in rb:
                for c in cb:
                    cell = outer[r - rb[0], c - cb[0]]
                    if (r, c) not in a and (r, c) not in b:
                        continue
                    if side:                                                   # the same membrane of the two runs, adjacent
                        inner = cell.subgridspec(1, 2, wspace=0.3)
                        axes = [fig.add_subplot(inner[0]), fig.add_subplot(inner[1])]
                    else:
                        axes = [fig.add_subplot(cell)] * 2
                    for ax, data, color in zip(axes, (a, b), COLORS):
                        if (r, c) in data:
                            ax.plot(base.LAM_GRID, data[(r, c)][1], color=color, lw=LINEWIDTH)
                        style(ax)
                    ids = [d[(r, c)][0] for d in (a, b) if (r, c) in d]
                    labels = [(ax, str(d[(r, c)][0])) for ax, d in zip(axes, (a, b)) if (r, c) in d] if side else \
                        [(axes[0], " | ".join(str(i) for i in dict.fromkeys(ids)))]                        # overlay: "A | B" only if the indices differ
                    for ax, text in labels:
                        ax.text(0.04, 0.93, text, transform=ax.transAxes, fontsize=FONT, va="top", fontweight="bold",
                                bbox={"facecolor": "white", "alpha": 0.85, "pad": 1.5, "lw": 0})                  # point index, readable even over the curve
            fig.supxlabel("wavelength (nm)", fontsize=FONT + 2, y=0.01)
            fig.supylabel(base.signal_name(), fontsize=FONT + 2, x=0.01)
            where = "" if total == 1 else f"   (block {n} of {total}: rows {rb[0] + 1}-{rb[-1] + 1}, columns {cb[0] + 1}-{cb[-1] + 1} of {n_rows} x {n_columns})"
            fig.suptitle(f"{names[0]} vs {names[1]}: {base.signal_name()}, {base.LAMBDA_0:g}-{base.LAMBDA_1:g} nm{where}", fontsize=FONT + 2, y=0.995)
            handles = [Line2D([], [], color=color, lw=3) for color in COLORS]
            fig.legend(handles, names, loc="upper center", ncol=2, fontsize=FONT + 2, frameon=False, bbox_to_anchor=(0.5, 1 - 0.42 / height))
            name = "mosaic_confronto" + base.suffix() + ("" if total == 1 else f"_part{n}of{total}") + ".png"
            fig.savefig(out / name, dpi=200)
            plt.close(fig)
            print("Saved", out / name)


def main():
    global RUN_A, RUN_B, REFERENCE_A, REFERENCE_B, POLARIZATION, MODE, BLOCKS, RECOMPUTE, OUT
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_a", nargs="?", default=RUN_A)
    p.add_argument("run_b", nargs="?", default=RUN_B)
    p.add_argument("--ref-a", default=REFERENCE_A, help="reference folder of run A (with --pol)")
    p.add_argument("--ref-b", default=REFERENCE_B, help="reference folder of run B (with --pol)")
    p.add_argument("--pol", default=POLARIZATION, help='"co" or "cross": divide by the transmission reference; "" = volts')
    p.add_argument("--mode", default=MODE, choices=["side", "overlay"])
    p.add_argument("--blocks", type=int, nargs=2, default=BLOCKS, metavar=("ROWS", "COLUMNS"))
    p.add_argument("--recompute", action="store_true", default=RECOMPUTE, help="read all the points again")
    p.add_argument("--out", default=OUT, help="folder of the pngs")
    args = p.parse_args()
    if not args.run_a or not args.run_b:
        raise SystemExit("Nothing to do: set RUN_A and RUN_B in WHAT TO COMPARE at the top of compare_runs.py, or write the two runs after the command")
    REFERENCE_A, REFERENCE_B = args.ref_a, args.ref_b
    POLARIZATION, MODE, BLOCKS, RECOMPUTE, OUT = args.pol, args.mode, tuple(args.blocks), args.recompute, args.out
    base.POLARIZATION = POLARIZATION                                           # the spectra of analysis.py are computed (and cached) with this
    compare_mosaic(args.run_a, args.run_b, RECOMPUTE)


if __name__ == "__main__":
    main()
