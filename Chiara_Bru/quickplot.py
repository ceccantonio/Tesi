import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# ---- EDIT ONLY THESE IF YOUR SWEEP SETTINGS CHANGE ----
F_START = 490e3   # Hz
F_STOP  = 550e3   # Hz
BIDIRECTIONAL = True   # sweep goes start->stop->start
# ----------------------------------------------------------

# Zurich-style colors: signal 1 = blue, signal 2 = green, signal 3 = red
COLORS = ["blue", "green", "red"]

folder = r"C:\Users\Labo\Documents\Zurich Instruments\LabOne\WebServer\session_20260826_170822_01\I3_bouncingSpectrum_ARharm123_500mV_002"
all_files = glob.glob(os.path.join(folder, "*.csv"))
data_files = [f for f in all_files if "header" not in os.path.basename(f).lower()]

def get_row_values(df, field):
    rows = df[df["fieldname"] == field]
    if rows.empty:
        return None
    best = rows.iloc[(~rows.iloc[:, 4:].isna()).sum(axis=1).values.argmax()]
    return best.iloc[4:].to_numpy(dtype=float)

def build_bidirectional_freq(n_points, f_start, f_stop):
    n_up = n_points // 2
    n_down = n_points - n_up
    freq_up = np.linspace(f_start, f_stop, n_up, endpoint=False)
    freq_down = np.linspace(f_stop, f_start, n_down)
    return np.concatenate([freq_up, freq_down])

plt.figure(figsize=(10, 6))

for i, f in enumerate(data_files):
    df = pd.read_csv(f, sep=None, engine="python")

    r = get_row_values(df, "r")
    if r is None:
        x = get_row_values(df, "x")
        y = get_row_values(df, "y")
        if x is not None and y is not None:
            r = np.sqrt(x**2 + y**2)

    if r is None:
        print(f"Skipping {os.path.basename(f)}: no r/x/y found")
        continue

    n = len(r)
    freq = build_bidirectional_freq(n, F_START, F_STOP) if BIDIRECTIONAL else np.linspace(F_START, F_STOP, n)

    label = os.path.splitext(os.path.basename(f))[0]
    color = COLORS[i % len(COLORS)]   # cycles if more than 3 files
    plt.plot(freq, r, lw=1.0, label=label, color=color)
    print(f"{os.path.basename(f)}: {n} points, color={color}")

plt.xlabel("Frequency (Hz)")
plt.ylabel("R (V)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()