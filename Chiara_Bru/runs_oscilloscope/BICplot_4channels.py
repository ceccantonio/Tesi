import os
import numpy as np
import matplotlib.pyplot as plt

#
date_str = "2026-08-28"   # cartella data
run_number = 2            # numero run

# costruisci il percorso (relativo alla cartella dello script, non alla cwd del processo)
base_folder = os.path.dirname(os.path.abspath(__file__))
run_folder = os.path.join(base_folder, date_str, f"run_{run_number:03d}")
csv_path = os.path.join(run_folder, "oscilloscope_sweep.csv")

# controllo esistenza file
if not os.path.exists(csv_path):
    raise FileNotFoundError(f"File non trovato:\n{csv_path}")

print(f"\nLoading data from:\n{csv_path}")

# carica dati
data = np.loadtxt(csv_path, delimiter=",", skiprows=1)

ch1 = data[:, 0]
ch2 = data[:, 1]
ch3 = data[:, 2]
ch4 = data[:, 3]

#PLOT - 4 subplot impilati, ognuno con la propria scala di tensione
colors = {"CH1": "tab:red", "CH2": "tab:blue", "CH3": "tab:green", "CH4": "tab:orange"}
labels = {
    "CH1": "CH1 (laser ramp)",
    "CH2": "CH2 (R spectrum)",
    "CH3": "CH3",
    "CH4": "CH4",
}
channels_data = {"CH1": ch1, "CH2": ch2, "CH3": ch3, "CH4": ch4}

fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

for ax, ch in zip(axes, ["CH1", "CH2", "CH3", "CH4"]):
    ax.plot(channels_data[ch], color=colors[ch], lw=1.0, label=labels[ch])
    ax.set_ylabel("Voltage (V)", color=colors[ch])
    ax.tick_params(axis='y', labelcolor=colors[ch])
    ax.legend(loc="best")

axes[-1].set_xlabel("Sample index")
fig.suptitle(f"Run {run_number:03d} - {date_str}")
plt.tight_layout()
plt.show()