import os
import numpy as np
import matplotlib.pyplot as plt

#
date_str = "2026-04-23"   # cartella data
run_number = 17            # numero run

# costruisci il percorso
run_folder = os.path.join(os.getcwd(), date_str, f"run_{run_number:03d}")
csv_path = os.path.join(run_folder, "oscilloscope_sweep.csv")

# controllo esistenza file
if not os.path.exists(csv_path):
    raise FileNotFoundError(f"File non trovato:\n{csv_path}")

print(f"\nLoading data from:\n{csv_path}")

# carica dati
data = np.loadtxt(csv_path, delimiter=",", skiprows=1)

ch1 = data[:, 0]
ch2 = data[:, 1]

#PLOT
fig, ax1 = plt.subplots(figsize=(10, 5))

# CH1
line1, = ax1.plot(ch1, color='tab:red', label="CH1 (laser ramp)", lw=1.0)
ax1.set_xlabel("Sample index")
ax1.set_ylabel("CH1 Voltage (V)", color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')

# CH2
ax2 = ax1.twinx()
line2, = ax2.plot(ch2, color='tab:blue', label="CH2 (T spectrum)", lw=1.0)
ax2.set_ylabel("CH2 Voltage (V)", color='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

# legenda
lines = [line1, line2]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc="best")

plt.title(f"Run {run_number:03d} - {date_str}")
plt.tight_layout()
plt.show()