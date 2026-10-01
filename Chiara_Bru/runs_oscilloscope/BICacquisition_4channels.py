# BIC Acquisition - 4 channel version

import os
import numpy as np
import datetime
import pyvisa
import matplotlib.pyplot as plt

#Data folder
today_str = datetime.date.today().strftime("%Y-%m-%d")
#base_folder = os.path.join(os.getcwd(), today_str)
base_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), today_str)
os.makedirs(base_folder, exist_ok=True)

existing_runs = [d for d in os.listdir(base_folder) if d.startswith("run_")]
run_number = len(existing_runs) + 1
run_folder = os.path.join(base_folder, f"run_{run_number:03d}")
os.makedirs(run_folder)

print(f"\nData will be saved in:\n{os.path.abspath(run_folder)}")


#Connect to oscilloscope
OSC_VISA = 'USB::0x0699::0x0408::C028813::INSTR'

print("\nConnecting to Tektronix MDO3104...")
rm = pyvisa.ResourceManager()
osc = rm.open_resource(OSC_VISA)
osc.timeout = 5000
print("Oscilloscope connected.")


#Configure full-screen waveform transfer
# Force manual horizontal mode (so record length is fixed)
osc.write("HORizontal:MODE MANUAL")

# Request whatever record length is currently displayed on screen
rec_len = int(osc.query("HORizontal:RECORDLENGTH?"))

# Read out entire record exactly as displayed
osc.write(f"DATA:START 1")
osc.write(f"DATA:STOP {rec_len}")

# Make sure formatting is correct
osc.write("DATA:WIDTH 1")
osc.write("DATA:ENCdg ASCii")


# List of channels to acquire
CHANNELS = ["CH1", "CH2", "CH3", "CH4"]

# Channels that need to be set to 50 ohm termination (edit as needed).
# If a channel is already set manually on the scope front panel, you can
# leave this list empty and skip the set_termination() call below.
CHANNELS_50OHM = ["CH3", "CH4"]


#Optional: set channel termination/impedance
def set_termination(channel, ohm50=True):
    """
    Set the input termination of a channel to 50 ohm or 1 Mohm.
    On Tektronix MDO3000-series scopes this is the CH<x>:TERmination command.
    """
    value = "FIFty" if ohm50 else "MEG"
    osc.write(f"{channel}:TERmination {value}")
    readback = osc.query(f"{channel}:TERmination?").strip()
    print(f"{channel}: termination set to {value} (scope reports: {readback})")


#Read waveform function
def acquire_channel(source="CH1"):
    """
    Acquire the exact waveform currently visible on the oscilloscope screen.
    """
    osc.write(f"DATA:SOURCE {source}")

    ymult = float(osc.query("WFMPRE:YMULT?"))
    yoff = float(osc.query("WFMPRE:YOFF?"))
    yzero = float(osc.query("WFMPRE:YZERO?"))
    print(f"{source}: YMULT={ymult}, YOFF={yoff}, YZERO={yzero}")

    raw = osc.query("CURVE?")
    ydata = np.array(raw.split(','), dtype=float)

    volts = (ydata-yoff)*ymult + yzero
    return volts


# Wait for triggered acquisition
def wait_for_acquisition():
    """
    Wait until the scope completes a single triggered acquisition.
    """
    print("\nWaiting for the oscilloscope to finish the triggered acquisition...")

    osc.write("ACQ:STOPAFTER SEQUENCE")
    osc.write("ACQ:STATE RUN")

    while True:
        status = int(osc.query("ACQ:STATE?"))
        if status == 0:
            print("Acquisition complete!")
            break


#Set termination for the channels that need 50 ohm (comment out if already set manually)
for ch in CHANNELS_50OHM:
    set_termination(ch, ohm50=True)


#MAIN: Wait for trigger, then read all four channels once
wait_for_acquisition()

data = {ch: acquire_channel(ch) for ch in CHANNELS}

# Save all four channels in a single CSV file (truncated to the shortest one)
csv_path = os.path.join(run_folder, "oscilloscope_sweep.csv")

min_len = min(len(v) for v in data.values())
all_cols = np.column_stack([data[ch][:min_len] for ch in CHANNELS])

np.savetxt(csv_path, all_cols, delimiter=",", header=",".join(CHANNELS), comments="")


#Plot all four channels, each on its own subplot with its own voltage axis
colors = {"CH1": "tab:red", "CH2": "tab:blue", "CH3": "tab:green", "CH4": "tab:orange"}
labels = {
    "CH1": "CH1 (laser ramp)",
    "CH2": "CH2 (R spectrum)",
    "CH3": "CH3",
    "CH4": "CH4",
}

fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

for ax, ch in zip(axes, CHANNELS):
    ax.plot(data[ch], color=colors[ch], lw=1.0, label=labels[ch])
    ax.set_ylabel("Voltage (V)", color=colors[ch])
    ax.tick_params(axis='y', labelcolor=colors[ch])
    ax.legend(loc="best")

axes[-1].set_xlabel("Sample index")
fig.suptitle("Full Sweep Acquisition (CH1, CH2, CH3, CH4)")
plt.tight_layout()
plt.show()

print("\nAcquisition complete.")
#print("Saved file:", os.path.abspath(os.path.join(run_folder, 'oscilloscope_sweep.npy')))