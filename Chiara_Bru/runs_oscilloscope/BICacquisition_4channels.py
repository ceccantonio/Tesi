# BIC Acquisition - 4 channel version

import os
import sys
import time
import numpy as np
import datetime
import pyvisa
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "laser"))
from tlb6700 import TLB6700

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
osc.timeout = 30000
print("Oscilloscope connected.")


#Connect to laser (TLB-6700) and configure the scan
# DEVICE_KEY: scoprila eseguendo scripts/laser/trova_tlb6700.py (è specifica
# del numero di serie del tuo laser, es. "6700 SN1234").
LASER_DEVICE_KEY = "6700 SN23163"  # numero di serie letto dall'etichetta del laser

LASER_WAVELENGTH_START = 1550.0  # nm, <-- da impostare
LASER_WAVELENGTH_STOP = 1560.0   # nm, <-- da impostare
LASER_SCAN_SPEED_FORWARD = 0.1   # nm/s, <-- da impostare
LASER_SCAN_SPEED_BACKWARD = 5.0  # nm/s, <-- da impostare

print("\nConnecting to TLB-6700 laser...")
laser = TLB6700(LASER_DEVICE_KEY)
laser.open()
print("Laser connected:", laser.idn())

laser.set_scan_limits(LASER_WAVELENGTH_START, LASER_WAVELENGTH_STOP)
laser.set_scan_speeds(LASER_SCAN_SPEED_FORWARD, LASER_SCAN_SPEED_BACKWARD)
laser.set_output(True)


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


# Arm the scope, then wait for the triggered acquisition to complete
def arm_acquisition():
    """
    Arm the scope for a single triggered acquisition. Call this BEFORE
    starting the laser scan, so the scope is already waiting for the
    trigger (laser ramp on CH1) when the scan begins.
    """
    osc.write("ACQ:STOPAFTER SEQUENCE")
    osc.write("ACQ:STATE RUN")


def wait_for_acquisition():
    """
    Wait until the scope completes a single triggered acquisition.
    Call arm_acquisition() (and start the laser scan) before this.
    """
    print("\nWaiting for the oscilloscope to finish the triggered acquisition...")

    while True:
        status = int(osc.query("ACQ:STATE?"))
        if status == 0:
            print("Acquisition complete!")
            break
        time.sleep(0.3)  # evita di martellare lo scope con query continue


#Set termination for the channels that need 50 ohm (comment out if already set manually)
for ch in CHANNELS_50OHM:
    set_termination(ch, ohm50=True)


#MAIN: arm the scope, start the laser scan (triggers the scope via CH1's
# laser-ramp signal), wait for acquisition, then read all four channels once
arm_acquisition()
laser.start_scan()

wait_for_acquisition()

laser.stop_scan()
laser.set_output(False)
laser.close()
print("Laser scan stopped and connection closed.")

time.sleep(1.0)  # lascia che lo scope si assesti dopo lo stop dell'acquisizione

data = {}
for ch in CHANNELS:
    print(f"\nLeggo {ch}...")
    data[ch] = acquire_channel(ch)
    time.sleep(0.5)  # respiro tra un canale e l'altro (trasferimenti ASCII grandi)

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