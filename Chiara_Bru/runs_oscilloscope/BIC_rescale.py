# BIC Acquisition

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


#MAIN: Wait for trigger, then read CH1 & CH2 once
wait_for_acquisition()

ch1 = acquire_channel("CH1")
ch2 = acquire_channel("CH2")

# Save both channels in a simple two-column CSV file
csv_path = os.path.join(run_folder, "oscilloscope_sweep.csv")

# Make sure both channels have the same length
min_len = min(len(ch1), len(ch2))
two_col_data = np.column_stack((ch1[:min_len], ch2[:min_len]))

np.savetxt(csv_path, two_col_data, delimiter=",", header="CH1,CH2", comments="")


#Plot the two channels
#plt.figure(figsize=(10, 5))

fig,ax1 = plt.subplots(figsize=(10, 5))

#CH1
#ax1.plot(ch1, label="CH1 (laser ramp)", lw=1.0)
line1, = ax1.plot(ch1, color='tab:red', label="CH1 (laser ramp)", lw=1.0)
ax1.set_xlabel("Sample index")
ax1.set_ylabel("CH1 Voltage (V)", color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')

#CH2
ax2 = ax1.twinx()
#ax2.plot(ch2, label="CH2 (BIC resonances)", lw=1.0)
line2, = ax2.plot(ch2, color='tab:blue', label="CH2 (R spectrum)", lw=1.0)
ax2.set_ylabel("CH2 Voltage (V)", color='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

#plt.legend()
lines = [line1, line2]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc="best")

plt.title("Full Sweep Acquisition (CH1 & CH2)")
plt.tight_layout()
plt.show()
 
print("\nAcquisition complete.")
#print("Saved file:", os.path.abspath(os.path.join(run_folder, 'oscilloscope_sweep.npy')))