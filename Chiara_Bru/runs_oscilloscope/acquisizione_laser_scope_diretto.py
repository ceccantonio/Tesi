"""
Come acquisizione_laser_scope.py, ma la parte del laser è scritta in modo
diretto/letterale (stesso stile di matlab_diretto_trigger_scan.py: nessuna
classe TLB6700, chiamate dirette a OpenDevices/Query/Write), per isolare se
il problema di comunicazione viene dalla classe o dal fatto di aprire
scope+laser nello stesso processo Python.

Flusso identico a acquisizione_laser_scope.py, solo il laser è "diretto".
"""

import os
import sys
import time
import datetime

import numpy as np
import pyvisa
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- CONFIGURAZIONE
OSC_VISA = "USB::0x0699::0x0408::C028813::INSTR"
OSC_TIMEOUT_MS = 30000

CHANNELS = ["CH1", "CH2"]
LABELS = {"CH1": "CH1 (laser ramp)", "CH2": "CH2 (R spectrum)"}
COLORS = {"CH1": "tab:red", "CH2": "tab:blue"}

USBADDR = 1
DEVICE_ID = 0x100A  # hex2dec('100A')

LASER_WAVELENGTH_START = 1520.0  # nm, come trigger_scan.m
LASER_WAVELENGTH_STOP = 1570.0   # nm, come trigger_scan.m
# Velocità di scan: NON impostata qui, come in trigger_scan.m ("Scan velocity
# MUST be set by hand") - resta quella già impostata a mano sul laser.

POLL_INTERVAL_S = 0.3
CHANNEL_READ_PAUSE_S = 0.5
# ---------------------------------------------------------------------------


def crea_cartella_run():
    oggi = datetime.date.today().strftime("%Y-%m-%d")
    base_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), oggi)
    os.makedirs(base_folder, exist_ok=True)
    run_esistenti = [d for d in os.listdir(base_folder) if d.startswith("run_")]
    run_numero = len(run_esistenti) + 1
    run_folder = os.path.join(base_folder, f"run_{run_numero:03d}")
    os.makedirs(run_folder)
    return run_folder


def connetti_oscilloscopio():
    print("\nConnessione all'oscilloscopio Tektronix MDO3104...")
    rm = pyvisa.ResourceManager()
    osc = rm.open_resource(OSC_VISA)
    osc.timeout = OSC_TIMEOUT_MS
    print("Oscilloscopio connesso.")
    return osc


def configura_scope(osc):
    osc.write("HORizontal:MODE MANUAL")
    rec_len = int(osc.query("HORizontal:RECORDLENGTH?"))
    osc.write("DATA:START 1")
    osc.write(f"DATA:STOP {rec_len}")
    osc.write("DATA:WIDTH 1")
    osc.write("DATA:ENCdg ASCii")


def arm_acquisition(osc):
    osc.write("ACQ:STOPAFTER SEQUENCE")
    osc.write("ACQ:STATE RUN")
    # In ARMED lo scope sta acquisendo il pre-trigger (PrTrig) e ignora i
    # trigger: lo scan parte solo quando lo stato non è più ARMED (READY).
    print("Attendo che lo scope finisca il pre-trigger (ARMED -> READY)...")
    while osc.query("TRIGger:STATE?").strip().upper() == "ARMED":
        time.sleep(POLL_INTERVAL_S)
    print("Scope pronto ad accettare il trigger.")


def wait_for_acquisition(osc):
    print("\nAttendo che l'oscilloscopio completi l'acquisizione triggerata...")
    while True:
        status = int(osc.query("ACQ:STATE?"))
        if status == 0:
            print("Acquisizione completata!")
            break
        time.sleep(POLL_INTERVAL_S)


def acquire_channel(osc, source="CH1"):
    osc.write(f"DATA:SOURCE {source}")
    ymult = float(osc.query("WFMPRE:YMULT?"))
    yoff = float(osc.query("WFMPRE:YOFF?"))
    yzero = float(osc.query("WFMPRE:YZERO?"))
    print(f"{source}: YMULT={ymult}, YOFF={yoff}, YZERO={yzero}")
    raw = osc.query("CURVE?")
    ydata = np.array(raw.split(","), dtype=float)
    return (ydata - yoff) * ymult + yzero


def connetti_laser():
    """Stesso codice, stesso ordine, di matlab_diretto_trigger_scan.py."""
    dll_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "laser", "UsbDllWrap.dll"
    )
    dll_path = os.path.abspath(dll_path)
    if not os.path.isfile(dll_path):
        raise RuntimeError(f"UsbDllWrap.dll non trovata in {dll_path}")

    print(f"DLL trovata: {dll_path}")
    sys.path.append(os.path.dirname(dll_path))

    import clr

    clr.AddReference("UsbDllWrap")
    import Newport
    from System.Text import StringBuilder

    np_usb = Newport.USBComm.USB()

    esito = np_usb.OpenDevices(DEVICE_ID)
    print(f"OpenDevices({hex(DEVICE_ID)}) -> {esito}")

    querydata = StringBuilder(64)
    return np_usb, querydata


def laser_query(np_usb, querydata, comando):
    esito = np_usb.Query(USBADDR, comando, querydata)
    risposta = querydata.ToString()
    print(f"Query {comando!r} -> esito={esito}, risposta={risposta!r}")
    return risposta


def laser_write(np_usb, comando):
    esito = np_usb.Write(USBADDR, comando)
    print(f"Write {comando!r} -> esito={esito}")


def main():
    run_folder = crea_cartella_run()
    print(f"\nI dati saranno salvati in:\n{os.path.abspath(run_folder)}")

    osc = connetti_oscilloscopio()
    configura_scope(osc)

    print("\nConnessione al laser TLB-6700 (stile diretto)...")
    np_usb, querydata = connetti_laser()
    laser_query(np_usb, querydata, "*IDN?")
    laser_query(np_usb, querydata, "OUTP:STAT?")

    laser_write(np_usb, f"SOUR:WAVE:START {LASER_WAVELENGTH_START}")
    laser_write(np_usb, f"SOUR:WAVE:STOP {LASER_WAVELENGTH_STOP}")

    arm_acquisition(osc)
    laser_write(np_usb, "OUTP:SCAN:START")

    wait_for_acquisition(osc)

    # Nessun CloseDevices(): i .m del laboratorio non lo chiamano mai.
    time.sleep(1.0)

    data = {}
    for ch in CHANNELS:
        print(f"\nLeggo {ch}...")
        data[ch] = acquire_channel(osc, ch)
        time.sleep(CHANNEL_READ_PAUSE_S)

    csv_path = os.path.join(run_folder, "oscilloscope_sweep.csv")
    min_len = min(len(v) for v in data.values())
    all_cols = np.column_stack([data[ch][:min_len] for ch in CHANNELS])
    np.savetxt(csv_path, all_cols, delimiter=",", header=",".join(CHANNELS), comments="")
    print(f"\nDati salvati in: {csv_path}")

    fig, axes = plt.subplots(len(CHANNELS), 1, figsize=(10, 6), sharex=True)
    for ax, ch in zip(axes, CHANNELS):
        ax.plot(data[ch], color=COLORS[ch], lw=1.0, label=LABELS[ch])
        ax.set_ylabel("Voltage (V)", color=COLORS[ch])
        ax.tick_params(axis="y", labelcolor=COLORS[ch])
        ax.legend(loc="best")
    axes[-1].set_xlabel("Sample index")
    fig.suptitle(f"Acquisizione sincronizzata con scan del laser ({', '.join(CHANNELS)})")
    plt.tight_layout()
    plt.show()

    print("\nAcquisizione completata.")


if __name__ == "__main__":
    main()
