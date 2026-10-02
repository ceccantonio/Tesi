"""
Acquisizione sincronizzata: scan del laser TLB-6700 + oscilloscopio Tektronix
MDO3104 (4 canali).

Flusso:
1. Connessione a oscilloscopio e laser
2. Configurazione scan del laser (limiti, velocità) e dei canali dello scope
3. Si arma lo scope per un'acquisizione singola (ACQ:STOPAFTER SEQUENCE)
4. Si avvia lo scan del laser: la rampa di scan è cablata su CH1 dello
   scope, che quindi si triggera da sola quando la rampa parte
5. Si aspetta che l'acquisizione finisca, si legge CH1-CH4, si salva un CSV
   e si plottano i 4 canali

Sostituisce BICacquisition_4channels.py (che resta come riferimento/storico,
ma questo è lo script da usare per le acquisizioni con scan del laser).
"""

import os
import sys
import time
import datetime

import numpy as np
import pyvisa
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "laser"))
from tlb6700 import TLB6700

# ---------------------------------------------------------------- CONFIGURAZIONE
OSC_VISA = "USB::0x0699::0x0408::C028813::INSTR"
OSC_TIMEOUT_MS = 30000

CHANNELS = ["CH1", "CH2"]  # solo i canali effettivamente collegati
CHANNELS_50OHM = []  # canali da mettere a 50 ohm (vuoto se già impostati a mano)

LABELS = {
    "CH1": "CH1 (laser ramp)",
    "CH2": "CH2 (R spectrum)",
}
COLORS = {"CH1": "tab:red", "CH2": "tab:blue"}

LASER_WAVELENGTH_START = 1520.0  # nm, come trigger_scan.m (range completo)
LASER_WAVELENGTH_STOP = 1570.0  # nm, come trigger_scan.m (range completo)
LASER_SCAN_SPEED_FORWARD = 0.1  # nm/s, <-- da impostare
LASER_SCAN_SPEED_BACKWARD = 5.0  # nm/s, <-- da impostare

POLL_INTERVAL_S = 0.3  # pausa tra una query di stato scope e la successiva
CHANNEL_READ_PAUSE_S = 0.5  # pausa tra la lettura di un canale e il successivo
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

    for ch in CHANNELS_50OHM:
        valore = "FIFty"
        osc.write(f"{ch}:TERmination {valore}")
        readback = osc.query(f"{ch}:TERmination?").strip()
        print(f"{ch}: termination impostata a {valore} (scope risponde: {readback})")


def arm_acquisition(osc):
    """Arma lo scope per una singola acquisizione triggerata. Chiamare PRIMA
    di avviare lo scan del laser, così lo scope è già in attesa del trigger
    (la rampa del laser su CH1) quando lo scan parte."""
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


def main():
    run_folder = crea_cartella_run()
    print(f"\nI dati saranno salvati in:\n{os.path.abspath(run_folder)}")

    osc = connetti_oscilloscopio()
    configura_scope(osc)

    print("\nConnessione al laser TLB-6700...")
    laser = TLB6700()
    laser.open()
    try:
        print("Laser connesso:", laser.idn())

        laser.setup_sweep(
            start_nm=LASER_WAVELENGTH_START,
            stop_nm=LASER_WAVELENGTH_STOP,
            speed=LASER_SCAN_SPEED_FORWARD,
            return_speed=LASER_SCAN_SPEED_BACKWARD,
        )
        laser.on()

        # Arma lo scope, poi avvia lo scan (la rampa su CH1 fa scattare il trigger)
        arm_acquisition(osc)
        laser.start_sweep()

        wait_for_acquisition(osc)

        laser.stop_sweep()
        laser.off()
    finally:
        # Chiude sempre la connessione al laser, anche se qualcosa sopra
        # fallisce: altrimenti il dispositivo resta "agganciato" e il
        # tentativo successivo fallisce con errori di comunicazione.
        laser.close()
        print("Connessione al laser chiusa.")

    time.sleep(1.0)  # lascia che lo scope si assesti dopo lo stop dell'acquisizione

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
