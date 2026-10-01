"""
Il TLB-6700 non usa VISA/porte seriali: usa il driver USB proprietario
Newport (Jungo) con una DLL .NET (UsbDllWrap.dll) installata insieme al
software di controllo del laser.

Questo script:
1. cerca automaticamente UsbDllWrap.dll nelle cartelle di installazione
   tipiche di Newport/New Focus
2. ci si connette e stampa la tabella dei dispositivi visti dal driver,
   per scoprire la DeviceKey esatta del tuo laser (è specifica per il
   numero di serie, es. "6700 SN1234")

Richiede:
    pip install pythonnet

Riferimento: esempio pubblico funzionante con lo stesso laser
https://github.com/marvync/Laser-automation-NewFocus
"""

import glob
import os
import sys

# ---------------------------------------------------------------- CONFIGURAZIONE
# ProductID USB standard dei controller TLB-6700 (dall'esempio di riferimento).
PRODUCT_ID = 4106

# Cartelle in cui cercare UsbDllWrap.dll (adatta/aggiungi percorsi se necessario)
CARTELLE_DA_CERCARE = [
    r"C:\Program Files\New Focus",
    r"C:\Program Files (x86)\New Focus",
    r"C:\Program Files\Newport",
    r"C:\Program Files (x86)\Newport",
]
# ---------------------------------------------------------------------------


def trova_dll():
    for cartella in CARTELLE_DA_CERCARE:
        if not os.path.isdir(cartella):
            continue
        risultati = glob.glob(
            os.path.join(cartella, "**", "UsbDllWrap.dll"), recursive=True
        )
        if risultati:
            return risultati[0]
    return None


dll_path = trova_dll()
if dll_path is None:
    print("UsbDllWrap.dll non trovata nelle cartelle note:")
    for c in CARTELLE_DA_CERCARE:
        print(f"  - {c}")
    print(
        "\nCercala manualmente (Esplora file -> cerca 'UsbDllWrap.dll' in C:\\) "
        "e aggiungi la cartella trovata a CARTELLE_DA_CERCARE sopra, oppure "
        "scrivimi il percorso completo che trovi."
    )
    sys.exit(1)

cartella_dll = os.path.dirname(dll_path)
print(f"Trovata DLL in: {dll_path}")

sys.path.append(cartella_dll)

import clr  # noqa: E402  (richiede pythonnet)

clr.AddReference("UsbDllWrap")

import Newport  # noqa: E402

tlb = Newport.USBComm.USB()

print(f"\nApro i dispositivi con ProductID={PRODUCT_ID}...")
tlb.OpenDevices(PRODUCT_ID, True)

tabella = tlb.GetDeviceTable()
print("\nDispositivi trovati (DeviceKey -> info):")
try:
    for voce in tabella:
        print(f"  {voce}")
except TypeError:
    # GetDeviceTable può restituire un oggetto .NET non iterabile direttamente
    print(f"  {tabella}")
    print(
        "\n(Se sopra non vedi le singole DeviceKey, dimmi cosa stampa esattamente "
        "cosi' capiamo come leggere la struttura restituita.)"
    )

tlb.CloseDevices()
print("\nFatto. Usa la DeviceKey qui sopra (es. '6700 SN1234') nello script successivo.")
