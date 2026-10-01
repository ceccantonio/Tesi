"""
Traduzione diretta, riga per riga, di NP_USB_connect_velocity.m
(scripts/laser/riferimento_matlab/NP_USB_connect_velocity.m).

Nessuna classe, nessun controllo di errori custom: fa esattamente quello
che fa lo script MATLAB del laboratorio, per isolare se un problema viene
dal protocollo di comunicazione di base o dal nostro wrapper (tlb6700.py).

Corrispondenza con il MATLAB:
    USBADDR = 1;                              -> USBADDR = 1
    deviceID = hex2dec('100A');               -> DEVICE_ID = 0x100A
    NET.addAssembly(...)                      -> clr.AddReference("UsbDllWrap")
    NPasm.AssemblyHandle.GetType(...)          -> import Newport
    System.Activator.CreateInstance(...)       -> Newport.USBComm.USB()
    NP_USB.OpenDevices(deviceID)                -> NP_USB.OpenDevices(DEVICE_ID)
    NP_USB.Query(USBADDR, '*IDN?', querydata)   -> NP_USB.Query(USBADDR, "*IDN?", querydata)
"""

import glob
import os
import sys

USBADDR = 1
DEVICE_ID = 0x100A  # hex2dec('100A')

CARTELLE_DA_CERCARE = [
    r"C:\Program Files\New Focus",
    r"C:\Program Files (x86)\New Focus",
    r"C:\Program Files\Newport",
    r"C:\Program Files (x86)\Newport",
]

dll_path = None
for cartella in CARTELLE_DA_CERCARE:
    if os.path.isdir(cartella):
        trovati = glob.glob(os.path.join(cartella, "**", "UsbDllWrap.dll"), recursive=True)
        if trovati:
            dll_path = trovati[0]
            break

if dll_path is None:
    raise RuntimeError(
        "UsbDllWrap.dll non trovata nelle cartelle note: "
        + ", ".join(CARTELLE_DA_CERCARE)
    )

print(f"DLL trovata: {dll_path}")
sys.path.append(os.path.dirname(dll_path))

import clr  # noqa: E402  (richiede pythonnet)

clr.AddReference("UsbDllWrap")
import Newport  # noqa: E402
from System.Text import StringBuilder  # noqa: E402

NP_USB = Newport.USBComm.USB()

# Open the USB device
esito = NP_USB.OpenDevices(DEVICE_ID)
print(f"\nOpenDevices({hex(DEVICE_ID)}) -> {esito}")

# The Query method sends the passed in command string to the specified
# device and reads the response data.
querydata = StringBuilder(64)

esito = NP_USB.Query(USBADDR, "*IDN?", querydata)
print(f"Query(USBADDR={USBADDR}, '*IDN?') -> esito={esito}, risposta={querydata.ToString()!r}")

esito = NP_USB.Query(USBADDR, "OUTP:STAT?", querydata)
print(f"Query(USBADDR={USBADDR}, 'OUTP:STAT?') -> esito={esito}, risposta={querydata.ToString()!r}")

NP_USB.CloseDevices()
print("\nDispositivo chiuso.")
