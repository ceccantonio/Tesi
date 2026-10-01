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

import os
import sys

USBADDR = 1
DEVICE_ID = 0x100A  # hex2dec('100A')

# Solo la DLL nel repo (confermata funzionante col MATLAB), nessun fallback
# su altre DLL eventualmente installate altrove sul PC.
dll_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "UsbDllWrap.dll")
if not os.path.isfile(dll_path):
    raise RuntimeError(f"UsbDllWrap.dll non trovata in {dll_path}")

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
