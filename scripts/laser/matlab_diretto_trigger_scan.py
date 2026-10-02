"""
Traduzione diretta, riga per riga, di NP_USB_connect_velocity.m seguito da
trigger_scan.m (scripts/laser/riferimento_matlab/) - la stessa identica
sequenza che hai già confermato funzionare nel MATLAB del laboratorio:
connessione + query di stato + avvio scan + pausa fissa.

Nessuna classe, nessuna logica aggiunta: stesse chiamate, stesso ordine.
"""

import os
import sys
import time

USBADDR = 1
DEVICE_ID = 0x100A  # hex2dec('100A')

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

# --- NP_USB_connect_velocity.m ---
esito = NP_USB.OpenDevices(DEVICE_ID)
print(f"\nOpenDevices({hex(DEVICE_ID)}) -> {esito}")

querydata = StringBuilder(64)

esito = NP_USB.Query(USBADDR, "*IDN?", querydata)
print(f"Query '*IDN?' -> esito={esito}, risposta={querydata.ToString()!r}")

esito = NP_USB.Query(USBADDR, "OUTP:STAT?", querydata)
print(f"Query 'OUTP:STAT?' -> esito={esito}, risposta={querydata.ToString()!r}")

# --- trigger_scan.m ---
pausetime = 10  # Pause during scan (stesso valore del .m)

print("\n*************************")
print("Waiting for laser scan")
print("*************************")

esito = NP_USB.Write(USBADDR, "OUTP:SCAN:START")
print(f"Write 'OUTP:SCAN:START' -> esito={esito}")

time.sleep(pausetime)

# Nessun CloseDevices(): i .m del laboratorio non lo chiamano mai.
