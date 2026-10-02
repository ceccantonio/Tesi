"""
Stampa i valori delle costanti di errore esposte dalla DLL Newport
(m_kn...) e prova una query *IDN?, per capire a quale errore corrisponde
il codice -1 quando il laser smette di rispondere.

Stesse chiamate di matlab_diretto.py, più la lettura delle costanti.
Lancialo quando il laser è nello stato bloccato (-1).
"""

import os
import sys

USBADDR = 1
DEVICE_ID = 0x100A

dll_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "UsbDllWrap.dll")
if not os.path.isfile(dll_path):
    raise RuntimeError(f"UsbDllWrap.dll non trovata in {dll_path}")

sys.path.append(os.path.dirname(dll_path))

import clr  # noqa: E402

clr.AddReference("UsbDllWrap")
import Newport  # noqa: E402
from System.Text import StringBuilder  # noqa: E402

NP_USB = Newport.USBComm.USB()

print("Costanti di errore della DLL:")
for nome in dir(NP_USB):
    if nome.startswith("m_kn"):
        try:
            print(f"  {nome} = {getattr(NP_USB, nome)}")
        except Exception as e:
            print(f"  {nome} -> ERRORE lettura: {e}")

esito = NP_USB.OpenDevices(DEVICE_ID)
print(f"\nOpenDevices({hex(DEVICE_ID)}) -> {esito}")

querydata = StringBuilder(64)
esito = NP_USB.Query(USBADDR, "*IDN?", querydata)
print(f"Query '*IDN?' -> esito={esito}, risposta={querydata.ToString()!r}")

try:
    print(f"NumProductsConnected = {NP_USB.NumProductsConnected}")
except Exception as e:
    print(f"NumProductsConnected -> ERRORE: {e}")
