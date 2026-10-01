# Codice MATLAB di riferimento

Codice del laboratorio (non scritto da noi) per pilotare il TLB-6700 da
MATLAB tramite la stessa DLL .NET (`UsbDllWrap.dll`) usata da `tlb6700.py`.
Tenuto qui come riferimento/documentazione del protocollo di comunicazione
provato funzionante su questo laser specifico — `tlb6700.py` è la
conversione Python di questa logica (indirizzo USB numerico, non DeviceKey
testuale; nessun comando `SYSTem:MCONtrol`).

- `NP_USB_connect_velocity.m` — apertura dispositivo, query di base
- `TLB6700_init.m` — inizializzazione (apertura + check stato)
- `trigger_scan.m` — avvio scan e attesa a tempo fisso
- `NP_USB_reperror.m` / `NP_USB_reperror_v2.m` — due versioni della
  funzione di interpretazione dei codici di errore del firmware USB
  (leggera inconsistenza tra le due nel codice originale: la v1 tratta
  sia 0 che 1 come successo, la v2 tratta solo 0 come successo; in
  `tlb6700.py` abbiamo seguito la v1, più permissiva)
