"""
Identifica quale risorsa VISA (porta seriale ASRL, USB, ecc.) corrisponde
a un certo strumento, interrogandolo con *IDN?.

Utile per distinguere, ad esempio, laser (TLB-6700) e oscilloscopio quando
compaiono entrambi nell'elenco di pyvisa.ResourceManager().list_resources()
e non è ovvio quale sia quale.

Uso:
1. Esegui prima questo per vedere tutte le risorse disponibili:
     python -c "import pyvisa; print(pyvisa.ResourceManager().list_resources())"
2. Metti le risorse da testare in RISORSE sotto ed esegui questo script.
"""

import pyvisa

# ---------------------------------------------------------------- CONFIGURAZIONE
RISORSE = ["ASRL1::INSTR", "ASRL5::INSTR"]
TIMEOUT_MS = 2000
# ---------------------------------------------------------------------------

rm = pyvisa.ResourceManager()

for res in RISORSE:
    inst = None
    try:
        inst = rm.open_resource(res)
        inst.timeout = TIMEOUT_MS
        inst.read_termination = "\r\n"
        inst.write_termination = "\r\n"
        idn = inst.query("*IDN?")
        print(f"{res} -> {idn.strip()}")
    except Exception as e:
        print(f"{res} -> ERRORE: {e}")
    finally:
        if inst is not None:
            try:
                inst.close()
            except Exception:
                pass
