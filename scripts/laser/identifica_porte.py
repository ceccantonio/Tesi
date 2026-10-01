"""
Identifica quale risorsa VISA (porta seriale ASRL, USB, ecc.) corrisponde
a un certo strumento, interrogandolo con *IDN? a diverse combinazioni di
baud rate / terminatore (su una porta seriale, se le impostazioni non sono
quelle giuste lo strumento non risponde e si ottiene solo un timeout).

Utile per distinguere, ad esempio, laser (TLB-6700) e oscilloscopio quando
compaiono entrambi nell'elenco di pyvisa.ResourceManager().list_resources()
e non è ovvio quale sia quale, o per capire con che baud rate risponde
uno strumento seriale sconosciuto.

Uso:
1. Esegui prima questo per vedere tutte le risorse disponibili:
     python -c "import pyvisa; print(pyvisa.ResourceManager().list_resources())"
2. Metti le risorse da testare in RISORSE sotto ed esegui questo script.
3. Se tutte le combinazioni danno timeout su una risorsa, probabilmente
   quella risorsa non è uno strumento che risponde a *IDN? (es. è lo
   scope già identificato altrove, o una porta non collegata a nulla).
"""

import pyvisa

# ---------------------------------------------------------------- CONFIGURAZIONE
RISORSE = ["ASRL1::INSTR", "ASRL4::INSTR", "ASRL5::INSTR"]
TIMEOUT_MS = 800  # basso apposta: dobbiamo provare molte combinazioni

BAUD_RATES = [9600, 19200, 38400, 57600, 115200]
TERMINATORI = ["\r\n", "\r", "\n"]
# ---------------------------------------------------------------------------

rm = pyvisa.ResourceManager()

for res in RISORSE:
    print(f"\n=== {res} ===")
    trovato = False

    for baud in BAUD_RATES:
        for term in TERMINATORI:
            inst = None
            try:
                inst = rm.open_resource(res)
                inst.baud_rate = baud
                inst.timeout = TIMEOUT_MS
                inst.read_termination = term
                inst.write_termination = term
                idn = inst.query("*IDN?")
                print(f"  baud={baud}, term={repr(term)} -> RISPOSTA: {idn.strip()}")
                trovato = True
            except Exception as e:
                # Timeout/errore: combinazione sbagliata, si prova la prossima
                pass
            finally:
                if inst is not None:
                    try:
                        inst.close()
                    except Exception:
                        pass

    if not trovato:
        print("  Nessuna combinazione ha dato risposta (sempre timeout/errore).")
