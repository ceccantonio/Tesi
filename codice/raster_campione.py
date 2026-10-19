"""Raster scan di tutto il campione: N_X x N_Y membrane, una scansione del laser su ognuna, percorso a meandro.

    python raster_campione.py

Il campione e' inclinato: i centri stanno su un reticolo ruotato, con l'angolo ricavato da PRIMO_Y e Y_FINE_PRIMA_RIGA.
Chiudere PI MicroMove prima di partire (la seriale puo' essere usata da un solo programma).
"""
import numpy as np

from strumenti import LIMITI_MM, init, scansiona_punti

PRIMO_X = None              # <-- X del centro della prima membrana (in alto a sinistra), in mm
PRIMO_Y = 133.978           # Y del centro della prima membrana
Y_FINE_PRIMA_RIGA = 134.078 # Y del centro dell'ultima membrana della prima riga (da qui l'inclinazione)
N_X, N_Y = 11, 11           # membrane lungo X e lungo Y
PASSO = 0.30                # mm tra membrane vicine
DESTRA_E_X_CRESCENTE = True # "destra" = X crescente
GIU_E_Y_DECRESCENTE = True  # "scendere" = Y decrescente


def percorso():
    """Centri delle membrane a meandro: prima riga verso destra, poi si scende e si torna indietro, ecc."""
    sx = 1.0 if DESTRA_E_X_CRESCENTE else -1.0
    sy = -1.0 if GIU_E_Y_DECRESCENTE else 1.0
    uy = (Y_FINE_PRIMA_RIGA - PRIMO_Y) / ((N_X - 1) * PASSO)       # seno dell'inclinazione
    ux = sx * np.sqrt(1 - uy ** 2)                                  # versore lungo la riga
    vx, vy = uy, -ux                                                # versore verso la riga successiva
    if np.sign(vy) != sy:
        vx, vy = -uy, ux
    punti = []
    for j in range(N_Y):
        for i in (range(N_X) if j % 2 == 0 else range(N_X - 1, -1, -1)):
            punti.append((PRIMO_X + PASSO * (i * ux + j * vx), PRIMO_Y + PASSO * (i * uy + j * vy)))
    return punti, np.degrees(np.arcsin(uy))


def main():
    if PRIMO_X is None:
        print("Scrivi PRIMO_X (X della prima membrana) in cima a questo file")
        return
    punti, angolo = percorso()
    x, y = np.array(punti).T
    print(f"{N_X} x {N_Y} membrane, passo {PASSO} mm, inclinazione {angolo:.3f} gradi")
    print(f"Si parte da X={x[0]:.4f} Y={y[0]:.4f} e si finisce a X={x[-1]:.4f} Y={y[-1]:.4f}")
    print(f"X da {x.min():.3f} a {x.max():.3f}, Y da {y.min():.3f} a {y.max():.3f} (limiti {LIMITI_MM[0]}..{LIMITI_MM[1]} mm)")
    print(f"Spazio su disco: circa {len(punti) * 0.24:.1f} GB")
    if min(x.min(), y.min()) < LIMITI_MM[0] or max(x.max(), y.max()) > LIMITI_MM[1]:
        print("Punti fuori dai limiti dei motori: controlla PRIMO_X, PRIMO_Y e il passo")
        return
    if input(f"Muovere i motori ed eseguire {len(punti)} scansioni? scrivi 'si': ").strip().lower() != "si":
        return
    laser, osc, mot = init()
    print("Dati in", scansiona_punti(punti, laser, osc, mot))


if __name__ == "__main__":
    main()
