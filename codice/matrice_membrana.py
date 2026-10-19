"""Matrice di punti su UNA membrana: N_X x N_Y punti dal bordo al bordo, una scansione del laser su ognuno, a meandro.

    python matrice_membrana.py

Chiudere PI MicroMove prima di partire (la seriale puo' essere usata da un solo programma).
"""
import numpy as np

from strumenti import LIMITI_MM, init, scansiona_punti

CENTRO_X, CENTRO_Y = 91.814, 146.508   # mm, centro della membrana
LATO_UM = 500.0                        # i punti vanno da -lato/2 a +lato/2 (bordi inclusi)
N_X, N_Y = 8, 8
DESTRA_E_X_CRESCENTE = True            # "destra" = X crescente
GIU_E_Y_DECRESCENTE = True             # "scendere" = Y decrescente; si parte dalla riga in alto


def percorso():
    """Punti a meandro: riga in alto verso destra, poi si scende e si torna indietro, ecc."""
    xs = CENTRO_X + np.linspace(-LATO_UM / 2, LATO_UM / 2, N_X) / 1000
    ys = CENTRO_Y + np.linspace(-LATO_UM / 2, LATO_UM / 2, N_Y) / 1000
    righe = ys[::-1] if GIU_E_Y_DECRESCENTE else ys
    verso_destra = xs if DESTRA_E_X_CRESCENTE else xs[::-1]
    return [(float(x), float(y)) for j, y in enumerate(righe) for x in (verso_destra if j % 2 == 0 else verso_destra[::-1])]


def main():
    punti = percorso()
    x, y = np.array(punti).T
    print(f"Matrice {N_X} x {N_Y}, passo X = {LATO_UM / (N_X - 1):.2f} um, passo Y = {LATO_UM / (N_Y - 1):.2f} um")
    print(f"X da {x.min():.4f} a {x.max():.4f}, Y da {y.min():.4f} a {y.max():.4f} (limiti {LIMITI_MM[0]}..{LIMITI_MM[1]} mm)")
    print(f"Spazio su disco: circa {len(punti) * 0.24:.1f} GB")
    if min(x.min(), y.min()) < LIMITI_MM[0] or max(x.max(), y.max()) > LIMITI_MM[1]:
        print("Punti fuori dai limiti dei motori: controlla il centro")
        return
    if input(f"Muovere i motori ed eseguire {len(punti)} scansioni? scrivi 'si': ").strip().lower() != "si":
        return
    laser, osc, mot = init()
    print("Dati in", scansiona_punti(punti, laser, osc, mot))


if __name__ == "__main__":
    main()
