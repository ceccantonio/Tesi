"""
Template generico per caricare dati da file di testo e fare un plot.

Formato atteso del file dati (separato da spazi/tab, righe che iniziano
con # sono ignorate):

    x    y    sigma_y
    0.0  0.98 0.10
    ...

Modifica solo la sezione "CONFIGURAZIONE" per adattarlo ai tuoi dati.
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- CONFIGURAZIONE
FILE_DATI = "data/esempio.txt"

# colonne nel file: indice a partire da 0
COL_X, COL_Y, COL_SIGMA_Y = 0, 1, 2

LABEL_X = "x [unità]"
LABEL_Y = "y [unità]"
TITOLO = "Titolo del grafico"
# ---------------------------------------------------------------------------

x, y, sigma_y = np.loadtxt(FILE_DATI, unpack=True, usecols=(COL_X, COL_Y, COL_SIGMA_Y))

fig, ax = plt.subplots(figsize=(7, 5))
ax.errorbar(x, y, yerr=sigma_y, fmt="o", markersize=4, capsize=3, label="dati")

ax.set_xlabel(LABEL_X)
ax.set_ylabel(LABEL_Y)
ax.set_title(TITOLO)
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
# fig.savefig("plot.png", dpi=300)   # decommenta per salvare
plt.show()
