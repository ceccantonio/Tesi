"""
Template generico per fit di dati con scipy.optimize.curve_fit.

Carica x, y, sigma_y da un file di testo (loadtxt), esegue il fit,
stampa i parametri con incertezza, il chi-quadro ridotto, e disegna
dati + curva di fit + pannello dei residui.

Modifica:
  - la sezione CONFIGURAZIONE (file, colonne, etichette)
  - la funzione modello() con la forma funzionale che ti serve
  - p0 (stima iniziale dei parametri) e NOMI_PARAMETRI
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ---------------------------------------------------------------- CONFIGURAZIONE
FILE_DATI = "data/esempio.txt"
COL_X, COL_Y, COL_SIGMA_Y = 0, 1, 2

LABEL_X = "x [unità]"
LABEL_Y = "y [unità]"
TITOLO = "Fit dei dati"

NOMI_PARAMETRI = ["m", "q"]
P0 = [1.0, 0.0]  # stima iniziale dei parametri, nell'ordine di modello()
# ---------------------------------------------------------------------------


def modello(x, m, q):
    """Sostituisci con la funzione da fittare, es. una retta m*x + q."""
    return m * x + q


def main():
    x, y, sigma_y = np.loadtxt(
        FILE_DATI, unpack=True, usecols=(COL_X, COL_Y, COL_SIGMA_Y)
    )

    popt, pcov = curve_fit(modello, x, y, p0=P0, sigma=sigma_y, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))

    residui = y - modello(x, *popt)
    chi2 = np.sum((residui / sigma_y) ** 2)
    ndof = len(x) - len(popt)
    chi2_rid = chi2 / ndof

    print("Parametri del fit:")
    for nome, val, err in zip(NOMI_PARAMETRI, popt, perr):
        print(f"  {nome} = {val:.6g} +/- {err:.2g}")
    print(f"chi2 = {chi2:.3f}, ndof = {ndof}, chi2/ndof = {chi2_rid:.3f}")

    # ------------------------------------------------------------ PLOT
    fig, (ax_fit, ax_res) = plt.subplots(
        2, 1, figsize=(7, 7), sharex=True, height_ratios=[3, 1]
    )

    x_fit = np.linspace(x.min(), x.max(), 500)
    ax_fit.errorbar(x, y, yerr=sigma_y, fmt="o", markersize=4, capsize=3, label="dati")
    ax_fit.plot(x_fit, modello(x_fit, *popt), "-", label="fit")
    ax_fit.set_ylabel(LABEL_Y)
    ax_fit.set_title(TITOLO)
    ax_fit.grid(True, alpha=0.3)
    ax_fit.legend()

    ax_res.errorbar(x, residui, yerr=sigma_y, fmt="o", markersize=4, capsize=3)
    ax_res.axhline(0, color="k", linewidth=1)
    ax_res.set_xlabel(LABEL_X)
    ax_res.set_ylabel("residui")
    ax_res.grid(True, alpha=0.3)

    fig.tight_layout()
    # fig.savefig("fit.png", dpi=300)   # decommenta per salvare
    plt.show()


if __name__ == "__main__":
    main()
