"""Analisi dei dati: ampiezza di modulazione, quality factor, heatmap e confronto di misure singole.

    python analisi.py run 2026-10-06/run_007       # modulazione e Q di ogni punto -> analisi.csv + due heatmap
    python analisi.py punto 2026-10-06/run_007 12  # grafico di controllo di un punto
    python analisi.py misure                       # elenco delle misure singole
    python analisi.py confronto 1 2 3 --rif 1      # confronta misure singole e le divide per la misura 1

Metodo (come movemin.py): media mobile di CH2 (np.convolve), dentro la rampa di CH1.
 - modulazione = max - min della media mobile nella rampa;
 - Q di ogni campana = altezza (prominenza) / larghezza a meta' altezza, in volt per campione.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks, peak_widths

DATI = Path(__file__).resolve().parent.parent / "Chiara_Bru" / "runs_oscilloscope"

FINESTRA = 10000            # lunghezza della media mobile, in campioni
V_INIZIO = 0.9              # la rampa di CH1 inizia quando CH1 supera questo valore...
TOLLERANZA_FINE = 0.3       # ...e finisce quando CH1 e' a meno di questo dal suo massimo
SOGLIA_PROMINENZA = 0.1     # una campana conta se alta almeno questa frazione dell'escursione (modulazione)...
SNR_MIN = 10                # ...e almeno questo multiplo del rumore residuo della media mobile
SOGLIA_RUN = 0.1            # ...e almeno questa frazione della campana piu' alta di tutto il run (senza risonanza = n.d.)
UNITA_Q = 1e-6              # Q mostrata in unita' di V per campione / UNITA_Q
LIVELLO_BUIO = 0.0          # livello di CH2 senza luce (confronto di misure singole)


# ---------------------------------------------------------------- lettura e metodo di base

def percorso_run(cartella):
    p = Path(cartella)
    return p if p.is_absolute() or p.exists() else DATI / p


def leggi_csv(percorso):
    df = pd.read_csv(percorso, dtype=np.float32)
    return df["CH1"].to_numpy(float), df["CH2"].to_numpy(float)


def media_mobile(y, n=FINESTRA):
    """Identica a movemin.py: finestra resa dispari, np.convolve in modo 'same'."""
    n = int(n) | 1
    return np.convolve(y, np.ones(n) / n, mode="same")


def rampa(ch1):
    """(inizio, fine) della rampa di CH1 in campioni, o None."""
    sopra = np.flatnonzero(ch1 > V_INIZIO)
    if sopra.size == 0:
        return None
    i0 = sopra[0]
    i1 = i0 + np.flatnonzero(ch1[i0:] >= ch1.max() - TOLLERANZA_FINE)[0]
    return (i0, i1) if i1 > i0 else None


def zona_valida(ch1):
    """La rampa senza FINESTRA/2 campioni ai bordi del record, dove la media mobile usa zeri."""
    r = rampa(ch1)
    if r is None:
        return None
    mezza = (int(FINESTRA) | 1) // 2
    a, b = max(r[0], mezza), min(r[1], len(ch1) - 1 - mezza)
    return (a, b) if b > a else None


def analizza_punto(ch1, ch2):
    """(modulazione in V, campane) di un punto. Ogni campana: dict con indice, altezza, larghezza, Q, sx, dx, livello."""
    zona = zona_valida(ch1)
    if zona is None:
        return np.nan, []
    a, b = zona
    curva = media_mobile(ch2)[a:b + 1]
    modulazione = curva.max() - curva.min()
    d = np.diff(ch2[a:b + 1])
    sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.sqrt(FINESTRA)   # rumore residuo della media mobile
    picchi, prop = find_peaks(curva, prominence=max(SOGLIA_PROMINENZA * modulazione, SNR_MIN * sigma), distance=FINESTRA)
    if len(picchi) == 0:
        return modulazione, []
    larg, livello, sx, dx = peak_widths(curva, picchi, rel_height=0.5,
                                        prominence_data=(prop["prominences"], prop["left_bases"], prop["right_bases"]))
    return modulazione, [{"indice": a + int(p), "altezza": float(h), "larghezza": float(w), "Q": float(h / w),
                          "sx": a + float(s), "dx": a + float(d_), "livello": float(l)}
                         for p, h, w, l, s, d_ in zip(picchi, prop["prominences"], larg, livello, sx, dx)]


# ---------------------------------------------------------------- analisi di un run

def analizza_run(cartella, ricalcola=False):
    """Modulazione e Q di ogni punto del run (cartella con griglia.csv). Tabella salvata in analisi.csv."""
    cartella = percorso_run(cartella)
    if (cartella / "analisi.csv").exists() and not ricalcola:
        return pd.read_csv(cartella / "analisi.csv")
    punti = []
    for r in pd.read_csv(cartella / "griglia.csv").itertuples():
        if (cartella / r.file).exists():
            modulazione, campane = analizza_punto(*leggi_csv(cartella / r.file))
            print(f"punto {r.indice}: modulazione {modulazione:.4f} V, {len(campane)} campane")
            punti.append((r, modulazione, campane))
    soglia = SOGLIA_RUN * max([c["altezza"] for _, _, campane in punti for c in campane], default=0)
    righe = []
    for r, modulazione, campane in punti:
        c = max([c for c in campane if c["altezza"] >= soglia], key=lambda c: c["altezza"], default=None)   # la campana piu' alta
        righe.append({"indice": r.indice, "x_mm": r.x_mm, "y_mm": r.y_mm, "modulazione_V": modulazione,
                      "Q": c["Q"] if c else np.nan, "altezza_V": c["altezza"] if c else np.nan,
                      "larghezza_campioni": c["larghezza"] if c else np.nan})
    tabella = pd.DataFrame(righe)
    tabella.to_csv(cartella / "analisi.csv", index=False)
    return tabella


def reticolo(x, y):
    """(riga, colonna) di ogni punto di un percorso a meandro; None se non e' ricostruibile."""
    salti = np.abs(np.diff(y))
    if len(x) < 4 or salti.max() == 0:
        return None
    righe = np.concatenate([[0], np.cumsum(salti > 0.5 * salti.max())])         # il passaggio di riga e' un salto grande in Y
    lunghezze = np.bincount(righe)
    if len(set(lunghezze)) != 1:
        return None
    colonne = np.concatenate([np.arange(n) if r % 2 == 0 else np.arange(n)[::-1] for r, n in enumerate(lunghezze)])
    return righe, colonne


def mappa(tabella, colonna, etichetta, file=None, scala=1.0, formato="{:.2f}"):
    """Heatmap con il valore scritto su ogni punto: a reticolo (righe/colonne del meandro) o alle posizioni X, Y reali."""
    x, y = tabella["x_mm"].to_numpy(), tabella["y_mm"].to_numpy()
    z = tabella[colonna].to_numpy() / scala
    valido = np.isfinite(z)
    if not valido.any():
        print(f"{colonna}: nessun valore valido, mappa non disegnata")
        return
    norma, cmap = Normalize(z[valido].min(), z[valido].max()), plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(9, 7))
    rc = reticolo(x, y)
    if rc is not None:
        righe, colonne = rc
        m = np.full((righe.max() + 1, colonne.max() + 1), np.nan)
        m[righe, colonne] = z
        ax.imshow(np.ma.masked_invalid(m), cmap=cmap, norm=norma, origin="upper")
        cx, cy = colonne, righe
        ax.set_xlabel("colonna (nel verso della prima riga)")
        ax.set_ylabel("riga (dall'alto)")
    else:
        d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
        passo = d[d > 0].min() if len(x) > 1 else 0.3
        for xi, yi, zi in zip(x, y, z):
            ax.add_patch(Rectangle((xi - passo / 2, yi - passo / 2), passo, passo, edgecolor="white",
                                   facecolor=cmap(norma(zi)) if np.isfinite(zi) else "lightgray"))
        ax.set_xlim(x.min() - passo, x.max() + passo)
        ax.set_ylim(y.min() - passo, y.max() + passo)
        ax.set_aspect("equal")
        cx, cy = x, y
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
    for xi, yi, zi in zip(cx, cy, z):
        ax.text(xi, yi, formato.format(zi) if np.isfinite(zi) else "n.d.", ha="center", va="center", fontsize=9,
                color="gray" if not np.isfinite(zi) else "black" if norma(zi) > 0.55 else "white")
    fig.colorbar(plt.cm.ScalarMappable(norm=norma, cmap=cmap), ax=ax, label=etichetta)
    ax.set_title(etichetta)
    plt.tight_layout()
    if file:
        fig.savefig(file, dpi=200)
        print("Salvata", file)
    plt.show()


def heatmap_run(cartella, ricalcola=False):
    """Heatmap della modulazione e del quality factor di un run."""
    cartella = percorso_run(cartella)
    t = analizza_run(cartella, ricalcola)
    mappa(t, "modulazione_V", "ampiezza di modulazione: max - min della media mobile di CH2 nella rampa (V)",
          cartella / "heatmap_modulazione.png", formato="{:.3f}")
    mappa(t, "Q", f"quality factor Q = altezza / larghezza  [{UNITA_Q:g} V per campione]",
          cartella / "heatmap_quality_factor.png", scala=UNITA_Q)


def grafico_punto(cartella, indice):
    """CH2 grezzo, media mobile, zona di misura e campane trovate per un punto."""
    cartella = percorso_run(cartella)
    r = pd.read_csv(cartella / "griglia.csv").query("indice == @indice").iloc[0]
    ch1, ch2 = leggi_csv(cartella / r.file)
    modulazione, campane = analizza_punto(ch1, ch2)
    passo = max(1, len(ch2) // 20000)
    asse = np.arange(len(ch2))[::passo]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(asse, ch2[::passo], color="gray", lw=0.6, alpha=0.6, label="CH2 grezzo")
    liscia = media_mobile(ch2)
    ax.plot(asse, liscia[::passo], color="red", lw=1.5, label=f"media mobile (N={FINESTRA})")
    zona = zona_valida(ch1)
    if zona:
        ax.axvspan(*zona, color="orange", alpha=0.12, label="zona di misura (rampa di CH1)")
    for n, c in enumerate(campane, start=1):
        ax.hlines(c["livello"], c["sx"], c["dx"], color="blue", lw=2)
        ax.annotate(f"{n}: Q={c['Q'] / UNITA_Q:.3g}", (c["indice"], liscia[c["indice"]]), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=9)
    ax.set_xlabel("campione")
    ax.set_ylabel("CH2 (V)")
    ax.set_title(f"Punto {indice} (X={r.x_mm}, Y={r.y_mm}): modulazione {modulazione:.4f} V, {len(campane)} campane")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------- misure singole

def misure():
    return pd.read_csv(DATI / "misure_singole" / "misure.csv", keep_default_na=False)


def spettro(ch1, ch2, n_punti=5000):
    """CH2 nella rampa di CH1 come (tensione di CH1, mediana a blocchi di CH2)."""
    r = rampa(ch1)
    i0, i1 = r
    m, q = np.polyfit(np.arange(i0, i1 + 1), ch1[i0:i1 + 1], 1)               # CH1 e' quantizzato: si usa la retta fittata
    n = max(1, (i1 - i0 + 1) // n_punti)
    k = (i1 - i0 + 1) // n
    return m * (i0 + (np.arange(k) + 0.5) * n - 0.5) + q, np.median(ch2[i0:i0 + k * n].reshape(k, n), axis=1) - LIVELLO_BUIO


def confronta(ids=None, rif=None):
    """Confronta misure singole (stesso tipo del riferimento) e le divide per la misura rif (default: la prima)."""
    m = misure()
    m = m[m["id"].isin(ids)] if ids else m
    rif = rif or int(m["id"].iloc[0])
    con_scansione = str(m[m["id"] == rif]["con_scansione"].iloc[0]) == "True"
    m = m[m["con_scansione"].astype(str) == str(con_scansione)]
    dati = {r.id: leggi_csv(DATI / "misure_singole" / r.file) for r in m.itertuples()}
    etichette = {r.id: f"{r.id} {r.nome}" for r in m.itertuples()}
    if con_scansione:
        curve = {i: spettro(*d) for i, d in dati.items()}
        v = np.linspace(max(c[0][0] for c in curve.values()), min(c[0][-1] for c in curve.values()), 2000)
        s = {i: np.interp(v, *c) for i, c in curve.items()}
        for i in s:
            plt.plot(v, s[i], lw=2.2 if i == rif else 1.2, label=etichette[i])
        plt.xlabel("CH1 (V)")
        plt.ylabel("CH2 (V)")
        plt.legend()
        plt.tight_layout()
        valori = {i: s[i].mean() for i in s}
        rapporti = {i: (s[i] / s[rif]).mean() for i in s}
        errori = {i: np.nan for i in s}
    else:
        liv, err = {}, {}
        for i, (_, ch2) in dati.items():
            blocchi = np.median(ch2[:len(ch2) // 1000 * 1000].reshape(-1, 1000), axis=1)
            liv[i], err[i] = np.median(ch2) - LIVELLO_BUIO, blocchi.std(ddof=1) / np.sqrt(len(blocchi))
        valori = liv
        rapporti = {i: liv[i] / liv[rif] for i in liv}
        errori = {i: abs(rapporti[i]) * np.hypot(err[i] / liv[i], err[rif] / liv[rif]) for i in liv}
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
        x = np.arange(len(liv))
        a1.bar(x, [liv[i] for i in liv], yerr=[err[i] for i in liv], capsize=3)
        a2.bar(x, [rapporti[i] for i in liv], yerr=[errori[i] for i in liv], capsize=3)
        for a, t in ((a1, "livello di CH2 (V)"), (a2, f"rapporto con la misura {rif}")):
            a.set_xticks(x)
            a.set_xticklabels([etichette[i] for i in liv], rotation=30, ha="right")
            a.set_ylabel(t)
        a2.axhline(1, color="gray", lw=0.8, ls=":")
        plt.tight_layout()
    tabella = pd.DataFrame({"id": list(valori), "nome": [etichette[i] for i in valori], "valore_V": list(valori.values()),
                            f"rapporto_con_{rif}": list(rapporti.values()), "errore_rapporto": list(errori.values()),
                            "note": [m[m["id"] == i]["note"].iloc[0] for i in valori]})
    print(tabella.to_string(index=False))
    plt.show()
    return tabella


# ---------------------------------------------------------------- da terminale

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="comando", required=True)
    a = sub.add_parser("run", help="modulazione, Q e heatmap di un run")
    a.add_argument("cartella")
    a.add_argument("--ricalcola", action="store_true", help="ignora analisi.csv e rilegge tutti i CSV")
    a = sub.add_parser("punto", help="grafico di controllo di un punto del run")
    a.add_argument("cartella")
    a.add_argument("indice", type=int)
    sub.add_parser("misure", help="elenco delle misure singole")
    a = sub.add_parser("confronto", help="confronto di misure singole")
    a.add_argument("ids", nargs="*", type=int, help="id delle misure (vuoto = tutte)")
    a.add_argument("--rif", type=int, help="id della misura di riferimento per i rapporti")
    args = p.parse_args()
    if args.comando == "run":
        heatmap_run(args.cartella, args.ricalcola)
    elif args.comando == "punto":
        grafico_punto(args.cartella, args.indice)
    elif args.comando == "misure":
        print(misure().to_string(index=False))
    else:
        confronta(args.ids, args.rif)


if __name__ == "__main__":
    main()
