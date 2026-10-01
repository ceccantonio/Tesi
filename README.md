# Tesi

- `logbook/logbook.ipynb` — diario di laboratorio da compilare in VSCode, in formato notebook. Contiene una cella template da copiare per ogni nuova giornata; puoi aggiungere celle di codice sotto una voce per un controllo rapido o un plot legato a quella giornata.
- `scripts/plot_template.py` / `scripts/plot_template.ipynb` — template generico per caricare dati con `np.loadtxt` e fare un plot con barre d'errore (stesso contenuto, in versione script e notebook).
- `scripts/fit_template.py` / `scripts/fit_template.ipynb` — template generico per fit con `scipy.optimize.curve_fit`: carica i dati, fitta, stampa parametri/incertezze/chi-quadro ridotto, plotta dati+fit+residui.
- `scripts/data/esempio.txt` — dati di esempio (x, y, sigma_y) per testare subito script e notebook senza modifiche.

## Setup

```bash
pip install -r requirements.txt
```

`ipykernel` serve per eseguire i notebook `.ipynb` dentro VSCode (estensione Jupyter): apri il file e seleziona come kernel l'interprete Python in cui hai fatto l'install.

## Uso

Per ogni nuovo set di dati:
1. copia il tuo file dati in `scripts/data/`
2. duplica il template (`.py` o `.ipynb`, quello che preferisci) con un nome descrittivo (es. `fit_michelson.ipynb`)
3. modifica solo la sezione/cella `CONFIGURAZIONE` in cima (percorso file, colonne, etichette) e, per il fit, la funzione `modello()` e `P0`
