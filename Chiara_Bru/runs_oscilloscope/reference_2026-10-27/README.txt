Transmission references for the NEW sample mounted on 2026-10-27 (measured ONCE and used every time).

Put here the two scans WITHOUT the sample, acquired like any other point (CSV with CH1, CH2 and the laser ramp):

    copol.csv     reference for the co-pol runs
    crosspol.csv  reference for the cross-pol runs

How to acquire them:  python single_measurements.py  (SCAN = True), one measurement for the co-pol set-up and one for the cross-pol set-up;
then copy the two files from Chiara_Bru/runs_oscilloscope/misure_singole/ into this folder and rename them copol.csv and crosspol.csv.

How to use them: in codice/analysis.py set
    REFERENCE_SET = "reference_2026-10-27"
    POLARIZATION  = "co"   (or "cross", according to the run you are looking at)
Every measured spectrum is then divided by the reference of that polarization, at the same wavelength, and you get the transmission T.
The first time a reference is used, its spectrum is saved next to it (copol.npz, crosspol.npz): the big CSV is not read again.

The references of the previous sample stay in the folder "reference" (REFERENCE_SET = "reference").

Note: the .csv files are not uploaded to GitHub (too big, see .gitignore); keep a copy of them elsewhere.
