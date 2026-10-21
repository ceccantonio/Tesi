Transmission references, measured ONCE and used every time (analysis.py and inverse_design.py).

Put here the two scans WITHOUT the sample, acquired like any other point (CSV with CH1, CH2 and the laser ramp):

    copol.csv     reference for the co-pol runs
    crosspol.csv  reference for the cross-pol runs

In codice/analysis.py set  POLARIZATION = "co"  or  "cross"  (according to the run you are looking at): every measured
spectrum is then divided by the reference of that polarization, at the same wavelength, and you get the transmission T
like in the simulation. With POLARIZATION = ""  nothing is divided.

The first time a reference is used, its spectrum is saved next to it (copol.npz, crosspol.npz): the next times the big CSV
is not read again. If you measure a new reference, replace the CSV: the .npz is recomputed by itself.

Note: the .csv files are not uploaded to GitHub (too big, see .gitignore); keep a copy of them elsewhere.
