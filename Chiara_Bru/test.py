import pandas as pd

f = r"C:\Users\Labo\Documents\Zurich Instruments\LabOne\WebServer\session_20260826_170822_01\I3_bouncingSpectrum_ARharm123_500mV_000\dev2233_demods_0_sample_00000.csv"

df = pd.read_csv(f, sep=None, engine="python")
print(df.columns.tolist())
print(df.dtypes)
print(df.head(10))