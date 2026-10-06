import numpy as np
import matplotlib

import matplotlib.pyplot as plt
import pandas as pd

csv_path = 'Chiara_Bru\\runs_oscilloscope\\2026-10-06\\run_007\\punto_011_x91.818_y146.526.csv'
#out_path = 'C:\\Users\\Labo\\Desktop\\Tesi\\Chiara_Bru\\runs_oscilloscope\\2026-10-06\\run_007\\punto_015_x91.789_y146.526.png'

def moving_average(values, window_length):
    """Media mobile con finestra di lunghezza variabile."""
    if window_length <= 0:
        raise ValueError('window_length deve essere > 0')

    window_length = int(window_length)
    if window_length % 2 == 0:
        window_length += 1

    kernel = np.ones(window_length) / window_length
    return np.convolve(values, kernel, mode='same')


df = pd.read_csv(csv_path)
y = df['CH2'].to_numpy(dtype=float)
window_length = 10000  # puoi cambiare questa lunghezza
smoothed = moving_average(y, window_length)
curve=smoothed[5000:60000]
modulation = np.max(curve) - np.min(curve)

plt.figure(figsize=(10, 5))
plt.plot(y, linewidth=1.0, color='gray', alpha=0.6, label='raw')
plt.plot(smoothed, linewidth=1.8, color='red', label=f'media mobile (N={window_length})')
plt.title('CH2 signal')
plt.xlabel('Sample index')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
print(modulation)
plt.show()
