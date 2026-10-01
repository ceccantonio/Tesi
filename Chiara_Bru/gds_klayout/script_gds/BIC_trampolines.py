print(">>> This macro is running!")

import gdstk
import math
import os
import numpy as np

# Library
lib = gdstk.Library()
top_cell = lib.new_cell("trampoline")

# Layers
LAYER_FRAME = (100, 0)
LAYER_TRENCH = (1, 0)
LAYER_HOLE = (3, 0)

# Parameters
chipS = 10000 / 2
Ns = 50

Nx = 2
Ny = 3
pitchx = 3500
pitchy = 1800

F = 800 / 2
L = 150
R_pad = 80
tw = [9, 10, 11, 9, 10, 11]
R_tet = 10

sf = [0.9, 0.95, 1, 1, 1.05, 1.1]
rad = 0.15
theta = np.linspace(0, 2 * np.pi, Ns)
a = 1
NSR = [80, 80, 80, 80, 80, 80]

# Chip outline
chip = [
    (-chipS, -chipS),
    (-chipS, chipS),
    (chipS, chipS),
    (chipS, -chipS),
]
top_cell.add(gdstk.Polygon(chip, layer=LAYER_FRAME[0], datatype=LAYER_FRAME[1]))

# Loop devices
for k in range(Nx):
    for kk in range(Ny):

        m_ind = kk + Ny * k
        x0 = -pitchx * (Nx - 1) / 2 + pitchx * k
        y0 = -pitchy * (Ny - 1) / 2 + pitchy * kk

        alpha = math.radians(45)

        the1 = np.linspace(math.pi, math.pi - alpha, Ns)
        the2 = -np.flip(the1)
        bet1 = np.linspace(math.pi - alpha, 0, Ns)
        bet2 = -np.flip(bet1)

        b = tw[m_ind] / (2 * math.sin(math.pi / 2 - alpha))

        # --- geometry construction (ONE continuous polyline) ---
        poly = []

        # start
        poly.append((L / 2, 0))

        # pad arc 1
        xp1 = L / 2 + R_pad
        yp1 = xp1 + R_pad * math.cos(math.pi - alpha) - R_pad * math.sin(math.pi - alpha) - b

        for th in the1:
            poly.append((xp1 + R_pad * math.cos(th),
                         yp1 + R_pad * math.sin(th)))

        # transition point
        xc1 = F - R_tet
        yc1 = xc1 + R_tet * math.cos(math.pi - alpha) - R_tet * math.sin(math.pi - alpha) - b

        FF = math.hypot(
            xc1 + R_tet * math.cos(math.pi - alpha),
            yc1 + R_tet * math.sin(math.pi - alpha),
        )

        xt = FF * math.cos(alpha)
        yt = FF * math.sin(alpha) - b
        poly.append((xt, yt))

        # inner arc 1
        for th in bet1:
            poly.append((xc1 + R_tet * math.cos(th),
                         yc1 + R_tet * math.sin(th)))

        # bottom symmetry connector
        poly.append((F, -poly[-1][1]))

        xc2 = xc1
        yc2 = poly[-1][1]

        # inner arc 2
        for th in bet2:
            poly.append((xc2 + R_tet * math.cos(th),
                         yc2 + R_tet * math.sin(th)))

        xp2 = xp1
        yp2 = -yp1

        # pad arc 2
        for th in the2:
            poly.append((xp2 + R_pad * math.cos(th),
                         yp2 + R_pad * math.sin(th)))

        # rotate 4 arms
        coords = np.array(poly)

        for j in range(4):
            angle = j * math.pi / 2
            rot = []

            for x, y in coords:
                xr = x * math.cos(angle) - y * math.sin(angle)
                yr = x * math.sin(angle) + y * math.cos(angle)
                rot.append((x0 + xr, y0 + yr))

            top_cell.add(
                gdstk.Polygon(rot, layer=LAYER_TRENCH[0], datatype=LAYER_TRENCH[1])
            )

        # BIC holes only (unchanged)
        for jSR in range(NSR[m_ind]):
            xi = x0 - ((NSR[m_ind] - 1) / 2 * a - a * jSR)

            for jjSR in range(NSR[m_ind]):
                yi = y0 - ((NSR[m_ind] - 1) / 2 * a - a * jjSR)

                xc = xi + sf[m_ind] * rad * np.cos(theta)
                yc = yi + sf[m_ind] * rad * np.sin(theta)

                circle = list(zip(xc, yc))
                top_cell.add(
                    gdstk.Polygon(circle, layer=LAYER_HOLE[0], datatype=LAYER_HOLE[1])
                )

# Save in local folder
output_path = "BIC_membrane.gds"
lib.write_gds(output_path)

print(">>> Saved:", os.path.abspath(output_path))