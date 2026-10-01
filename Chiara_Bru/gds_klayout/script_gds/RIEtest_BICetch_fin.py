"""
EBL GDS Generator — gdstk version
==================================
Single column of 5 submatrices, 80×80 holes each, pitch 1 µm.
Radii centered on 150 nm with ±5% and ±10% steps.
Order from bottom (closest to marker) to top: largest → smallest radius.

Writefield layout (FIELD_SIZE = 100 µm):
  - Field 0 (  0→100 µm): marker ring only, centred at (50, 50)
  - Field 1 (100→200 µm): R=165.0 nm, centred at (50, 150)
  - Field 2 (200→300 µm): R=157.5 nm, centred at (50, 250)
  - Field 3 (300→400 µm): R=150.0 nm, centred at (50, 350)
  - Field 4 (400→500 µm): R=142.5 nm, centred at (50, 450)
  - Field 5 (500→600 µm): R=135.0 nm, centred at (50, 550)

Total working area bottom-left at (0, 0).
All units in µm. Requires: gdstk  →  pip install gdstk
"""

import math
import gdstk

# ── Parameters ───────────────────────────────────────────────────────────────

BASE_R = 0.150   # µm — base radius (150 nm)
STEPS  = [+0.10, +0.05, 0.0, -0.05, -0.10]   # k=0 (bottom) → k=4 (top)

SUBMATRICES = [
    (BASE_R * (1 + s), 1.0, 80)
    for s in STEPS
]

N_VERTS        = 64
LAYER_HOLES    = 1
LAYER_MARKER   = 2
DATATYPE       = 0

FIELD_SIZE     = 100.0   # µm — EBL writefield size
MARKER_R_IN    = 15.0    # µm
MARKER_R_OUT   = 20.0    # µm

OUTPUT_FILE    = "ebl_holes_150nm.gds"

# ── Helpers ───────────────────────────────────────────────────────────────────

def circle_points(cx, cy, radius, n_verts):
    return [
        (cx + radius * math.cos(2 * math.pi * i / n_verts),
         cy + radius * math.sin(2 * math.pi * i / n_verts))
        for i in range(n_verts)
    ]

def submatrix_extent(pitch, grid_n):
    return (grid_n - 1) * pitch

# ── Build GDS ─────────────────────────────────────────────────────────────────

lib  = gdstk.Library()
cell = lib.new_cell("EBL_150NM")

# All submatrices have the same extent (same pitch and grid_n)
extent   = submatrix_extent(1.0, 80)   # = 79 µm
centre_x = FIELD_SIZE / 2.0 + 100           # = 50 µm — shared X centre for all fields

# Place each submatrix centred in its own writefield (k+1) * FIELD_SIZE
for k, (radius, pitch, grid_n) in enumerate(SUBMATRICES):
    field_centre_y = (k + 1) * FIELD_SIZE + FIELD_SIZE / 2.0  # 150, 250, 350, 450, 550
    origin_x = centre_x - extent / 2.0
    origin_y = field_centre_y - extent / 2.0

    for row in range(grid_n):
        for col in range(grid_n):
            cx = origin_x + col * pitch
            cy = origin_y + row * pitch
            pts = circle_points(cx, cy, radius, N_VERTS)
            cell.add(gdstk.Polygon(pts, layer=LAYER_HOLES, datatype=DATATYPE))

# Marker ring centred in field 0 at (50, 50)
marker_cx = FIELD_SIZE / 2.0
marker_cy = FIELD_SIZE / 2.0
outer_pts  = circle_points(marker_cx, marker_cy, MARKER_R_OUT, N_VERTS * 2)
inner_pts  = circle_points(marker_cx, marker_cy, MARKER_R_IN,  N_VERTS * 2)
outer_poly = gdstk.Polygon(outer_pts, layer=LAYER_MARKER, datatype=DATATYPE)
inner_poly = gdstk.Polygon(inner_pts, layer=LAYER_MARKER, datatype=DATATYPE)
ring = gdstk.boolean(outer_poly, inner_poly, "not", layer=LAYER_MARKER, datatype=DATATYPE)
for poly in ring:
    cell.add(poly)

# ── Write file ────────────────────────────────────────────────────────────────
lib.write_gds(OUTPUT_FILE)
total = sum(n * n for _, _, n in SUBMATRICES)
print(f"✓  GDS file written : {OUTPUT_FILE}")
print(f"   Total holes      : {total}")
print(f"   Marker           : inner={MARKER_R_IN} µm  outer={MARKER_R_OUT} µm  at ({marker_cx},{marker_cy})")
print(f"   Field size       : {FIELD_SIZE} µm")
print()
print(f"   {'Field':>6}  {'Y centre':>10}  {'Radius':>10}  {'Grid':>8}")
print(f"   {'0':>6}  {'50':>9} µm  {'marker':>10}")
for k, (r, p, n) in enumerate(SUBMATRICES):
    fc = (k + 1) * FIELD_SIZE + FIELD_SIZE / 2.0
    print(f"   {k+1:>6}  {fc:>9.4g} µm  {r*1000:>8.4g} nm  {n:>4}×{n:<4}")