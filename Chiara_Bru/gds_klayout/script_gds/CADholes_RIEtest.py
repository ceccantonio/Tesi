"""
EBL GDS Generator — gdstk version

Single column of 5 submatrices, stacked vertically, centres aligned.
No rotation. Marker ring (inner r=15µm, outer r=20µm) at (20,20).
Pattern bounding box bottom-left offset from origin.

Submatrices in decreasing radius order:
  - R = 2.0 µm  →  pitch = 6 µm  →  10 × 10  holes
  - R = 1.0 µm  →  pitch = 5 µm  →  10 × 10  holes
  - R = 0.5 µm  →  pitch = 2 µm  →  20 × 20  holes
  - R = 0.3 µm  →  pitch = 1 µm  →  50 × 50  holes
  - R = 0.15 µm →  pitch = 1 µm  →  50 × 50  holes

Submatrix centres are spaced at multiples of FIELD_SIZE (100 µm)
to avoid EBL stitching artefacts at pattern boundaries.
"""

import math
import gdstk

# Parameters
#
#   Each row: (radius_µm, pitch_µm, grid_N)
#   Order: largest radius first (top submatrix) → smallest last (bottom)
#
SUBMATRICES = [
    (2.0,  6.0, 10),   # R=2.0 µm,  pitch=6 µm,  10×10
    (1.0,  5.0, 10),   # R=1.0 µm,  pitch=5 µm,  10×10
    (0.5,  2.0, 20),   # R=0.5 µm,  pitch=2 µm,  20×20
    (0.3,  1.0, 50),   # R=0.3 µm,  pitch=1 µm,  50×50
    (0.2, 1.0, 50),   # R=0.15 µm, pitch=1 µm,  50×50
    (0.1, 0.5, 70),   # R=0.15 µm, pitch=1 µm,  50×50
]

N_VERTS         = 64      # vertices per circle polygon
LAYER_HOLES     = 1
LAYER_MARKER    = 2
DATATYPE        = 0

FIELD_SIZE      = 100.0  # µm — EBL field size; submatrix centres spaced at multiples of this
MARKER_R_IN     = 15.0   # µm — inner radius of marker ring
MARKER_R_OUT    = 20.0   # µm — outer radius of marker ring
PATTERN_OFFSET  = 100.0   # µm — distance from origin to pattern bbox corner

OUTPUT_FILE     = "ebl_holes.gds"

# Functions

def circle_points(cx, cy, radius, n_verts):
    return [
        (cx + radius * math.cos(2 * math.pi * i / n_verts),
         cy + radius * math.sin(2 * math.pi * i / n_verts))
        for i in range(n_verts)
    ]

def submatrix_extent(pitch, grid_n):
    return (grid_n - 1) * pitch

# Build GDS

lib  = gdstk.Library()
cell = lib.new_cell("EBL_PATTERN")

# Centre X: align all submatrices on the same vertical axis
max_extent = max(submatrix_extent(p, n) for _, p, n in SUBMATRICES)
centre_x   = max_extent / 2.0

# Stack submatrices vertically with centres at k * FIELD_SIZE
all_polys = []

for k, (radius, pitch, grid_n) in enumerate(SUBMATRICES):
    extent   = submatrix_extent(pitch, grid_n)
    origin_x = centre_x - extent / 2.0        # centre-align on X
    origin_y = k * FIELD_SIZE - extent / 2.0  # centre at k * FIELD_SIZE

    for row in range(grid_n):
        for col in range(grid_n):
            cx = origin_x + col * pitch
            cy = origin_y + row * pitch
            pts = circle_points(cx, cy, radius, N_VERTS)
            all_polys.append(pts)

# Shift pattern so bbox bottom-left is at (PATTERN_OFFSET, PATTERN_OFFSET)
all_x = [x for pts in all_polys for x, y in pts]
all_y = [y for pts in all_polys for x, y in pts]
shift_x = PATTERN_OFFSET - min(all_x)
shift_y = PATTERN_OFFSET - min(all_y)

for pts in all_polys:
    shifted = [(x + shift_x, y + shift_y) for x, y in pts]
    cell.add(gdstk.Polygon(shifted, layer=LAYER_HOLES, datatype=DATATYPE))

# Marker ring at (20,20)
outer_pts  = circle_points(40, 40, MARKER_R_OUT, N_VERTS * 2)
inner_pts  = circle_points(40, 40, MARKER_R_IN,  N_VERTS * 2)
outer_poly = gdstk.Polygon(outer_pts, layer=LAYER_MARKER, datatype=DATATYPE)
inner_poly = gdstk.Polygon(inner_pts, layer=LAYER_MARKER, datatype=DATATYPE)
ring = gdstk.boolean(outer_poly, inner_poly, "not", layer=LAYER_MARKER, datatype=DATATYPE)
for poly in ring:
    cell.add(poly)

# Write file
lib.write_gds(OUTPUT_FILE)
total = sum(n * n for _, _, n in SUBMATRICES)
print(f"✓  GDS file written : {OUTPUT_FILE}")
print(f"   Total holes      : {total}")
print(f"   Marker           : inner={MARKER_R_IN} µm  outer={MARKER_R_OUT} µm  at (20,20)")
print(f"   Field size       : {FIELD_SIZE} µm  (submatrix centres at 0, 100, 200, 300, 400 µm)")
print()
print(f"   {'Radius':>10}  {'Pitch':>8}  {'Grid':>8}  {'Extent':>10}  {'Centre Y':>10}")
for k, (r, p, n) in enumerate(SUBMATRICES):
    print(f"   {r:>8.3g} µm  {p:>6.3g} µm  {n:>4}×{n:<4}  {submatrix_extent(p,n):>8.3g} µm  {k*FIELD_SIZE:>8.3g} µm")