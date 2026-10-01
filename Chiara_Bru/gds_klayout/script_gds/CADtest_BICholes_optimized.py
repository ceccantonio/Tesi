import gdstk
import math

# parameters (microns)
chip_size = 900.0  # chip size (µm)

matrix_sub = 5
submatrix_size = 100.0  # 100x100 µm writefield

spacing = 200  # distance between submatrix centers

# photonic patterns
radius_start = 0.1
radius_step = 0.025

pitch_start = 0.8
pitch_step = 0.1

# matrix offset from chip edge
matrix_offset_x = 0
matrix_offset_y = 0

# FIXED number of holes per direction
n_holes = 50

# library and single cell
lib = gdstk.Library()
top = lib.new_cell("TOP")

matrix_layer = 1
marker_layer = 10
chip_layer = 20

# chip area
top.add(
    gdstk.rectangle((0, 0), (chip_size, chip_size), layer=chip_layer)
)

# circle using ellipse with fine tolerance
def circle(cx, cy, r):
    return gdstk.ellipse(
        (cx, cy),
        r,
        layer=matrix_layer,
        tolerance=0.001
    )

# parameter arrays
radius_values = []
pitch_values = []

for i in range(matrix_sub):
    radius_values.append(radius_start + i * radius_step)
    pitch_values.append(pitch_start + i * pitch_step)

# big matrix directly in TOP
for row in range(matrix_sub):

    pitch = pitch_values[row]

    for col in range(matrix_sub):

        radius = radius_values[col]

        # total size occupied by 50 holes
        total_size = (n_holes - 1) * pitch

        # safety check (optional but recommended)
        if total_size > submatrix_size:
            raise ValueError(
                f"Pitch too large for 50 holes: pitch={pitch}, total_size={total_size}"
            )

        # center the pattern in the 100x100 µm writefield
        offset = total_size / 2

        base_x = matrix_offset_x + col * spacing
        base_y = matrix_offset_y + row * spacing

        for i in range(n_holes):
            for j in range(n_holes):

                x = base_x + (i * pitch - offset + submatrix_size / 2)
                y = base_y + (j * pitch - offset + submatrix_size / 2)

                top.add(circle(x, y, radius))

'''
# bottom-left marker (optional)
def create_marker(x0, y0, size=100):

    half = size / 2

    r1 = gdstk.rectangle(
        (x0, y0),
        (x0 + half, y0 + half),
        layer=marker_layer,
    )

    r2 = gdstk.rectangle(
        (x0 + half, y0 + half),
        (x0 + size, y0 + size),
        layer=marker_layer,
    )

    top.add(r1, r2)

create_marker(500, 500, 100)
'''

# save GDS
lib.write_gds(r"C:\Users\Labo\Documents\BIG_MATRIX_50x50holes.gds")

print("GDS generated")