import gdstk
import math

# parameters (microns)
chip_size = 5000.0  # 5x5 mm

matrix_sub = 7
submatrix_size = 100.0  # WF 100x100um

holes = 50

spacing = 300  # submatrix centers distance

# photonic patterns
radius_start = 0.06
radius_step = 0.03

pitch_start = 0.7
pitch_step = 0.1

circle_pts = 32  # now used as tolerance control

# matrix offset from chip edge
matrix_offset_x = 500
matrix_offset_y = 500


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


# circle using ellipse
def circle(cx, cy, r):

    return gdstk.ellipse(
        (cx, cy),
        r,
        layer=matrix_layer,
        tolerance=0.0001  #0.001 controls smoothness (smaller = smoother)
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

        offset = (holes - 1) * pitch / 2

        base_x = matrix_offset_x + col * spacing
        base_y = matrix_offset_y + row * spacing

        for i in range(holes):
            for j in range(holes):

                x = base_x + (i * pitch - offset + submatrix_size / 2)
                y = base_y + (j * pitch - offset + submatrix_size / 2)

                top.add(circle(x, y, radius))


# bottom-left marker
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


create_marker(300, 300, 100)


# save GDS
lib.write_gds(r"C:\Users\Labo\Documents\BIG_MATRIX_ellipse.gds")

print("GDS generated")