import gdstk
import math

# parameters (microns)
chip_size = 1700.0  # working area

matrix_sub = 5
submatrix_size = 100.0  # 100x100 um writefield

spacing = 300  # submatrix centers distance

# photonic patterns
radius_start = 0.1
radius_step = 0.025

pitch_start = 0.8
pitch_step = 0.1

# matrix offset from chip edge
matrix_offset_x = 200
matrix_offset_y = 200

# library and TOP cell
lib = gdstk.Library()
top = lib.new_cell("TOP")

matrix_layer = 1
marker_layer = 10
chip_layer = 20

# chip area
top.add(
    gdstk.rectangle((0, 0), (chip_size, chip_size), layer=chip_layer)
)


# !CACHE (avoids duplicated identical cells)
hole_cache = {}
submatrix_cache = {}

# circle (unchanged)
def get_hole_cell(radius):
    key = round(radius, 6)

    if key not in hole_cache:
        cell = lib.new_cell(f"HOLE_R_{radius:.3f}")
        cell.add(
            gdstk.ellipse(
                (0, 0),
                radius,
                layer=matrix_layer,
                tolerance=0.001
            )
        )
        hole_cache[key] = cell

    return hole_cache[key]


def get_submatrix_cell(radius, pitch):
    key = (round(radius, 6), round(pitch, 6))

    if key in submatrix_cache:
        return submatrix_cache[key]

    name = f"SUB_R_{radius:.3f}_P_{pitch:.3f}"
    cell = lib.new_cell(name)

    hole_cell = get_hole_cell(radius)

    # same math as your code
    n_holes = int((submatrix_size - 2 * radius) / pitch) + 1
    total_size = (n_holes - 1) * pitch

    if total_size > submatrix_size:
        n_holes -= 1
        total_size = (n_holes - 1) * pitch

    offset = total_size / 2

    # !Repetition (instead of loops)
    ref = gdstk.Reference(
        hole_cell,
        (submatrix_size / 2 - offset, submatrix_size / 2 - offset),
    )

    ref.repetition = gdstk.Repetition(
        columns=n_holes,
        rows=n_holes,
        spacing=(pitch, pitch),
    )

    cell.add(ref)
    submatrix_cache[key] = cell

    return cell


# parameter arrays
radius_values = []
pitch_values = []

for i in range(matrix_sub):
    radius_values.append(radius_start + i * radius_step)
    pitch_values.append(pitch_start + i * pitch_step)



# Place submatrices in TOP
for row in range(matrix_sub):

    pitch = pitch_values[row]

    for col in range(matrix_sub):

        radius = radius_values[col]

        sub_cell = get_submatrix_cell(radius, pitch)

        base_x = matrix_offset_x + col * spacing
        base_y = matrix_offset_y + row * spacing

        top.add(gdstk.Reference(sub_cell, (base_x, base_y)))


# save GDS
lib.write_gds(r"C:\Users\Labo\Documents\BIG_MATRIX_hierarchical.gds")

print("GDS generated")