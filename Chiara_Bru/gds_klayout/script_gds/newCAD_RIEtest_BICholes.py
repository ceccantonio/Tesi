import gdstk
import numpy as np


#Parameters

#holes with r from 100nm to 2um
radii = np.array([
    2.00, 1.50, 1.00, 0.75,
    0.50, 0.20, 0.15, 0.10
])

num_columns = 550

pitch = 4.0
column_spacing = 4.0

#pattern angle
theta = np.deg2rad(45)

ellipse_tolerance = 0.01

#gap between mirrored sequences
mirror_gap = 1.0

#number of repeated sequences
num_repeats = 35


#Setup

lib = gdstk.Library()
cell = lib.new_cell("ALTERNATING_MIRROR_PATTERN")

vx, vy = np.cos(theta), np.sin(theta)
px, py = -np.sin(theta), np.cos(theta)


#Generate

for col in range(num_columns):

    base_x = col * column_spacing * px
    base_y = col * column_spacing * py

   
    #alternating column orientation

    if col % 2 == 0:
        first_half = radii                # big to small
        second_half = radii[::-1]         # small to big
    else:
        first_half = radii[::-1]          # small to big
        second_half = radii               # big to small

    half_len = len(first_half)

    #total size of one repeated pair
    motif_length = (
        half_len * pitch
        + mirror_gap
        + half_len * pitch
        + mirror_gap
    )


    #Repeat

    for repeat in range(num_repeats):

        repeat_offset = repeat * motif_length


        #First sequence

        for row, r in enumerate(first_half):

            d = repeat_offset + row * pitch

            x = base_x + d * vx
            y = base_y + d * vy

            hole = gdstk.ellipse(
                center=(x, y),
                radius=r,
                tolerance=ellipse_tolerance,
                layer=1,
                datatype=0,
            )

            cell.add(hole)

        
        #Second sequence

        second_offset = (
            repeat_offset
            + half_len * pitch
            + mirror_gap
        )

        for row, r in enumerate(second_half):

            d = second_offset + row * pitch

            x = base_x + d * vx
            y = base_y + d * vy

            hole = gdstk.ellipse(
                center=(x, y),
                radius=r,
                tolerance=ellipse_tolerance,
                layer=1,
                datatype=0,
            )

            cell.add(hole)



#Marker ring (bottom left)

marker_x = -1540
marker_y = -40

ring_radius = 50.0     # center radius
ring_width = 5.0       # thickness

ring = gdstk.RobustPath(
    (marker_x + ring_radius, marker_y),
    ring_width,
    layer=2,
    datatype=0,
)

ring.arc(
    ring_radius,
    0,
    2 * np.pi,
)

cell.add(ring)


#CAD shift

dx = 1600
dy = 100

for polygon in cell.polygons:
    polygon.translate(dx, dy)

for path in cell.paths:
    path.translate(dx, dy)

for label in cell.labels:
    label.translate(dx, dy)

for reference in cell.references:
    reference.translate(dx, dy)


#Save

lib.write_gds("alternating_mirror_pattern.gds")

print("Generated: alternating_mirror_pattern.gds")