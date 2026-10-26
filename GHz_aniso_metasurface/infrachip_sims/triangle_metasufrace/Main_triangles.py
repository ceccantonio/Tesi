# Import the modules
import pya                        # KLayout python API
import numpy as np                # numpy module
from pathlib import Path
import math

from Python_triangles import create_triangle  # Import the function

# Parameters of the triangles 
cx, cy = 0, 0 # is in the middle of the bottom side or O of the paper
a1, a2 = 388e-3, 388e-3 # lenght of the bottom and right side
theta1 = 60 #angle between a1 and a2 it has to be less thatn 90!

# Find the center of the triangles 
paper_tet2 = 96.5
paper_tet1 = 120

alpha2 = paper_tet2 - 90
alpha1 = 180 - paper_tet1 - alpha2 

d1 = 70e-3
d2 = 50e-3
L = 265e-3

A = pya.DPoint(cx - a1 / 2 ,cy - a2 / 2 * math.sin(math.radians(theta1))) #center of 1 triangle 
B = pya.DPoint(cx - d2 , cy + d1) #new origin O' of the paper
slope = (B.y - A.y) / (B.x - A.x)
    
slope_rad = math.atan(slope)
slope_deg = math.degrees(slope_rad)

gamma1 = 90 - slope_deg
gamma1_rad = math.radians(gamma1)

alpha1 = paper_tet2 - 90
alpha1_rad = math.radians(alpha1)

alpha2 = 180 - paper_tet1 - alpha1
alpha2_rad = math.radians(alpha2)

#centers of the 3 triangles 

C1 = [-d2 -(L * math.sin(gamma1_rad)), d1 -(L * math.cos(gamma1_rad))]
C2 = [-d2 + L * math.cos(alpha1_rad), d1 + L * math.sin(alpha1_rad)]
C3 = [-d2 -(L * math.cos(alpha2_rad)), d1 + L * math.sin(alpha2_rad)]

# Target file name
fn = 'Infrachip_triangle_test.gds'

dbu = 1e-3 #data base unit

# Create the basic layout objects
layout = pya.Layout()
layout.dbu = dbu


# Create a single cell to hold all designs
cell = layout.create_cell('TrianglesDesign')
layer = layout.layer(1, 0)

# Loop through different centers and generate multiple geometries
centers = [(C1[0], C1[1]), (C2[0], C2[1]), (C3[0], C3[1])]   


for cx, cy in centers:
    create_triangle(cell, layer, cx, cy, a1, a2, theta1)

# Save the layout to a GDS file
layout.write(fn)
print(f"Design saved to {fn}")


###########################################################################################################
# Function to determine desktop folder for saving
def desktop_folder():
    folder = Path.home() / 'Desktop'
    if not folder.exists():
        for item in Path.home().iterdir():
            if item.is_dir() and item.name.startswith('OneDrive '):
                folder = item / 'Desktop'
                break
    return folder

# Save the layout to a file
filename = str(desktop_folder() / 'università' / 'assegno_di_ricerca' / 'Infrachip' / 'CAD' / fn).replace('\\','\\\\')
