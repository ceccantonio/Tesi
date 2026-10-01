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
chipS = 3000

Ns = 64


# 6x6 trampoline matrix
Nx = 6
Ny = 6


# spacing adjusted to fit 500um writefield
pitchx = 500
pitchy = 500


F = 300 / 2
L = 100

R_pad = 80
R_tet = 10


tw = np.ones((Nx*Ny))*10


# photonic parameters
a_values = [0.9, 1.0, 1.1]        # um
r_values = [0.145, 0.150, 0.165] # um


NSR = 50


sf = 1
rad0 = 0.15


theta = np.linspace(0, 2*np.pi, Ns, endpoint=False)


# Chip outline
chip = [
    (0 , 0),
    (0 , chipS),
    (chipS, chipS),
    (chipS, 0)
]

top_cell.add(
    gdstk.Polygon(
        chip,
        layer=LAYER_FRAME[0],
        datatype=LAYER_FRAME[1]
    )
)



# LOOP 6x6 TRAMPOLINES

for k in range(Nx):

    for kk in range(Ny):


        # position
        x0 = chipS/2 - pitchx*(Nx-1)/2 + pitchx*k
        y0 = chipS/2 - pitchy*(Ny-1)/2 + pitchy*kk



        # 2x2 block parameter selection

        block_col = k//2     # pitch variation
        block_row = kk//2    # radius variation


        a = a_values[block_col]
        rad = r_values[block_row]



        alpha = math.radians(45)



        the1 = np.linspace(
            math.pi,
            math.pi-alpha,
            Ns
        )

        the2 = -np.flip(the1)


        bet1 = np.linspace(
            math.pi-alpha,
            0,
            Ns
        )

        bet2 = -np.flip(bet1)



        b = 15/(2*math.sin(math.pi/2-alpha))



        # trampoline geometry


        poly=[]


        poly.append((L/2,0))



        xp1=L/2+R_pad

        yp1 = (
            xp1
            +R_pad*math.cos(math.pi-alpha)
            -R_pad*math.sin(math.pi-alpha)
            -b
        )


        for th in the1:

            poly.append(
                (
                xp1+R_pad*math.cos(th),
                yp1+R_pad*math.sin(th)
                )
            )



        xc1=F-R_tet


        yc1 = (
            xc1
            +R_tet*math.cos(math.pi-alpha)
            -R_tet*math.sin(math.pi-alpha)
            -b
        )


        FF=math.hypot(
            xc1+R_tet*math.cos(math.pi-alpha),
            yc1+R_tet*math.sin(math.pi-alpha)
        )


        xt=FF*math.cos(alpha)
        yt=FF*math.sin(alpha)-b


        poly.append((xt,yt))



        for th in bet1:

            poly.append(
                (
                xc1+R_tet*math.cos(th),
                yc1+R_tet*math.sin(th)
                )
            )



        poly.append(
            (
            F,
            -poly[-1][1]
            )
        )
        



        xc2=xc1
        yc2=poly[-1][1]
        



        for th in bet2:

            poly.append(
                (
                xc2+R_tet*math.cos(th),
                yc2+R_tet*math.sin(th)
                )
            )



        xp2=xp1
        yp2=-yp1



        for th in the2:

            poly.append(
                (
                xp2+R_pad*math.cos(th),
                yp2+R_pad*math.sin(th)
                )
            )




        coords=np.array(poly)



        # four arms

        for j in range(4):

            angle=j*math.pi/2

            rot=[]


            for x,y in coords:

                xr=x*math.cos(angle)-y*math.sin(angle)

                yr=x*math.sin(angle)+y*math.cos(angle)


                rot.append(
                    (
                    x0+xr,
                    y0+yr
                    )
                )



            top_cell.add(
                gdstk.Polygon(
                    rot,
                    layer=LAYER_TRENCH[0],
                    datatype=LAYER_TRENCH[1]
                )
            )



        # BIC holes with variable a,r


        for jSR in range(NSR):


            xi = (
                x0
                -((NSR-1)/2*a - a*jSR)
            )


            for jjSR in range(NSR):


                yi = (
                    y0
                    -((NSR-1)/2*a - a*jjSR)
                )



                xc = xi + rad*np.cos(theta)

                yc = yi + rad*np.sin(theta)



                circle=list(zip(xc,yc))


                top_cell.add(
                    gdstk.Polygon(
                        circle,
                        layer=LAYER_HOLE[0],
                        datatype=LAYER_HOLE[1]
                    )
                )




# Save

output_path="BICtrampolines_6x6.gds"

lib.write_gds(output_path)
print(top_cell.bounding_box())

print(">>> Saved:", os.path.abspath(output_path))