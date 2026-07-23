#!/usr/bin/env python

import sys, os, math

if (len(sys.argv) < 6):
    print "Pas assez d'arguments, mettez au moins 6 arguments"
    line = raw_input()

# on prend en argument le rayon de la sphere
# et les deux facteurs de raffinement
param_r = float(sys.argv[1])
param_a = float(sys.argv[2])
param_L = float(sys.argv[3])
nx = int(sys.argv[4])
ny = int(sys.argv[5])


def GenerateCell(a, b, x0, y0, z0, ref_sphere, ref_pave, first_cell = True):
    fid = open("don", "w")
    fid.write("1\n4\n3\n" + str(a) + "\n" + str(b) + "\n "+ str(x0) + " " + str(y0) + " " + str(z0) + "\n")
    fid.write(str(nx) + " " + str(ny) + "\n 2 \n 2 1\n "+ str(ref_sphere) + "\n " + str(ref_pave[0])+ " " + str(ref_pave[1])+ " " + str(ref_pave[2])+ " " + str(ref_pave[3])+ " " + str(ref_pave[4])+ " " + str(ref_pave[5]) + "\n")
    fid.write("6\n7\n")
    fid.close()
    
    if (first_cell):
        os.system("manipule.x < don > sort && mv toto.mesh cell_sphere.mesh")
    else:
        os.system("manipule.x cell_sphere.mesh < don > sort && mv toto.mesh output.mesh")

def AddLayer(x0, xN, y0, yN, z0, zN, nx, ny, nz, ref_pave):
    fid = open("don", "w")
    fid.write("1\n4\n1\n"+ str(x0) + " " + str(xN) + " " + str(y0) + " " + str(yN) + " " + str(z0) + " " + str(zN) + "\n")
    fid.write(str(nx) + " " + str(ny) + " " + str(nz) + "\n")
    fid.write("1\n"+ str(ref_pave[0])+ " " + str(ref_pave[1])+ " " + str(ref_pave[2])+ " " + str(ref_pave[3])+ " " + str(ref_pave[4])+ " " + str(ref_pave[5])+ "\n2\n6\n7\n")
    fid.close()
    os.system("manipule.x output.mesh < don > sort && mv toto.mesh output.mesh")


# generation de la premiere cellule
ref_pave = [3]*6
ref_pave[5] = 4;
GenerateCell(param_r, param_a, 0, 0, 0, 1, ref_pave)

# generation de la deuxieme cellule
ref_pave = [3]*6
ref_pave[0] = 4;
GenerateCell(param_r, param_a, 2*param_a + param_L, 0, 0, 2, ref_pave, False)

# on remplit le trou
if (param_L>0):
    ref_pave = [3]*6
    ref_pave[0] = 4; ref_pave[5] = 4;
    nb_points_interval = int(math.ceil(param_L / param_a * nx))
    print "Nombre d'intervalles pour boucher le trou =", nb_points_interval
    AddLayer(param_a, param_a+param_L, -param_a, param_a, -param_a, param_a, nb_points_interval+1, 2*nx+1, 2*nx+1, ref_pave)

# on enleve la reference 4
fid = open("don", "w")
fid.write("2\n4\n4\n6\n7\n")
fid.close()
os.system("manipule.x output.mesh < don > sort && mv toto.mesh output.mesh")





