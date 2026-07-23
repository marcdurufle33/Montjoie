#!/usr/bin/env python

import sys, os, math

# on prend en argument le rayon de la sphere
# et les deux facteurs de raffinement
a = float(sys.argv[1])
nx = int(sys.argv[2])
ny = int(sys.argv[3])

# taille du cube exterieur
b = 1.0

def GenerateCell():
    fid = open("don", "w")
    fid.write("1\n4\n3\n" + str(a) + "\n" + str(b) + "\n 0 0 0\n")
    fid.write(str(nx) + " " + str(ny) + "\n2 \n2 1\n2\n3 3 3 3 3 3\n")
    fid.write("6\n7\n")
    fid.close()
    
    os.system("./manipule.x < don > sort && mv toto.mesh cell_sphere.mesh")

# generation de la cellule
GenerateCell()

# on periodise cette cellule
ix0 = -3; jy0 = -3; kz0 = -1;
ixN = 3; jyN = 3; kzN = 1;

commande_per = "./per_mesh.x cell_sphere.mesh output.mesh 1 x " + str(ix0) + " " + str(ixN) + " y " + str(jy0) + " " + str(jyN) + " z " + str(kz0) + " " + str(kzN) + " r 2"
os.system(commande_per)

# on change la reference de la premiere inclusion
fid = open("don", "w")
fid.write("3\n6\n2\n1\n6\n7\n")
fid.close()
os.system("./manipule.x output.mesh < don > sort && mv toto.mesh output.mesh")


def AddLayer(x0, xN, y0, yN, z0, zN, nx, ny, nz):
    fid = open("don", "w")
    fid.write("1\n4\n1\n"+ str(x0) + " " + str(xN) + " " + str(y0) + " " + str(yN) + " " + str(z0) + " " + str(zN) + "\n")
    fid.write(str(nx) + " " + str(ny) + " " + str(nz) + "\n")
    fid.write("1\n3 3 3 3 3 3\n2\n6\n7\n")
    fid.close()
    os.system("./manipule.x output.mesh < don > sort && mv toto.mesh output.mesh")

# boite exterieure
xmin = -b + 2*ix0*b;
xmax = b + 2*ixN*b;
ymin = -b + 2*jy0*b;
ymax = b + 2*jyN*b;
zmin = -b + 2*kz0*b;
zmax = b + 2*kzN*b;
Nx = 2*nx*(ixN-ix0+1)
Ny = 2*nx*(jyN-jy0+1)
Nz = 2*nx*(kzN-kz0+1)

# on rajoute des zones rectangulaires autour
# attention : dx et dy doivent etre des multiples de b
dx = 3.0;
dy = 3.0;
Npx = nx*int(math.ceil(dx/b))
Npy = nx*int(math.ceil(dy/b))

AddLayer(xmin-dx, xmin, ymin, ymax, zmin, zmax, Npx+1, Ny+1, Nz+1)
AddLayer(xmax, xmax+dx, ymin, ymax, zmin, zmax, Npx+1, Ny+1, Nz+1)
AddLayer(xmin-dx, xmax+dx, ymin-dx, ymin, zmin, zmax, Nx+2*Npx+1, Npy+1, Nz+1)
AddLayer(xmin-dx, xmax+dx, ymax, ymax+dx, zmin, zmax, Nx+2*Npx+1, Npy+1, Nz+1)

# on change les references
fid = open("don", "w")
fid.write("2\n4\n3\n6\n7\n")
fid.close()
os.system("./manipule.x output.mesh < don > sort && mv toto.mesh output.mesh")

fid = open("don", "w")
fid.write("3\n18\n" + str(xmin-dx) + " " + str(xmax+dx) + " " + str(ymin-dx) + " " + str(ymax+dx) + " " + str(zmin) + " " + str(zmin) +"\n2\n")
fid.write("3\n18\n" + str(xmin-dx) + " " + str(xmax+dx) + " " + str(ymin-dx) + " " + str(ymax+dx) + " " + str(zmax) + " " + str(zmax) +"\n3\n")
fid.write("3\n18\n" + str(xmin-dx) + " " + str(xmin-dx) + " " + str(ymin-dx) + " " + str(ymax+dx) + " " + str(zmin) + " " + str(zmax) +"\n3\n")
fid.write("3\n18\n" + str(xmax+dx) + " " + str(xmax+dx) + " " + str(ymin-dx) + " " + str(ymax+dx) + " " + str(zmin) + " " + str(zmax) +"\n3\n")
fid.write("3\n18\n" + str(xmin-dx) + " " + str(xmax+dx) + " " + str(ymin-dx) + " " + str(ymin-dx) + " " + str(zmin) + " " + str(zmax) +"\n3\n")
fid.write("3\n18\n" + str(xmin-dx) + " " + str(xmax+dx) + " " + str(ymax+dx) + " " + str(ymax+dx) + " " + str(zmin) + " " + str(zmax) +"\n3\n")
fid.write("6\n7\n")
fid.close()
os.system("./manipule.x output.mesh < don > sort && mv toto.mesh final.mesh")
