from pylab import *

import os

Nf = 101
Nx = 301

choix_ordre = "7 3"

gam = linspace(0.262, 0.280, 19);
print "gam = ", gam
xmin = 1.0; xmax = 4.2
ymin = -36.45; ymax = -24.8;
dx_min = 0.01; dx_max = 0.06;
dy_min = -0.36; dy_max = 0.084;

def LaunchSdirk(gami, x0, x1, y0, y1, Np):
    print "gami = ", gami
    command = "./sdirk.x " + choix_ordre + " " + str(gami) + " " + str(gami) + " " + str(x0) + " " + str(x1) + " " + str(y0) + " " + str(y1) + " 1 " + str(Np) + " " + str(Np) + " 1"
    print "command = ", command
    os.system(command)
    cfl = reshape(loadtxt('CflR.dat'), (Np, Np))
    return cfl

def FindIntervalCFL(cfl):
    sum_x = sum(cfl, 1)
    sum_y = sum(cfl, 0)
    i0 = find(sum_x>1e100).min()
    i1 = find(sum_x>1e100).max()
    j0 = find(sum_y>1e100).min()
    j1 = find(sum_y>1e100).max()
    return i0, i1, j0, j1

Ng = len(gam) - 1
for i in range(len(gam)):
    gami = gam[i]
    N = int(round(gami*1000))
    print "N = ", N, gami
    cfl = LaunchSdirk(gami, xmin, xmax, ymin, ymax, Nf)    
    if (cfl.max() < 1e100):
        print "on est perdus : ", gam[Ng-i]
        line = raw_input()

    i0, i1, j0, j1 = FindIntervalCFL(cfl)

    x = linspace(xmin, xmax, Nf)
    y = linspace(ymin, ymax, Nf)
    
    print "i, j", i0, i1, j0, j1

    while ((i0 <2) or (i1 >= Nf-2) or (j0 < 2) or (j1 >= Nf-2)):
        Lx = x[i1] - x[i0]
        Ly = y[j1] - y[j0]
        
        if (i0 < 2):
            xmin -= Lx
        else:
            xmin = x[i0-2]

        if (i1 >= Nf-2):
            xmax += Lx
        else:
            xmax = x[i1+2]
            
        if (j0 < 2):
            ymin -= Ly
        else:
            ymin = y[j0-2]

        if (j1 >= Nf-2):
            ymax += Ly
        else:
            ymax = y[j1+2]
        
        cfl = LaunchSdirk(gami, xmin, xmax, ymin, ymax, Nf)
        i0, i1, j0, j1 = FindIntervalCFL(cfl)
        x = linspace(xmin, xmax, Nf)
        y = linspace(ymin, ymax, Nf)
        
        print "i, j", i0, i1, j0, j1, Lx, Ly

    xmin = x[i0-2]; xmax = x[i1+2]
    ymin = y[j0-2]; ymax = y[j1+2]
    Lx = x[i1] - x[i0]
    Ly = y[j1] - y[j0]
    print "xmin, xmax", xmin, xmax, ymin, ymax
    cfl = LaunchSdirk(gami, xmin, xmax, ymin, ymax, Nx)
    i0, i1, j0, j1 = FindIntervalCFL(cfl)
    print "i, j", i0, i1, j0, j1

    command = "./stab.x CflR.dat " + str(Nx) + " " + str(Nx) + " 1 " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " " + str(gami) + " " + str(gami)
    os.system(command)

    os.system("./sdirk.x " + choix_ordre + " test.mesh " + str(gami))
    os.system("mv contour.mesh contourN" + str(N) + ".mesh")

    if (i > 0):
        dx_min = xmin - xold_min
        dx_max = xmax - xold_max
        dy_min = ymin - yold_min
        dy_max = ymax - yold_max

    xold_min = xmin
    xold_max = xmax
    yold_min = ymin
    yold_max = ymax
    xmin += dx_min; xmax += dx_max; ymin += dy_min; ymax += dy_max
    xmin -= 0.03*Lx; xmax += 0.03*Lx; ymin -= 0.03*Lx; ymax += 0.03*Ly
    print "nouvelle boite = ", xmin, xmax, ymin, ymax
    print "dx, dy", dx_min, dx_max, dy_min, dy_max
    
