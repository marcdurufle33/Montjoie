from visuND import *
from pade_fit import *
from scipy.optimize import *
from numpy import *

Ni = 26; Nfreq = 3; Nmat = 5

A = zeros([Nfreq, Nmat]) + 1j*zeros([Nfreq, Nmat])
for j in range(Nfreq):
    for k in range(Nmat):
        coef = zeros(Ni) + 1j*zeros(Ni)
        kx_r = zeros(Ni)
        for i in range(Ni):
            kx = i*0.2; kx_r[i] = kx
            num = k*Nfreq*Ni + j*Ni + i
            if (kx < 0.5):
                n = 5; y = 5.0;
            elif (kx < 1.0):
                n = 4; y = 2.0;
            elif (kx < 2.0):
                n = 3; y = 1.0;
            elif (kx < 4.0):
                n = 2; y = 0.5;
            else:
                n = 1; y = 0.2;
            
            X, Y, Z, coor, V = loadND('../diffracLine' + str(num) + '_U0.dat', n)
            val = V*exp(-1j*kx*X)
            coef[i] = sum(val)/(200*exp(-kx*y));

        f = lambda x : index_frac(x, kx_r, coef);
        sol = least_squares(f, array([-0.3, 0.3]));
        print(sol.x); print(sol.optimality);
        A[j, k] = sol.x[0] + 1j*sol.x[1]

savetxt("coef.dat", A)
