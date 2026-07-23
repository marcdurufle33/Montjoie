from pylab import *
import numpy as np
import mpmath as mp
import matrices as mat

mp.mp.dps = 200
pasNewton = 1e-5

def det2d(a,b,c,d):
    return a*c-b*d

def Jac(equation, x, y, n, f, eps, sigma, mu, rayon, m) :
    mp.mp.dps = 200
    J = mp.zeros(2)
    # schéma centré (ordre 2) pour les dérivées
    DxM1y = mat.ComputeDet(equation, mp.mpc(x - pasNewton, y), n, f, eps, sigma, mu, rayon, m)
    DxP1y = mat.ComputeDet(equation, mp.mpc(x + pasNewton, y), n, f, eps, sigma, mu, rayon, m)
    DxyM1 = mat.ComputeDet(equation, mp.mpc(x, y - pasNewton), n, f, eps, sigma, mu, rayon, m)
    DxyP1 = mat.ComputeDet(equation, mp.mpc(x, y + pasNewton), n, f, eps, sigma, mu, rayon, m)
    J[0, 0] = 0.5 * (mp.re(DxP1y) - mp.re(DxM1y)) / pasNewton
    J[1, 0] = 0.5 * (mp.im(DxP1y) - mp.im(DxM1y)) / pasNewton
    J[0, 1] = 0.5 * (mp.re(DxyP1) - mp.re(DxyM1)) / pasNewton
    J[1, 1] = 0.5 * (mp.im(DxyP1) - mp.im(DxyM1)) / pasNewton
    return J


def Newton(equation, B0, expTol, n, f, eps, sigma, mu, rayon, m) :
    mp.mp.dps = 200
    print("B0 : ", B0)
    x = mp.re(B0); y = mp.im(B0)
    X0 = mp.matrix([x, y])
    DComp = mat.ComputeDet(equation, mp.mpc(x, y), n, f, eps, sigma, mu, rayon, m)
    D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
    Jac0 = Jac(equation, x, y, n, f, eps, sigma, mu, rayon, m)
    i = 0
    vectDiff = mp.lu_solve(Jac0, D0)
    norme2 = (vectDiff.T * vectDiff)[0]
    X0 = X0 - vectDiff
    # mp.mp.dps = 15 # pour l'affichage
    print()
    print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
    while norme2 >= 10**(-expTol * 2) :
        i += 1
        DComp = mat.ComputeDet(equation, mp.mpc(X0[0], X0[1]), n, f, eps, sigma, mu, rayon, m)
        D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
        Jac0 = Jac(equation, X0[0], X0[1], n, f, eps, sigma, mu, rayon, m)
        vectDiff = mp.lu_solve(Jac0, D0)
        norme2 = (vectDiff.T * vectDiff)[0]
        X0 = X0 - vectDiff
        # mp.mp.dps = 15 # pour l'affichage
        print()
        print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
        mp.mp.dps = 100
    return X0[0] +1j * X0[1]
# print(X)
