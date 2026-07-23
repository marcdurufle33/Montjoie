from pylab import *
import numpy as np
import mpmath as mp
from enum import *

mp.mp.dps = 200

def createMatAMaxwellDirichlet(n, k, rayon, c) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = c[i]*k[i]*mp.besselj(1, k[i]*rayon[i]);
        A[2*i+1, 2*i] = c[i]*k[i]*mp.bessely(1, k[i]*rayon[i]);
        A[2*i+1, 2*i+1] = -c[i+1]*k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
        A[2*i+1, 2*i+2] = -c[i+1]*k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);
    # Condition de Dirichlet homogène assez loin
    if (n > 0) :
        A[2*n, 2*n-1] = mp.besselj(0, k[n]*rayon[n]);
        A[2*n, 2*n] = mp.bessely(0, k[n]*rayon[n])
    # Faire cas Hankel
    # Faire cas BGT2
    return A

def createMatAMaxwellOpen(n, k, rayon, c, m = 0) :
    A = mp.zeros(2*n)
    if (n > 1) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);
        for i in range(1, n-1) :
            # Continuité de Ez
            A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
            A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
            A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
            A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
            # Continuité de la dérivée
            A[2*i+1, 2*i-1] = c[i]*k[i]*mp.besselj(1, k[i]*rayon[i]);
            A[2*i+1, 2*i] = c[i]*k[i]*mp.bessely(1, k[i]*rayon[i]);
            A[2*i+1, 2*i+1] = -c[i+1]*k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
            A[2*i+1, 2*i+2] = -c[i+1]*k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);

        A[2*n-2, 2*n-3] = mp.besselj(0, k[n-1]*rayon[n-1]);
        A[2*n-2, 2*n-2] = mp.bessely(0, k[n-1]*rayon[n-1]);
        A[2*n-2, 2*n-1] = -mp.hankel1(0, k[n]*rayon[n-1]);
        # Continuité de la dérivée
        A[2*n-1, 2*n-3] = -c[n-1]*k[n-1]*mp.besselj(1, k[n-1]*rayon[n-1]);
        A[2*n-1, 2*n-2] = -c[n-1]*k[n-1]*mp.bessely(1, k[n-1]*rayon[n-1]);
        A[2*n-1, 2*n-1] = c[n]*k[n]*mp.hankel1(1, k[n]*rayon[n-1]);

    elif (n == 1) :
        A[0, 0] = mp.besselj(0, k[n-1]*rayon[n-1]);
        A[0, 1] = -mp.hankel1(0, k[n]*rayon[n-1]);
        # Continuité de la dérivée
        A[1, 0] = -c[n-1]*k[n-1]*mp.besselj(1, k[n-1]*rayon[n-1]);
        A[1, 1] = c[n]*k[n]*mp.hankel1(1, k[n]*rayon[n-1]);
    # Faire cas Hankel
    # Faire cas BGT2
    return A

def createMatAMaxwellBGT2(n, k, rayon, c) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = c[i]*k[i]*mp.besselj(1, k[i]*rayon[i]);
        A[2*i+1, 2*i] = c[i]*k[i]*mp.bessely(1, k[i]*rayon[i]);
        A[2*i+1, 2*i+1] = -c[i+1]*k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
        A[2*i+1, 2*i+2] = -c[i+1]*k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);
    # Condition de Dirichlet homogène assez loin
    if (n > 0) :
        sBeta = - 1j*k[n]
        R = rayon[n]
        num = 3/4 + 3*R*sBeta + 2*R*R*sBeta*sBeta
        denom = 2 * R * (1 + sBeta * R)

        X = num / denom

        A[2*n, 2*n-1] = -k[n]*mp.besselj(1, k[n]*R) + X * mp.besselj(0, k[n]*R)
        A[2*n, 2*n] = -k[n]*mp.bessely(1, k[n]*R) + X * mp.bessely(0, k[n]*R)

    # Faire cas Hankel
    # Faire cas BGT2
    return A



def createMatAMaxwellBGT2N3(n, k, rayon, c, appRacine) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = c[i]*k[i]*mp.besselj(1, k[i]*rayon[i]);
        A[2*i+1, 2*i] = c[i]*k[i]*mp.bessely(1, k[i]*rayon[i]);
        A[2*i+1, 2*i+1] = -c[i+1]*k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
        A[2*i+1, 2*i+2] = -c[i+1]*k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);
    # Condition de Dirichlet homogène assez loin
    if (n > 0) :
        sBeta = - 1j*k[n]
        R = rayon[n]

        num = - R * sBeta * sBeta + 9/(8*R)
        denom = 1 + appRacine * R

        A[2*n, 2*n-1] = -k[n]*mp.besselj(1, k[n]*R) + (3/(2*R) - num/denom) * mp.besselj(0, k[n]*R)
        A[2*n, 2*n] = -k[n]*mp.bessely(1, k[n]*R) + (3/(2*R) - num/denom) * mp.bessely(0, k[n]*R)

    # Faire cas Hankel
    # Faire cas BGT2
    return A

# matrice Helmholtz
def createMatAHelmholtzDirichlet(n, k, rayon, c) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -k[1]*mp.bessely(1, k[1]*rayon[0]);
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = k[i]*mp.besselj(1, k[i]*rayon[i]);
        A[2*i+1, 2*i] = k[i]*mp.bessely(1, k[i]*rayon[i]);
        A[2*i+1, 2*i+1] = -k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
        A[2*i+1, 2*i+2] = -k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);
    # Condition de Dirichlet homogène assez loin
    if (n > 0) :
        A[2*n, 2*n-1] = mp.besselj(0, k[n]*rayon[n]);
        A[2*n, 2*n] = mp.bessely(0, k[n]*rayon[n])
    # Faire cas Hankel
    # Faire cas BGT2
    return A


def createMatAHelmholtzOpen(n, k, rayon, c, m = 0) :
    A = mp.zeros(2*n)
    if (n > 1) :
        A[0, 0] = mp.besselj(m, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(m, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(m, k[1]*rayon[0]);
        A[1, 0] = -k[0]*mp.besselj(m+1, k[0]*rayon[0]) + (m/rayon[0]) * mp.besselj(m, k[0]*rayon[0]);
        A[1, 1] = k[1]*mp.besselj(m+1, k[1]*rayon[0]) - (m/rayon[0]) * mp.besselj(m, k[1]*rayon[0]);
        A[1, 2] = k[1]*mp.bessely(m+1, k[1]*rayon[0]) - (m/rayon[0]) * mp.bessely(m, k[1]*rayon[0]);
        for i in range(1, n-1) :
            # Continuité de Ez
            A[2*i, 2*i-1] = mp.besselj(m, k[i]*rayon[i]);
            A[2*i, 2*i] = mp.bessely(m, k[i]*rayon[i]);
            A[2*i, 2*i+1] = -mp.besselj(m, k[i+1]*rayon[i]);
            A[2*i, 2*i+2] = -mp.bessely(m, k[i+1]*rayon[i]);
            # Continuité de la dérivée
            A[2*i+1, 2*i-1] = -k[i]*mp.besselj(m+1, k[i]*rayon[i]) + (m/rayon[i]) * mp.besselj(m, k[i]*rayon[i])
            A[2*i+1, 2*i] = -k[i]*mp.bessely(m+1, k[i]*rayon[i]) + (m/rayon[i]) * mp.bessely(m, k[i]*rayon[i])
            A[2*i+1, 2*i+1] = k[i+1]*mp.besselj(m+1, k[i+1]*rayon[i]) - (m/rayon[i]) * mp.besselj(m, k[i+1]*rayon[i])
            A[2*i+1, 2*i+2] = k[i+1]*mp.bessely(m+1, k[i+1]*rayon[i]) - (m/rayon[i]) * mp.bessely(m, k[i+1]*rayon[i])

        A[2*n-2, 2*n-3] = mp.besselj(m, k[n-1]*rayon[n-1]);
        A[2*n-2, 2*n-2] = mp.bessely(m, k[n-1]*rayon[n-1]);
        A[2*n-2, 2*n-1] = -mp.hankel1(m, k[n]*rayon[n-1]);
        # Continuité de la dérivée
        A[2*n-1, 2*n-3] = -k[n-1]*mp.besselj(m+1, k[n-1]*rayon[n-1]) + (m/rayon[n-1]) * mp.besselj(m, k[n-1]*rayon[n-1])
        A[2*n-1, 2*n-2] = -k[n-1]*mp.bessely(m+1, k[n-1]*rayon[n-1]) + (m/rayon[n-1]) * mp.bessely(m, k[n-1]*rayon[n-1])
        A[2*n-1, 2*n-1] = k[n]*mp.hankel1(m+1, k[n]*rayon[n-1]) - (m/rayon[n-1]) * mp.hankel1(m, k[n]*rayon[n-1])
    elif (n == 1) :
        A[0, 0] = mp.besselj(m, k[0]*rayon[0]);
        A[0, 1] = -mp.hankel1(m, k[1]*rayon[0]);
        # Continuité de la dérivée
        A[1, 0] = -k[0]*mp.besselj(m+1, k[0]*rayon[0]) + (m/rayon[0]) * mp.besselj(m, k[0]*rayon[0]);
        A[1, 1] = k[1]*mp.hankel1(m+1, k[1]*rayon[0]) - (m/rayon[0]) * mp.hankel1(m, k[1]*rayon[0]);
    # Faire cas Hankel
    # Faire cas BGT2
    # print(A)
    return A


# matrice Helmholtz BGT2 (Newton complexe converge mal)
def createMatAHelmholtzBGT2(n, k, rayon, c, m = 0) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(m, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(m, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(m, k[1]*rayon[0]);
        A[1, 0] = k[0]*mp.besselj(m+1, k[0]*rayon[0]) - (m/rayon[0]) * mp.besselj(m, k[0]*rayon[0])
        A[1, 1] = -k[1]*mp.besselj(m+1, k[1]*rayon[0]) + (m/rayon[0]) * mp.besselj(m, k[1]*rayon[0])
        A[1, 2] = -k[1]*mp.bessely(m+1, k[1]*rayon[0]) + (m/rayon[0]) * mp.bessely(m, k[1]*rayon[0])
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(m, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(m, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(m, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(m, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = k[i]*mp.besselj(m+1, k[i]*rayon[i]) - (m/rayon[i]) * mp.besselj(m, k[i]*rayon[i])
        A[2*i+1, 2*i] = k[i]*mp.bessely(m+1, k[i]*rayon[i]) - (m/rayon[i]) * mp.bessely(m, k[i]*rayon[i])
        A[2*i+1, 2*i+1] = -k[i+1]*mp.besselj(m+1, k[i+1]*rayon[i]) + (m/rayon[i]) * mp.besselj(m, k[i+1]*rayon[i])
        A[2*i+1, 2*i+2] = -k[i+1]*mp.bessely(m+1, k[i+1]*rayon[i]) + (m/rayon[i]) * mp.bessely(m, k[i+1]*rayon[i])
    # Condition BGT2
    if (n > 0) :
        sBeta = - 1j*k[n]
        R = rayon[n]
        num = 3/4 + 3*R*sBeta + 2*R*R*sBeta*sBeta
        denom = 2 * R * (1 + sBeta * R)

        X = num / denom
        Y = 1.0 / denom

        A[2*n, 2*n-1] = -k[n]*mp.besselj(m+1, k[n]*R) + X * mp.besselj(m, k[n]*R) + (m/R) * mp.besselj(m, k[n]*R) + m*m*Y*mp.besselj(m, k[n]*R)
        A[2*n, 2*n]   = -k[n]*mp.bessely(m+1, k[n]*R) + X * mp.bessely(m, k[n]*R) + (m/R) * mp.bessely(m, k[n]*R) + m*m*Y*mp.bessely(m, k[n]*R)
    return A



# matrice Helmholtz BGT2N3
def createMatAHelmholtzBGT2N3(n, k, rayon, c, appRacine) :
    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(m, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(m, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(m, k[1]*rayon[0]);
        A[1, 0] = k[0]*mp.besselj(m+1, k[0]*rayon[0]) - (m/rayon[0]) * mp.besselj(m, k[0]*rayon[0])
        A[1, 1] = -k[1]*mp.besselj(m+1, k[1]*rayon[0]) + (m/rayon[0]) * mp.besselj(m, k[1]*rayon[0])
        A[1, 2] = -k[1]*mp.bessely(m+1, k[1]*rayon[0]) + (m/rayon[0]) * mp.bessely(m, k[1]*rayon[0])
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(m, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(m, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(m, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(m, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = k[i]*mp.besselj(m+1, k[i]*rayon[i]) - (m/rayon[i]) * mp.besselj(m, k[i]*rayon[i])
        A[2*i+1, 2*i] = k[i]*mp.bessely(m+1, k[i]*rayon[i]) - (m/rayon[i]) * mp.bessely(m, k[i]*rayon[i])
        A[2*i+1, 2*i+1] = -k[i+1]*mp.besselj(m+1, k[i+1]*rayon[i]) + (m/rayon[i]) * mp.besselj(m, k[i+1]*rayon[i])
        A[2*i+1, 2*i+2] = -k[i+1]*mp.bessely(m+1, k[i+1]*rayon[i]) + (m/rayon[i]) * mp.bessely(m, k[i+1]*rayon[i])


    A = mp.zeros(2*n+1)
    if (n > 0) :
        A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
        A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
        A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
        A[1, 0] = k[0]*mp.besselj(1, k[0]*rayon[0]);
        A[1, 1] = -k[1]*mp.besselj(1, k[1]*rayon[0]);
        A[1, 2] = -k[1]*mp.bessely(1, k[1]*rayon[0]);
    for i in range(1, n) :
        # Continuité de Ez
        A[2*i, 2*i-1] = mp.besselj(0, k[i]*rayon[i]);
        A[2*i, 2*i] = mp.bessely(0, k[i]*rayon[i]);
        A[2*i, 2*i+1] = -mp.besselj(0, k[i+1]*rayon[i]);
        A[2*i, 2*i+2] = -mp.bessely(0, k[i+1]*rayon[i]);
        # Continuité de la dérivée
        A[2*i+1, 2*i-1] = k[i]*mp.besselj(1, k[i]*rayon[i]);
        A[2*i+1, 2*i] = k[i]*mp.bessely(1, k[i]*rayon[i]);
        A[2*i+1, 2*i+1] = -k[i+1]*mp.besselj(1, k[i+1]*rayon[i]);
        A[2*i+1, 2*i+2] = -k[i+1]*mp.bessely(1, k[i+1]*rayon[i]);
    # Condition BGT2
    if (n > 0) :
        sBeta = - 1j*k[n]
        R = rayon[n]

        num = - R * sBeta * sBeta + 9/(8*R)
        denom = 1 + appRacine * R

        A[2*n, 2*n-1] = -k[n]*mp.besselj(1, k[n]*R) + (3/(2*R) - num/denom) * mp.besselj(0, k[n]*R)
        A[2*n, 2*n] = -k[n]*mp.bessely(1, k[n]*R) + (3/(2*R) - num/denom) * mp.bessely(0, k[n]*R)
    return A


""" Choix de la matrice """

class Equation(Enum):
    MAXWELL_DIRICHLET = auto()
    MAXWELL_OPEN = auto()
    MAXWELL_BGT2 = auto()
    MAXWELL_BGT2N3 = auto()
    HELMHOLTZ_DIRICHLET = auto()
    HELMHOLTZ_OPEN = auto()
    HELMHOLTZ_BGT2 = auto()
    HELMHOLTZ_BGT2N3 = auto()

def defineEquation(str) :
    if (str == "MAXWELL_DIRICHLET") :
        return Equation.MAXWELL_DIRICHLET
    elif (str == "MAXWELL_OPEN") :
        return Equation.MAXWELL_OPEN
    elif (str == "MAXWELL_BGT2") :
        return Equation.MAXWELL_BGT2
    elif (str == "MAXWELL_BGT2N3") :
        return Equation.MAXWELL_BGT2N3
    elif (str == "HELMHOLTZ_DIRICHLET") :
        return Equation.HELMHOLTZ_DIRICHLET
    elif (str == "HELMHOLTZ_OPEN") :
        return Equation.HELMHOLTZ_OPEN
    elif (str == "HELMHOLTZ_BGT2") :
        return Equation.HELMHOLTZ_BGT2
    elif (str == "HELMHOLTZ_BGT2N3") :
        return Equation.HELMHOLTZ_BGT2N3
    else :
        print("PB EQUATION :", str)
        exit()

# liste toutes les matrices faites pour assigner la bonne
def createMatA(equation, n, k, rayon, c, appRacine = None, m = 0) :
    if equation == Equation.MAXWELL_DIRICHLET :
        return createMatAMaxwellDirichlet(n, k, rayon, c)
    elif equation == Equation.MAXWELL_OPEN :
        return createMatAMaxwellOpen(n, k, rayon, c)
    elif equation == Equation.MAXWELL_BGT2 :
        return createMatAMaxwellBGT2(n, k, rayon, c)
    elif equation == Equation.MAXWELL_BGT2N3 :
        return createMatAMaxwellBGT2N3(n, k, rayon, c, appRacine)
    elif equation == Equation.HELMHOLTZ_DIRICHLET :
        return createMatAHelmholtzDirichlet(n, k, rayon, c)
    elif equation == Equation.HELMHOLTZ_OPEN :
        return createMatAHelmholtzOpen(n, k, rayon, c, m)
    elif equation == Equation.HELMHOLTZ_BGT2 :
        return createMatAHelmholtzBGT2(n, k, rayon, c, m)
    elif equation == Equation.HELMHOLTZ_BGT2N3 :
        return createMatAHelmholtzBGT2N3(n, k, rayon, c, appRacine)
    else :
        print("PB MATRICES !")
        exit()


""" Déterminants """

# toutes les valeurs doivent être en mp précision
def ComputeDet(equation, beta, n, f, eps, sigma, mu, rayon, m):
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    rho = mp.zeros(n+1, 1)
    k = mp.zeros(n+1, 1)
    c = mp.zeros(n+1, 1)

    rhoInfini = eps[n]*mu[n] + 1j*sigma[n]*mu[n]/omega
    rhoPlus = 0.0
    appRacine = 0.0
    # on va essayer de lui mettre la partie réelle maximale

    for i in range(n+1):
        rhoPlus = max(eps[i]*mu[i], rhoPlus)
        Delta= omega*omega*eps[i]*mu[i] + 1j*omega*sigma[i]*mu[i] - beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = -mp.sqrt(Delta)
        # print("D :", Delta, "sqrt :", k[i])
        # print(k[i])
        # print("k : ", k[i])
        c[i] = (omega*omega*eps[i] + 1j*omega*sigma[i]) / Delta
        # print("c : ", c[i])

    if equation == Equation.HELMHOLTZ_BGT2N3 or equation == Equation.MAXWELL_BGT2N3 :
        ww = omega
        alpha = mp.sqrt(rhoPlus - rhoInfini) * omega
        appRacine = alpha/8-(rhoInfini*ww*ww-beta*beta*ww*ww)/(8*alpha) - \
                    (rhoInfini*ww*ww-beta*beta*ww*ww)/\
                    (2.0*alpha - 2.0*(rhoInfini*ww*ww - beta*beta*ww*ww)/alpha) \
                     - (rhoInfini*ww*ww-beta*beta*ww*ww)/\
                     (alpha/2.0 - (rhoInfini*ww*ww - beta*beta*ww*ww)/\
                     (2.0*alpha) - 2.0*(rhoInfini*ww*ww - beta*beta*ww*ww)/\
                     (alpha - (rhoInfini*ww*ww - beta*beta*ww*ww)/alpha))

    A = createMatA(equation, n, k, rayon, c, appRacine, m)
    #print("k = ", k)
    return mp.det(A)


if __name__ == "__main__":
    for eq in Equation:
        print(eq)

    eq = defineEquation("MAXWELL_DIRICHLET")
    print(eq)
    eq2 = defineEquation("HELMHOLTZ_OPEN")
    print(eq2)
