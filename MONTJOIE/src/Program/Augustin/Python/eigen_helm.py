from pylab import *
from scipy.special import *
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

mp.mp.dps = 200

pasNewton = 1e-5
rayonCLA = 0.4

def createMatASommerfeld(k1, k2, a) :
    A = mp.zeros(2)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.hankel1(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.hankel1(1, k2*a);
    return A

def createMatATrueABC(k1, k2, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    A[2, 1] = -mp.besselj(1, k2*R) - 1j*mp.besselj(0, k2*R);
    A[2, 2] = -mp.bessely(1, k2*R) - 1j*mp.bessely(0, k2*R);
    return A

def createMatA_BGT1(k1, k2, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    A[2, 1] = -k2*mp.besselj(1, k2*R) + (1/(2*R) - 1j*k2) * mp.besselj(0, k2*R);
    A[2, 2] = -k2*mp.bessely(1, k2*R) + (1/(2*R) - 1j*k2) * mp.bessely(0, k2*R);
    return A

def createMatA_BGT2(k1, k2, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);

    num = 3/4 - 3*1j*R*k2 - 2*R*R*k2*k2
    denom = 2 * R * (1 - 1j * k2 * R)


    A[2, 0] = 0.0;
    A[2, 1] = -k2*mp.besselj(1, k2*R) + num/denom * mp.besselj(0, k2*R);
    A[2, 2] = -k2*mp.bessely(1, k2*R) + num/denom * mp.bessely(0, k2*R);
    return A

def createMatA_BGT2_Newton(k1, k2, appRacine, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);

    num = k2*k2*R + 9/(8*R)
    # num = -appRacine*appRacine*R + 9/(8*R)
    denom = 1 + appRacine * R
    # denom = 1 - 1j * k2 * R

    #
    # num = 3/4 - 2*1j*R*k2 - 2*R*R*k2*k2
    # denom = R * (1 - 2*1j*k2*R)

    A[2, 0] = 0.0;
    A[2, 1] = -k2*mp.besselj(1, k2*R) + (3/(2*R) - num/denom) * mp.besselj(0, k2*R);
    A[2, 2] = -k2*mp.bessely(1, k2*R) + (3/(2*R) - num/denom) * mp.bessely(0, k2*R);
    return A

def createMatAK0(k1, k2, kBeta, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    eps = 0.001
    q = mp.besselk(1, kBeta*R)/mp.besselk(0, kBeta*R)
    A[2, 1] = - k2 * mp.besselj(1, k2*R) + kBeta * q * mp.besselj(0, k2*R);
    A[2, 2] = - k2 * mp.bessely(1, k2*R) + kBeta * q * mp.bessely(0, k2*R);
    return A

def createMatA1ABC(k1, k2, kRob, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    A[2, 1] = k2 * mp.besselj(1, k2*R) + 1j*kRob*mp.besselj(0, k2*R);
    A[2, 2] = k2 * mp.bessely(1, k2*R) + 1j*kRob*mp.bessely(0, k2*R);
    return A

def createMatA_NewtonABC(k1, k2, appRacine, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = -k1*mp.besselj(1, k1*a);
    A[1, 1] = k2*mp.besselj(1, k2*a);
    A[1, 2] = k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    A[2, 1] = - k2 * mp.besselj(1, k2*R) + appRacine*mp.besselj(0, k2*R);
    A[2, 2] = - k2 * mp.bessely(1, k2*R) + appRacine*mp.bessely(0, k2*R);
    return A

def createMatADirichlet(k1, k2, a, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(0, k1*a);
    A[0, 1] = -mp.besselj(0, k2*a);
    A[0, 2] = -mp.bessely(0, k2*a);
    A[1, 0] = k1*mp.besselj(1, k1*a);
    A[1, 1] = -k2*mp.besselj(1, k2*a);
    A[1, 2] = -k2*mp.bessely(1, k2*a);
    A[2, 0] = 0.0;
    A[2, 1] = mp.besselj(0, k2*R);
    A[2, 2] = mp.bessely(0, k2*R);
    return A

def det2D(A) :
    return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]


def det3D(A) :
    return A[0, 0]*A[1, 1]*A[2, 2] + A[0, 1]*A[1, 2]*A[2, 0] + A[0, 2]*A[1, 0]*A[2, 1] - A[0, 2]*A[1, 1]*A[2,0] - A[0, 1]*A[1, 0]*A[2,2] - A[0, 0]*A[1, 2]*A[2,1]

# toutes les valeurs doivent être en mp précision
def ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta*beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta*beta*ww*ww)
    A = createMatASommerfeld(k1, k2, a)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det2D(A)

def ComputeDetTrueABC(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatATrueABC(k1, k2, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT1(k1, k2, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT2(k1, k2, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDet_BGT2_Newton(beta, ww, rhoMoins, rhoPlus, a, R, step):
    mp.mp.dps = 200
    alpha = np.sqrt(rhoPlus - rhoMoins) * ww
    if step == 0 :
        appRacine = alpha
    elif step == 1 :
        appRacine = alpha/2-(rhoMoins*ww*ww-beta*beta*ww*ww)/(2*alpha)
    elif step == 2 :
        appRacine = alpha/4-(rhoMoins*ww*ww-beta*beta*ww*ww)/(4*alpha) - (rhoMoins*ww*ww-beta*beta*ww*ww)/(alpha + (-rhoMoins*ww*ww + beta*beta*ww*ww)/alpha)
    elif step == 3 :
        appRacine = alpha/8-(rhoMoins*ww*ww-beta*beta*ww*ww)/(8*alpha) - (rhoMoins*ww*ww-beta*beta*ww*ww)/(2.0*alpha - 2.0*(rhoMoins*ww*ww - beta*beta*ww*ww)/alpha) \
                     - (rhoMoins*ww*ww-beta*beta*ww*ww)/(alpha/2.0 - (rhoMoins*ww*ww - beta*beta*ww*ww)/(2.0*alpha) - 2.0*(rhoMoins*ww*ww - beta*beta*ww*ww)/(alpha - (rhoMoins*ww*ww - beta*beta*ww*ww)/alpha))
    else :
        sys.exit("Racine non calculée")
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT2_Newton(k1, k2, appRacine, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDetK0(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    kBeta = mp.sqrt(beta * beta*ww*ww - rhoMoins*ww*ww)
    A = createMatAK0(k1, k2, kBeta, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)


def ComputeDet1ABC(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    kRob = mp.sqrt(rhoMoins*ww*ww - rhoPlus*ww*ww)
    A = createMatA1ABC(k1, k2, kRob, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDet_NewtonABC(beta, ww, rhoMoins, rhoPlus, a, R, step):
    mp.mp.dps = 200
    eps0 = - np.sqrt(rhoPlus - rhoMoins) * ww
    if step == 1 :
        appRacine = -eps0/2+(rhoMoins*ww*ww-beta*beta*ww*ww)/(2*eps0)
    elif step == 2 :
        appRacine = -eps0/4+(rhoMoins*ww*ww-beta*beta*ww*ww)/(4*eps0) + (rhoMoins*ww*ww-beta*beta*ww*ww)/(eps0 + (-rhoMoins*ww*ww + beta*beta*ww*ww)/eps0)
    else :
        sys.exit("Racine non calculée")
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_NewtonABC(k1, k2, appRacine, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)

def ComputeDetDirichlet(beta, ww, rhoMoins, rhoPlus, a, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatADirichlet(k1, k2, a, R)
    # return - mp.besselj(0, k1*a) * k2*mp.hankel1(1, k2*a) + k1*mp.besselj(1, k1*a) * mp.hankel1(0, k2*a);
    return det3D(A)


def createMatA(k1, k2, a, ray, cond, kBord = 0) :
    if cond == 0 :
        return createMatASommerfeld(k1, k2, a)
    elif cond == 1 :
        return createMatATrueABC(k1, k2, a, ray)
    elif cond == 2 :
        return createMatA1ABC(k1, k2, kBord, a, ray)
    elif (cond == 3 or cond == 4) :
        return createMatA_NewtonABC(k1, k2, kBord, a, ray)
    elif cond == 5 :
        return createMatADirichlet(k1, k2, a, ray)
    elif cond == 6 :
        return createMatAK0(k1, k2, kBord, a, ray)
    elif cond == 7 :
        return createMatA_BGT1(k1, k2, a, ray)
    elif cond == 8 :
        return createMatA_BGT2(k1, k2, a, ray)
    elif (cond == 9 or cond == 10 or cond == 11 or cond == 12) :
        return createMatA_BGT2_Newton(k1, k2, kBord, a, ray)

def ComputeDet(beta, ww, rhoMoins, rhoPlus, a, ray, cond) :
    if cond == 0 :
        return ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a)
    elif cond == 1 :
        return ComputeDetTrueABC(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 2 :
        return ComputeDet1ABC(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 3 :
        return ComputeDet_NewtonABC(beta, ww, rhoMoins, rhoPlus, a, ray, 1)
    elif cond == 4 :
        return ComputeDet_NewtonABC(beta, ww, rhoMoins, rhoPlus, a, ray, 2)
    elif cond == 5 :
        return ComputeDetDirichlet(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 6 :
        return ComputeDetK0(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 7 :
        return ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 8 :
        return ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, ray)
    elif cond == 9 :
        return ComputeDet_BGT2_Newton(beta, ww, rhoMoins, rhoPlus, a, ray, 0)
    elif cond == 10 :
        return ComputeDet_BGT2_Newton(beta, ww, rhoMoins, rhoPlus, a, ray, 1)
    elif cond == 11 :
        return ComputeDet_BGT2_Newton(beta, ww, rhoMoins, rhoPlus, a, ray, 2)
    elif cond == 12 :
        return ComputeDet_BGT2_Newton(beta, ww, rhoMoins, rhoPlus, a, ray, 3)

def ComputeFLambda(ww, rhoMoins, rhoPlus, a, ray, cond) :
    if cond == 0 :
        return lambda x : np.float64(imag(ComputeDetSommerfeld(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a)))
    elif cond == 1 :
        return lambda x : np.float64(imag(ComputeDetTrueABC(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 2 :
        return lambda x : np.float64(real(ComputeDet1ABC(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 3 :
        return lambda x : np.float64(real(ComputeDet_NewtonABC(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 1)))
    elif cond == 4 :
        return lambda x : np.float64(real(ComputeDet_NewtonABC(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 2)))
    elif cond == 5 :
        return lambda x : np.float64(real(ComputeDetDirichlet(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 6 :
        return lambda x : np.float64(real(ComputeDetK0(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 7 :
        return lambda x : np.float64(real(ComputeDet_BGT1(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 8 :
        return lambda x : np.float64(real(ComputeDet_BGT2(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray)))
    elif cond == 9 :
        return lambda x : np.float64(real(ComputeDet_BGT2_Newton(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 0)))
    elif cond == 10 :
        return lambda x : np.float64(real(ComputeDet_BGT2_Newton(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 1)))
    elif cond == 11 :
        return lambda x : np.float64(real(ComputeDet_BGT2_Newton(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 2)))
    elif cond == 12 :
        return lambda x : np.float64(real(ComputeDet_BGT2_Newton(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, ray, 3)))

def Jac(x, y, ww, rhoMoins, rhoPlus, a, ray, cond) :
    mp.mp.dps = 200
    J = mp.zeros(2)
    # schéma centré (ordre 2) pour les dérivées
    DxM1y = ComputeDet(mp.mpc(x - pasNewton, y), ww, rhoMoins, rhoPlus, a, ray, cond)
    DxP1y = ComputeDet(mp.mpc(x + pasNewton, y), ww, rhoMoins, rhoPlus, a, ray, cond)
    DxyM1 = ComputeDet(mp.mpc(x, y - pasNewton), ww, rhoMoins, rhoPlus, a, ray, cond)
    DxyP1 = ComputeDet(mp.mpc(x, y + pasNewton), ww, rhoMoins, rhoPlus, a, ray, cond)
    J[0, 0] = 0.5 * (mp.re(DxP1y) - mp.re(DxM1y)) / pasNewton
    J[1, 0] = 0.5 * (mp.im(DxP1y) - mp.im(DxM1y)) / pasNewton
    J[0, 1] = 0.5 * (mp.re(DxyP1) - mp.re(DxyM1)) / pasNewton
    J[1, 1] = 0.5 * (mp.im(DxyP1) - mp.im(DxyM1)) / pasNewton
    return J


def Newton(B0, expTol, ww, rhoMoins, rhoPlus, a, ray, cond) :
    mp.mp.dps = 200
    # print('Newton')
    # print('Type de rhoPlus :', type(rhoPlus))
    # print('Type de rhoMoins :', type(rhoMoins))
    # print('Type de a :', type(a))
    # print('Type de ww :', type(ww))
    print("B0 : ", B0)
    x = mp.re(B0); y = mp.im(B0)
    X0 = mp.matrix([x, y])
    DComp = ComputeDet(mp.mpc(x, y), ww, rhoMoins, rhoPlus, a, ray, cond)
    D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
    Jac0 = Jac(x, y, ww, rhoMoins, rhoPlus, a, ray, cond)
    i = 0
    vectDiff = mp.lu_solve(Jac0, D0)
    norme2 = (vectDiff.T * vectDiff)[0]
    X0 = X0 - vectDiff
    # mp.mp.dps = 15 # pour l'affichage
    print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
    while norme2 >= 10**(-expTol * 2) :
        i += 1
        DComp = ComputeDet(mp.mpc(X0[0], X0[1]), ww, rhoMoins, rhoPlus, a, ray, cond)
        D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
        Jac0 = Jac(X0[0], X0[1], ww, rhoMoins, rhoPlus, a, ray, cond)
        vectDiff = mp.lu_solve(Jac0, D0)
        norme2 = (vectDiff.T * vectDiff)[0]
        X0 = X0 - vectDiff
        # mp.mp.dps = 15 # pour l'affichage
        print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
        mp.mp.dps = 100
    return X0[0] +1j * X0[1]
    # print(X)

def ComputeParam(beta, ww, rhoMoins, rhoPlus, a, ray, cond) :
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    kBord = 0
    if cond == 2 :
        kBord = mp.sqrt(rhoMoins*ww*ww - rhoPlus*ww*ww)
    elif cond == 3 :
        eps0 = - np.sqrt(rhoPlus - rhoMoins) * ww
        kBord = -eps0/2+(rhoMoins*ww*ww-beta*beta*ww*ww)/(2*eps0)
    elif cond == 4 :
        eps0 = - np.sqrt(rhoPlus - rhoMoins) * ww
        kBord = -eps0/4+(rhoMoins*ww*ww-beta*beta*ww*ww)/(4*eps0) + (rhoMoins*ww*ww-beta*beta*ww*ww)/(eps0 + (-rhoMoins*ww*ww + beta*beta*ww*ww)/eps0)
    elif cond == 6 :
        kBord = mp.sqrt(beta*beta*ww*ww - rhoMoins*ww*ww)
    elif cond == 9 :
        alpha = ww * np.sqrt(rhoPlus - rhoMoins)
        kBord = alpha
    elif cond == 10 :
        alpha = ww * mp.sqrt(rhoPlus - rhoMoins)
        kBord = alpha/2-(rhoMoins*ww*ww-beta*beta*ww*ww)/(2*alpha)
    elif cond == 11 :
        alpha = ww * mp.sqrt(rhoPlus - rhoMoins)
        kBord = alpha/4-(rhoMoins*ww*ww-beta*beta*ww*ww)/(4*alpha) - (rhoMoins*ww*ww-beta*beta*ww*ww)/(alpha + (-rhoMoins*ww*ww + beta*beta*ww*ww)/alpha)
    elif cond == 12 :
        alpha = ww * mp.sqrt(rhoPlus - rhoMoins)
        kBord = alpha/8-(rhoMoins*ww*ww-beta*beta*ww*ww)/(8*alpha) - (rhoMoins*ww*ww-beta*beta*ww*ww)/(2.0*alpha - 2.0*(rhoMoins*ww*ww - beta*beta*ww*ww)/alpha) \
                     - (rhoMoins*ww*ww-beta*beta*ww*ww)/(alpha/2.0 - (rhoMoins*ww*ww - beta*beta*ww*ww)/(2.0*alpha) - 2.0*(rhoMoins*ww*ww - beta*beta*ww*ww)/(alpha - (rhoMoins*ww*ww - beta*beta*ww*ww)/alpha))
    A = createMatA(k1, k2, a, ray, cond, kBord)
    val, EL, vect = mp.eig(A, left = True, right = True)
    print("Valeurs :", val)
    if(abs(val[0])<abs(val[1])) :
        if (len(val) == 3) and (abs(val[2])<abs(val[0])) :
            print("valeur propre nulle = ", val[2])
            alphas = vect[:,2]
        else :
            print("valeur propre nulle = ", val[0])
            alphas = vect[:,0]
    else :
        if (len(val) == 3) and (abs(val[2])<abs(val[1])) :
            print("valeur propre nulle = ", val[2])
            alphas = vect[:,2]
        else :
            print("valeur propre nulle = ", val[1])
            alphas = vect[:,1]
    return k1, k2, alphas

def drawSolution(rMin, rMax, nb, alphas, k1, k2, a) :
    rr = np.linspace(rMin, rMax, nb)
    coef = 1/(alphas[0]*mp.besselj(0, 0))
    if len(alphas) == 2 :
        uu = np.array([real(mp.besselj(0, k1*r) if r <= a else alphas[1]/alphas[0]*mp.hankel1(0, k2*r)) for r in rr])
    else :
        uu = np.array([real(mp.besselj(0, k1*r) if r <= a else alphas[1]/alphas[0]*mp.besselj(0, k2*r) + alphas[2]/alphas[0]*mp.bessely(0, k2*r)) for r in rr])
    return rr, uu

def ComputeModesParam(ww, rhoMoins, rhoPlus, a, cond) :
    # print('ComputeModes')
    # print('Type de rhoPlus :', type(rhoPlus))
    # print('Type de rhoMoins :', type(rhoMoins))
    # print('Type de a :', type(a))
    # print('Type de ww :', type(ww))
    mp.mp.dps = 200
    ptInit = mp.mpc(1.4,0.0)
    Beta = Newton(ptInit, 50, ww, rhoMoins, rhoPlus, a, cond)
    mp.mp.dps = 15
    print("Résultat : ", Beta, ", Det = ", ComputeDet(Beta, ww, rhoMoins, rhoPlus, a, ray, cond))
    mp.mp.dps = 200
    k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, cond)
    return real(Beta), k1, k2, alphas

# programme principal (le faire tourner donne le mode propagatif)
def ComputeModes(ww, rhoMoins, rhoPlus, a, ray, rMin, rMax, nb, cond) :
    mp.mp.dps = 200
    ptInit = mp.mpc(2.9,0.0)
    Beta = Newton(ptInit, 50, ww, rhoMoins, rhoPlus, a, ray , cond)
    mp.mp.dps = 15
    print("Résultat : ", Beta, ", Det = ", ComputeDet(Beta, ww, rhoMoins, rhoPlus, a, ray, cond))
    mp.mp.dps = 200
    k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, ray, cond)
    print("Alpha, k = ", alphas, k1, k2)


    return drawSolution(rMin, rMax, nb, alphas, k1, k2, a)

def ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, rMin, rMax, nb, cond) :
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, ray, cond)
    mp.mp.dps = 200
    xMin = np.float64(1.0)
    xMax = np.float64(10.0)
    Betas = zer0s(fCLA, xMin, xMax)
    mp.mp.dps = 15


    if len(Betas) > 0 :
        # print("Mode le plus propagatif")
        Beta = Betas[-1]
        print("Résultat : ", Beta)
        mp.mp.dps = 200
        k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, ray, cond)
        print("Alpha, k = ", alphas, k1, k2)
        return drawSolution(rMin, rMax, nb, alphas, k1, k2, a)
    else :
        print("Aucun mode trouvé...")
        return np.linspace(rMin, rMax, nb), 0.0 * np.linspace(rMin, rMax, nb)


def dichotomie(f, a, b, tol) :
    croissant = (f(b) > f(a))
    u = (a+b)/2
    while abs(f(u)) > tol :
        # print(real(a), real(b), real(u), f(u))
        # print("    ", abs(f(u)))
        # print("    Type a : ", type(a))
        # print("    Type b : ", type(b))
        # print("    Type u : ", type(u))
        # print(real(a), real(b), real(u), f(u))
        if croissant :
            if f(u) > 0 :
                b = u
            else :
                a = u
        else :
            if f(u) < 0 :
                b = u
            else :
                a = u
        u = (a+b)/2
    # return np.float64(real(u))
    return(real(u))

def zer0s(f, a, b) :
    # print(type(a))
    U = np.linspace(a, b, 100) # evaluate function at 100 different points
    U = [mp.mpc(x, 0.0) for x in np.linspace(a, b, 100)]
    c = [f(mp.mpc(x, 0.0)) for x in np.linspace(a, b, 100)]
    s = np.sign(c)
    zer0 = []
    for i in range(100-1):
        if s[i] + s[i+1] == 0: # opposite signs
            u = dichotomie(f, U[i], U[i+1], 1e-36)
            zer0.append(u)
    return zer0


if __name__ == "__main__":
    print("1) Calcul du mode")
    print("2) Étude du déterminant")
    print("3) Brefs calculs d'erreurs relatives à rayon de CLA constant")
    print("4) Étude d'erreurs relatives (plusieurs CLA)")
    print("5) Étude d'erreurs relatives (plusieurs fréquences)")
    choix = int(input("Votre choix : "))
    if choix == 1 :
        ww = mp.mpf(10.0)
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)
        cond = 8

        rr, uu = ComputeModes(ww, rhoMoins, rhoPlus, a, 0.4, 0, 0.4, 5000, cond)
        # rr1, uu1 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 0.7, 0, 0.7, 5000, cond)
        # rr2, uu2 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 1.0, 0, 1.0, 5000, cond)
        # rr3, uu3 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 1.5, 0, 1.5, 5000, cond)
        # rrS, uuSomm = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 0.4, 0, 0.4, 5000, 0)
        matplotlib.rcParams.update({'font.size': 20})
        # plt.plot(rr, uuSomm, "b-")
        plt.plot(rr, uu, "r--")
        # plt.plot(rr1, uu1, "g--")
        # plt.plot(rr2, uu2, "b--")
        # plt.plot(rr3, uu3, "r-.")
        plt.xlabel("Rayon")
        plt.ylabel("Mode u")
        plt.show()
    elif choix == 2 :
        ww = 2
        xMin = 1
        xMax = np.sqrt(10)
        R = 0.4
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)
        B = np.linspace(xMin, xMax, 10000)
        print("Open")
        fSom = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 0)
        Somi = [fSom(x) for x in B]

        # zer0_Open = np.array(zer0s(fSom, xMin, xMax))
        # print(zer0_Open)
        # CLAi = [imag(ComputeDetTrueABC(mp.mpc(x, 0.0), ww, rhoMoins, rhoPlus, a, 1.0)) for x in B]

        # print("True CLA")
        #
        # # fCLA = lambda x : np.float64(imag(ComputeDetTrueABC(x, ww, rhoMoins, rhoPlus, a, R)))
        # fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 1)
        # zer0TrueABC = np.array(zer0s(fCLA, xMin, xMax))
        # print(zer0TrueABC)
        #
        # # print("BGT 1")
        # #
        # # fBGT1 = lambda x : np.float64(real(ComputeDet_BGT1(x, ww, rhoMoins, rhoPlus, a, R)))
        # # zer0_BGT1 = np.array(zer0s(fBGT1, xMin, xMax))
        # # print(zer0_BGT1)
        #
        print("BGT 2")

        fBGT2 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 8)
        # BGT2 = [fBGT2(x) for x in B]
        zer0_BGT2 = np.array(zer0s(fBGT2, xMin, xMax))
        print(zer0_BGT2)
        #
        print("BGT 2 Newton 0")

        fBGT2_0 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 9)
        # BGT20 = [fBGT2_0(x) for x in B]
        zer0_BGT2_0 = np.array(zer0s(fBGT2_0, xMin, xMax))
        print(zer0_BGT2_0)

        print("BGT 2 Newton 1")

        fBGT2_1 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 10)
        zer0_BGT2_1 = np.array(zer0s(fBGT2_1, xMin, xMax))
        print(zer0_BGT2_1)
        # BGT21 = [fBGT2_1(x) for x in B]
        #
        print("BGT 2 Newton 2")
        #
        fBGT2_2 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 11)
        # BGT22 = [fBGT2_2(x) for x in B]
        zer0_BGT2_2 = np.array(zer0s(fBGT2_2, xMin, xMax))
        print(zer0_BGT2_2)
        #
        print("BGT 2 Newton 3")
        #
        fBGT2_3 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, R, 12)
        # BGT22 = [fBGT2_2(x) for x in B]
        zer0_BGT2_3 = np.array(zer0s(fBGT2_3, xMin, xMax))
        print(zer0_BGT2_3)

        # fCLA1 = lambda x : np.float64(real(ComputeDet1ABC(x, ww, rhoMoins, rhoPlus, a, R)))
        # zer0ABC1 = np.array(zer0s(fCLA1, xMin, xMax))
        # print(zer0ABC1)
        # fCLA2 = lambda x : np.float64(real(ComputeDet_NewtonABC(x, ww, rhoMoins, rhoPlus, a, R, 1)))
        # zer0ABC2 = np.array(zer0s(fCLA2, xMin, xMax))
        # print(zer0ABC2)
        # fCLA3 = lambda x : np.float64(real(ComputeDet_NewtonABC(x, ww, rhoMoins, rhoPlus, a, R, 2)))
        # zer0ABC3 = np.array(zer0s(fCLA3, xMin, xMax))
        # print(zer0ABC3)
        # fDir = lambda x : np.float64(real(ComputeDetDirichlet(x, ww, rhoMoins, rhoPlus, a, R)))
        # zer0Dir = np.array(zer0s(fDir, xMin, xMax))
        # print(zer0Dir)
        # fK0 = lambda x : np.float64(real(ComputeDetK0(x, ww, rhoMoins, rhoPlus, a, R)))
        # zer0K0 = np.array(zer0s(fK0, xMin, xMax))
        # print(zer0K0)

        matplotlib.rcParams.update({'font.size': 20})
        plt.plot(B*ww, Somi, label="Open", linewidth=2)
        # plt.plot(B*ww,[fK0(x) for x in B], label="K0 (R = {})".format(R), linewidth=2)
        plt.plot(B*ww, 0.0*B*ww, "r--", linewidth=2)

        # plt.plot(zer0TrueABC*ww,0.0*zer0TrueABC,'x', color="red", label="True ABC (R = {})".format(R), markersize=16, markeredgewidth = 3, marker="s")
        # plt.plot(zer0_BGT1*ww,0.0*zer0_BGT1,'x', color="red", label="BGT 1 (R = {})".format(R), markersize=16, markeredgewidth = 3)




        plt.plot(zer0_BGT2*ww,0.0*zer0_BGT2,'o', color="cyan", label="BGT 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        plt.plot(zer0_BGT2_0*ww,0.0*zer0_BGT2_0,'^', color="black", label="BGT 2 Newton 0 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        plt.plot(zer0_BGT2_1*ww,0.0*zer0_BGT2_1,'s', color="blue", label="BGT 2 Newton 1 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        plt.plot(zer0_BGT2_2*ww,0.0*zer0_BGT2_2,'v', color="magenta", label="BGT 2 Newton 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        plt.plot(zer0_BGT2_3*ww,0.0*zer0_BGT2_3,'x', color="green", label="BGT 2 Newton 3 (R = {})".format(R), markersize=16, markeredgewidth = 3)




        # plt.plot(zer0ABC1*ww, 0.0*zer0ABC1,'x', color="magenta", label="ABC Robin (R = {})".format(R), markersize=16, markeredgewidth = 3, marker="v")
        # plt.plot(zer0ABC2*ww, 0.0*zer0ABC2, 'x', color="cyan", label="ABC Newton 1 (R = {})".format(R), markersize=16, markeredgewidth = 3, marker="^")
        # plt.plot(zer0ABC3*ww, 0.0*zer0ABC3, 'x', color="green", label="ABC Newton 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(zer0Dir*ww, 0.0*zer0Dir, 'x', color="black", label="Dirichlet (R = {})".format(R), markersize=16, markeredgewidth = 3, marker="d")
        # plt.plot(zer0K0*ww, 0.0*zer0K0, 'x', label="K0 (R = {})".format(R), markersize=16, markeredgewidth = 3)


        # plt.plot(B*ww,[fCLA(x) for x in B], label="True ABC (R = {})".format(R))
        # plt.plot(B*ww,[fCLA1(x) for x in B], label="ABC Robin (R = {})".format(R))
        # plt.plot(B*ww, 0.0*zer0ABC1,'x', label="ABC Robin (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(B*ww, 0.0*zer0ABC2, 'x', label="ABC Newton 1 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(B*ww, 0.0*zer0ABC3, 'x', label="ABC Newton 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(B*ww, 0.0*zer0Dir, 'x', label="Dirichlet (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.rcParams.update({
        #     "text.usetex": True,
        #     "font.family": "monospace",
        #     "font.monospace": 'Computer Modern Typewriter',
        #     'figure.figsize' : (10,6)
        # })
        plt.rcParams["text.usetex"] = True
        plt.xlabel(r'$\beta$')
        plt.ylabel(r'$det(A(\beta))$')
        plt.legend()
        plt.show()
    elif choix == 3 :
        ww = mp.mpf(2.0)
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)

        print("Sommerfeld")
        print()
        # rr, uu = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 0)
        rr, uu = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 0)
        print()
        print("Non-linear")
        print()
        # rr, uu1 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 1)
        rr, uu1 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 1)
        # print()
        # print("Robin")
        # print()
        # # rr, uu2 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 2)
        # rr, uu2 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 2)
        # print()
        # print("Newton 1")
        # print()
        # # rr, uu3 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 3)
        # rr, uu3 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 3)
        # print()
        # print("Newton 2")
        # print()
        # # rr, uu4 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 4)
        # rr, uu4 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 4)
        # print()
        # print("Dirichlet")
        # print()
        # # rr, uu5 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 5)
        # rr, uu5 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 5)
        print()
        print("K0")
        print()
        # rr, uu5 = ComputeModes(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 5)
        rr, uu6 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 5000, 6)

        print()
        print("Erreurs relatives")
        diff = uu1 - uu
        eAbs = np.dot((diff),np.conj(diff))
        eAn = np.dot((uu),np.conj(uu))
        eRel = np.sqrt(eAbs / eAn)
        print("Non-linear ABC :", np.format_float_scientific(np.float64(eRel), precision = 1, exp_digits=3))

        # diff = uu2 - uu
        # eAbs = np.dot((diff),np.conj(diff))
        # eRel = np.sqrt(eAbs / eAn)
        # print("ABC Robin :", np.format_float_scientific(np.float64(eRel), precision = 1, exp_digits=3))
        #
        # diff = uu3 - uu
        # eAbs = np.dot((diff),np.conj(diff))
        # eRel = np.sqrt(eAbs / eAn)
        # print("ABC Newton 1 :", np.format_float_scientific(np.float64(eRel), precision = 1, exp_digits=3))
        #
        # diff = uu4 - uu
        # eAbs = np.dot((diff),np.conj(diff))
        # eRel = np.sqrt(eAbs / eAn)
        # print("ABC Newton 2 :", np.format_float_scientific(np.float64(eRel), precision = 1, exp_digits=3))
        #
        # diff = uu5 - uu
        # eAbs = np.dot((diff),np.conj(diff))
        # eRel = np.sqrt(eAbs / eAn)
        # print("Dirichlet :", np.format_float_scientific(np.float64(eRel), precision = 3))

        diff = uu6 - uu
        eAbs = np.dot((diff),np.conj(diff))
        eRel = np.sqrt(eAbs / eAn)
        print("K0 :", np.format_float_scientific(np.float64(eRel), precision = 3))

        # matplotlib.rcParams.update({'font.size': 20})
        # plt.plot(rr, uu, "r-")
        # plt.plot(rr, uu1)
        # plt.xlabel("Rayon")
        # plt.ylabel("Mode u")
        # # plt.vlines(0.2, 0, 0.01, linestyles = 'dashed')
        # plt.show()
        #
        # plt.plot(rr, uu, "r-")
        # plt.plot(rr, uu2)
        # plt.xlabel("Rayon")
        # plt.ylabel("Mode u")
        # # plt.vlines(0.2, 0, 0.01, linestyles = 'dashed')
        # plt.show()
        #
        # plt.plot(rr, uu, "r-")
        # plt.plot(rr, uu3)
        # plt.xlabel("Rayon")
        # plt.ylabel("Mode u")
        # # plt.vlines(0.2, 0, 0.01, linestyles = 'dashed')
        # plt.show()
        #
        # plt.plot(rr, uu, "r-")
        # plt.plot(rr, uu4)
        # plt.xlabel("Rayon")
        # plt.ylabel("Mode u")
        # # plt.vlines(0.2, 0, 0.01, linestyles = 'dashed')
        # plt.show()
        #
        # plt.plot(rr, uu, "r-")
        # plt.plot(rr, uu5)
        # plt.xlabel("Rayon")
        # plt.ylabel("Mode u")
        # # plt.vlines(0.2, 0, 0.01, linestyles = 'dashed')
        # plt.show()

    elif choix == 4 :
        ww = mp.mpf(2.0)
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)
        # rays = [0.4, 0.7, 1.0, 1.5, 2.0]
        rays = np.linspace(0.4, 1.0, 30)
        eRelTrueABC = []
        eRelBGT2 = []
        eRel1 = []
        eRel2 = []
        eRel3 = []
        eRel4 = []

        for ray in rays :
            print()
            print(ray)
            print()
            rr, uu = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 0)
            # rr, uu0 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 1)

            print()
            print("BGT2, rayon :", ray)
            print()
            rr, uu1 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 8)

            diff = uu1 - uu
            eAbs = np.dot((diff),np.conj(diff))
            eAn = np.dot((uu),np.conj(uu))
            eRelBGT2.append(np.sqrt(eAbs / eAn))

            print()
            print("Erreur BGT2, rayon", ray, ":", np.sqrt(eAbs / eAn))
            print()

            print()
            print("Approx BGT2N0, rayon :", ray)
            print()
            rr, uu2 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 9)


            diff = uu2 - uu
            eAbs = np.dot((diff),np.conj(diff))
            eRel1.append(np.sqrt(eAbs / eAn))
            print()
            print("Erreur BGT2N0, rayon", ray, ":", np.sqrt(eAbs / eAn))
            print()

            print()
            print("Approx BGT2N1, rayon :", ray)
            print()
            rr, uu3 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 10)

            diff = uu3 - uu
            eAbs = np.dot((diff),np.conj(diff))
            eRel2.append(np.sqrt(eAbs / eAn))
            print()
            print("Erreur BGT2N1, rayon", ray, ":", np.sqrt(eAbs / eAn))
            print()

            print()
            print("Approx BGT2N2, rayon :", ray)
            print()
            rr, uu4 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 11)

            diff = uu4 - uu
            eAbs = np.dot((diff),np.conj(diff))
            eRel3.append(np.sqrt(eAbs / eAn))
            print()
            print("Erreur BGT2N2, rayon", ray, ":", np.sqrt(eAbs / eAn))
            print()

            print()
            print("Approx BGT2N3, rayon :", ray)
            print()
            rr, uu5 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 12)

            diff = uu5 - uu
            eAbs = np.dot((diff),np.conj(diff))
            eRel4.append(np.sqrt(eAbs / eAn))
            print()
            print("Erreur BGT2N3, rayon", ray, ":", np.sqrt(eAbs / eAn))
            print()
            # diff = uu0 - uu
            # eAbs = np.dot((diff),np.conj(diff))
            # eAn = np.dot((uu),np.conj(uu))
            # eRelTrueABC.append(np.sqrt(eAbs / eAn))



            #


        matplotlib.rcParams.update({'font.size': 20})
        # plt.semilogy(rays, eRelTrueABC, 'o-', label='Fixed Sommerfeld')
        plt.semilogy(rays, eRelBGT2, 'o-', label='ABC BGT2')
        plt.semilogy(rays, eRel1, 'o-', label='ABC BGT2 Newton 0')
        plt.semilogy(rays, eRel2, 'o-', label='ABC BGT2 Newton 1')
        plt.semilogy(rays, eRel3, 'o-', label='ABC BGT2 Newton 2')
        plt.semilogy(rays, eRel4, 'o-', label='ABC BGT2 Newton 3')
        plt.xlabel("Boundary radius")
        plt.ylabel("Relative error")
        plt.legend()
        plt.show()
    elif choix == 5 :
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)
        freq = [2.0, 3.0, 5.0, 10.0]
        # freq = [mp.mpf(10.0)]
        rays = [0.21, 0.3, 0.7, 1.0, 1.5, 2.0]
        eRays = []
        eRels = []
        for ww in freq :
            print("Fréquence :", ww)
            rr, uu = ComputeModes1D(ww, rhoMoins, rhoPlus, a, rayonCLA, 0, rayonCLA, 100, 0)
            if uu[0] == 0 :
                print("Pas de mode...")
                exit()
            else :
                print("OK")
                # eRays.append([])
                # eRels.append([])
                eRays = []
                eRels = []
                for ray in rays :
                    print("Rayon :", ray)
                    rr, uu = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 0)
                    rr, uuCLA = ComputeModes1D(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 6)
                    if uuCLA[0] == 0 :
                        print("Pas de mode approché...")
                    else :
                        diff = uuCLA - uu
                        eAbs = np.dot((diff),np.conj(diff))
                        eAn = np.dot((uu),np.conj(uu))
                        eRays.append(ray)
                        eRels.append(np.sqrt(eAbs / eAn))
                        # eRays[-1].append(ray)
                        # eRels[-1].append(np.sqrt(eAbs / eAn))
                plt.semilogy(eRays, eRels, 'o-', label='omega = {}'.format(ww))
        plt.xlabel("Rayon de la frontière")
        plt.ylabel("Erreur relative")
        plt.legend()
        plt.show()
        # for ray in rays :
        #     rr, uu = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 0)
        #     rr, uu1 = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 1)
        #     rr, uu2 = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 2)
        #     rr, uu3 = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 3)
        #     rr, uu4 = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 4)
        #     rr, uu5 = ComputeModes(ww, rhoMoins, rhoPlus, a, ray, 0, ray, 5000, 5)
        #
        #     diff = uu1 - uu
        #     eAbs = np.dot((diff),np.conj(diff))
        #     eAn = np.dot((uu),np.conj(uu))
        #     eRelTrueABC.append(np.sqrt(eAbs / eAn))
        #
        #     diff = uu2 - uu
        #     eAbs = np.dot((diff),np.conj(diff))
        #     eRel1.append(np.sqrt(eAbs / eAn))
        #
        #     diff = uu3 - uu
        #     eAbs = np.dot((diff),np.conj(diff))
        #     eRel2.append(np.sqrt(eAbs / eAn))
        #
        #     diff = uu4 - uu
        #     eAbs = np.dot((diff),np.conj(diff))
        #     eRel3.append(np.sqrt(eAbs / eAn))
        #
        #     diff = uu5 - uu
        #     eAbs = np.dot((diff),np.conj(diff))
        #     eRelDir.append(np.sqrt(eAbs / eAn))
        #
        # plt.semilogy(rays, eRelTrueABC, 'o-', label='Non-linear ABC')
        # plt.semilogy(rays, eRel1, 'o-', label='ABC Robin')
        # plt.semilogy(rays, eRel2, 'o-', label='ABC Newton 1')
        # plt.semilogy(rays, eRel3, 'o-', label='ABC Newton 2')
        # plt.semilogy(rays, eRelDir, 'o-', label='Dirichlet')
        # plt.xlabel("Rayon de la frontière")
        # plt.ylabel("Erreur relative")
        # plt.legend()
        # plt.show()

    # Test
    # k1, k2, alphas = ComputeParam(mp.mpc(2.04, 0.0), ww, rhoMoins, rhoPlus, a)
    # rr, uu = drawSolution(0.0, 2.0, 200, alphas, k1, k2, a)

    # xMin = 1.0
    # xMax = 3.0
    # yMin = -1.0
    # yMax = 1.0
    # X = np.linspace(xMin, xMax, 10)
    # Y = np.linspace(yMin, yMax, 10)
    # xx,yy = np.meshgrid(X,Y)
    # print(xx)
    # print(yy)
    # bb = np.zeros(xx.shape)
    # # print(bb)
    # i = 0
    # for x in X :
    #     print("i = ", i)
    #     j = 0
    #     for y in Y :
    #         bb[i,j] = abs(ComputeDet(mp.mpc(x, y), ww, rhoMoins, rhoPlus, a))
    #         j+=1
    #     i+=1
    # print(bb)-1]
    # # print(abs(ComputeDet(mp.mpc(2.0, 1.0), ww, rhoMoins, rhoPlus, a)))
    # plt.imshow(bb.T, extent=(xMin, xMax,yMax,yMin))
    # plt.colorbar()
    # plt.show()

    # diff = uNum - uAn
    # eAbs = np.dot((MMref * diff),np.conj(diff))
    # eAn = np.dot((MMref * uAn),np.conj(uAn))
    # eRel = eAbs / eAn



    # ww = 10
    # xMin = 1.0
    # xMax = np.sqrt(10.0)
    # rhoMoins = mp.mpf(1.0)
    # rhoPlus = mp.mpf(10.0)
    # a = mp.mpf(0.2)
    # B = np.linspace(xMin, xMax, 200)
    # Jrr = [real(mp.bessely(1, mp.mpc(0.0, x)) + 1j*mp.bessely(0, mp.mpc(0.0, x))) for x in B]
    # Jri = [imag(mp.bessely(1, mp.mpc(0.0, x)) + 1j*mp.bessely(0, mp.mpc(0.0, x))) for x in B]
    # # Jrr = [real(mp.bessely(1, mp.mpc(x, 0))) for x in B]
    # # Jri = [imag(mp.bessely(1, mp.mpc(x, 0))) for x in B]
    # # print(abs(ComputeDet(mp.mpc(2.0, 1.0), ww, rhoMoins, rhoPlus, a)))
    # plt.plot(B, Jrr)
    # plt.plot(B, Jri)
    # plt.show()
