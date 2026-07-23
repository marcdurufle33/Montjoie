from pylab import *
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp


pasNewton = 1e-5
rayonCLA = 0.4

def createMatASommerfeld(k1, k2, a, m) :
    A = np.zeros((2,2), dtype = 'complex_')
    A[0, 0] = sp.jn(m, k1*a);
    A[0, 1] = -sp.hankel1(m, k2*a);
    A[1, 0] = -k1*sp.jn(m+1, k1*a) + (m/a) * sp.jn(m, k1*a);
    A[1, 1] = k2*sp.hankel1(m+1, k2*a) - (m/a) * sp.hankel1(m, k2*a);
    return A

def createMatA_BGT1(k1, k2, a, m, R) :
    A = np.zeros((3,3), dtype = 'complex_')
    A[0, 0] = sp.jn(m, k1*a);
    A[0, 1] = -sp.jn(m, k2*a);
    A[0, 2] = -sp.yv(m, k2*a);
    A[1, 0] = - k1*sp.jn(m+1, k1*a) + (m/a) * sp.jn(m, k1*a);
    A[1, 1] = k2*sp.jn(m+1, k2*a) - (m/a) * sp.jn(m, k2*a);
    A[1, 2] = k2*sp.yv(m+1, k2*a) - (m/a) * sp.yv(m, k2*a);

    X = (1/(2*R) - 1j*k2)

    A[2, 0] = 0.0;
    A[2, 1] = -k2*sp.jn(m+1, k2*R) + (m/R) * sp.jn(m, k2*R) + X * sp.jn(m, k2*R);
    A[2, 2] = -k2*sp.yv(m+1, k2*R) + (m/R) * sp.yv(m, k2*R) + X * sp.yv(m, k2*R);
    return A

def createMatA_BGT2(k1, k2, a, m, R) :
    A = np.zeros((3,3), dtype = 'complex_')
    A[0, 0] = sp.jn(m, k1*a);
    A[0, 1] = -sp.jn(m, k2*a);
    A[0, 2] = -sp.yv(m, k2*a);
    A[1, 0] = - k1*sp.jn(m+1, k1*a) + (m/a) * sp.jn(m, k1*a);
    A[1, 1] = k2*sp.jn(m+1, k2*a) - (m/a) * sp.jn(m, k2*a);
    A[1, 2] = k2*sp.yv(m+1, k2*a) - (m/a) * sp.yv(m, k2*a);

    num = 3/4 - 3*1j*R*k2 - 2*R*R*k2*k2
    denom = 2 * R * (1 - 1j * k2 * R)

    X = num / denom
    Y = 1.0 / denom

    A[2, 0] = 0.0;
    A[2, 1] = -k2*sp.jn(m+1, k2*R) + (m/R) * sp.jn(m, k2*R) + X * sp.jn(m, k2*R) + m*m*Y*sp.jn(m, k2*R);
    A[2, 2] = -k2*sp.yv(m+1, k2*R) + (m/R) * sp.yv(m, k2*R) + X * sp.yv(m, k2*R) + m*m*Y*sp.yv(m, k2*R);
    return A

def det2D(A) :
    return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]


def det3D(A) :
    return A[0, 0]*A[1, 1]*A[2, 2] + A[0, 1]*A[1, 2]*A[2, 0] + A[0, 2]*A[1, 0]*A[2, 1] - A[0, 2]*A[1, 1]*A[2,0] - A[0, 1]*A[1, 0]*A[2,2] - A[0, 0]*A[1, 2]*A[2,1]

# toutes les valeurs doivent être en mp précision
def ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a, m):
    k1 = np.emath.sqrt(rhoPlus*ww*ww - beta*beta*ww*ww)
    k2 = np.emath.sqrt(rhoMoins*ww*ww - beta*beta*ww*ww)
    A = createMatASommerfeld(k1, k2, a, m)
    return det2D(A)

def ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, m, R):
    k1 = np.emath.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = np.emath.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT1(k1, k2, a, m, R)
    return det3D(A)

def ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, m, R):
    k1 = np.emath.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = np.emath.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT2(k1, k2, a, m, R)
    return det3D(A)


def createMatA(k1, k2, a, m, ray, cond, kBord = 0) :
    if cond == 0 :
        return createMatASommerfeld(k1, k2, a, m)
    elif cond == 1 :
        return createMatA_BGT1(k1, k2, a, m, ray)
    elif cond == 2 :
        return createMatA_BGT2(k1, k2, a, m, ray)

def ComputeDet(beta, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    if cond == 0 :
        return ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a, m)
    elif cond == 1 :
        return ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, m, ray)
    elif cond == 2 :
        return ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, m, ray)

def ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    if cond == 0 :
        return lambda x : np.float64(real(ComputeDetSommerfeld(x, ww, rhoMoins, rhoPlus, a, m) * ((-1j)**(m+1))))
    elif cond == 1 :
        return lambda x : np.float64(real(ComputeDet_BGT1(x, ww, rhoMoins, rhoPlus, a, m, ray)))
    elif cond == 2 :
        return lambda x : np.float64(real(ComputeDet_BGT2(x, ww, rhoMoins, rhoPlus, a, m, ray)))

# def Jac(x, y, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
#     J = np.zeros((2,2), dtype = 'complex_')
#     # schéma centré (ordre 2) pour les dérivées
#     DxM1y = ComputeDet(x - pasNewton + 1j * y, ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     DxP1y = ComputeDet(mp.mpc(x + pasNewton, y), ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     DxyM1 = ComputeDet(mp.mpc(x, y - pasNewton), ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     DxyP1 = ComputeDet(mp.mpc(x, y + pasNewton), ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     J[0, 0] = 0.5 * (mp.re(DxP1y) - mp.re(DxM1y)) / pasNewton
#     J[1, 0] = 0.5 * (mp.im(DxP1y) - mp.im(DxM1y)) / pasNewton
#     J[0, 1] = 0.5 * (mp.re(DxyP1) - mp.re(DxyM1)) / pasNewton
#     J[1, 1] = 0.5 * (mp.im(DxyP1) - mp.im(DxyM1)) / pasNewton
#     return J
#
#
# def Newton(B0, expTol, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
#     mp.mp.dps = 200
#     # print('Newton')
#     # print('Type de rhoPlus :', type(rhoPlus))
#     # print('Type de rhoMoins :', type(rhoMoins))
#     # print('Type de a :', type(a))
#     # print('Type de ww :', type(ww))
#     print("B0 : ", B0)
#     x = mp.re(B0); y = mp.im(B0)
#     X0 = mp.matrix([x, y])
#     DComp = ComputeDet(mp.mpc(x, y), ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
#     Jac0 = Jac(x, y, ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     i = 0
#     vectDiff = mp.lu_solve(Jac0, D0)
#     norme2 = (vectDiff.T * vectDiff)[0]
#     X0 = X0 - vectDiff
#     # mp.mp.dps = 15 # pour l'affichage
#     print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
#     while norme2 >= 10**(-expTol * 2) :
#         i += 1
#         DComp = ComputeDet(mp.mpc(X0[0], X0[1]), ww, rhoMoins, rhoPlus, a, m, ray, cond)
#         D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
#         Jac0 = Jac(X0[0], X0[1], ww, rhoMoins, rhoPlus, a, m, ray, cond)
#         vectDiff = mp.lu_solve(Jac0, D0)
#         norme2 = (vectDiff.T * vectDiff)[0]
#         X0 = X0 - vectDiff
#         # mp.mp.dps = 15 # pour l'affichage
#         print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
#         mp.mp.dps = 100
#     return X0[0] +1j * X0[1]
    # print(X)

def ComputeParam(beta, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    k1 = np.emath.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = np.emath.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    kBord = 0
    # if cond == 2 :
    #     kBord = mp.sqrt(rhoMoins*ww*ww - rhoPlus*ww*ww)
    A = createMatA(k1, k2, a, m, ray, cond, kBord)
    print(A)
    val, vect = np.linalg.eig(A)
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

def drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m) :
    rr = np.linspace(rMin, rMax, nb)
    if len(alphas) == 2 :
        uu = np.array([sp.jn(m, k1*r) if r <= a else alphas[1]/alphas[0]*sp.hankel1(m, k2*r) for r in rr])
    else :
        uu = np.array([sp.jn(m, k1*r) if r <= a else alphas[1]/alphas[0]*sp.jn(m, k2*r) + alphas[2]/alphas[0]*sp.yv(m, k2*r) for r in rr])
    return rr, uu

def drawSolution3D(rMax, nb, alphas, k1, k2, a, m) :
    rr = np.linspace(-rMax, rMax, nb)
    XX, YY = np.meshgrid(rr, rr)
    RR = np.sqrt(XX**2 + YY**2)
    # TT = np.arctan2(YY, XX)
    # imshow(TT)
    # show()
    TT = -np.arccos(XX/RR) * (np.sign(YY)+(YY == 0.0))
    TT[nb//2, nb//2] = 0.0
    computeAnalyCable = lambda r : sp.jn(m, k1*r)
    A = computeAnalyCable(RR) * np.exp(1j * m * TT)
    computeAnalyDehors = lambda r :  alphas[1]/alphas[0]*sp.hankel1(m,k2*r)
    if len(alphas) == 3 :
        computeAnalyDehors = lambda r :  alphas[1]/alphas[0]*sp.jn(m,k2*r) + alphas[2]/alphas[0]*sp.yv(m,k2*r)
    B = computeAnalyDehors(RR) * np.exp(1j * m * TT)
    B[nb//2, nb//2] = 0.0
    PetitRayon = (RR <= a)
    GrandRayon = (RR > a) * (RR <= rMax)
    UU = (A * PetitRayon + B * GrandRayon)

    # imshow(RR <= rMax)
    # show()
    return XX, YY, UU

# def ComputeModesParam(ww, rhoMoins, rhoPlus, a, m, cond) :
#     # print('ComputeModes')
#     # print('Type de rhoPlus :', type(rhoPlus))
#     # print('Type de rhoMoins :', type(rhoMoins))
#     # print('Type de a :', type(a))
#     # print('Type de ww :', type(ww))
#     ptInit = mp.mpc(1.4,0.0)
#     Beta = Newton(ptInit, 50, ww, rhoMoins, rhoPlus, a, m, cond)
#     print("Résultat : ", Beta, ", Det = ", ComputeDet(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond))
#     k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, cond)
#     return real(Beta), k1, k2, alphas

# programme principal (le faire tourner donne le mode propagatif)
# def ComputeModes(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond) :
#     ptInit = mp.mpc(2.6,0.0)
#     Beta = Newton(ptInit, 50, ww, rhoMoins, rhoPlus, a, m, ray , cond)
#     print("Résultat : ", Beta, ", Det = ", ComputeDet(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond))
#     k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond)
#     print("Alpha, k = ", alphas, k1, k2)
#
#
#     return drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m)

def ComputeModes1D(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond) :
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond)
    xMin = np.float64(1.0)
    xMax = np.float64(np.sqrt(10.0))
    Betas = zer0s(fCLA, xMin, xMax)
    print(Betas)


    if len(Betas) > 0 :
        # print("Mode le plus propagatif")
        Beta = Betas[0]
        print("Résultat : ", Beta)
        k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond)
        print("Alpha, k = ", alphas, k1, k2)
        return drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m)
    else :
        print("Aucun mode trouvé...")
        return np.linspace(rMin, rMax, nb), 0.0 * np.linspace(rMin, rMax, nb)

def ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond, _2D = True) :
    # print("Compute")
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond)
    # print("Init")
    xMin = np.float64(1.01)
    xMax = np.float64(np.sqrt(10.0))
    # print("Calcul betas")
    Betas = zer0s(fCLA, xMin, xMax)
    print(Betas)

    print("Avant boucle")

    if len(Betas) > 0 :
        eigenfunctions = []
        for Beta in Betas :
            print()
            print("Beta = ", Beta)
            k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond)
            print("Alpha, k = ", alphas, k1, k2)
            print()
            if _2D :
                eigenfunctions.append(drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m))
            else :
                eigenfunctions.append(drawSolution3D(rMax, nb, alphas, k1, k2, a, m))
        return Betas, eigenfunctions
    else :
        print("Aucun mode trouvé...")
        return np.linspace(rMin, rMax, nb), 0.0 * np.linspace(rMin, rMax, nb)


def dichotomie(f, a, b, tol) :
    # a = np.longdouble(a)
    # b = np.longdouble(b)
    mf = max(abs(f(a)),abs(f(b)))
    fnorm = lambda x : f(x)/mf
    croissant = (f(b) > f(a))
    u = (a+b)/2
    compteur = 0
    while (abs(fnorm(u)) > tol) and (compteur <= 200) :
        # print(real(a), real(b), real(u), f(u))
        # print("    ", abs(f(u)))
        # print(type(a))
        # print(type(b))
        # print(type(u))
        if croissant :
            if fnorm(u) > 0 :
                b = u
            else :
                a = u
        else :
            if fnorm(u) < 0 :
                b = u
            else :
                a = u
        u = (a+b)/2
        compteur += 1

    if tol < abs(fnorm(u)) :
        print("Attention, pas de convergence à ", tol, " : ", abs(f(u)))
    # return np.float64(real(u))
    return(real(u))

def zer0s(f, a, b) :
    N = 10000
    U = np.linspace(a, b, N) # evaluate function at 100 different points
    # U = [x for x in np.linspace(a, b, 100)]
    c = [f(x) for x in U]
    s = np.sign(c)
    zer0 = []
    for i in range(N-1):
        # print(i, U[i])
        if s[i] + s[i+1] == 0: # opposite signs
            u = dichotomie(f, U[i], U[i+1], 1e-16)
            zer0.append(u)
    return zer0


if __name__ == "__main__":
    print("1) Calcul du mode")
    print("2) Étude du déterminant")
    print("3) Étude erreur selon k")
    choix = int(input("Votre choix : "))
    if choix == 1 :
        ww = 10.0
        rhoMoins = 1.0
        rhoPlus = 10.0+10j
        a = 0.2
        cond = 0
        m = -1
        _2D = False

        # rr, uu = ComputeModes1D(ww, rhoMoins, rhoPlus, a, m, 0.4, 0, rayonCLA, 5000, cond)
        _, eigenfunctions = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, 0.4, 0, rayonCLA, 501, cond, _2D)
        # print("Eig : ", eigenfunctions)
        # rr1, uu1 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 0.7, 0, 0.7, 5000, cond)
        # rr2, uu2 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 1.0, 0, 1.0, 5000, cond)
        # rr3, uu3 = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 1.5, 0, 1.5, 5000, cond)
        # rrS, uuSomm = ComputeModes1D(ww, rhoMoins, rhoPlus, a, 2.0, 0, 2.0, 5000, 0)
        matplotlib.rcParams.update({'font.size': 20})
        # plt.plot(rrS, uuSomm, "b-")
        for paire in eigenfunctions :
            if _2D :
                rr = paire[0]
                uu = paire[1]
                plt.plot(rr, real(uu), "r-")
                # plt.plot(rr1, uu1, "g--")
                # plt.plot(rr2, uu2, "b--")
                # plt.plot(rr3, uu3, "r-.")
                plt.xlabel("Rayon")
                plt.ylabel("Mode u")
                plt.show()
            else :
                XX = paire[0]
                YY = paire[1]
                UU = paire[2]
                plt.imshow(real(UU))
                plt.colorbar()
                plt.show()
                plt.imshow(imag(UU))
                plt.colorbar()
                plt.show()
    elif choix == 2 :
        ww = 10
        xMin = 1.01
        xMax = np.sqrt(10)
        # xMin = 1.3
        # xMax = 1.35
        R = 1.0
        m = -1
        rhoMoins = 1.0
        rhoPlus = 10.0+10j
        a = 0.2
        B = np.linspace(xMin, xMax, 100000)
        print("Sommerfeld")
        fSom = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 0)
        Somi = [fSom(x) for x in B]
        zer0 = np.array(zer0s(fSom, xMin, xMax))
        print(zer0)

        # print("BGT1")
        # fBGT1 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 1)
        # BGT1 = [fBGT1(x) for x in B]
        # zer0_BGT1 = np.array(zer0s(fBGT1, xMin, xMax))
        # print(zer0_BGT1)

        # print("BGT2")
        # fBGT2 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 2)
        # BGT2 = [fBGT2(x) for x in B]
        # zer0_BGT2 = np.array(zer0s(fBGT2, xMin, xMax))
        # print(zer0_BGT2)

        matplotlib.rcParams.update({'font.size': 20})
        plt.plot(B*ww, Somi, label="Sommerfeld", linewidth=2)
        # plt.plot(B*ww, BGT1, label="BGT1", linewidth=2)
        # plt.plot(B*ww, BGT2, label="BGT2", linewidth=2)
        # plt.plot(zer0_BGT1*ww,0.0*zer0_BGT1,'x', color="red", label="BGT 1 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(zer0_BGT2*ww,0.0*zer0_BGT2,'x', color="green", label="BGT 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(B*ww,[fK0(x) for x in B], label="K0 (R = {})".format(R), linewidth=2)
        plt.plot(B*ww, 0.0*B*ww, "r--", linewidth=2)

        plt.xlabel("Beta")
        plt.ylabel("det(A(Beta))")
        plt.legend()
        plt.show()
    elif choix == 3 :
        mm = [0, 1, 2, 3]
        rr = [0.4, 0.6, 0.8]
        mark = [".", "x", "s"]
        R = 0.4
        ww = 10.0
        rhoMoins = 1.0
        rhoPlus = 10.0
        a = 0.2

        matplotlib.rcParams.update({'font.size': 20})
        figErr, axErr = plt.subplots(figsize = (9, 6))
        j = 0
        for R in rr :
            kk = []
            ee1 = []
            ee2 = []
            for m in mm :
                betas, eigenfunctionsOpen = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, 501, 0, False)
                _, eigenfunctionsBGT1 = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, 501, 1, False)
                _, eigenfunctionsBGT2 = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, 501, 2, False)

                print(len(eigenfunctionsOpen), " modes trouvés")
                print()

                nbModes = len(eigenfunctionsOpen)

                print("Visualisation")
                for i in range(nbModes) :
                    try :
                        uOpen = eigenfunctionsOpen[i][2]
                        uBGT1 = eigenfunctionsBGT1[i][2]
                        uBGT2 = eigenfunctionsBGT2[i][2]
                        coef1 = 1
                        # if m == 0 :
                        #     coef1 = uOpen[250, 250] / uBGT1[250, 250]
                        # else :
                        #     coef1 = uOpen[250, 300] / uBGT1[250, 300]
                        diff = uOpen - coef1 * uBGT1
                        eAn = np.sum(uOpen*np.conj(uOpen))
                        eAbs = np.sum(diff*np.conj(diff))
                        eRel1 = np.emath.sqrt(eAbs / eAn)
                        print()
                        print("Erreur Approx BGT1, rayon", R, ":", np.format_float_scientific(real(eRel1), precision = 3, exp_digits=1))
                        print()


                        # if m == 0 :
                        #     coef2 = uOpen[250, 250] / uBGT2[250, 250]
                        # else :
                        #     coef2 = uOpen[250, 300] / uBGT2[250, 300]
                        coef2 = 1
                        diff = uOpen - coef2 * uBGT2
                        eAbs = np.sum(diff*np.conj(diff))
                        eRel2 = np.emath.sqrt(eAbs / eAn)
                        print()
                        print("Erreur Approx BGT2, rayon", R, ":", np.format_float_scientific(real(eRel2), precision = 3, exp_digits=1))
                        # print(type(eRel2))
                        print()

                        # fig, axs = plt.subplots(2, 3)
                        #
                        # axs[0, 0].set_title('Open')
                        # axs[0,0].imshow(real(uOpen))
                        # axs[1,0].imshow(imag(uOpen))
                        #
                        # axs[0, 1].set_title('BGT1')
                        # axs[0,1].imshow(real(uBGT1))
                        # axs[1,1].imshow(imag(uBGT1))
                        #
                        # axs[0, 2].set_title('BGT2')
                        # axs[0,2].imshow(real(uBGT2))
                        # axs[1,2].imshow(imag(uBGT2))
                        # plt.colorbar()

                        # plt.pause(0.1)
                        bonMode = int(input("Bon mode ? (oui : 1 / non : 0) : "))
                        # plt.close()
                        if bonMode :
                            kk.append(np.sqrt(betas[i]*betas[i] - rhoMoins))
                            # kk.append(betas[i])
                            ee1.append(real(eRel1))
                            ee2.append(real(eRel2))

                    except IndexError as err :
                        print("Le mode n'existe pas en fait... (singularité Bessel)")

            # K = [ (kk[i],i) for i in range(len(kk)) ]
            # K.sort()
            # sorted_kk,permutation = zip(*K)
            # print(sorted_kk)
            # print(permutation)

            print(kk)
            print(ee1)
            print(ee2)
            axErr.scatter(kk, ee1, marker = mark[j], label = "BGT1, R = {}".format(R), color="red")
            axErr.scatter(kk, ee2, marker = mark[j], label = "BGT2, R = {}".format(R), color="blue")
            # figErr.update()
            j+=1
        # ax.set_xscale("log");
        axErr.set_yscale("log");
        axErr.legend()
        plt.show()
