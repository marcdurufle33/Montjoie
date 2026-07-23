from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

mp.mp.dps = 200

rayonCLA = 0.4

def createMatASommerfeld(k1, k2, a, m) :
    A = mp.zeros(2)
    A[0, 0] = mp.besselj(m, k1*a);
    A[0, 1] = -mp.hankel1(m, k2*a);
    A[1, 0] = -k1*mp.besselj(m+1, k1*a) + (m/a) * mp.besselj(m, k1*a);
    A[1, 1] = k2*mp.hankel1(m+1, k2*a) - (m/a) * mp.hankel1(m, k2*a);
    return A

def createMatA_BGT1(k1, k2, a, m, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(m, k1*a);
    A[0, 1] = -mp.besselj(m, k2*a);
    A[0, 2] = -mp.bessely(m, k2*a);
    A[1, 0] = - k1*mp.besselj(m+1, k1*a) + (m/a) * mp.besselj(m, k1*a);
    A[1, 1] = k2*mp.besselj(m+1, k2*a) - (m/a) * mp.besselj(m, k2*a);
    A[1, 2] = k2*mp.bessely(m+1, k2*a) - (m/a) * mp.bessely(m, k2*a);

    X = (1/(2*R) - 1j*k2)

    A[2, 0] = 0.0;
    A[2, 1] = -k2*mp.besselj(m+1, k2*R) + (m/R) * mp.besselj(m, k2*R) + X * mp.besselj(m, k2*R);
    A[2, 2] = -k2*mp.bessely(m+1, k2*R) + (m/R) * mp.bessely(m, k2*R) + X * mp.bessely(m, k2*R);
    return A

def createMatA_BGT2(k1, k2, a, m, R) :
    A = mp.zeros(3)
    A[0, 0] = mp.besselj(m, k1*a);
    A[0, 1] = -mp.besselj(m, k2*a);
    A[0, 2] = -mp.bessely(m, k2*a);
    A[1, 0] = - k1*mp.besselj(m+1, k1*a) + (m/a) * mp.besselj(m, k1*a);
    A[1, 1] = k2*mp.besselj(m+1, k2*a) - (m/a) * mp.besselj(m, k2*a);
    A[1, 2] = k2*mp.bessely(m+1, k2*a) - (m/a) * mp.bessely(m, k2*a);

    num = 3/4 - 3*1j*R*k2 - 2*R*R*k2*k2
    denom = 2 * R * (1 - 1j * k2 * R)

    X = num / denom
    Y = 1.0 / denom

    A[2, 0] = 0.0;
    A[2, 1] = -k2*mp.besselj(m+1, k2*R) + (m/R) * mp.besselj(m, k2*R) + X * mp.besselj(m, k2*R) + m*m*Y*mp.besselj(m, k2*R);
    A[2, 2] = -k2*mp.bessely(m+1, k2*R) + (m/R) * mp.bessely(m, k2*R) + X * mp.bessely(m, k2*R) + m*m*Y*mp.bessely(m, k2*R);
    return A

def det2D(A) :
    return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]


def det3D(A) :
    return A[0, 0]*A[1, 1]*A[2, 2] + A[0, 1]*A[1, 2]*A[2, 0] + A[0, 2]*A[1, 0]*A[2, 1] - A[0, 2]*A[1, 1]*A[2,0] - A[0, 1]*A[1, 0]*A[2,2] - A[0, 0]*A[1, 2]*A[2,1]

# toutes les valeurs doivent être en mp précision
def ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a, m):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta*beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta*beta*ww*ww)
    A = createMatASommerfeld(k1, k2, a, m)
    return det2D(A)

def ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, m, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT1(k1, k2, a, m, R)
    return det3D(A)

def ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, m, R):
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    A = createMatA_BGT2(k1, k2, a, m, R)
    return det3D(A)


def createMatA(k1, k2, a, m, ray, cond, kBord = 0) :
    if cond == 0 :
        return createMatASommerfeld(k1, k2, a, m)
    elif cond == 1 :
        return createMatA_BGT1(k1, k2, a, m, ray)
    elif cond == 2 :
        return createMatA_BGT2(k1, k2, a, m, ray)

# def ComputeDet(beta, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
#     if cond == 0 :
#         return ComputeDetSommerfeld(beta, ww, rhoMoins, rhoPlus, a, m)
#     elif cond == 1 :
#         return ComputeDet_BGT1(beta, ww, rhoMoins, rhoPlus, a, m, ray)
#     elif cond == 2 :
#         return ComputeDet_BGT2(beta, ww, rhoMoins, rhoPlus, a, m, ray)

def ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    if cond == 0 :
        return lambda x : np.float128(real(ComputeDetSommerfeld(x, ww, rhoMoins, rhoPlus, a, m) * ((-1j)**(m+1))))
    elif cond == 1 :
        return lambda x : np.float128(real(ComputeDet_BGT1(x, ww, rhoMoins, rhoPlus, a, m, ray)))
    elif cond == 2 :
        return lambda x : np.float128(real(ComputeDet_BGT2(x, ww, rhoMoins, rhoPlus, a, m, ray)))

def ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    if cond == 0 :
        return lambda x : real(ComputeDetSommerfeld(x, ww, rhoMoins, rhoPlus, a, m) * ((-1j)**(m+1)))
    elif cond == 1 :
        return lambda x : real(ComputeDet_BGT1(x, ww, rhoMoins, rhoPlus, a, m, ray))
    elif cond == 2 :
        return lambda x : real(ComputeDet_BGT2(x, ww, rhoMoins, rhoPlus, a, m, ray))


def ComputeParam(beta, ww, rhoMoins, rhoPlus, a, m, ray, cond) :
    mp.mp.dps = 200
    k1 = mp.sqrt(rhoPlus*ww*ww - beta * beta*ww*ww)
    k2 = mp.sqrt(rhoMoins*ww*ww - beta * beta*ww*ww)
    kBord = 0
    # if cond == 2 :
    #     kBord = mp.sqrt(rhoMoins*ww*ww - rhoPlus*ww*ww)
    A = createMatA(k1, k2, a, m, ray, cond, kBord)
    print(A)
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

def drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m) :
    rr = np.linspace(rMin, rMax, nb)
    if len(alphas) == 2 :
        uu = np.array([real(mp.besselj(m, k1*r) if r <= a else alphas[1]/alphas[0]*mp.hankel1(m, k2*r)) for r in rr])
    else :
        uu = np.array([real(mp.besselj(m, k1*r) if r <= a else alphas[1]/alphas[0]*mp.besselj(m, k2*r) + alphas[2]/alphas[0]*mp.bessely(m, k2*r)) for r in rr])
    return rr, uu

# def drawSolution3D(rMax, nb, alphas, k1, k2, a, m) :
#     rr = np.linspace(-rMax, rMax, nb)
#     XX, YY = np.meshgrid(rr, rr)
#     RR = np.sqrt(XX**2 + YY**2)
#     # TT = np.arctan2(YY, XX)
#     # imshow(TT)
#     # show()
#     UU = mp.mpf(0.0) * XX
#     for i in range(nb) :
#         print(i)
#         for j in range(nb) :
#             x = XX[i, j]
#             y = YY[i, j]
#             r = mp.sqrt(x**2 + y**2)
#             t = mp.mpf(np.arctan2(y, x))
#             if r <= a :
#                 UU[i, j] = mp.besselj(m, k1*r) * mp.exp(1j * m * t)
#             elif (r > a) and (r <= rMax) :
#                 if len(alphas) == 2 :
#                     UU[i, j] = alphas[1]/alphas[0]*mp.hankel1(m,k2*r)
#                 elif len(alphas) == 3 :
#                     UU[i, j] = alphas[1]/alphas[0]*mp.besselj(m,k2*r) + alphas[2]/alphas[0]*mp.bessely(m,k2*r)
#                 UU[i,j] *= mp.exp(mp.mpc(0.0, 1.0) * m * t)
#     # imshow(RR <= rMax)
#     # show()
#     return XX, YY, UU

def drawSolution3D(rMax, nb, alphas, k1, k2, a, m) :
    rr = np.linspace(-rMax, rMax, nb)
    XX, YY = np.meshgrid(rr, rr)
    RR = np.sqrt(XX**2 + YY**2)
    # TT = np.arctan2(YY, XX)
    # imshow(TT)
    # show()
    UU = np.complex256(0.0) * XX
    for i in range(nb) :
        print(i)
        for j in range(nb) :
            x = XX[i, j]
            y = YY[i, j]
            r = mp.sqrt(x**2 + y**2)
            t = mp.mpf(np.arctan2(y, x))
            if r <= a :
                UU[i, j] = np.complex256(mp.besselj(m, k1*r) * mp.exp(1j * m * t))
            elif (r > a) and (r <= rMax) :
                if len(alphas) == 2 :
                    UU[i, j] = np.complex256(alphas[1]/alphas[0]*mp.hankel1(m,k2*r))
                elif len(alphas) == 3 :
                    UU[i, j] = np.complex256(alphas[1]/alphas[0]*mp.besselj(m,k2*r) + alphas[2]/alphas[0]*mp.bessely(m,k2*r))
                UU[i,j] *= np.complex256(mp.exp(mp.mpc(0.0, 1.0) * m * t))
                # UU[i,j] = np.float64(real(UU[i,j])) + 1j * np.float64(imag(UU[i,j]))
                # print(type(UU[i,j]))
    # imshow(RR <= rMax)
    # show()
    return XX, YY, UU

def ComputeModes1D(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond) :
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond)
    mp.mp.dps = 200
    xMin = np.float64(1.0)
    xMax = np.float64(np.sqrt(10.0))
    Betas = zer0s(fCLA, xMin, xMax)
    mp.mp.dps = 15
    print(Betas)


    if len(Betas) > 0 :
        # print("Mode le plus propagatif")
        Beta = Betas[0]
        print("Résultat : ", Beta)
        mp.mp.dps = 200
        k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond)
        print("Alpha, k = ", alphas, k1, k2)
        return drawSolutionTheta0(rMin, rMax, nb, alphas, k1, k2, a, m)
    else :
        print("Aucun mode trouvé...")
        return np.linspace(rMin, rMax, nb), 0.0 * np.linspace(rMin, rMax, nb)

def ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond, _2D = True) :
    # print("Compute")
    mp.mp.dps = 200
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond)
    # print("Init")
    xMin = np.float64(1.01)
    xMax = np.float64(np.sqrt(10.0))
    # print("Calcul betas")
    Betas = zer0s(fCLA, xMin, xMax)
    mp.mp.dps = 15
    print(Betas)

    print("Avant boucle")

    if len(Betas) > 0 :
        eigenfunctions = []
        for Beta in Betas :
            print()
            mp.mp.dps = 15
            print("Beta = ", Beta)
            mp.mp.dps = 200
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

def ComputeParam1DAll(ww, rhoMoins, rhoPlus, a, m, ray, rMin, rMax, nb, cond) :
    # print("Compute")
    mp.mp.dps = 200
    fCLA = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, ray, cond)
    # print("Init")
    xMin = np.float64(1.01)
    xMax = np.float64(np.sqrt(10.0))
    # print("Calcul betas")
    Betas = zer0s(fCLA, xMin, xMax)
    mp.mp.dps = 15
    print(Betas)

    print("Avant boucle")

    if len(Betas) > 0 :
        AlphasTab = []
        for Beta in Betas :
            print()
            mp.mp.dps = 15
            print("Beta = ", Beta)
            mp.mp.dps = 200
            k1, k2, alphas = ComputeParam(Beta, ww, rhoMoins, rhoPlus, a, m, ray, cond)
            print("Alpha, k = ", alphas, k1, k2)
            print()
            AlphasTab.append(alphas)
        return Betas, AlphasTab
    else :
        print("Aucun mode trouvé...")
        return np.linspace(rMin, rMax, nb), 0.0 * np.linspace(rMin, rMax, nb)


def dichotomie(f, a, b, tol) :
    # a = np.longdouble(a)
    # b = np.longdouble(b)
    croissant = (f(b) > f(a))
    u = (a+b)/2
    compteur = 0
    while (abs(f(u)) > tol) and (compteur <= 200) :
        # print(real(a), real(b), real(u), f(u))
        # print("    ", abs(f(u)))
        # print("    Type a : ", type(a))
        # print("    Type b : ", type(b))
        # print("    Type u : ", type(u))
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
        compteur += 1

    if tol < abs(f(u)) :
        print("Attention, pas de convergence à ", tol, " : ", abs(f(u)))
    # return np.float64(real(u))
    return(real(u))

def zer0s(f, a, b) :
    # print(type(a))
    N = 100
    # U = np.linspace(np.float128(a), np.float128(b), N) # evaluate function at 100 different points
    U = [mp.mpc(x, 0.0) for x in np.linspace(a, b, 100)]
    c = [f(mp.mpc(x, 0.0)) for x in U]
    s = np.sign(c)
    zer0 = []
    for i in range(N-1):
        # print(i, U[i])
        if s[i] + s[i+1] == 0: # opposite signs
            u = dichotomie(f, U[i], U[i+1], 1e-64)
            zer0.append(u)
    return zer0


if __name__ == "__main__":
    print("1) Calcul du mode")
    print("2) Étude du déterminant")
    print("3) Étude erreur selon k")
    print("4) Étude erreur selon k (distance param)")
    choix = int(input("Votre choix : "))
    if choix == 1 :
        ww = mp.mpf(10.0)
        rhoMoins = mp.mpf(1.0)
        rhoPlus = mp.mpf(10.0)
        a = mp.mpf(0.2)
        cond = 0
        m = 0
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
                print(UU)
                plt.imshow(real(UU))
                plt.colorbar()
                plt.show()
                plt.imshow(imag(UU))
                plt.colorbar()
                plt.show()
    elif choix == 2 :
        ww = 10
        xMin = 0.01
        xMax = np.sqrt(10)
        # xMin = 1.3
        # xMax = 1.35
        R = 1.0
        m = 0
        rhoMoins = 1.0
        rhoPlus = 10.0
        a = 0.2
        B = np.linspace(xMin, xMax, 200)
        print("Sommerfeld")
        fSom = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 0)
        Somi = [fSom(x) for x in B]

        # print("BGT1")
        # fBGT1 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 1)
        # BGT1 = [fBGT1(x) for x in B]
        # zer0_BGT1 = np.array(zer0s(fBGT1, xMin, xMax))
        # print(zer0_BGT1)

        print("BGT2")
        fBGT2 = ComputeFLambda(ww, rhoMoins, rhoPlus, a, m, R, 2)
        # BGT2 = [fBGT2(x) for x in B]
        zer0_BGT2 = np.array(zer0s(fBGT2, xMin, xMax))
        print(zer0_BGT2)

        matplotlib.rcParams.update({'font.size': 20})
        plt.plot(B*ww, Somi, label="Sommerfeld", linewidth=2)
        # plt.plot(B*ww, BGT1, label="BGT1", linewidth=2)
        # plt.plot(B*ww, BGT2, label="BGT2", linewidth=2)
        # plt.plot(zer0_BGT1*ww,0.0*zer0_BGT1,'x', color="red", label="BGT 1 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        plt.plot(zer0_BGT2*ww,0.0*zer0_BGT2,'x', color="green", label="BGT 2 (R = {})".format(R), markersize=16, markeredgewidth = 3)
        # plt.plot(B*ww,[fK0(x) for x in B], label="K0 (R = {})".format(R), linewidth=2)
        plt.plot(B*ww, 0.0*B*ww, "r--", linewidth=2)

        plt.xlabel("Beta")
        plt.ylabel("det(A(Beta))")
        plt.legend()
        plt.show()
    elif choix == 3 :
        mm = [0, 1, 2, 3]
        # mm = [0]
        rr = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
        mark = [".", "x", "s"]
        R = 0.4
        ww = 10.0
        rhoMoins = 1.0
        rhoPlus = 10.0
        a = 0.2
        N = 5

        matplotlib.rcParams.update({'font.size': 20})
        figErr, axErr = plt.subplots(figsize = (9, 6))
        j = 0
        for R in rr :
            print("R = ", R)
            kk = []
            ee1 = []
            ee2 = []
            for m in mm :
                print("m = ", m)
                print("    Open")
                betas, eigenfunctionsOpen = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 0, False)
                # print("    BGT1")
                # _, eigenfunctionsBGT1 = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 1, False)
                print("    BGT2")
                _, eigenfunctionsBGT2 = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 2, False)

                print(len(eigenfunctionsOpen), " modes trouvés")
                print()

                nbModes = len(eigenfunctionsOpen)

                print("Visualisation")
                for i in range(nbModes) :
                    try :
                        uOpen = eigenfunctionsOpen[i][2]
                        # uBGT1 = eigenfunctionsBGT1[i][2]
                        uBGT2 = eigenfunctionsBGT2[i][2]
                        # coef1 = mp.sqrt(np.sum(uOpen*np.conj(uOpen)) / np.sum(uBGT1*np.conj(uBGT1)))
                        # # if m == 0 :
                        # #     coef1 = uOpen[250, 250] / uBGT1[250, 250]
                        # # else :
                        # #     coef1 = uOpen[250, 300] / uBGT1[250, 300]
                        # diff = uOpen - coef1 * uBGT1
                        eAn = np.sum(uOpen*np.conj(uOpen))
                        # eAbs = np.sum(diff*np.conj(diff))
                        # eRel1 = mp.sqrt(eAbs / eAn)
                        # print()
                        # print("VP = ", betas[i])
                        # print("Erreur Approx BGT1, rayon", R, ":", np.format_float_scientific(real(eRel1), precision = 3, exp_digits=1))
                        # print()


                        # if m == 0 :
                        #     coef2 = uOpen[250, 250] / uBGT2[250, 250]
                        # else :
                        #     coef2 = uOpen[250, 300] / uBGT2[250, 300]
                        coef2 = mp.sqrt(np.sum(uOpen*np.conj(uOpen)) / np.sum(uBGT2*np.conj(uBGT2)))
                        diff = uOpen - coef2 * uBGT2
                        eAbs = np.sum(diff*np.conj(diff))
                        eRel2 = mp.sqrt(eAbs / eAn)
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
                            # ee1.append(real(eRel1))
                            ee2.append(real(eRel2))

                    except IndexError as err :
                        print("Le mode n'existe pas en fait... (singularité Bessel)")

            K = [ (kk[i],i) for i in range(len(kk)) ]
            K.sort()
            sorted_kk,permutation = zip(*K)
            sorted_ee2 = [ee2[i] for i in permutation]
            print(sorted_kk)
            print(sorted_ee2)
            print(permutation)

            print(kk)
            # print(ee1)
            print(ee2)
            # axErr.scatter(kk, ee1, marker = mark[j], label = "BGT1, R = {}".format(R), color="red")
            axErr.plot(sorted_kk, sorted_ee2, label = "BGT2, R = {}".format(R))
            # figErr.update()
            j+=1
        # axErr.set_xscale("log");
        axErr.set_yscale("log");
        axErr.legend()
        plt.show()
    elif choix == 4 :
        mm = [0, 1, 2, 3]
        # mm = [0]
        # rr = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
        rr = [0.4]
        mark = [".", "x", "s"]
        ww = 10.0
        rhoMoins = 1.0
        rhoPlus = 10.0
        a = 0.2
        N = 5

        matplotlib.rcParams.update({'font.size': 20})
        figErr, axErr = plt.subplots(figsize = (9, 6))
        j = 0
        for R in rr :
            print("R = ", R)
            kk = []
            nn = []
            ee1 = []
            ee2 = []
            eeA2 = []
            eeB2 = []
            for m in mm :
                print("m = ", m)
                print("    Open")
                betasOpen, alphasOpen = ComputeParam1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 0)
                # print("    BGT1")
                # _, eigenfunctionsBGT1 = ComputeModes1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 1, False)
                print("    BGT2")
                betasBGT2, alphasBGT2 = ComputeParam1DAll(ww, rhoMoins, rhoPlus, a, m, R, 0, R, N, 2)

                print(len(alphasOpen), " modes trouvés")
                print()

                nbModes = len(alphasOpen)

                print("Visualisation")
                for i in range(nbModes) :
                    try :
                        aOpen = alphasOpen[i]
                        bOpen = betasOpen[i]
                        print("VP Ouvert : ", bOpen)
                        # uBGT1 = eigenfunctionsBGT1[i][2]
                        aBGT2 = alphasBGT2[i]
                        bBGT2 = betasBGT2[i]
                        print("VP BGT2 : ", bBGT2)

                        distAlpha2 = (aBGT2[1]/aBGT2[0] - aOpen[1]/aOpen[0]) * mp.conj(aBGT2[1]/aBGT2[0] - aOpen[1]/aOpen[0]) \
                                + (aBGT2[2]/aBGT2[0] - 1j * aOpen[1]/aOpen[0]) * mp.conj(aBGT2[2]/aBGT2[0] - 1j * aOpen[1]/aOpen[0])
                        distBeta2 = (bBGT2 - bOpen) * mp.conj(bBGT2 - bOpen)
                        dist2 = distAlpha2 + distBeta2

                        normAlpha2Open = (aOpen[1]/aOpen[0]) * mp.conj(aOpen[1]/aOpen[0]) \
                                    + (1j * aOpen[1]/aOpen[0]) * mp.conj(1j * aOpen[1]/aOpen[0])
                        normBeta2Open = (bOpen) * mp.conj(bOpen)
                        norm2Open = normAlpha2Open + normBeta2Open

                        distAlpha = real(mp.sqrt(distAlpha2 / normAlpha2Open))
                        distBeta = real(mp.sqrt(distBeta2 / normBeta2Open))
                        dist = real(mp.sqrt(dist2 / norm2Open))
                        print()
                        print("Erreur Approx BGT2, rayon", R, ":", np.format_float_scientific(dist, precision = 3, exp_digits=1))
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
                            sBeta = np.sqrt(betasOpen[i]*betasOpen[i] - rhoMoins)
                            normdOmega = mp.sqrt(2.0 * mp.pi * (1+m*m)) * mp.fabs(aOpen[1] * mp.hankel1(m, 1j * R * sBeta))
                            kk.append(sBeta)
                            nn.append(normdOmega)
                            # kk.append(betas[i])
                            # ee1.append(real(eRel1))
                            ee2.append(dist)
                            eeA2.append(distAlpha)
                            eeB2.append(distBeta)

                    except IndexError as err :
                        print("Le mode n'existe pas en fait... (singularité Bessel)")

            K = [ (kk[i],i) for i in range(len(kk)) ]
            K.sort()
            sorted_kk,permutation = zip(*K)
            sorted_ee2 = [ee2[i] for i in permutation]
            sorted_eeA2 = [eeA2[i] for i in permutation]
            sorted_eeB2 = [eeB2[i] for i in permutation]
            print(sorted_kk)
            print(sorted_ee2)
            print(permutation)

            print(kk)
            # print(ee1)
            print(ee2)
            # axErr.scatter(kk, ee1, marker = mark[j], label = "BGT1, R = {}".format(R), color="red")
            # axErr.plot(nn, ee2, label = "BGT2, R = {}".format(R), color="red")
            axErr.plot(sorted_kk, sorted_ee2, label = "BGT2, R = {}".format(R))
            # figErr.update()
            j+=1
        # axErr.set_xscale("log");
        axErr.set_yscale("log");
        axErr.legend()
        plt.show()
