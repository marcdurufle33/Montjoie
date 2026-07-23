from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
import matrices as mat
import algo
from visuND import *

mp.mp.dps = 200

# pour récupérer récupérer la fréquence et les paramètres de chaque couche
def getVariables(file_name) :
    # Forme fichiers données
    # n (sans compter le conducteur du milieu)
    # f
    # eps mu sigma rayon (sur n+1 lignes)
    variables = []
    eps = []
    sigma = []
    mu = []
    rayon = []
    filin = open(file_name, "r")
    lignes = filin.readlines()
    equation = mat.defineEquation(lignes[0].split()[0])
    ptInit = mp.mpc(mp.mpf(lignes[1].split()[0]), mp.mpf(lignes[1].split()[1]))
    n = int(lignes[2].split()[0])
    m = int(lignes[2].split()[1])
    f = mp.mpf(lignes[3])


    for l in lignes[4:] :
        variables = [mp.mpf(j) for j in l.split()]
        eps.append(variables[0])
        mu.append(variables[1])
        sigma.append(variables[2])
        rayon.append(variables[3])
    eps = mp.matrix(eps)
    mu = mp.matrix(mu)
    sigma = mp.matrix(sigma)
    rayon = mp.matrix(rayon)
    print("Équation : ", equation)
    print("Couches dans câble : ", n)
    print("m : ", m)
    print("Eps : ", eps)
    print("Mu : ", mu)
    print("Sigma : ", sigma)
    print("Rayons : ", rayon)
    return equation, ptInit, n, m, f, eps, sigma, mu, rayon

def ComputeParam(equation, beta, n, f, eps, sigma, mu, rayon, m)     :
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    rho = mp.zeros(n+1, 1)
    k = mp.zeros(n+1, 1)
    c = mp.zeros(n+1, 1)

    rhoInfini = eps[n]*mu[n] + 1j*sigma[n]*mu[n]/omega
    rhoPlus = 0.0
    appRacine = 0.0

    for i in range(n+1):
        rhoPlus = max(eps[i]*mu[i], rhoPlus)
        Delta= omega*omega*eps[i]*mu[i] + 1j*omega*sigma[i]*mu[i] - beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = mp.sqrt(Delta)
        c[i] = (omega*omega*eps[i] + 1j*omega*sigma[i]) / Delta
    # print("k = ", k)


    if equation == mat.Equation.HELMHOLTZ_BGT2N3 or equation == mat.Equation.MAXWELL_BGT2N3 :
        ww = omega
        alpha = mp.sqrt(rhoPlus - rhoInfini) * omega
        appRacine = alpha/8-(rhoInfini*ww*ww-beta*beta*ww*ww)/(8*alpha) - \
                    (rhoInfini*ww*ww-beta*beta*ww*ww)/\
                    (2.0*alpha - 2.0*(rhoInfini*ww*ww - beta*beta*ww*ww)/alpha) \
                     - (rhoInfini*ww*ww-beta*beta*ww*ww)/\
                     (alpha/2.0 - (rhoInfini*ww*ww - beta*beta*ww*ww)/\
                     (2.0*alpha) - 2.0*(rhoInfini*ww*ww - beta*beta*ww*ww)/\
                     (alpha - (rhoInfini*ww*ww - beta*beta*ww*ww)/alpha))

    A = mat.createMatA(equation, n, k, rayon, c, appRacine, m)
    val, EL, vect = mp.eig(A, left = True, right = True)
    val, EL, vect = mp.eig_sort(val, EL, vect, f = lambda x: mp.fabs(x))
    print()
    print(val)
    print()
    print("valeur propre nulle = ", val[0])
    alphas = vect[:,0]
    return k, alphas

# calcul des champs en un point r
def ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu, m) :
    mp.mp.dps = 200
    # Recherche de la couche
    omega = mp.mpf(f)*2*pi;
    n = len(rayon) - 1
    indCouche = len(rayon)
    # print("r = ", r)
    try :
        for i, ray in enumerate(rayon) :
            if i == 0 :
                if r < ray :
                    indCouche = 0
            else :
                if r >= rayon[i-1] and r <= ray :
                    indCouche = i
        eps_tild = eps[indCouche] + 1j * sigma[indCouche]/omega
        Delta = - omega*omega*eps[indCouche]*mu[indCouche] - 1j*omega*sigma[indCouche]*mu[indCouche] + beta*beta*omega*omega
        Ez = 0.0
        dEz = 0.0
        if indCouche == 0 :
            Ez = mp.besselj(m, k[0]*r)
            # dEz = - k[0] * mp.besselj(m+1, k[0]*r) + (m/r) * mp.besselj(m, k[0]*r)
        elif (indCouche == n) and (len(alphas) % 2 == 0) :
            print("Cas Hankel, m = ", m)
            Ez = alphas[2*indCouche-1] / alphas[0] * mp.hankel1(m, k[indCouche]*r)
            # dEz = alphas[2*indCouche-1] / alphas[0] * (- k[indCouche] * mp.hankel1(m+1, k[indCouche]*r) + (m/r) * mp.hankel1(m, k[indCouche]*r))
        else :
            #print("alpha", alphas[2*indCouche-1], alphas[2*indCouche], mp.besselj(0, k[indCouche]*r), mp.bessely(0, k[indCouche]*r), k[indCouche])
            Ez = alphas[2*indCouche-1] / alphas[0] * mp.besselj(m, k[indCouche]*r) + alphas[2*indCouche] / alphas[0] * mp.bessely(m, k[indCouche]*r)
            # dEz = alphas[2*indCouche-1] / alphas[0] * (- k[indCouche] * mp.besselj(m+1, k[indCouche]*r) + (m/r) * mp.besselj(m, k[indCouche]*r))\
                # + alphas[2*indCouche] / alphas[0] * (- k[indCouche] * mp.bessely(m+1, k[indCouche]*r) + (m/r) * mp.bessely(m, k[indCouche]*r))
        # Er = -1j*beta * omega * dEz / Delta
        # Ht = -1j * omega * eps_tild * dEz / Delta
        # return Ez, dEz, Er, Ht
        return Ez, 0, 0, 0
    except :
        # print("PB COMPUTECHAMPS")
        return 0.0, 0.0, 0.0, 0.0

# calcul des champs sur un intervalle de r
def ComputeChampsR(rMin, rMax, pas, f, k, alphas, beta, rayon, eps, sigma, mu, m) :
    taille_1 = int((rMax - rMin) / pas)
    print(rMax - rMin, pas, taille_1)
    rr = mp.linspace(rMin, rMax, taille_1)
    EEz = mp.zeros(taille_1, 1)
    EEr = mp.zeros(taille_1, 1)
    HHt = mp.zeros(taille_1, 1)
    print("Calcul des champs")
    champs = np.array([ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu, m) for r in rr])
    EEz = champs[:,0]
    EEr = champs[:,2]
    HHt = champs[:,3]
    print("Fin du calcul")
    return rr, EEz, EEr, HHt

def drawSolution3D(rMax, nb, f, k, alphas, beta, rayon, eps, sigma, mu, m) :
    rr = np.linspace(-rMax, rMax, nb)
    XX, YY = np.meshgrid(rr, rr)
    RR = np.sqrt(XX**2 + YY**2)
    # TT = np.arctan2(YY, XX)
    # imshow(TT)
    # show()
    UU = np.complex256(0.0) * XX
    for i in range(nb) :
        print(i+1, "/", nb)
        for j in range(nb) :
            x = XX[i, j]
            y = YY[i, j]
            r = mp.sqrt(x**2 + y**2)
            t = mp.mpf(np.arctan2(y, x))
            UU[i, j], _, _, _ = ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu, m)
            UU[i,j] *= np.complex256(mp.exp(mp.mpc(0.0, 1.0) * m * t))
                # UU[i,j] = np.float64(real(UU[i,j])) + 1j * np.float64(imag(UU[i,j]))
                # print(type(UU[i,j]))
    # imshow(RR <= rMax)
    # show()
    return XX, YY, UU

# programme principal (le faire tourner donne le mode propagatif)
def ComputeModes(ficName, rMin, rMax, pas, _2D = True) :
    mp.mp.dps = 200
    equation, ptInit, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)
    Beta = algo.Newton(equation, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
    mp.mp.dps = 15
    print()
    print("Résultat : ", Beta, ", Det = ", mat.ComputeDet(equation, Beta, n, f, eps, sigma, mu, rayon, m))
    mp.mp.dps = 200
    k, alphas = ComputeParam(equation, Beta, n, f, eps, sigma, mu, rayon, m)
    print("Alpha, k = ", alphas, k)
    if _2D :
        rr, EEz, EEr, HHt = ComputeChampsR(rMin, min(rMax,rayon[n]), pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)
        return rr, EEz, EEr, HHt
    else :
        XX, YY, UU = drawSolution3D(min(rMax,np.float64(rayon[n])), int(2 * min(rMax,rayon[n]) / pas), f, k, alphas, Beta, rayon, eps, sigma, mu, m)
        return XX, YY, UU

def displayDet(ficName, rMin, rMax, iMin, iMax, prec) :
    mp.mp.dps = 200
    equation, _, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)
    nbr = int((rMax - rMin) / prec)
    nbi = int((iMax - iMin) / prec)
    rr = np.linspace(rMin, rMax, nbr)
    ii = np.linspace(iMin, iMax, nbi)
    RR, II = np.meshgrid(rr, ii)
    UU = np.complex256(0.0) * RR
    for r in range(nbr) :
        print(r, "/", nbr)
        for i in range(nbi) :
            print("   ", i, "/", nbi)
            x = RR[i, r]
            y = II[i, r]
            beta = mp.mpc(x,y)
            UU[i, r] = mat.ComputeDet(equation, beta, n, f, eps, sigma, mu, rayon, m)
    plot2dinst(RR, II, log(real(UU)))
    plot2dinst(RR, II, log(imag(UU)))
    plot2dinst(RR, II, log(abs(UU)))
    plt.show()

if __name__ == "__main__":
    print("1) Calcul du mode")
    print("2) Affichage du déterminant")
    print("3) Erreurs par rayons")
    print("4) Erreurs par rayons CLA dans le câble")
    choix = int(input("Votre choix : "))
    if choix == 1 :
        ficName = "donneesCable60kHz"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        rr, EEz, EEr, HHt = ComputeModes(ficName, 0, 20.0, 0.001, _2D = True)
        # XX, YY, UU = ComputeModes(ficName, 0, 20.0, 0.01, _2D = False)
        plt.semilogy(rr, abs(EEz), "r-")
        plt.show()
        # plt.imshow(real(UU))
        # plt.colorbar()
        # plt.show()
        # plt.imshow(imag(UU))
        # plt.colorbar()
        # plt.show()
    elif choix == 2 :
        ficName = "exempleHelm"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        displayDet(ficName, 3.0, 3.5, 0.8, 1.3, 0.01)
    elif choix == 3 :
        mp.mp.dps = 200
        rays = np.linspace(0.3, 2.0, 25)
        ficName = "donneesComp"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        equation, ptInit, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)

        eRelBGT2 = []
        eRelBGT2N3 = []
        # eRelDir = []
        Beta_Open = algo.Newton(mat.Equation.HELMHOLTZ_OPEN, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
        k_Open, alphas_Open = ComputeParam(mat.Equation.HELMHOLTZ_OPEN, Beta_Open, n, f, eps, sigma, mu, rayon, m)
        for ray in rays :
            rayon[n] = ray
            N = 1000
            pas = ray / N
            print("Rayon : ", ray)
            print("Open")
            rr, EEz_Open, _, _ = ComputeChampsR(0, ray, pas, f, k_Open, alphas_Open, Beta_Open, rayon, eps, sigma, mu, m)

            # Beta = algo.Newton(mat.Equation.MAXWELL_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon)
            # k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2, Beta, n, f, eps, sigma, mu, rayon)
            # _, EEz_BGT2, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu)

            print("BGT2N3")
            Beta = algo.Newton(mat.Equation.HELMHOLTZ_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.HELMHOLTZ_BGT2, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_BGT2N3, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)

            # print("Dirichlet")
            # Beta = algo.Newton(mat.Equation.MAXWELL_DIRICHLET, ptInit, 50, n, f, eps, sigma, mu, rayon)
            # k, alphas = ComputeParam(mat.Equation.MAXWELL_DIRICHLET, Beta, n, f, eps, sigma, mu, rayon)
            # _, EEz_Dir, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu)

            eAn = np.dot((EEz_Open),np.conj(EEz_Open))

            # diff = EEz_BGT2 - EEz_Open
            # eAbs = np.dot((diff),np.conj(diff))
            # eRelBGT2.append(real(mp.sqrt(eAbs / eAn)))

            diff = EEz_BGT2N3 - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRelBGT2N3.append(real(mp.sqrt(eAbs / eAn)))

            # diff = EEz_Dir - EEz_Open
            # eAbs = np.dot((diff),np.conj(diff))
            # eRelDir.append(real(mp.sqrt(eAbs / eAn)))

        matplotlib.rcParams.update({'font.size': 20})
        # plt.semilogy(rays, eRelTrueABC, 'o-', label='Fixed Sommerfeld')
        # plt.semilogy(rays, eRelBGT2, 'o-', label='ABC BGT2')
        plt.semilogy(rays, eRelBGT2N3, 'o-', label='ABC BGT2 Newton 3')
        # plt.semilogy(rays, eRelDir, 'o--', label='Dirichlet')
        plt.xlabel("Boundary radius")
        plt.ylabel("Relative error")
        plt.legend()
        plt.show()

    elif choix == 4 :
        mp.mp.dps = 200
        rays = np.linspace(1.0, 5.0, 15)
        ficName = "donneesCable60kHz"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        equation, ptInit, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)

        eRelBGT2 = []
        eRelBGT2N3 = []
        eRelDir = []

        N = 100
        rayCable = rayon[n-1]
        rayDebut = rayCable + 0
        print(rayCable)
        pas = rayCable / N

        print("Open")
        Beta_Open = algo.Newton(mat.Equation.MAXWELL_OPEN, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
        k_Open, alphas_Open = ComputeParam(mat.Equation.MAXWELL_OPEN, Beta_Open, n, f, eps, sigma, mu, rayon, m)

        debut = 0
        for ray in rays :
            rayon[n] = ray
            print("Rayon : ", ray)
            if ray < rayDebut :
                debut += 1
                continue
            #
            # print("BGT2")
            # Beta = algo.Newton(mat.Equation.MAXWELL_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon)
            # k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2, Beta, n, f, eps, sigma, mu, rayon)
            # _, EEz_BGT2, _, _ = ComputeChampsR(0, rayCable, pas, f, k, alphas, Beta, rayon, eps, sigma, mu)

            pas = (ray - rayCable) / N

            print("Open")
            rr, EEz_Open, _, _ = ComputeChampsR(rayCable, ray, pas, f, k_Open, alphas_Open, Beta_Open, rayon, eps, sigma, mu, m)
            eAn = np.dot((EEz_Open),np.conj(EEz_Open))

            print("BGT2N3")
            Beta = algo.Newton(mat.Equation.MAXWELL_BGT2N3, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2N3, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_BGT2N3, _, _ = ComputeChampsR(rayCable, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)

            print("Dirichlet")
            Beta = algo.Newton(mat.Equation.MAXWELL_DIRICHLET, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.MAXWELL_DIRICHLET, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_Dir, _, _ = ComputeChampsR(rayCable, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)

            # diff = EEz_BGT2 - EEz_Open
            # eAbs = np.dot((diff),np.conj(diff))
            # eRelBGT2.append(real(mp.sqrt(eAbs / eAn)))

            diff = EEz_BGT2N3 - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRelBGT2N3.append(real(mp.sqrt(eAbs / eAn)))

            diff = EEz_Dir - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRelDir.append(real(mp.sqrt(eAbs / eAn)))

        matplotlib.rcParams.update({'font.size': 20})
        # plt.semilogy(rays, eRelTrueABC, 'o-', label='Fixed Sommerfeld')
        # plt.semilogy(rays, eRelBGT2, 'o-', label='ABC BGT2')
        plt.semilogy(rays[debut:], eRelDir, 'o-', label='Dirichlet')
        plt.semilogy(rays[debut:], eRelBGT2N3, 'o--', label='ABC BGT2 Newton 3')
        plt.xlabel("Boundary radius")
        plt.ylabel("Relative error")
        plt.legend()
        plt.show()

    elif choix == 5 :
        mp.mp.dps = 200
        nbCLA = 10
        nbDir = 10
        rayCLA = np.linspace(0.5, 5.0, nbCLA)
        rayDir = np.linspace(25.0, 100.0, nbDir)
        ficName = "donneesCable60kHz"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        equation, ptInit, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)

        eRelBGT2 = []
        eRelDir = []
        Beta_Open = algo.Newton(mat.Equation.MAXWELL_OPEN, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
        k_Open, alphas_Open = ComputeParam(mat.Equation.MAXWELL_OPEN, Beta_Open, n, f, eps, sigma, mu, rayon, m)

        Beta_Dir = []
        k_Dir = []
        alphas_Dir = []
        for ray in rayDir :
            rayon[n] = ray
            bD = algo.Newton(mat.Equation.MAXWELL_DIRICHLET, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            kD, aD = ComputeParam(mat.Equation.MAXWELL_DIRICHLET, bD, n, f, eps, sigma, mu, rayon, m)

            Beta_Dir.append(bD)
            k_Dir.append(kD)
            alphas_Dir.append(aD)

        for ray in rayCLA :
            rayon[n] = ray
            N = 1000
            pas = ray / N
            print("Rayon : ", ray)
            print("Open")
            rr, EEz_Open, _, _ = ComputeChampsR(0, ray, pas, f, k_Open, alphas_Open, Beta_Open, rayon, eps, sigma, mu, m)
            eAn = np.dot((EEz_Open),np.conj(EEz_Open))

            # Beta = algo.Newton(mat.Equation.MAXWELL_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon)
            # k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2, Beta, n, f, eps, sigma, mu, rayon)
            # _, EEz_BGT2, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu)

            for i in range(len(Beta_Dir)) :
                print("Rayon : ", ray)
                print("Dirichlet", rayDir[i])
                _, EEz_Dir, _, _ = ComputeChampsR(0, ray, pas, f, k_Dir[i], alphas_Dir[i], Beta_Dir[i], rayon, eps, sigma, mu, m)

                diff = EEz_Dir - EEz_Open
                eAbs = np.dot((diff),np.conj(diff))
                eRelDir.append(real(mp.sqrt(eAbs / eAn)))

            print("Rayon : ", ray)
            print("BGT2")
            Beta = algo.Newton(mat.Equation.MAXWELL_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_BGT2, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)
            # plt.semilogy(rr, abs(EEz_Open))
            # plt.semilogy(rr, abs(EEz_BGT2), "r-")
            # plt.show()
            diff = EEz_BGT2 - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRelBGT2.append(real(mp.sqrt(eAbs / eAn)))

        print()
        print("Rayons CLA")
        print(rayCLA)
        print("Erreurs BGT2")
        print(eRelBGT2)
        print()
        print("Rayon Dirichlet")
        print(rayDir)
        for i in range(nbCLA) : # 10 dans rayCLA
            print()
            print("Pour un domaine tronqué au rayon des CLA ", rayCLA[i]," pour le pb de Dirichlet ")
            print(eRelDir[i*nbCLA:(i+1)*nbCLA])

    elif choix == 6 :
        mp.mp.dps = 200
        nbCLA = 10
        nbDir = 10
        rayCLA = np.linspace(0.5, 5.0, nbCLA)
        tol = 1e-5
        ficName = "donneesCable60kHz"
        if len(sys.argv) == 2 :
            ficName = sys.argv[1]
        equation, ptInit, n, m, f, eps, sigma, mu, rayon = getVariables(ficName)

        rayDir = []
        Beta_Open = algo.Newton(mat.Equation.MAXWELL_OPEN, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
        k_Open, alphas_Open = ComputeParam(mat.Equation.MAXWELL_OPEN, Beta_Open, n, f, eps, sigma, mu, rayon, m)


        for ray in rayCLA :
            rayon[n] = ray
            N = 1000
            pas = ray / N
            print("Rayon : ", ray)
            print("Open")
            rr, EEz_Open, _, _ = ComputeChampsR(0, ray, pas, f, k_Open, alphas_Open, Beta_Open, rayon, eps, sigma, mu, m)
            eAn = np.dot((EEz_Open),np.conj(EEz_Open))

            print("Rayon : ", ray)
            print("BGT2")
            Beta = algo.Newton(mat.Equation.MAXWELL_BGT2, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.MAXWELL_BGT2, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_BGT2, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)
            diff = EEz_BGT2 - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRel_BGT2 = real(mp.sqrt(eAbs / eAn))

            ray1 = ray
            ray2 = 100

            print("Dirichlet it 0")
            rayD = 50 + ray/2
            rayon[n] = rayD
            Beta = algo.Newton(mat.Equation.MAXWELL_DIRICHLET, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
            k, alphas = ComputeParam(mat.Equation.MAXWELL_DIRICHLET, Beta, n, f, eps, sigma, mu, rayon, m)
            _, EEz_Dir, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)
            diff = EEz_Dir - EEz_Open
            eAbs = np.dot((diff),np.conj(diff))
            eRelDir = real(mp.sqrt(eAbs / eAn))

            while eRelDir - eRel_BGT2 > tol :
                if eRelDir > eRel_BGT2 :
                    ray1 = rayD
                else :
                    ray2 = rayD
                rayD = (ray1 + ray2) /2
                print("Dirichlet ray ", rayD)
                rayon[n] = rayD
                Beta = algo.Newton(mat.Equation.MAXWELL_DIRICHLET, ptInit, 50, n, f, eps, sigma, mu, rayon, m)
                k, alphas = ComputeParam(mat.Equation.MAXWELL_DIRICHLET, Beta, n, f, eps, sigma, mu, rayon, m)
                _, EEz_Dir, _, _ = ComputeChampsR(0, ray, pas, f, k, alphas, Beta, rayon, eps, sigma, mu, m)
                diff = EEz_Dir - EEz_Open
                eAbs = np.dot((diff),np.conj(diff))
                eRelDir = real(mp.sqrt(eAbs / eAn))
            rayDir.append(rayD)
        matplotlib.rcParams.update({'font.size': 20})
        # plt.semilogy(rays, eRelTrueABC, 'o-', label='Fixed Sommerfeld')
        # plt.semilogy(rays, eRelBGT2, 'o-', label='ABC BGT2')
        plt.plot(rayCLA, rayDir, 'o-')
        plt.xlabel("Rayon de CLA")
        plt.ylabel("Rayon de Dirichlet")
        plt.show()
