from pylab import *
from scipy.special import *
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

mp.mp.dps = 200

pasNewton = 1e-5

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
    n = int(lignes[0])
    f = mp.mpf(lignes[1])


    for l in lignes[2:] :
        variables = [mp.mpf(j) for j in l.split()]
        eps.append(variables[0])
        mu.append(variables[1])
        sigma.append(variables[2])
        rayon.append(variables[3])
    eps = mp.matrix(eps)
    mu = mp.matrix(mu)
    sigma = mp.matrix(sigma)
    rayon = mp.matrix(rayon)
    print("Eps : ", eps)
    print("Mu : ", mu)
    print("Sigma : ", sigma)
    print("Rayons : ", rayon)
    return n, f, eps, sigma, mu, rayon

def createMatA(k, rayon, c) :
    A = mp.zeros(2*3+1)
    # condition conducteur <-> isolant
    A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
    A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
    A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
    A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
    A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
    A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);

    # condition isolant <-> blindage
    A[2, 1] = mp.besselj(0, k[1]*rayon[1]);
    A[2, 2] = mp.bessely(0, k[1]*rayon[1]);
    A[2, 3] = -mp.besselj(0, k[2]*rayon[1]);
    A[2, 4] = -mp.bessely(0, k[2]*rayon[1]);
    # Continuité de la dérivée
    A[3, 1] = c[1]*k[1]*mp.besselj(1, k[1]*rayon[1]);
    A[3, 2] = c[1]*k[1]*mp.bessely(1, k[1]*rayon[1]);
    A[3, 3] = -c[2]*k[2]*mp.besselj(1, k[2]*rayon[1]);
    A[3, 4] = -c[2]*k[2]*mp.bessely(1, k[2]*rayon[1]);

    # condition blindage <-> eau
    A[4, 3] = mp.besselj(0, k[2]*rayon[2]);
    A[4, 4] = mp.bessely(0, k[2]*rayon[2]);
    A[4, 5] = -mp.besselj(0, k[3]*rayon[2]);
    A[4, 6] = -mp.bessely(0, k[3]*rayon[2]);
    # Continuité de la dérivée
    A[5, 3] = c[2]*k[2]*mp.besselj(1, k[2]*rayon[2]);
    A[5, 4] = c[2]*k[2]*mp.bessely(1, k[2]*rayon[2]);
    A[5, 5] = -c[3]*k[3]*mp.besselj(1, k[3]*rayon[2]);
    A[5, 6] = -c[3]*k[3]*mp.bessely(1, k[3]*rayon[2]);

    # Condition de Dirichlet homogène assez loin
    A[6, 5] = mp.besselj(0, k[3]*rayon[3]);
    A[6, 6] = mp.bessely(0, k[3]*rayon[3])
    return A


def createMatA_GITC(beta, omega, k, kappa2, mu, rayon, c) :
    A = mp.zeros(2*2+1)
    # condition conducteur <-> isolant
    A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
    A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
    A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
    A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
    A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
    A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);

    # condition GITC
    a = 0.5
    b = 0.5
    d = rayon[2] - rayon[1]
    Ab = mu[2] / kappa2[2] - b * mu[1] / kappa2[1] - (1-b) * mu[3] / kappa2[3]
    Bb = mu[2] - b * mu[1] - (1-b) * mu[3]
    # Cb = 1/mu[2] - b / mu[1] - (1-b) / mu[3]
    Db = kappa2[2] / (mu[2]) - b * kappa2[1]/mu[1] - (1-b) * kappa2[3]/mu[3]
    # Ab = mu[2] / (k[2]*k[2])
    # Bb = 0
    # # Cb = 1/mu[2] - b / mu[1] - (1-b) / mu[3]
    # Db = 0
    r_Gamma = rayon[1] + b * d
    # print("r_Gamma =", r_Gamma)
    A[2, 1] = mp.besselj(0, k[1]*r_Gamma) - d*(Bb - Ab*beta*beta*omega*omega)*(1-a)*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma)
    A[2, 2] = mp.bessely(0, k[1]*r_Gamma) - d*(Bb - Ab*beta*beta*omega*omega)*(1-a)*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma)
    A[2, 3] = -mp.besselj(0, k[3]*r_Gamma) - d*(Bb - Ab*beta*beta*omega*omega)*a*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma)
    A[2, 4] = -mp.bessely(0, k[3]*r_Gamma) - d*(Bb - Ab*beta*beta*omega*omega)*a*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma)

    A[3, 1] = -k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma) - d*Db*a*mp.besselj(0, k[1]*r_Gamma)
    A[3, 2] = -k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma) - d*Db*a*mp.bessely(0, k[1]*r_Gamma)
    A[3, 3] = k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma) - d*Db*(1-a)*mp.besselj(0, k[3]*r_Gamma)
    A[3, 4] = k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma) - d*Db*(1-a)*mp.bessely(0, k[3]*r_Gamma)

    # Condition de Dirichlet homogène assez loin
    # A[4, 3] = mp.besselj(0, k[3]*rayon[3]);
    # A[4, 4] = mp.bessely(0, k[3]*rayon[3])

    # Condition limite absorbante
    A[4, 3] = k[3]*c[3]*mu[3]*mp.besselj(1, k[3]*rayon[3]) + 1j * mp.sqrt(kappa2[3]) * mp.besselj(0, k[3]*rayon[3])
    A[4, 4] = k[3]*c[3]*mu[3]*mp.bessely(1, k[3]*rayon[3]) + 1j * mp.sqrt(kappa2[3]) * mp.bessely(0, k[3]*rayon[3])
    return A

def createMatA_ITC(beta, omega, eps, mu, sigma, rayon, c, k, kappa2) :
    print("sigma = ", sigma)
    A = mp.zeros(2*2+1)
    # condition conducteur <-> isolant
    A[0, 0] = mp.besselj(0, k[0]*rayon[0]);
    A[0, 1] = -mp.besselj(0, k[1]*rayon[0]);
    A[0, 2] = -mp.bessely(0, k[1]*rayon[0]);
    A[1, 0] = c[0]*k[0]*mp.besselj(1, k[0]*rayon[0]);
    A[1, 1] = -c[1]*k[1]*mp.besselj(1, k[1]*rayon[0]);
    A[1, 2] = -c[1]*k[1]*mp.bessely(1, k[1]*rayon[0]);

    # condition ITC
    d = rayon[2] - rayon[1]
    sigmaTilde = sigma[2] * d * d
    gamma = mp.exp(3*1j*mp.pi/4)*mp.sqrt(omega * mu[2] * sigmaTilde)
    A1 = - 1/(omega*omega)* 0.5 *(1/(eps[3] + 1j * sigma[3]/omega) + 1/(eps[1] + 1j * sigma[1]/omega))
    B1 = 2 * mu[2] / gamma * mp.tanh(gamma/2) - 0.5 * (mu[1] + mu[3])
    A2 = A1 / 4
    B2 = mu[2] / (2*gamma) * mp.coth(gamma/2) - (mu[1] + mu[3]) / 8
    A3 = - 1/(4 * omega*omega) * (1/(eps[3] + 1j * sigma[3]/omega) - 1/(eps[1] + 1j * sigma[1]/omega))
    B3 = - (mu[3] - mu[1]) / 4

    L1 = A1 * beta*beta*omega*omega - B1
    L2 = A2 * beta*beta*omega*omega - B2
    L3 = A3 * beta*beta*omega*omega - B3

    r_Gamma = rayon[1] + 0.5 * d
    # print("r_Gamma =", r_Gamma)
    A[2, 1] = -mp.besselj(0, k[1]*r_Gamma) - d*(0.5*L1-L3)*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma)
    A[2, 2] = -mp.bessely(0, k[1]*r_Gamma) - d*(0.5*L1-L3)*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma)
    A[2, 3] = mp.besselj(0, k[3]*r_Gamma) - d*(0.5*L1+L3)*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma)
    A[2, 4] = mp.bessely(0, k[3]*r_Gamma) - d*(0.5*L1+L3)*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma)

    A[3, 1] = 0.5*mp.besselj(0, k[1]*r_Gamma) - d*(0.5*L3-L2)*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma)
    A[3, 2] = 0.5*mp.bessely(0, k[1]*r_Gamma) - d*(0.5*L3-L2)*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma)
    A[3, 3] = 0.5*mp.besselj(0, k[3]*r_Gamma) - d*(0.5*L3+L2)*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma)
    A[3, 4] = 0.5*mp.bessely(0, k[3]*r_Gamma) - d*(0.5*L3+L2)*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma)

    # Condition de Dirichlet homogène assez loin
    A[4, 3] = mp.besselj(0, k[3]*rayon[3]);
    A[4, 4] = mp.bessely(0, k[3]*rayon[3])

    # Condition limite absorbante
    # A[4, 3] = k[3]*c[3]*mu[3]*mp.besselj(1, k[3]*rayon[3]) + 1j * mp.sqrt(kappa2[3]) * mp.besselj(0, k[3]*rayon[3])
    # A[4, 4] = k[3]*c[3]*mu[3]*mp.bessely(1, k[3]*rayon[3]) + 1j * mp.sqrt(kappa2[3]) * mp.bessely(0, k[3]*rayon[3])
    return A

# toutes les valeurs doivent être en mp précision
def ComputeDet(beta, f, eps, sigma, mu, rayon, gitc):
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    k = mp.zeros(4, 1)
    c = mp.zeros(4, 1)
    kappa2 = mp.zeros(4, 1)
    for i in range(4):
        Delta= -omega*omega*eps[i]*mu[i] - 1j*omega*sigma[i]*mu[i] + beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = mp.sqrt(-Delta)
        # print("k : ", k[i])
        c[i] = (- omega*omega*eps[i] - 1j*omega*sigma[i]) / Delta
        # print("c : ", c[i])
        kappa2[i] = omega*omega*mu[i]*eps[i] + 1j*omega*sigma[i]*mu[i]
    if gitc :
        # A = createMatA_GITC(beta, omega, k, kappa2, mu, rayon, c)
        A = createMatA_ITC(beta, omega, eps, mu, sigma, rayon, c, k, kappa2)
    else :
        A = createMatA(k, rayon, c)
    #print("k = ", k)
    return mp.det(A)


def det2d(a,b,c,d):
    return a*c-b*d

def Jac(x, y, f, eps, sigma, mu, rayon, gitc) :
    mp.mp.dps = 200
    J = mp.zeros(2)
    # schéma centré (ordre 2) pour les dérivées
    DxM1y = ComputeDet(mp.mpc(x - pasNewton, y), f, eps, sigma, mu, rayon, gitc)
    DxP1y = ComputeDet(mp.mpc(x + pasNewton, y), f, eps, sigma, mu, rayon, gitc)
    DxyM1 = ComputeDet(mp.mpc(x, y - pasNewton), f, eps, sigma, mu, rayon, gitc)
    DxyP1 = ComputeDet(mp.mpc(x, y + pasNewton), f, eps, sigma, mu, rayon, gitc)
    J[0, 0] = 0.5 * (mp.re(DxP1y) - mp.re(DxM1y)) / pasNewton
    J[1, 0] = 0.5 * (mp.im(DxP1y) - mp.im(DxM1y)) / pasNewton
    J[0, 1] = 0.5 * (mp.re(DxyP1) - mp.re(DxyM1)) / pasNewton
    J[1, 1] = 0.5 * (mp.im(DxyP1) - mp.im(DxyM1)) / pasNewton
    return J


def Newton(B0, expTol, f, eps, sigma, mu, rayon, gitc=False) :
    mp.mp.dps = 200
    print("B0 : ", B0)
    x = mp.re(B0); y = mp.im(B0)
    X0 = mp.matrix([x, y])
    DComp = ComputeDet(mp.mpc(x, y), f, eps, sigma, mu, rayon, gitc)
    D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
    Jac0 = Jac(x, y, f, eps, sigma, mu, rayon, gitc)
    i = 0
    vectDiff = mp.lu_solve(Jac0, D0)
    norme2 = (vectDiff.T * vectDiff)[0]
    X0 = X0 - vectDiff
    # mp.mp.dps = 15 # pour l'affichage
    # print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
    while norme2 >= 10**(-expTol * 2) :
        i += 1
        DComp = ComputeDet(mp.mpc(X0[0], X0[1]), f, eps, sigma, mu, rayon, gitc)
        D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
        Jac0 = Jac(X0[0], X0[1], f, eps, sigma, mu, rayon, gitc)
        vectDiff = mp.lu_solve(Jac0, D0)
        norme2 = (vectDiff.T * vectDiff)[0]
        X0 = X0 - vectDiff
        # mp.mp.dps = 15 # pour l'affichage
        # print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
        mp.mp.dps = 100
    return X0[0] +1j * X0[1]
    # print(X)

def ComputeParam(beta, f, eps, sigma, mu, rayon, gitc = False) :
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    k = mp.zeros(4, 1)
    c = mp.zeros(4, 1)
    kappa2 = mp.zeros(4, 1)
    for i in range(4):
        Delta= - omega*omega*eps[i]*mu[i] - 1j*omega*sigma[i]*mu[i] + beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = mp.sqrt( - Delta)
        c[i] = (- omega*omega*eps[i] - 1j*omega*sigma[i]) / Delta
        kappa2[i] = omega*omega*mu[i]*eps[i] + 1j*omega*sigma[i]*mu[i]
    # print("k = ", k)
    if gitc :
        A = createMatA_ITC(beta, omega, eps, mu, sigma, rayon, c, k, kappa2)
    else :
        A = createMatA(k, rayon, c)
    val, EL, vect = mp.eig(A, left = True, right = True)
    val, EL, vect = mp.eig_sort(val, EL, vect, f = lambda x: mp.fabs(x))
    # print("valeur propre nulle = ", val[0])
    alphas = vect[:,0]
    return k, alphas

# calcul des champs en un point r
def ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu) :
    mp.mp.dps = 200
    # Recherche de la couche
    omega = mp.mpf(f)*2*pi;
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
            Ez = alphas[0] * mp.besselj(0, k[0]*r)
            dEz = - alphas[0] * k[0] * mp.besselj(1, k[0]*r)
        else :
            #print("alpha", alphas[2*indCouche-1], alphas[2*indCouche], mp.besselj(0, k[indCouche]*r), mp.bessely(0, k[indCouche]*r), k[indCouche])
            Ez = alphas[2*indCouche-1] * mp.besselj(0, k[indCouche]*r) + alphas[2*indCouche] * mp.bessely(0, k[indCouche]*r)
            dEz = - alphas[2*indCouche-1] * k[indCouche] * mp.besselj(1, k[indCouche]*r) - alphas[2*indCouche] * k[indCouche] * mp.bessely(1, k[indCouche]*r)
        Er = -1j*beta * omega * dEz / Delta
        Ht = -1j * omega * eps_tild * dEz / Delta
        return Ez, dEz, Er, Ht
    except :
        return 0.0, 0.0, 0.0, 0.0

def ComputeChamps_GITC(r, f, k, alphas, beta, rayon, eps, sigma, mu) :
    mp.mp.dps = 200
    # Recherche de la couche
    omega = mp.mpf(f)*2*pi;
    indCouche = len(rayon)
    # print("r = ", r)
    for i, ray in enumerate(rayon) :
        if i == 0 :
            if r < ray :
                indCouche = 0
        else :
            if r >= rayon[i-1] and r <= ray :
                indCouche = i
    if indCouche == 2 :
        d = rayon[2] - rayon[1]
        b = 0.5
        r_Gamma = rayon[1] + b * d
        if r < r_Gamma :
            indCouche = 1
        else :
            indCouche = 3
    iAlpha = min(2, indCouche)
    eps_tild = eps[indCouche] + 1j * sigma[indCouche]/omega
    Delta = - omega*omega*eps[indCouche]*mu[indCouche] - 1j*omega*sigma[indCouche]*mu[indCouche] + beta*beta*omega*omega
    Ez = 0.0
    dEz = 0.0
    if indCouche == 0 :
        Ez = alphas[0] * mp.besselj(0, k[0]*r)
        dEz = - alphas[0] * k[0] * mp.besselj(1, k[0]*r)
    else :
        #print("alpha", alphas[2*indCouche-1], alphas[2*indCouche], mp.besselj(0, k[indCouche]*r), mp.bessely(0, k[indCouche]*r), k[indCouche])
        Ez = alphas[2*iAlpha-1] * mp.besselj(0, k[indCouche]*r) + alphas[2*iAlpha] * mp.bessely(0, k[indCouche]*r)
        dEz = - alphas[2*iAlpha-1] * k[indCouche] * mp.besselj(1, k[indCouche]*r) - alphas[2*iAlpha] * k[indCouche] * mp.bessely(1, k[indCouche]*r)
    Er = -1j*beta * omega * dEz / Delta
    Ht = -1j * omega * eps_tild * dEz / Delta
    return Ez, dEz, Er, Ht
    # except :
    #     return 0.0, 0.0, 0.0, 0.0



# calcul des champs sur un intervalle de r
def ComputeChampsR(rMin, rMax, pas, f, k, k_GITC, alphas, alphas_GITC, beta, beta_GITC, rayon, eps, sigma, mu) :
    taille_1 = int((rMax - rMin) / pas)
    rr = mp.linspace(rMin, rMax, taille_1)
    EEz = mp.zeros(taille_1, 1)
    EEr = mp.zeros(taille_1, 1)
    HHt = mp.zeros(taille_1, 1)
    print("Calcul des champs")
    champs = np.array([ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu) for r in rr])
    champs_GITC = np.array([ComputeChamps_GITC(r, f, k_GITC, alphas_GITC, beta_GITC, rayon, eps, sigma, mu) for r in rr])
    EEz = champs[:,0]
    EEr = champs[:,2]
    HHt = champs[:,3]
    EEz_GITC = champs_GITC[:,0]
    EEr_GITC = champs_GITC[:,2]
    HHt_GITC = champs_GITC[:,3]
    print("Fin du calcul")
    return rr, EEz, EEr, HHt, EEz_GITC, EEr_GITC, HHt_GITC

# programme principal (le faire tourner donne le mode propagatif)
def ComputeModes(ficName, rMin, rMax, pas) :
    mp.mp.dps = 200
    ptInit = mp.mpc(2.0,1.0)
    n, f, eps, sigma, mu, rayon = getVariables(ficName)
    Beta = Newton(ptInit, 60, f, eps, sigma, mu, rayon)
    mp.mp.dps = 15
    print("VP : ", Beta)
    mp.mp.dps = 200
    k, alphas = ComputeParam(Beta, f, eps, sigma, mu, rayon)
    # print("Alpha, k = ", alphas, k)

    # autre calcul, condition GITC
    print("___Calcul GITC___")
    Beta_GITC = Newton(ptInit, 60, f, eps, sigma, mu, rayon, True)
    mp.mp.dps = 15
    print("VP : ", Beta_GITC)
    mp.mp.dps = 200
    k_GITC, alphas_GITC = ComputeParam(Beta_GITC, f, eps, sigma, mu, rayon, True)
    # print("Alpha, k = ", alphas_GITC, k_GITC)
    rr, EEz, EEr, HHt, EEz_GITC, EEr_GITC, HHt_GITC = ComputeChampsR(rMin, rMax, pas, f, k, k_GITC, alphas, alphas_GITC, Beta, Beta_GITC, rayon, eps, sigma, mu)
    return rr, EEz, EEr, HHt, EEz_GITC, EEr_GITC, HHt_GITC


if __name__ == "__main__":
    rr, EEz, EEr, HHt, EEz_GITC, EEr_GITC, HHt_GITC = ComputeModes("donneesCable2", 0, 0.3, 0.001)
    coef = EEr_GITC[91]/EEr[91]
    plt.figure(1)
    plt.semilogy(rr, abs(coef*EEr), "r-")
    plt.semilogy(rr, abs(EEr_GITC), "b--")
    coef2 = HHt_GITC[91]/HHt[91]
    plt.figure(2)
    plt.semilogy(rr, abs(coef2*HHt), "r-")
    plt.semilogy(rr, abs(HHt_GITC), "b--")
    plt.show()
    plt.figure(3)
    plt.semilogy(rr, abs(coef2*EEz), "r-")
    plt.semilogy(rr, abs(EEz_GITC), "b--")
    plt.show()
