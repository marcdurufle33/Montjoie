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

def createMatA(n, k, rayon, c) :
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
    return A

# toutes les valeurs doivent être en mp précision
def ComputeDet(beta, n, f, eps, sigma, mu, rayon):
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    k = mp.zeros(n+1, 1)
    c = mp.zeros(n+1, 1)
    for i in range(n+1):
        Delta= omega*omega*eps[i]*mu[i] + 1j*omega*sigma[i]*mu[i]-beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = mp.sqrt(Delta)
        # print("k : ", k[i])
        c[i] = (omega*omega*eps[i] + 1j*omega*sigma[i]) / Delta
        # print("c : ", c[i])
    A = createMatA(n, k, rayon, c)
    #print("k = ", k)
    return mp.det(A)


def det2d(a,b,c,d):
    return a*c-b*d

def Jac(x, y, n, f, eps, sigma, mu, rayon) :
    mp.mp.dps = 200
    J = mp.zeros(2)
    # schéma centré (ordre 2) pour les dérivées
    DxM1y = ComputeDet(mp.mpc(x - pasNewton, y), n, f, eps, sigma, mu, rayon)
    DxP1y = ComputeDet(mp.mpc(x + pasNewton, y), n, f, eps, sigma, mu, rayon)
    DxyM1 = ComputeDet(mp.mpc(x, y - pasNewton), n, f, eps, sigma, mu, rayon)
    DxyP1 = ComputeDet(mp.mpc(x, y + pasNewton), n, f, eps, sigma, mu, rayon)
    J[0, 0] = 0.5 * (mp.re(DxP1y) - mp.re(DxM1y)) / pasNewton
    J[1, 0] = 0.5 * (mp.im(DxP1y) - mp.im(DxM1y)) / pasNewton
    J[0, 1] = 0.5 * (mp.re(DxyP1) - mp.re(DxyM1)) / pasNewton
    J[1, 1] = 0.5 * (mp.im(DxyP1) - mp.im(DxyM1)) / pasNewton
    return J


def Newton(B0, expTol, n, f, eps, sigma, mu, rayon) :
    mp.mp.dps = 200
    print("B0 : ", B0)
    x = mp.re(B0); y = mp.im(B0)
    X0 = mp.matrix([x, y])
    DComp = ComputeDet(mp.mpc(x, y), n, f, eps, sigma, mu, rayon)
    D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
    Jac0 = Jac(x, y, n, f, eps, sigma, mu, rayon)
    i = 0
    vectDiff = mp.lu_solve(Jac0, D0)
    norme2 = (vectDiff.T * vectDiff)[0]
    X0 = X0 - vectDiff
    # mp.mp.dps = 15 # pour l'affichage
    print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
    while norme2 >= 10**(-expTol * 2) :
        i += 1
        DComp = ComputeDet(mp.mpc(X0[0], X0[1]), n, f, eps, sigma, mu, rayon)
        D0 = mp.matrix([mp.re(DComp), mp.im(DComp)])
        Jac0 = Jac(X0[0], X0[1], n, f, eps, sigma, mu, rayon)
        vectDiff = mp.lu_solve(Jac0, D0)
        norme2 = (vectDiff.T * vectDiff)[0]
        X0 = X0 - vectDiff
        # mp.mp.dps = 15 # pour l'affichage
        print("  Étape ", i, " : beta = ", mp.mpc(X0[0], X0[1]), " erreur = ", norme2)
        mp.mp.dps = 100
    return X0[0] +1j * X0[1]
    # print(X)

def ComputeParam(beta, n, f, eps, sigma, mu, rayon) :
    mp.mp.dps = 200
    omega = mp.mpf(f)*2*pi;
    k = mp.zeros(n+1, 1)
    c = mp.zeros(n+1, 1)
    for i in range(n+1):
        Delta= - omega*omega*eps[i]*mu[i] - 1j*omega*sigma[i]*mu[i] + beta*beta*omega*omega
        # print("Delta : ", Delta)
        k[i] = mp.sqrt( - Delta)
        c[i] = (- omega*omega*eps[i] - 1j*omega*sigma[i]) / Delta
    # print("k = ", k)
    A = createMatA(n, k, rayon, c)
    val, EL, vect = mp.eig(A, left = True, right = True)
    val, EL, vect = mp.eig_sort(val, EL, vect, f = lambda x: mp.fabs(x))
    print("valeur propre nulle = ", val[0])
    alphas = vect[:,0]
    return k, alphas

# calcul des champs en un point r
def ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu) :
    mp.mp.dps = 200
    # Recherche de la couche
    omega = mp.mpf(f)*2*pi;
    indCouche = len(rayon)
    print("r = ", r)
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

# calcul des champs sur un intervalle de r
def ComputeChampsR(rMin, rMax, pas, f, k, alphas, beta, rayon, eps, sigma, mu) :
    taille_1 = int((rMax - rMin) / pas)
    rr = mp.linspace(rMin, rMax, taille_1)
    EEz = mp.zeros(taille_1, 1)
    EEr = mp.zeros(taille_1, 1)
    HHt = mp.zeros(taille_1, 1)
    print("Calcul des champs")
    champs = np.array([ComputeChamps(r, f, k, alphas, beta, rayon, eps, sigma, mu) for r in rr])
    EEz = champs[:,0]
    EEr = champs[:,2]
    HHt = champs[:,3]
    print("Fin du calcul")
    return rr, EEz, EEr, HHt

# programme principal (le faire tourner donne le mode propagatif)
def ComputeModes(ficName, rMin, rMax, pas) :
    mp.mp.dps = 200
    ptInit = mp.mpc(2.0,1.0)
    n, f, eps, sigma, mu, rayon = getVariables(ficName)
    Beta = Newton(ptInit, 50, n, f, eps, sigma, mu, rayon)
    mp.mp.dps = 15
    print("Résultat : ", Beta, ", Det = ", ComputeDet(Beta, n, f, eps, sigma, mu, rayon))
    mp.mp.dps = 200
    k, alphas = ComputeParam(Beta, n, f, eps, sigma, mu, rayon)
    print("Alpha, k = ", alphas, k)
    rr, EEz, EEr, HHt = ComputeChampsR(rMin, rMax, pas, f, k, alphas, Beta, rayon, eps, sigma, mu)
    return rr, EEz, EEr, HHt


if __name__ == "__main__":
    rr, EEz, EEr, HHt = ComputeModes("donneesCable60kHz", 0, 2.0, 0.001)
    plt.semilogy(rr, abs(EEz), "r-")
    plt.show()
