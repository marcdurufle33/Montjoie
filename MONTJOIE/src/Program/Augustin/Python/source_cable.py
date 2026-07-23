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

def createMatAdir(k, rayon, c) :
    A = mp.zeros(6)
    # condition Dirichlet pour le premier rayon
    A[0, 0] = mp.besselj(0, k[1]*rayon[0]);
    A[0, 1] = mp.bessely(0, k[1]*rayon[0]);

    # condition isolant <-> blindage
    A[1, 0] = mp.besselj(0, k[1]*rayon[1]);
    A[1, 1] = mp.bessely(0, k[1]*rayon[1]);
    A[1, 2] = -mp.besselj(0, k[2]*rayon[1]);
    A[1, 3] = -mp.bessely(0, k[2]*rayon[1]);
    # Continuité de la dérivée avec le coefficient c = 1/mu
    A[2, 0] = c[1]*k[1]*mp.besselj(1, k[1]*rayon[1]);
    A[2, 1] = c[1]*k[1]*mp.bessely(1, k[1]*rayon[1]);
    A[2, 2] = -c[2]*k[2]*mp.besselj(1, k[2]*rayon[1]);
    A[2, 3] = -c[2]*k[2]*mp.bessely(1, k[2]*rayon[1]);

    # condition blindage <-> eau
    A[3, 2] = mp.besselj(0, k[2]*rayon[2]);
    A[3, 3] = mp.bessely(0, k[2]*rayon[2]);
    A[3, 4] = -mp.besselj(0, k[3]*rayon[2]);
    A[3, 5] = -mp.bessely(0, k[3]*rayon[2]);
    # Continuité de la dérivée
    A[4, 2] = c[2]*k[2]*mp.besselj(1, k[2]*rayon[2]);
    A[4, 3] = c[2]*k[2]*mp.bessely(1, k[2]*rayon[2]);
    A[4, 4] = -c[3]*k[3]*mp.besselj(1, k[3]*rayon[2]);
    A[4, 5] = -c[3]*k[3]*mp.bessely(1, k[3]*rayon[2]);

    # Condition de Dirichlet homogène assez loin
    A[5, 4] = mp.besselj(0, k[3]*rayon[3]);
    A[5, 5] = mp.bessely(0, k[3]*rayon[3])
    return A

def createMatAdir_GITC(k, omega, rayon, c, mu, sigma, type_gitc = 0) :
    A = mp.zeros(4)
    # condition Dirichlet pour le premier rayon
    A[0, 0] = mp.besselj(0, k[1]*rayon[0]);
    A[0, 1] = mp.bessely(0, k[1]*rayon[0]);

    d = rayon[2] - rayon[1]
    r_Gamma = rayon[1] + 0.5*d
    if (type_gitc == 0):
        # condition [ n \times E] = delta B Lambda
        Ab = mu[2] / k[2]**2 - 0.5 * mu[1] / k[1]**2 - 0.5 * mu[3] / k[3]**2
        Bb = mu[2] - 0.5 * mu[1] - 0.5 * mu[3]
        Dd = k[2]**2 / mu[2] - 0.5 * k[1]**2/mu[1] - 0.5 * k[3]**2/mu[3]
        A[1, 0] = mp.besselj(0, k[1]*r_Gamma) - 0.5*d*Bb*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma);
        A[1, 1] = mp.bessely(0, k[1]*r_Gamma) - 0.5*d*Bb*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma);
        A[1, 2] = -mp.besselj(0, k[3]*r_Gamma)- 0.5*d*Bb*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma);
        A[1, 3] = -mp.bessely(0, k[3]*r_Gamma)- 0.5*d*Bb*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma);
        # condition [n/mu \times rot E] = delta D <E_T>
        A[2, 0] = -c[1]*k[1]*mp.besselj(1, k[1]*r_Gamma) - 0.5*d*Dd*mp.besselj(0, k[1]*r_Gamma)
        A[2, 1] = -c[1]*k[1]*mp.bessely(1, k[1]*r_Gamma) - 0.5*d*Dd*mp.bessely(0, k[1]*r_Gamma)
        A[2, 2] = c[3]*k[3]*mp.besselj(1, k[3]*r_Gamma) - 0.5*d*Dd*mp.besselj(0, k[3]*r_Gamma)
        A[2, 3] = c[3]*k[3]*mp.bessely(1, k[3]*r_Gamma) - 0.5*d*Dd*mp.bessely(0, k[3]*r_Gamma)
    else:
        # cas fortement conducteur
        sigma_tilde = sigma[2]*d*d
        gam = exp(3j*pi/4)*sqrt(omega*mu[2]*sigma_tilde)
        B1 = 2.0*mu[2] / gam*mp.tanh(0.5*gam) - 0.5*(mu[1] + mu[3])
        B2 = 0.5*mu[2] / (gam*mp.tanh(0.5*gam)) - 0.125*(mu[1] + mu[3])
        B3 = -0.25*(mu[3] - mu[1]);
        L1 = -B1; L2 = -B2; L3 = -B3;
        A[1, 0] = -mp.besselj(0, k[1]*r_Gamma) - d*(L1*0.5 - L3)*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma);
        A[1, 1] = -mp.bessely(0, k[1]*r_Gamma) - d*(L1*0.5 - L3)*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma);
        A[1, 2] = mp.besselj(0, k[3]*r_Gamma) - d*(L1*0.5 + L3)*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma);
        A[1, 3] = mp.bessely(0, k[3]*r_Gamma) - d*(L1*0.5 + L3)*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma);

        A[2, 0] = 0.5*mp.besselj(0, k[1]*r_Gamma) - d*(L3*0.5 - L2)*k[1]*c[1]*mp.besselj(1, k[1]*r_Gamma);
        A[2, 1] = 0.5*mp.bessely(0, k[1]*r_Gamma) - d*(L3*0.5 - L2)*k[1]*c[1]*mp.bessely(1, k[1]*r_Gamma);
        A[2, 2] = 0.5*mp.besselj(0, k[3]*r_Gamma) - d*(L3*0.5 + L2)*k[3]*c[3]*mp.besselj(1, k[3]*r_Gamma);
        A[2, 3] = 0.5*mp.bessely(0, k[3]*r_Gamma) - d*(L3*0.5 + L2)*k[3]*c[3]*mp.bessely(1, k[3]*r_Gamma);
    
    # Condition de Dirichlet homogène assez loin
    A[3, 2] = mp.besselj(0, k[3]*rayon[3]);
    A[3, 3] = mp.bessely(0, k[3]*rayon[3])
    return A


def ComputeChampsR(r, f, k, alpha, alpha_G, rayon, eps, sigma, mu) :
    Ez = mp.zeros(len(r), 1)
    Ht = mp.zeros(len(r), 1)
    Ez_G = mp.zeros(len(r), 1)
    Ht_G = mp.zeros(len(r), 1)
    # pour le premier rayon, on cherche la couche
    num = 0
    for i in range(len(r)):
        if (r[0] < rayon[i]):
            num = i
            break

    i0 = 0
    if (abs(r[0] - rayon[0]) < 1e-12):
        num = 1
    
    # cas particulier ou r[0] < rayon, Ez = 0 pour les premiers points
    if (num == 0):
        num = 1
        for i in range(len(r)):
            if (r[i] >= rayon[0]):
                i0 = i
                break

    # on boucle sur les points a partir de j0
    lastR = rayon[len(rayon)-1]
    for i in range(i0, len(r)):
        # on change de couche si necessaire
        if (r[i] > lastR):
            break;

        while (r[i] > rayon[num]):
            num += 1

        evalJ0 = mp.besselj(0, k[num]*r[i])
        evalJ1 = mp.besselj(1, k[num]*r[i])
        evalY0 = mp.bessely(0, k[num]*r[i])
        evalY1 = mp.bessely(1, k[num]*r[i])
        Ez[i] = alpha[2*num-2]*evalJ0 + alpha[2*num-1] * evalY0
        Ht[i] = alpha[2*num-2]*evalJ1 + alpha[2*num-1] * evalY1
        Ht[i] *= mp.mpf('1.0') *k[num] / mu[num]

        Ez_G[i] = alpha_G[2*num-2]*evalJ0 + alpha_G[2*num-1] * evalY0
        Ht_G[i] = alpha_G[2*num-2]*evalJ1 + alpha_G[2*num-1] * evalY1
        Ht_G[i] *= mp.mpf('1.0') *k[num] / mu[num]
    
    return Ez, Ht, Ez_G, Ht_G

def SolveSourceProblem(f, eps, sigma, mu, rayon, type_gitc = 0):
    omega = 2*pi*f
    k = mp.zeros(len(eps), 1)
    c = mp.zeros(len(eps), 1)
    for i in range(len(eps)):
        Delta = omega*omega*eps[i]*mu[i] + 1j*omega*sigma[i]*mu[i]
        k[i] = mp.sqrt(Delta)
        c[i] = mp.mpf('1.0')/mu[i]

    A = createMatAdir(k, rayon, c)
    rhs = mp.zeros(6, 1)
    rhs[0] = mp.mpf('1.0')
    alpha = mp.lu_solve(A, rhs)

    A = createMatAdir_GITC(k, omega, rayon, c, mu, sigma, type_gitc)
    rhs = mp.zeros(4, 1)
    rhs[0] = mp.mpf('1.0')
    aG = mp.lu_solve(A, rhs)
    alpha_G = mp.zeros(6, 1)
    alpha_G[0:2] = aG[0:2]; alpha_G[4:6] = aG[2:4]
    return k, alpha, alpha_G

# programme principal (le faire tourner donne le mode propagatif)
def ComputeSolution(ficName, rMin, rMax, Nr, type_gitc = 0) :
    r = mp.linspace(rMin, rMax, Nr)
    n, f, eps, sigma, mu, rayon = getVariables(ficName)
    k, alpha, alpha_G = SolveSourceProblem(f, eps, sigma, mu, rayon, type_gitc)
    print("Résultat : ", array(alpha, dtype='complex128'))
    
    Ez, Ht, Ez_G, Ht_G = ComputeChampsR(r, f, k, alpha, alpha_G, rayon, eps, sigma, mu)
    return np.array(r,dtype='float64'), np.array(k, dtype='complex128'), np.array(Ez, dtype='complex128'), np.array(Ht, dtype='complex128'), np.array(Ez_G, dtype='complex128'), np.array(Ht_G, dtype='complex128')
