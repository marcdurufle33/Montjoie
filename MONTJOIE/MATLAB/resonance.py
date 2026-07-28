from scipy.special import *
from pylab import *
from visuND import *
import mpmath as mp

def EigenDisque(r, nb_bessel, nb_zeros, nb_eig):
    val = array([])
    num = []
    for m in range(nb_bessel):
        L = jn_zeros(m, nb_zeros)**2/r**2
        val = append(val, L)
        num = num + [m]*nb_zeros
        if (m != 0):
            val = append(val, L)
            num = num + [m]*nb_zeros
    
    ind = argsort(val)
    all_eig = val[ind]
    num_mode = array(num)[ind]
    return num_mode[0:nb_eig], all_eig[0:nb_eig]

def EigenSphere(r, nb_bessel, nb_zeros, nb_eig):
    val = array([])
    num = []
    for m in range(nb_bessel):
        for k in range(nb_zeros):
            L = mp.besseljzero(m+0.5, k+1)**2/r**2
            val = append(val, float(L))
            num = append(num, m)
    
    ind = argsort(val)
    all_eig = val[ind]
    num_mode = array(num)[ind]
    return num_mode[0:nb_eig], all_eig[0:nb_eig]
                

def AfficheNumericalBessel(num, L, Lexact, m, racine_fichier):
    if (m == 0):
        pas = 1
    else:
        pas = 2
    
    for i in range(0, len(num), pas):
        [X, Y, Z, coor, V1] = loadND(racine_fichier + EntierToString(num[i]) + ".dat")
        if (m == 0):
            Z = V1[100:,100]
        else:
            [X, Y, Z, coor, V2] = loadND(racine_fichier + EntierToString(num[i+1]) + ".dat")
            Z1 = V1[100:,100]
            Z2 = V2[100:,100]
            if (norm(Z2) > norm(Z1)):
                Z = V2[100:,100]
            else:
                Z = V1[100:,100]
        
        R = Y[100:]
        Zref = jn(m, sqrt(Lexact[num[i]])*R)
        coef = abs(Zref).max() / abs(Z).max()
        err_p = norm(coef*Z-Zref)
        err_m = norm(coef*Z+Zref)
        print "Erreur = ", min(err_m, err_p)/norm(Zref)
        clf();
        if (err_m < err_p):
            plot(R, -coef*Z); 
        else:
            plot(R, coef*Z);
        
        plot(R, Zref)
        pause(1.0)

