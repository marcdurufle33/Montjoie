from pylab import *

def GetStability(zTab, A, B):
    """ Calcul du domaine de stabilite pour une methode de Runge-Kutta
    zTab : valeurs propres de dt A pour lesquelles on veut savoir si on a un algo stable 
    A, B : coefficients de la methode de Runge-Kutta
    renvoie un tableau qui contient des 1 (instables) et des 0 (stables) """
    stab = ones(zTab.shape)
    dim = len(zTab.shape)
    k = zeros(A.shape[0]) + 1j*zeros(A.shape[0])
    for i1 in range(zTab.shape[0]):
        N = 1
        if (dim == 2):
            N = zTab.shape[1]
        
        for i2 in range(N):
            if (dim == 2):
                z = zTab[i1, i2]
            else:
                z = zTab[i1]
                
            G = 1.0+0j
            for i in range(A.shape[0]):
                k[i] = z+0j
                for j in range(i):
                    k[i] += A[i, j]*z*k[j]
                
                k[i] = k[i]/(1.0-A[i, i]*z+0j)
                G += B[i]*k[i]
            
            if (abs(G) <= 1+1e-12):
                if (dim == 2):
                    stab[i1, i2] = 0
                else:
                    stab[i1] = 0

    return stab

def GetStabilityMultistep(zTab, A, B):
    G = zeros(zTab.shape)
    for i in range(len(zTab)):
        P = poly1d(A) + zTab[i]*poly1d(B)
        G[i] = abs(P.r).max()
        
    return G

def GetStabilityNystrom(zTab, Abar, Bbar, B, C):
    """ Calcul du domaine de stabilite pour une methode de Runge-Kutta Nystrom
    z : valeurs propres de dt^2 A pour lesquelles on veut savoir si on a un algo stable
    Abar, Bbar, B, C : coefficients de Runge-Kutta 
    renvoie un tableau qui contient des 1 (instables) et des 0 (stables) """
    stab = ones(zTab.shape)
    G = zeros(zTab.shape)
    G1 = zeros(zTab.shape)+1j*zeros(zTab.shape)
    G2 = zeros(zTab.shape)+1j*zeros(zTab.shape)
    # Runge-Kutta Nystrom :
    # k_i = A ( y^n + c_i dt yprime^n + dt^2 \sum_j Abar_{i, j} k_j)
    # y^{n+1} = y^n + dt yprime^n + dt^2 \sum_i Bbar_i k_i
    # yprime^{n+1} = yprime^n + dt \sum_i B_i k_i
    for i1 in range(zTab.shape[0]):
        z = zTab[i1]
        kw = zeros(Abar.shape[0])
        ky = zeros(Abar.shape[0])
        # on calcule la matrice D telle que :
        # y^{n+1}  = D00 y^n + D01 w^n
        # w^{n+1} = D10 y^n + D11 w^n
        # On note w^n = yprime^n / (dt A) et z = dt^2 A
        # on a alors les relations :
        # dt^2 k_i = z y^n + c_i z^2 w^n + z \sum_j Abar_{i, j} dt^2 k_j
        # y^{n+1} = y + z w^n + \sum_i Bbar_i dt^2 k_i
        # w^{n+1} = w_n + 1/z \sum_i B_i dt^2 k_i
        # on utilise ces relations pour construire la matrice D
        D = array([[1.0, z],[0,1.0]])
        for i in range(Abar.shape[0]):
            kw[i] = C[i]*z*z
            ky[i] = z
            for j in range(i):
                kw[i] += Abar[i, j]*z*kw[j]
                ky[i] += Abar[i, j]*z*ky[j]
            
            D[0, 0] += Bbar[i]*ky[i]
            D[0, 1] += Bbar[i]*kw[i]
            D[1, 0] += 1.0/z*B[i]*ky[i]
            D[1, 1] += 1.0/z*B[i]*kw[i]
            
        # si les vps de D sont de module inferieur a 1, la methode est stable
        if (z == 0):
            G[i1] = 1.0
            stab[i1] = 0.0
        else :
            L, V = eig(D)
            G[i1] = abs(L).max()
            G1[i1] = L[0]
            G2[i1] = L[1]
            if (G[i1] <= 1):
                stab[i1] = 0.0
    
    return G, stab, G1, G2


def GetStabilityMatrix(Abar, Bbar, B, C):
    # Runge-Kutta Nystrom :
    # k_i = A ( y^n + c_i dt yprime^n + dt^2 \sum_j Abar_{i, j} k_j)
    # y^{n+1} = y^n + dt yprime^n + dt^2 \sum_i Bbar_i k_i
    # yprime^{n+1} = yprime^n + dt \sum_i B_i k_i
   
    # on calcule la matrice D telle que :
    # y^{n+1}  = D00 y^n + D01 w^n
    # w^{n+1} = D10 y^n + D11 w^n
    # On note w^n = yprime^n / (dt A) et z = dt^2 A
    # on a alors les relations :
    # dt^2 k_i = z y^n + c_i z^2 w^n + z \sum_j Abar_{i, j} dt^2 k_j
    # y^{n+1} = y + z w^n + \sum_i Bbar_i dt^2 k_i
    # w^{n+1} = w_n + 1/z \sum_i B_i dt^2 k_i
    # on utilise ces relations pour construire la matrice D
    # on stocke directement les polynomes en z
    D = [[poly1d([1]), poly1d([1, 0])], [poly1d([0]), poly1d([1])]]
    kw = [0]*Abar.shape[0]
    ky = [0]*Abar.shape[0]
    for i in range(Abar.shape[0]):
        kw[i] = poly1d([C[i], 0, 0])
        ky[i] = poly1d([1, 0])
        for j in range(i):
            kw[i] += float(Abar[i, j])*poly1d([1,0])*kw[j]
            ky[i] += float(Abar[i, j])*poly1d([1,0])*ky[j]
            
        ky_div_z = poly1d(ky[i].coeffs[0:ky[i].order])
        kw_div_z = poly1d(kw[i].coeffs[0:kw[i].order])
        
        D[0][0] += Bbar[i]*ky[i]
        D[0][1] += Bbar[i]*kw[i]
        D[1][0] += B[i]*ky_div_z
        D[1][1] += B[i]*kw_div_z

    return D
