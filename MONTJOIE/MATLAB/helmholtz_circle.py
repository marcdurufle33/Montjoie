from mpmath import *

mp.dps = 100;

N_max = 80;

rho = matrix([mpf('2.0'), mpf('3.5'), mpf('1.0')])
mu = matrix([mpf('1.5'), mpf('2.0'), mpf('1.0')])
sigma = matrix([mpf('0.5'), mpf('50.0'), mpf('0.0')])

radius = matrix([mpf('0.0'), mpf('0.95'), mpf('1.05'), mpf('1.5')])

impedance_condition = True
beta_impedance = mpf('0')

freq = mpf('0.5')
omega = 2.0*pi*freq

first_order_abc = True

kwave = rho.copy()
for i in range(rho.rows):
    eps = rho[i]
    eps += 1j*sigma[i]/omega
    kwave[i] = omega*sqrt(eps/mu[i])

def GetHnPrime(z, n):
    jn = besselj(n, z)
    yn = bessely(n, z)
    hn = jn + 1j*yn
    
    if (n == 0):
        jn_next = besselj(n+1, z)
        yn_next = bessely(n+1, z)
        hn_next = jn_next + 1j*yn_next

        jn_prime = -jn_next
        hn_prime = -hn_next
    else:
        jn_prev = besselj(n-1, z)
        yn_prev = bessely(n-1, z)
        hn_prev = jn_prev + 1j*yn_prev
        
        jn_next = besselj(n+1, z)
        yn_next = bessely(n+1, z)
        hn_next = jn_next + 1j*yn_next
        
        jn_prime = 0.5*(jn_prev - jn_next)
        hn_prime = 0.5*(hn_prev - hn_next)
        
    hn2 = -hn + 2.0*jn
    hn2_prime = -hn_prime + 2.0*jn_prime
    
    return hn, hn2, hn_prime, hn2_prime, jn, jn_prime

# boucle sur tous les modes
N = rho.rows
alpha = matrix(N_max, N)
beta = matrix(N_max, N)
for n in range(0, N_max):
    sys_alpha = zeros(2*N)
    rhs_alpha = zeros(2*N, 1)
    
    # premiere equation donnee avec la condition pour r = r0
    if (radius[0] == mpf(0)):
        # cas d'une disque dielectrique
        # alpha_n = beta_n
        sys_alpha[0, 0] = 1.0
        sys_alpha[0, 1] = -1.0
    else:
        hn, hn2, hn_prime, hn2_prime, jn, jn_prime = GetHnPrime(kwave[0]*radius[0], n)
        
        hn_prime = hn_prime*kwave[0]
        hn2_prime = hn2_prime*kwave[0]
        jn_prime = jn_prime*kwave[0]

        if (impedance_condition):
            # condition de Robin -du/dn + beta u = 0
            sys_alpha[0, 0] = -beta_impedance*hn + hn_prime
            sys_alpha[0, 1] = -beta_impedance*hn2 + hn2_prime
            if (N == 1):
                rhs_alpha[0] = -(-beta_impedance*jn + jn_prime)
        else:
            # condition de Dirichlet
            sys_alpha[0, 0] = hn
            sys_alpha[0, 1] = hn2
            if (N == 1):
                rhs_alpha[0] = -jn
        
    
    # boucle sur les interfaces
    num = 1
    for i in range(1, radius.rows-1):
        hn, hn2, hn_prime, hn2_prime, jn, jn_prime = GetHnPrime(kwave[i-1]*radius[i], n)
        hn_next, hn2_next, hn_prime_next, hn2_prime_next, jn_next, jn_prime_next = GetHnPrime(kwave[i]*radius[i], n)
        
        hn_prime = hn_prime*kwave[i-1]
        hn2_prime = hn2_prime*kwave[i-1]
        jn_prime = jn_prime*kwave[i-1]
        hn_prime_next = hn_prime_next*kwave[i]
        hn2_prime_next = hn2_prime_next*kwave[i]
        jn_prime_next = jn_prime_next*kwave[i]

        # continuite de la solution
        sys_alpha[num, 2*i-2] = hn
        sys_alpha[num, 2*i-1] = hn2
        sys_alpha[num, 2*i] = -hn_next
        sys_alpha[num, 2*i+1] = -hn2_next
        if (i == N-1):
            rhs_alpha[num,0] = jn_next
          
        num += 1
          
        # continuity of mu du_dn
        sys_alpha[num, 2*i-2] = hn_prime * mu[i-1]
        sys_alpha[num, 2*i-1] = hn2_prime * mu[i-1]
        sys_alpha[num, 2*i] = -hn_prime_next * mu[i]
        sys_alpha[num, 2*i+1] = -hn2_prime_next * mu[i]
        if (i == N-1):
            rhs_alpha[num,0] = jn_prime_next * mu[i]
            
        num += 1
        
    # last condition due to Sommerfeld condition (at finite distance or not)
    if (first_order_abc):
        hn, hn2, hn_prime, hn2_prime, jn, jn_prime = GetHnPrime(kwave[N-1]*radius[N], n)
        sys_alpha[num, 2*N-2] = hn_prime - 1j*hn
        sys_alpha[num, 2*N-1] = hn2_prime - 1j*hn2
    else:
        # exact Sommerfeld condition => no component in Hn_n^{(2)}
        sys_alpha[num, 2*N-1] = 1.0
                
    sol_alpha = lu_solve(sys_alpha, rhs_alpha)
    
    print "Norm rhs = ", float(norm(rhs_alpha))
    for i in range(N):
        alpha[n, i] = sol_alpha[2*i]
        beta[n, i] = sol_alpha[2*i+1]
        

def WriteText(A, nom):
    file_out = open(nom, 'w')
    for i in range(A.rows):
        for j in range(A.cols):
            file_out.write("(" + str(A[i, j].real) + "," + str(A[i, j].imag) + ") ")
        
        file_out.write("\n")
    
    file_out.close()

# on imprime les coefs dans un fichier
WriteText(alpha, "alpha.dat")
WriteText(beta, "beta.dat")
