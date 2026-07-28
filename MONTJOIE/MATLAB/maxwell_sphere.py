from mpmath import *

mp.dps = 200;

N_max = 30;

epsilon = matrix([mpf('2.0'), mpf('3.5'), mpf('1.0')])
mu = matrix([mpf('1.5'), mpf('2.0'), mpf('1.0')])
sigma = matrix([mpf('0.5'), mpf('50.0'), mpf('0.0')])

radius = matrix([mpf('0.0'), mpf('0.95'), mpf('1.05'), mpf('1.5')])

impedance_condition = True
beta_impedance = mpf('0')

freq = mpf('0.5')
omega = 2.0*pi*freq

first_order_abc = True

kwave = epsilon.copy()
Z0 = epsilon.copy()
for i in range(epsilon.rows):
    eps = epsilon[i]
    eps += 1j*sigma[i]/omega
    Z0[i] = sqrt(mu[i]/eps)
    kwave[i] = omega*sqrt(mu[i]*eps)

def GetXiPrime(z, n):
    jn = besselj(n+0.5, z)
    yn = bessely(n+0.5, z)
    hn = jn + 1j*yn
    
    if (n == 0):
        jn_next = besselj(n+1.5, z)
        yn_next = bessely(n+1.5, z)
        hn_next = jn_next + 1j*yn_next

        jn_prime = -jn_next + 0.5*jn/z
        hn_prime = -hn_next + 0.5*hn/z
    else:
        jn_prev = besselj(n-0.5, z)
        yn_prev = bessely(n-0.5, z)
        hn_prev = jn_prev + 1j*yn_prev
        
        jn_next = besselj(n+1.5, z)
        yn_next = bessely(n+1.5, z)
        hn_next = jn_next + 1j*yn_next
        
        jn_prime = 0.5*(jn_prev - jn_next)
        hn_prime = 0.5*(hn_prev - hn_next)
        
    jn = sqrt(0.5*pi/z) * jn
    psi_n = z*jn
    hn = sqrt(0.5*pi/z) * hn
    xi_n = z*hn
    psi_prime = 0.5*jn + sqrt(0.5*pi*z) * jn_prime
    xi_prime = 0.5*hn + sqrt(0.5*pi*z) * hn_prime
    jn_prime = -0.5*jn/z + sqrt(0.5*pi/z) * jn_prime
    hn_prime = -0.5*hn/z + sqrt(0.5*pi/z) * hn_prime
    
    hn2 = -hn + 2.0*jn
    hn2_prime = -hn_prime + 2.0*jn_prime
    xi2_n = -xi_n + 2.0*psi_n
    xi2_prime = -xi_prime + 2.0*psi_prime
    
    return xi_n, xi2_n, xi_prime, xi2_prime, psi_n, psi_prime

# boucle sur tous les modes
N = epsilon.rows
alpha = matrix(N_max, N)
gamma = matrix(N_max, N)
beta = matrix(N_max, N)
delta = matrix(N_max, N)
for n in range(1, N_max):
    sys_alpha = zeros(2*N)
    sys_beta = zeros(2*N)
    rhs_alpha = zeros(2*N, 1)
    rhs_beta = zeros(2*N, 1)
    
    # premiere equation donnee avec la condition pour r = r0
    if (radius[0] == mpf(0)):
        # cas d'une sphere dielectrique
        # alpha_n = gamma_n, beta_n = delta_n
        sys_alpha[0, 0] = 1.0
        sys_alpha[0, 1] = -1.0
        sys_beta[0, 0] = 1.0
        sys_beta[0, 1] = -1.0
    else:
        xi_n, xi2_n, xi_prime, xi2_prime, psi_n, psi_prime = GetXiPrime(kwave[0]*radius[0], n)
        if (impedance_condition):
            coef = omega / Z0[0]
            sys_alpha[0, 0] = beta_impedance*xi_prime + coef*xi_n
            sys_alpha[0, 1] = beta_impedance*xi2_prime + coef*xi2_n
            if (N == 1):
                rhs_alpha[0] = -beta_impedance*psi_prime - coef*psi_n
                
            sys_beta[0, 0] = -coef*xi_prime + beta_impedance*xi_n
            sys_beta[0, 1] = -coef*xi2_prime + beta_impedance*xi2_n
            if (N == 1):
                rhs_beta[0] = coef*psi_prime - beta_impedance*psi_n
        else:
            # condition de Dirichlet
            sys_alpha[0, 0] = 1j*xi_prime
            sys_alpha[0, 1] = 1j*xi2_prime
            if (N == 1):
                rhs_alpha[0] = -1j*psi_prime
                
            sys_beta[0, 0] = xi_n
            sys_beta[0, 1] = xi2_n
            if (N == 1):
                rhs_beta[0] = -psi_n
        
    # boucle sur les interfaces
    num = 1
    for i in range(1, radius.rows-1):
        xi_n, xi2_n, xi_prime, xi2_prime, psi_n, psi_prime = GetXiPrime(kwave[i-1]*radius[i], n)
        xi_next, xi2_next, xi_prime_next, xi2_prime_next, psi_next, psi_prime_next = GetXiPrime(kwave[i]*radius[i], n)

        # continuite de la trace tangentielle du champ electrique
        sys_alpha[num, 2*i-2] = xi_prime / kwave[i-1]
        sys_alpha[num, 2*i-1] = xi2_prime / kwave[i-1]
        sys_alpha[num, 2*i] = -xi_prime_next / kwave[i]
        sys_alpha[num, 2*i+1] = -xi2_prime_next / kwave[i]
        if (i == N-1):
            rhs_alpha[num,0] = psi_prime_next / kwave[i]
        
        sys_beta[num, 2*i-2] = xi_n / kwave[i-1]
        sys_beta[num, 2*i-1] = xi2_n / kwave[i-1]
        sys_beta[num, 2*i] = -xi_next / kwave[i]
        sys_beta[num, 2*i+1] = -xi2_next / kwave[i];
        if (i == N-1):
            rhs_beta[num,0] = psi_next / kwave[i]
          
        num += 1
          
        # continuity of the tangential trace of magnetic field
        sys_alpha[num, 2*i-2] = xi_n / (kwave[i-1]*Z0[i-1])
        sys_alpha[num, 2*i-1] = xi2_n / (kwave[i-1]*Z0[i-1])
        sys_alpha[num, 2*i] = -xi_next / (kwave[i]*Z0[i])
        sys_alpha[num, 2*i+1] = -xi2_next / (kwave[i]*Z0[i])
        if (i == N-1):
            rhs_alpha[num,0] = psi_next / (kwave[i]*Z0[i])

        sys_beta[num, 2*i-2] = xi_prime / (kwave[i-1]*Z0[i-1])
        sys_beta[num, 2*i-1] = xi2_prime / (kwave[i-1]*Z0[i-1])
        sys_beta[num, 2*i] = -xi_prime_next / (kwave[i]*Z0[i])
        sys_beta[num, 2*i+1] = -xi2_prime_next / (kwave[i]*Z0[i])
        if (i == N-1):
            rhs_beta[num,0] = psi_prime_next / (kwave[i]*Z0[i])
            
        num += 1
        
    # last condition due to Sommerfeld condition (at finite distance or not)
    if (first_order_abc):
        xi_n, xi2_n, xi_prime, xi2_prime, psi_n, psi_prime = GetXiPrime(kwave[N-1]*radius[N], n)
        sys_alpha[num, 2*N-2] = xi_prime - 1j*xi_n
        sys_alpha[num, 2*N-1] = xi2_prime - 1j*xi2_n
        
        sys_beta[num, 2*N-2] = xi_prime - 1j*xi_n
        sys_beta[num, 2*N-1] = xi2_prime - 1j*xi2_n
    else:
        # exact Sommerfeld condition => no component in Xi_n^{(2)}
        sys_alpha[num, 2*N-1] = 1.0
        sys_beta[num, 2*N-1] = 1.0
        
    #print "sys_alpha", sys_alpha
    #print "sys_beta", sys_beta
    #print "rhs_alpha", rhs_alpha
    #print "rhs_beta", rhs_beta
    
    #sol_alpha = sys_alpha**-1*rhs_alpha
    #sol_beta = sys_beta**-1*rhs_beta
    sol_alpha = lu_solve(sys_alpha, rhs_alpha)
    sol_beta = lu_solve(sys_beta, rhs_beta)
    
    print("Norm rhs = ", float(norm(rhs_alpha)), float(norm(rhs_beta)))
    #print "n = ", n
    #print "sol_alpha", sol_alpha
    #print "sol_beta", sol_beta
    for i in range(N):
        alpha[n, i] = sol_alpha[2*i]
        gamma[n, i] = sol_alpha[2*i+1]
        beta[n, i] = sol_beta[2*i]
        delta[n, i] = sol_beta[2*i+1]


def WriteText(A, nom):
    file_out = open(nom, 'w')
    for i in range(A.rows):
        for j in range(A.cols):
            file_out.write("(" + str(A[i, j].real) + "," + str(A[i, j].imag) + ") ")
        
        file_out.write("\n")
    
    file_out.close()

# on imprime les coefs dans un fichier
WriteText(alpha, "alpha.dat")
WriteText(gamma, "gamma.dat")
WriteText(beta, "beta.dat")
WriteText(delta, "delta.dat")
