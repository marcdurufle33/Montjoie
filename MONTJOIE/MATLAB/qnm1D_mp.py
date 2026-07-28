from mpmath import *

mp.dps = 100;
threshold = 1e-40

def qnm_circle_exact_mode(n, kmax, rho, mu, kimag = mpf('-0.5')):
    k0 = arange(0.3, kmax, 0.5)
    omega = []
    def my_f(x):
        k = x[0] + 1j*x[1]
        ki = sqrt(rho/mu)*k
        jnm1 = besselj(n-1, ki)
        jn = besselj(n, ki)
        jnp1 = besselj(n+1, ki)
        hnm1 = besselj(n-1, k) + 1j*bessely(n-1, k)
        hn = besselj(n, k) + 1j*bessely(n, k)
        hnp1 = besselj(n+1, k) + 1j*bessely(n+1, k)
        w = k*(hnm1 - hnp1)*jn - mu*ki*(jnm1 - jnp1)*hn
        return matrix([w.real, w.imag])
    
    for i in range(len(k0)):
        x = matrix([k0[i], kimag])
        test_loop = True
        num = 0; h = mpf('1e-12')
        #print "kinit = ", x
        while (test_loop):
            fn = my_f(x)
            dfx = (my_f(matrix([x[0]+h, x[1]])) - my_f(matrix([x[0]-h, x[1]]))) / (2.0*h)
            dfy = (my_f(matrix([x[0], x[1]+h])) - my_f(matrix([x[0], x[1]-h]))) / (2.0*h)
            dfn = matrix([[dfx[0], dfy[0]], [dfx[1], dfy[1]]])
            x = x - lu_solve(dfn, fn)
            eps = norm(my_f(x))
            num = num+1
            print "norme de eps = ", eps
            if ((num > 30) or (eps < threshold)):
                test_loop = False

        if (eps < threshold):
            k = x[0] + 1j*x[1]
            #print "solution = ", k
            sol_already = False
            for j in range(len(omega)):
                if (abs(omega[j] -k) < threshold):
                    sol_already = True
                    
            if (not sol_already):
                omega.append(k)
    
    return omega

def qnm_circle_exact(kmax, rho, mu, inc = 1, omegaI = mpf('2.0')):
    Nmax = 1000
    omega = []
    mode = []
    n = 0
    while (n < Nmax):
        omega_tmp = sort(qnm_circle_exact_mode(n, kmax, rho, mu))
        nb_elt = 0
        for i in range(len(omega_tmp)):
            if ((omega_tmp[i].real >= 0) and (omega_tmp[i].real <= kmax) and (omega_tmp[i].imag >= -omegaI)):
                nb_elt = nb_elt+1
                omega.append(omega_tmp[i])

        if (nb_elt == 0):
            Nmax = n

        n = n + inc

    return sort(omega)
