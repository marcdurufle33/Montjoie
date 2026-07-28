from pylab import *

# pour recuperer omega_num et calcul d'erreur
#Lr2 = loadtxt('../eigenval_real.dat')
#Li2 = loadtxt('../eigenval_imag.dat')
#num = find(Li2 >= 0)
#permut = argsort(Li2[num] -1j*Lr2[num])
#omega_num = Li2[num[permut]] - 1j*Lr2[num[permut]]
#for i in range(45):
#        X, Y, Z, coor, E = loadND('../ModeE' + str(num[permut[i]]) + '.dat')
#        erreur1[i] = norm(Eexact[:, i]/Eexact[750, i] - E/E[750])/norm(E/E[750])

# pour projeter la solution sur les modes numeriques
#numb = find(Li2 <= 0); permut_b = argsort(-Li2[numb] +1j*Lr2[numb])
#alpha_num = load1D('../alpha1.dat')[num[permut]]; beta_num = load1D('../alpha1.dat')[numb[permut_b]];
#f = zeros(len(X)) + 1j*zeros(len(X))
#for i in range(45):
#        X, Y, Z, coor, E = loadND('../ModeE' + str(num[permut[i]]) + '.dat'); f = f + alpha_num[i]*E
#for i in range(1, 45):
#        X, Y, Z, coor, E = loadND('../ModeE' + str(numb[permut_b[i]]) + '.dat'); f = f + beta_num[i]*E


#freq = linspace(0.01, 2.0, 3)
#X, Y, Z, coor, sol_num = loadND('../Un1.dat')
#omega = 2*pi*freq[1]
#ol_qnm = qnm_proj_incident(omega, n0, N, X)


def qnm_exact_1D(n0, N, X = array([])):
    xsi = 0.5*log((n0+1)/(n0-1))
    omega_exact = arange(0, N)*pi/2/n0 - 1j*xsi/n0
    if (len(X) == 0):
        return omega_exact
    
    i0 = find(X == -1.0)[0]
    i1 = find(X == 1.0)[0]
    Eexact = zeros([len(X), N]) + 1j*zeros([len(X), N])
    for j in range(N):
        omegaj = omega_exact[j]
        if (j%2 == 0):
            Eexact[i0:i1, j] = cos(n0*omegaj*X[i0:i1])
            Eexact[i1:, j] = cos(omegaj*(X[i1:]-1.0)+n0*omegaj) + (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0, j] = cos(omegaj*(X[0:i0]+1.0)-n0*omegaj) - (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))
        else:
            Eexact[i0:i1, j] = sin(n0*omegaj*X[i0:i1])
            Eexact[i1:, j] = sin(omegaj*(X[i1:]-1.0)+n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0, j] = sin(omegaj*(X[0:i0]+1.0)-n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))
    
    return omega_exact, Eexact

from scipy.integrate import *

def qnm_proj_incident(omega, n0, N, X = array([]), f = None):
    alpha = zeros(N) + 1j*zeros(N)
    beta = zeros(N) + 1j*zeros(N)
    omega_exact = qnm_exact_1D(n0, N)
        
    for j in range(N):
        omegaj = omega_exact[j]
        if (j%2 == 0):
            if (f != None):
                alpha[j] = -1j/(omega-omegaj)*(quad(lambda x : real(cos(n0*omegaj*x)*f(x)), -1.0, 1.0)[0] + 1j*quad(lambda x : imag(cos(n0*omegaj*x)*f(x)), -1.0, 1.0)[0])
                beta[j] = -1j/(omega+conj(omegaj))*(quad(lambda x : real(cos(-n0*conj(omegaj)*x)*f(x)), -1.0, 1.0)[0] + 1j*quad(lambda x : imag(cos(-n0*conj(omegaj)*x)*f(x)), -1.0, 1.0)[0])
            else:
                alpha[j] = (n0*n0-1.0)*omega/(omega-omegaj)*(sin(n0*omegaj+omega)/(n0*omegaj+omega) + sin(-n0*omegaj+omega)/(-n0*omegaj+omega))
                beta[j] = (n0*n0-1.0)*omega/(omega+conj(omegaj))*(sin(-n0*conj(omegaj)+omega)/(-n0*conj(omegaj)+omega) + sin(n0*conj(omegaj)+omega)/(n0*conj(omegaj)+omega))
        else:
            if (f != None):
                alpha[j] = -1j/(omega-omegaj)*(quad(lambda x : real(sin(n0*omegaj*x)*f(x)), -1.0, 1.0)[0] + 1j*quad(lambda x : imag(sin(n0*omegaj*x)*f(x)), -1.0, 1.0)[0])
                beta[j] = -1j/(omega+conj(omegaj))*(quad(lambda x : real(sin(-n0*conj(omegaj)*x)*f(x)), -1.0, 1.0)[0] + 1j*quad(lambda x : imag(sin(-n0*conj(omegaj)*x)*f(x)), -1.0, 1.0)[0])
            else:
                alpha[j] = (n0*n0-1.0)*omega/(omega-omegaj)*(sin(n0*omegaj+omega)/(1j*n0*omegaj+1j*omega) - sin(-n0*omegaj+omega)/(-1j*n0*omegaj+1j*omega))
                beta[j] = (n0*n0-1.0)*omega/(omega+conj(omegaj))*(sin(-n0*conj(omegaj)+omega)/(-1j*n0*conj(omegaj)+1j*omega) - sin(n0*conj(omegaj)+omega)/(1j*n0*conj(omegaj)+1j*omega))

    # sans normaliser les modes QNM (c'est directement sin et cos)
    alpha *= -1.0/(n0*sqrt(2.0))**2
    beta *= -1.0/(n0*sqrt(2.0))**2

    if (len(X) == 0):
        return alpha, beta

    sol = zeros(len(X)) + 1j*zeros(len(X))
    Eexact = zeros(len(X)) + 1j*zeros(len(X))
    i0 = find(X == -1.0)[0]
    i1 = find(X == 1.0)[0]
    for j in range(N):
        omegaj = omega_exact[j]
        if (j%2 == 0):
            Eexact[i0:i1] = cos(n0*omegaj*X[i0:i1])
            Eexact[i1:] = cos(omegaj*(X[i1:]-1.0)+n0*omegaj) + (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0] = cos(omegaj*(X[0:i0]+1.0)-n0*omegaj) - (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))
        else:
            Eexact[i0:i1] = sin(n0*omegaj*X[i0:i1])
            Eexact[i1:] = sin(omegaj*(X[i1:]-1.0)+n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0] = sin(omegaj*(X[0:i0]+1.0)-n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))

        sol += alpha[j]*Eexact

    for j in range(1, N):
        omegaj = -conj(omega_exact[j])
        if (j%2 == 0):
            Eexact[i0:i1] = cos(n0*omegaj*X[i0:i1])
            Eexact[i1:] = cos(omegaj*(X[i1:]-1.0)+n0*omegaj) + (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0] = cos(omegaj*(X[0:i0]+1.0)-n0*omegaj) - (1.0-n0)*sin(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))
        else:
            Eexact[i0:i1] = sin(n0*omegaj*X[i0:i1])
            Eexact[i1:] = sin(omegaj*(X[i1:]-1.0)+n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[i1:]-1.0))
            Eexact[0:i0] = sin(omegaj*(X[0:i0]+1.0)-n0*omegaj) + (n0-1.0)*cos(n0*omegaj)*sin(omegaj*(X[0:i0]+1.0))
    
        sol += beta[j]*Eexact

    return sol

def layered_incident_wave(kx, ky_above, a, b, d):
    """
    Usage ky, A, B = layered_incident_wave(kx, ky_above, a, b, d)
    warning : a, b and d are ordered with increasing y
              kx, ky wave vector for the above part of the plane
    """
    signe = sign(ky_above); # signe des ky
    N = len(d)-1
    omega = sqrt(a[N+1] / b[N+1]*(kx**2 + ky_above**2))
    ky = zeros(N+2)+1j*zeros(N+2)
    for i in range(N+2):
        epsilon = b[i] / a[i]
        ky[i] = signe * sqrt(omega*omega*epsilon - kx*kx)   
    # the unknowns are ordered as follows
    # A : 2*i 
    # B : 2*i+1
    # A(0) and B(0) correspond to the lower infinite medium 
    # A(i) and B(i) correspond to the layer between d(i-1) and d(i)
    # A(N+1) and B(N+1) correspond to the upper infinite medium
    
    # the plane wave is supposed to come from above, we set 
    # A(N+1) = 1
    # and no reflected wave from low infty 
    # B(0) = 0

    M = zeros([2*N+4, 2*N+4]) + 1j*zeros([2*N+4, 2*N+4])
    for i in range(N+1):
        M[2*i+1,2*i]       = exp( 1j * ky[i] * d[i] ); # *A(i)
        M[2*i+1,2*i+1]     = exp( -1j * ky[i] * d[i] ); # *B(i)

        M[2*i+1,2*i+2]     = - exp(1j * ky[i+1] * d[i] ); # *A(i+1)
        M[2*i+1, 2*i+3]    = - exp(-1j * ky[i+1] * d[i] );	# *B(i+1)
	
        M[2*i+2, 2*i]    = 1j * ky[i] * a[i] * M[2*i+1,2*i] ; 
        M[2*i+2, 2*i+1]  = -1j * ky[i] * a[i] * M[2*i+1,2*i+1] ;

        M[2*i+2, 2*i+2]  = 1j * ky[i+1]  *a[i+1] * M[2*i+1,2*i+2] ;
        M[2*i+2, 2*i+3]  = -1j*ky[i+1]  *a[i+1] * M[2*i+1, 2*i+3] ;
      
    M[0,0] = 0.0; 
    M[0,1] = 1.0; # * B(0)

    M[2*N+3, 2*N+2] = 1.0; # *A(N+1)

    rhs = zeros(2*N+4)
    rhs[2*N+3] = 1.0;
    
    sol = solve(M, rhs)
    A = sol[0::2]
    B = sol[1::2]
    return ky, A, B

from scipy.special import *

def qnm_circle_exact_mode(n, kmax, rho, mu, kimag = -0.5):
    k0 = arange(0.3, kmax, 0.5)
    omega = []
    def my_f(x):
        k = x[0] + 1j*x[1]
        ki = sqrt(rho/mu)*k
        w = k*(hankel1(n-1, k) - hankel1(n+1, k))*jn(n, ki) - mu*ki*(jn(n-1, ki) - jn(n+1, ki))*hankel1(n, k)
        return array([real(w), imag(w)])
    
    for i in range(len(k0)):
        x = array([k0[i], kimag])
        test_loop = True
        num = 0; h = 1e-6
        #print "kinit = ", x
        while (test_loop):
            fn = my_f(x)
            dfx = (my_f(array([x[0]+h, x[1]])) - my_f(array([x[0]-h, x[1]]))) / (2.0*h)
            dfy = (my_f(array([x[0], x[1]+h])) - my_f(array([x[0], x[1]-h]))) / (2.0*h)
            dfn = array([[dfx[0], dfy[0]], [dfx[1], dfy[1]]])
            x = x - solve(dfn, fn)
            eps = norm(my_f(x))
            num = num+1
            #print "norme de eps = ", eps
            if ((num > 20) or (eps < 1e-12)):
                test_loop = False

        if (eps < 1e-12):
            k = x[0] + 1j*x[1]
            #print "solution = ", k
            sol_already = False
            for j in range(len(omega)):
                if (abs(omega[j] -k) < 1e-10):
                    sol_already = True
                    
            if (not sol_already):
                omega.append(k)
    
    return omega

def qnm_circle_exact_mode_drude(n, kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                                polaTE=True, kimag = -0.5, k0 = 0.3, k0_tab = []):
    if (len(k0_tab) == 0):
        k0 = arange(k0, kmax, 0.5)
    else:
        k0 = k0_tab.copy()
    
    omega = []
    def my_f(x):
        k = x[0] + 1j*x[1]
        rho = eps_inf
        mu = 1.0
        for ii in range(len(gamma)):
            rho -= (eps_omega_p2[ii] - 1j*k*eps_sigma[ii])/(k**2 - omega_02[ii] + 1j*k*gamma[ii])
        
        if (polaTE != True):
            mu = 1.0/rho
            rho = 1.0
        
        ki = sqrt(rho/mu)*k
        w = k*(hankel1(n-1, k) - hankel1(n+1, k))*jn(n, ki) - mu*ki*(jn(n-1, ki) - jn(n+1, ki))*hankel1(n, k)
        return array([real(w), imag(w)])
    
    for i in range(len(k0)):
        if (len(k0_tab) == 0):
            x = array([k0[i], kimag])
        else:
            x = array([real(k0[i]), imag(k0[i])])
        
        test_loop = True
        num = 0; h = 1e-6
        #print "kinit = ", x
        while (test_loop):
            fn = my_f(x)
            dfx = (my_f(array([x[0]+h, x[1]])) - my_f(array([x[0]-h, x[1]]))) / (2.0*h)
            dfy = (my_f(array([x[0], x[1]+h])) - my_f(array([x[0], x[1]-h]))) / (2.0*h)
            dfn = array([[dfx[0], dfy[0]], [dfx[1], dfy[1]]])
            if (isnan(cond(dfn))):
                print(x, n, dfn, fn)
                line = input("Press o to keep the mode")
            if (cond(dfn) > 1e16):
                break;
            
            x = x - solve(dfn, fn)
            eps = norm(my_f(x))
            num = num+1
            #print "norme de eps = ", eps
            if ((num > 50) or (eps < 1e-12)):
                test_loop = False

        if (eps < 1e-10):
            k = x[0] + 1j*x[1]
            #print "solution = ", k
            sol_already = False
            for j in range(len(omega)):
                if (abs(omega[j] -k) < 1e-10):
                    sol_already = True
                    
            if (not sol_already):
                omega.append(k)
    
    return omega


def qnm_circle_exact(kmax, rho, mu, inc = 1, omegaI = 2.0, X = array([]), Y = array([]), kimag = -0.5, k0_tab = []):
    Nmax = 100
    omega = []
    mode = []
    if (len(X) != 0):
        R = sqrt(X**2 + Y**2) + 1e-15
        Teta = arccos(X/R)
    n = 0
    while (n < Nmax):
        omega_tmp = sort(qnm_circle_exact_mode(n, kmax, rho, mu, kimag))
        nb_elt = 0
        for i in range(len(omega_tmp)):
            if ((real(omega_tmp[i]) >= 0) and (real(omega_tmp[i] <= kmax)) and (imag(omega_tmp[i]) >= -omegaI)):
                nb_elt = nb_elt+1
                omega.append(omega_tmp[i])
                if (len(X) != 0):
                    k = omega_tmp[i]
                    ki = sqrt(rho/mu)*k
                    beta = jn(n, ki) / hankel1(n, k)
                    mode_tmp = (R < 1)*jn(n, ki*R)*cos(n*Teta) + (R >= 1)*beta*hankel1(n, k*R)*cos(n*Teta)
                    mode.append(mode_tmp)


        if (nb_elt == 0):
            Nmax = n

        n = n + inc

    omega = array(omega)
    if (len(X) != 0):
        num=argsort(omega)
        omega2=omega[num]
        mode2=[]
        for i in range(len(num)):
            mode2.append(mode[num[i]])
        return omega2, mode2

    return sort(omega)


def qnm_circle_exact_drude(kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                           polaTE=True, inc = 1, omegaI = 2.0, X = array([]), Y = array([]), kimag = -0.5, k0 = 0.3, k0_tab = []):
    Nmax = 100
    omega = []
    mode = []
    if (len(X) != 0):
        R = sqrt(X**2 + Y**2) + 1e-15
        Teta = arccos(X/R)
    
    n = 0
    while (n < Nmax):
        omega_tmp = sort(qnm_circle_exact_mode_drude(n, kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma, polaTE, kimag, k0, k0_tab))
        nb_elt = 0
        for i in range(len(omega_tmp)):
            if ((real(omega_tmp[i]) >= 0) and (real(omega_tmp[i] <= kmax)) and (imag(omega_tmp[i]) >= -omegaI)):
                nb_elt = nb_elt+1
                omega.append(omega_tmp[i])
                if (len(X) != 0):
                    k = omega_tmp[i]
                    rho = epsiloninf
                    mu = 1.0
                    for ii in range(len(gamma)):
                        rho -= (eps_omega_p2[ii] - 1j*k*eps_sigma[ii])/(k**2 - omega_02[ii] + 1j*k*gamma[ii])
                    
                    if (polaTE != True):
                        mu = 1.0/rho
                        rho = 1.0
                    
                    ki = sqrt(rho/mu)*k
                    beta = jn(n, ki) / hankel1(n, k)
                    mode_tmp = (R < 1)*jn(n, ki*R)*cos(n*Teta) + (R >= 1)*beta*hankel1(n, k*R)*cos(n*Teta)
                    mode.append(mode_tmp)


        if (nb_elt == 0):
            Nmax = n

        n = n + inc

    omega = array(omega)
    if (len(X) != 0):
        num = argsort(omega)
        omega2 = omega[num]
        mode2 = []
        for i in range(len(num)):
            mode2.append(mode[num[i]])
        return omega2, mode2

    return sort(omega)


def find_qnm_circle_exact_drude(k, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                                polaTE = True, inc = 1, Nmax = 10, threshold=1e-3):
    if (type(k) == ndarray):
        tab_om = []
        for i in range(len(k)):
            if ((imag(k[i]) > -1.0) and (real(k[i])<10.0) and (abs(k[i])>1e-3)):
                print(i)
                om = find_qnm_circle_exact_drude(k[i], eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                                                 polaTE, inc, Nmax, threshold)

                if (om != 0.0):
                    tab_om.append(om)
        
        num = argsort(tab_om)
        res = [tab_om[num[0]]]
        om_prev = tab_om[num[0]]
        for i in range(1, len(tab_om)):
            if (tab_om[num[i]] != om_prev):
                res.append(tab_om[num[i]])
                om_prev = tab_om[num[i]]

        return array(res)
    
    n = 0
    while (n < Nmax):
        omega_tmp = qnm_circle_exact_mode_drude(n, 1.0, eps_inf, eps_omega_p2, omega_02, gamma,
                                                eps_sigma, polaTE, 0.0, 0.0, array([k]))

        if (len(omega_tmp) >= 1):
            err = abs(array(omega_tmp) - k)
            num = argsort(err)
            if (err[num[0]] < threshold):
                #print("Root found = ", omega_tmp[num[0]], n)
                return omega_tmp[num[0]]
        
        n = n + inc

    return 0.0


def min_QNM(omega,L):
    MIN=zeros(len(omega))
    pos=zeros(len(omega), dtype='int64')
    for j in range(len(omega)):
        pop=10000000
        it=0
        for i in range(len(L)):
            pop1=abs(omega[j]-L[i])
            if pop1<pop:
                pop=pop1
                it=i
        MIN[j]=pop
        pos[j]=it
    return MIN, pos
    
def writemyaxes(a):
    xlabel(r'$Re(\omega_m)$')
    ylabel(r'$Im(\omega_m)$')
    return


from visuND import *

def calcul_erreur_mode(pos,mode, racine ='../ModeU'):
    error=zeros(len(pos))
    for i in range(len(pos)):
        if pos[i]<1000:
            Nameoffile=racine+"0%dG0_U0.dat"%(pos[i])
            if pos[i]<100:
                Nameoffile=racine+"00%dG0_U0.dat"%(pos[i])
                if pos[i]<10: 
                    Nameoffile=racine + "000%dG0_U0.dat"%(pos[i])
        [X,Y,Z,coor,V]=loadND(Nameoffile)
        error[i]=erreurMode(conj(V), mode[i])
    return error

def toolazy(omega,mode):
    L = load1D("../Lambda.dat")
    L=-imag(L)-1j*real(L)
    figure(1)   
    plot(real(L),imag(L),'.')
    MIN, pos=min_QNM(omega,L)
    figure(2)
    semilogy(MIN)
    figure(3)
    error=calcul_erreur_mode(pos,mode)
    semilogy(error)
    return 

def animate_modes(n,racine ='../EigenvalCarre',end ='.dat',sauve_film = None,dt=0):
    sauve_fig = figure()
    ylim(-2,0)
    xlim(0,22)
    col=['b','r','g','y','c','m','k']
    symb=['o','s','v','^','<','>','d']
    SHIFT = load1D('../liste_shift.dat')
    SHIFT = -1j*SHIFT
    for i in range(n):
        Nameoffile_r=racine+"%d_real"%i+end
        Nameoffile_i=racine+"%d_imag"%i+end
        Lr = loadtxt(Nameoffile_r)
        Li = loadtxt(Nameoffile_i)
        w=int(floor(float(i)/len(col))%len(symb))
        p=int(i%len(col))
        stul=col[p]+symb[w]
        plot(-Li, -Lr, stul,fillstyle = 'none')
        plot(real(SHIFT[i]),imag(SHIFT[i]),'+k',fillstyle='full',markersize=10)
        draw()
        if (sauve_film != None):
            savefig('sauve_film'+EntierToString(i)+".png")
        os.system('sleep ' +str(dt))
    return

def qnm_sphere_exact_mode(n, kmax, eps, mu, kimag = -0.5): #fonction qui calcule les frequences propres des modes d'une sphere analytiquement 
    k0 = arange(0.3, kmax, 0.5)   
    omega = []
    def my_f(x):
        k = x[0] + 1j*x[1]
        ki = sqrt(eps*mu)*k
        Z0 = sqrt(mu/eps)
        jbessel = sqrt(pi/(2*ki))*jn(n+0.5,ki) # petit jn (fonction de Bessel spherique du nieme ordre)
        djbessel = sqrt(pi/(2*ki))*(jn(n-0.5,ki) - jn(n+1.5,ki))/2 - sqrt(pi/(2*ki))/(2*ki)*jn(n+0.5,ki) # derivee de petit jn (jbessel)
        Psy = ki*jbessel
        dPsy = ki*djbessel + jbessel # derivee de dPsy
        hankelsph = sqrt(pi/(2*k))*hankel1(n+0.5,k) # fonction de hankel spherique du premier type (petit hn)
        dhankelsph = sqrt(pi/(2*k))*(hankel1(n-0.5,k) - hankel1(n+1.5,k))/2 - sqrt(pi/(2*k))/(2*k)*hankel1(n+0.5,k) # derivee de petit hn         
        Ksi = k*hankelsph
        dKsi = k*dhankelsph + hankelsph # derivee de dKsi      
        #w = k*(hankel1(n-1, k) - hankel1(n+1, k))*jn(n, ki) - mu*ki*(jn(n-1, ki) - jn(n+1, ki))*hankel1(n, k) 
        w =  Z0*dPsy*Ksi - dKsi*Psy # Wronskien (1) dont on cherche les 0 afin de trouver les modes
        w# =  Z0*dKsi*Psy - dPsy*Ksi # Wronskien (2) dont on cherche les 0 afin de trouver les modes
        return array([real(w), imag(w)])
    
    for i in range(len(k0)):
        x = array([k0[i], kimag])
        test_loop = True
        num = 0; h = 1e-6
        #print "kinit = ", x
        while (test_loop):
            fn = my_f(x)
            dfx = (my_f(array([x[0]+h, x[1]])) - my_f(array([x[0]-h, x[1]]))) / (2.0*h)
            dfy = (my_f(array([x[0], x[1]+h])) - my_f(array([x[0], x[1]-h]))) / (2.0*h)
            dfn = array([[dfx[0], dfy[0]], [dfx[1], dfy[1]]])
            x = x - solve(dfn, fn)
            eps2 = norm(my_f(x))
            num = num+1
            print("norme de eps2 = ", eps2)
            if ((num > 20) or (eps2 < 1e-12)):
                test_loop = False

        if (eps2 < 1e-12):
            k = x[0] + 1j*x[1]
            #print "solution = ", k
            sol_already = False
            for j in range(len(omega)):
                if (abs(omega[j] -k) < 1e-10):
                    sol_already = True
                    
            if (not sol_already):
                omega.append(k)
    
    return omega




def qnm_sphere_exact(kmax, eps, mu, inc = 1, omegaI = 2.0, X = array([]), Y = array([]), kimag = -0.5):
    Nmax = 1000
    omega = []
    mode = []
    if (len(X) != 0):
        R = sqrt(X**2 + Y**2) + 1e-15
        Teta = arccos(X/R)
    n = 0
    while (n < Nmax):
        omega_tmp = sort(qnm_circle_exact_mode(n, kmax, eps, mu, kimag))
        nb_elt = 0
        for i in range(len(omega_tmp)):
            if ((real(omega_tmp[i]) >= 0) and (real(omega_tmp[i] <= kmax)) and (imag(omega_tmp[i]) >= -omegaI)):
                nb_elt = nb_elt+1
                omega.append(omega_tmp[i])
                if (len(X) != 0):
                    k = omega_tmp[i]
                    ki = sqrt(rho/mu)*k
                    beta = jn(n, ki) / hankel1(n, k)
                    mode_tmp = (R < 1)*jn(n, ki*R)*cos(n*Teta) + (R >= 1)*beta*hankel1(n, k*R)*cos(n*Teta)
                    mode.append(mode_tmp)


        if (nb_elt == 0):
            Nmax = n

        n = n + inc

    omega = array(omega)
    if (len(X) != 0):
        num=argsort(omega)
        omega2=omega[num]
        mode2=[]
        for i in range(len(num)):
            mode2.append(mode[num[i]])
        return omega2, mode2

    return sort(omega)

#fonction qui calcule les frequences propres des modes d'une sphere analytiquement 
def qnm_sphere_exact_mode_drude(n, kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma, mu, numw = 1, kimag = -0.5, k0_tab = [], display=False):
    if (len(k0_tab) == 0):
        k0 = arange(0.3, kmax, 0.5)
    else:
        k0 = k0_tab.copy()
    
    omega = []
    def my_f(x):
        k = x[0] + 1j*x[1]
        eps = eps_inf
        for ii in range(len(gamma)):
            eps -= (eps_omega_p2[ii] - 1j*k*eps_sigma[ii])/(k**2 - omega_02[ii] + 1j*k*gamma[ii])

        ki = sqrt(eps*mu)*k
        Z0 = sqrt(mu/eps)
        jbessel = sqrt(pi/(2*ki))*jn(n+0.5,ki) # petit jn (fonction de Bessel spherique du nieme ordre)
        djbessel = sqrt(pi/(2*ki))*(jn(n-0.5,ki) - jn(n+1.5,ki))/2 - sqrt(pi/(2*ki))/(2*ki)*jn(n+0.5,ki) # derivee de petit jn (jbessel)
        Psy = ki*jbessel
        dPsy = ki*djbessel + jbessel # derivee de dPsy
        hankelsph = sqrt(pi/(2*k))*hankel1(n+0.5,k) # fonction de hankel spherique du premier type (petit hn)
        dhankelsph = sqrt(pi/(2*k))*(hankel1(n-0.5,k) - hankel1(n+1.5,k))/2 - sqrt(pi/(2*k))/(2*k)*hankel1(n+0.5,k) # derivee de petit hn         
        Ksi = k*hankelsph
        dKsi = k*dhankelsph + hankelsph # derivee de dKsi      
        #w = k*(hankel1(n-1, k) - hankel1(n+1, k))*jn(n, ki) - mu*ki*(jn(n-1, ki) - jn(n+1, ki))*hankel1(n, k) 
        if (numw == 1):
            w =  Z0*dPsy*Ksi - dKsi*Psy # Wronskien (1) dont on cherche les 0 afin de trouver les modes
        elif (numw ==0):
            w =  Z0*dKsi*Psy - dPsy*Ksi # Wronskien (2) dont on cherche les 0 afin de trouver les modes
        else:
            w = Z0*(jbessel*dhankelsph)-hankelsph*djbessel
        return array([real(w), imag(w)])
    
    for i in range(len(k0)):
        if (len(k0_tab) == 0):
            x = array([k0[i], kimag])
        else:
            x = array([real(k0[i]), imag(k0[i])])
        test_loop = True
        num = 0; h = 1e-6
        if (display):
            print("kinit = ", x, n, numw)
        while (test_loop):
            fn = my_f(x)
            dfx = (my_f(array([x[0]+h, x[1]])) - my_f(array([x[0]-h, x[1]]))) / (2.0*h)
            dfy = (my_f(array([x[0], x[1]+h])) - my_f(array([x[0], x[1]-h]))) / (2.0*h)
            dfn = array([[dfx[0], dfy[0]], [dfx[1], dfy[1]]])
            x = x - solve(dfn, fn)
            eps2 = norm(my_f(x))
            num = num+1
            if (display):
                print("norme de eps2 = ", eps2)
            if ((num > 30) or (eps2 < 1e-12)):
                test_loop = False

        if (eps2 < 1e-8):
            k = x[0] + 1j*x[1]
            if (display):
                print("solution = ", k)
            sol_already = False
            for j in range(len(omega)):
                if (abs(omega[j] -k) < 1e-10):
                    sol_already = True
                    
            if (not sol_already):
                omega.append(k)
    
    return omega

def qnm_sphere_exact_drude(kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                           n0 = 0, omegaI = 2.0, X = array([]), Y = array([]), kimag = -0.5, k0_tab = []):
    Nmax = 1000
    omega = []
    mode = []
    if (len(X) != 0):
        R = sqrt(X**2 + Y**2) + 1e-15
        Teta = arccos(X/R)
    n = n0
    while (n < Nmax):
        omega_tmp = sort(qnm_sphere_exact_mode_drude(n//2, kmax, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma, 1.0, n%2, kimag, k0_tab))
        #print("mode trouve", n, omega_tmp)
        nb_elt = 0
        for i in range(len(omega_tmp)):
            if ((real(omega_tmp[i]) >= 0) and (real(omega_tmp[i] <= kmax)) and (imag(omega_tmp[i]) >= -omegaI)):
                nb_elt = nb_elt+1
                omega.append(omega_tmp[i])
        if (nb_elt == 0):
            Nmax = n
        n = n + 1
    omega = array(omega)
    if (len(X) != 0):
        num=argsort(omega)
        omega2=omega[num]
        mode2=[]
        for i in range(len(num)):
            mode2.append(mode[num[i]])
        return omega2, mode2

    return sort(omega)

def find_qnm_sphere_exact_drude(k, eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                                n0 = 1, Nmax = 10, threshold=1e-3):
    if (type(k) == ndarray):
        tab_om = []
        for i in range(len(k)):
            if ((imag(k[i]) > -1.0) and (real(k[i])<10.0) and (abs(k[i])>1e-3)):
                print(i)
                om = find_qnm_sphere_exact_drude(k[i], eps_inf, eps_omega_p2, omega_02, gamma, eps_sigma,
                                                 n0, Nmax, threshold)

                if (om != 0.0):
                    tab_om.append(om)
        
        num = argsort(tab_om)
        res = [tab_om[num[0]]]
        om_prev = tab_om[num[0]]
        for i in range(1, len(tab_om)):
            if (tab_om[num[i]] != om_prev):
                res.append(tab_om[num[i]])
                om_prev = tab_om[num[i]]

        return array(res)
    
    n = n0
    while (n < Nmax):
        omega_tmp = qnm_sphere_exact_mode_drude(n//2, 1.0, eps_inf, eps_omega_p2, omega_02, gamma,
                                                eps_sigma, 1.0, n%2, 0.0, array([k]))

        if (len(omega_tmp) >= 1):
            err = abs(array(omega_tmp) - k)
            num = argsort(err)
            if (err[num[0]] < threshold):
                #print("Root found = ", omega_tmp[num[0]], n)
                return omega_tmp[num[0]]
        
        n = n + 1

    return 0.0
