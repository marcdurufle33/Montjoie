# inclut les fonctions necessaires au TP:
from pylab import *

def euler(f, tn, xn, h):
    x_next = xn + h*f(tn, xn)
    return x_next

def euler_implicite(f, tn, xn, h):
    x = xn.copy()
    gx = -h*f(tn+h, xn)
    m = len(xn)
    alpha = 1e-6
    while (norm(gx) >= 1e-12):
        ## Computing numerically the Jacobian
        A = eye(m)
        for l in range(m) :
            xm = x.copy(); xp = x.copy()
            xm[l] = xm[l]-alpha
            xp[l] = xp[l]+alpha
            A[:,l] -= h*(f(tn+h, xp) - f(tn+h, xm)) / (2*alpha)

        print("A = ", A)

        ## Newton Iterate
        x = x - solve(A, gx)
        gx = x-xn - h*f(tn+h, x)
        print(norm(gx))
    
    return x

def ode_solve(f, t0, x0, tf, h, schema):
    # nombre d'iterations en temps
    nb_iter = int(ceil((tf-t0) / h))
    h = (tf - t0) / nb_iter
    
    # on initialise les vecteurs tn et xn
    tn = zeros(nb_iter+1)
    m = 1
    if (type(x0) == ndarray):
        m = len(x0)
    
    xn = zeros([nb_iter+1, m])
    
    # ecrire ici votre code qui calcule la suite des xn avec le schema d'Euler
    xn[0,:] = x0; tn[0] = t0;
    for i in range(nb_iter):
        xn[i+1,:] = schema(f, tn[i], xn[i,:], h)
        tn[i+1] = t0 + (i+1)*h
    
    # a la fin de la fonction il faut retourner t_n et x_n
    return tn, xn

def ode_order(tf, h0, ratio, level, schema):
    f = lambda t, x : -x
    err = zeros(level)
    h = zeros(level) 
    for i in range(level):
        h[i] = h0 / ratio**i
        tn, xn = ode_solve(f, 0.0, 1.0, tf, h[i], schema)
        err[i] = norm(reshape(xn, len(xn)) - exp(-tn)) / norm(xn)

    return h, err

def RK2(f, tn, xn, h):
    k1 = f(tn, xn)
    k2 = f(tn + 0.5*h, xn + 0.5*h*k1)
    return xn + h*k2

def RK3(f, tn, xn, h):
    k1 = f(tn, xn)
    k2 = f(tn + 0.5*h, xn + 0.5*h*k1)
    k3 = f(tn + h, xn - h*k1 + 2.0*h*k2)
    return xn + h/6*(k1 + 4*k2 + k3)

def RK4(f, tn, xn, h):
    k1 = f(tn, xn)
    k2 = f(tn + 0.5*h, xn + 0.5*h*k1)
    k3 = f(tn + 0.5*h, xn + 0.5*h*k2)
    k4 = f(tn + h, xn + h*k3)
    return xn + h/6*(k1 + 2*k2 + 2*k3 + k4)

def FuncGrav(t, x):
    f = zeros(4)
    f[0] = x[2]
    f[1] = x[3]
    r = sqrt(x[0]**2 + x[1]**2)
    f[2] = -x[0] / r**3
    f[3] = -x[1] / r**3
    return f

def FuncExo7(t, x):
    f = zeros(2)
    f[0] = x[1]
    f[1] = t-2.0/t*x[1] + 2.0/(t*t)*x[0]
    return f

def FuncExo10(t, x):
    f = zeros(2)
    f[0] = x[0]*(1.0-x[1])
    f[1] = x[1]*(x[0]-1.0)
    return f

def FuncExo14(t, x):
    f = zeros(2)
    f[0] = x[1]
    f[1] = -sin(x[0])
    return f

def FuncExo17(t, pt):
    f = zeros(2)
    x = pt[0]; y = pt[1]
    f[0] = -2.0*sin(x)**4*sin(y)*cos(y) + 2.0*y*cos(x)**2
    f[1] = -4.0*sin(x)**3*cos(x)*cos(y)**2 + 2.0*cos(x)*sin(x)*y*y
    return f

def ode_solveAB3(f, t0, x0, tf, h):
    # nombre d'iterations en temps
    nb_iter = int(ceil((tf-t0) / h))
    h = (tf - t0) / nb_iter
    
    # on initialise les vecteurs tn et xn
    tn = zeros(nb_iter+1)
    m = 1
    if (type(x0) == ndarray):
        m = len(x0)
    
    xn = zeros([nb_iter+1, m])
    
    # ecrire ici votre code qui calcule la suite des xn avec le schema de AB3
    xn[0,:] = x0; tn[0] = t0;
    for i in range(nb_iter):
        if (i == 0):
            fnm2 = f(t0, xn[0, :])
            xn[i+1,:] = RK3(f, tn[i], xn[i,:], h)
        elif (i == 1):
            fnm1 = f(tn[1], xn[1, :])
            xn[i+1,:] = RK3(f, tn[i], xn[i,:], h)
        else:
            fn = f(tn[i], xn[i, :])
            xn[i+1, :] = xn[i, :] + h/12*(23.0*fn - 16.0*fnm1 + 5.0*fnm2)
            fnm2 = fnm1.copy(); fnm1 = fn.copy();
            
        tn[i+1] = t0 + (i+1)*h
    
    # a la fin de la fonction il faut retourner t_n et x_n
    return tn, xn

