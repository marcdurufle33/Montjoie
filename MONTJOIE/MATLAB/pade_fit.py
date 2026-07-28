from numpy import *

def index_multi(x, omega, epsmes):
    L = len(x)//4;
    eps = x[0]*ones(len(omega))+x[1]*1j*ones(len(omega));
    for kk in range(L):
        omega_r = x[kk*4+2]
        omega_i = x[kk*4+3]
        sig_r = x[kk*4+4]
        sig_i = x[kk*4+5]
        eps = eps + (sig_r + 1j*sig_i) / (omega - (omega_r+1j*omega_i))
    
    fun = zeros(2*len(omega))
    for kk in range(len(omega)):
        fun[2*kk]=real(eps[kk]-epsmes[kk]);
        fun[2*kk+1]=imag(eps[kk]-epsmes[kk]);
    
    return fun

def index_frac(x, omega, epsmes):
    L = len(x)//4;
    eps = zeros(len(omega))+1j*zeros(len(omega));
    omega_r = x[0]
    omega_i = x[1]
    eps = eps + (omega_r + 1j*omega_i) / (omega - (omega_r+1j*omega_i))
    
    fun = zeros(2*len(omega))
    for kk in range(len(omega)):
        fun[2*kk]=real(eps[kk]-epsmes[kk]);
        fun[2*kk+1]=imag(eps[kk]-epsmes[kk]);
    
    return fun
