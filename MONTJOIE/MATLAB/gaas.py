#La = linspace(150e-9, 2000e-9, 1001);
from numpy import *

def gaas_index(La):
    eps_inf = -0.54651;
    Om_p = array([1.4377, 2.7229, 2.8922, 4.5222, 4.9278])
    Om2_p = array([-0.05948, -1.2972, -0.23992, -0.42072, -0.19972])
    Sig_p = array([0.01981, 7.8336, 2.706, 2.1137, -1.243])
    Sig_p2 = array([0.01122, 8.3274, 1.616, 4.6445, 1.4424])
    e = 1.602176634e-19
    h = 6.62607015e-34
    Om = (Om_p+1j*Om2_p)*e/(h/(2.0*pi))
    Sig = (Sig_p+1j*Sig_p2)*e/(h/(2.0*pi))
    c = 299792458.0;
    om = (2*pi*c)/La;
    xOm = h/(2*pi)*om/e;
    eps = ones(om.shape)*eps_inf + 1j*zeros(om.shape);
    for k in range(5):
        eps = eps + 1j*Sig[k] / (om - Om[k]) + 1j*conj(Sig[k])/(om + conj(Om[k]))

    return eps

def get_poly_index(eps_inf, Om_p, Om2_p, Sig_p, Sig_p2):
    Denom = poly1d([1.0])
    for i in range(len(Om_p)):
        Denom = Denom*poly1d([1.0, -(Om_p[i]+1j*Om2_p[i])])
        Denom = Denom*poly1d([1.0, +(Om_p[i]-1j*Om2_p[i])])

    Numer = eps_inf*Denom
    for i in range(len(Om_p)):
        Pol = 1j*poly1d([Sig_p[i]+1j*Sig_p2[i]])*poly1d([1.0, (Om_p[i]-1j*Om2_p[i])])
        Pol = Pol + 1j*poly1d([Sig_p[i]-1j*Sig_p2[i]])*poly1d([1.0, -(Om_p[i]+1j*Om2_p[i])])
        for j in range(len(Om_p)):
            if (j != i):
                Pol = Pol*poly1d([1.0, -(Om_p[j]+1j*Om2_p[j])])
                Pol = Pol*poly1d([1.0, (Om_p[j]-1j*Om2_p[j])])

        Numer = Numer + Pol

    return Numer, Denom
                
                
def gaas_first(x, omega, epsmes):
    eps_inf = x[0];
    omega1r = x[1];
    omega1i = x[2];
    omega1bisr = x[3];
    omega1bisi = x[4];
    omega1 = omega1r+1j*omega1i;
    omega1bis = omega1bisr+1j*omega1bisi;
    eps = eps_inf*(omega-omega1)/(omega-omega1bis) * (omega+conj(omega1))/(omega+conj(omega1bis));
    fun = zeros(2*len(omega))
    for kk in range(len(omega)):
        fun[2*kk] = real(eps[kk]-epsmes[kk]);
        fun[2*kk+1] = imag(eps[kk]-epsmes[kk]);
    return fun

def gaas_multi(x, omega, epsmes):
    L = len(x)//4;
    eps = x[0]*(ones(len(omega))+1j*zeros(len(omega)));
    for kk in range(L):
        omega1r = x[kk*4+1]
        omega1i = x[kk*4+2]
        omega1bisr = x[kk*4+3]
        omega1bisi = x[kk*4+4]
        omega1 = omega1r+1j*omega1i;
        omega1bis = omega1bisr+1j*omega1bisi;                                                                           
        eps = eps * (omega-omega1)/(omega-omega1bis)*(omega+conj(omega1))/(omega+conj(omega1bis))

    fun = zeros(2*len(omega))
    for kk in range(len(omega)):
        fun[2*kk]=real(eps[kk]-epsmes[kk]);
        fun[2*kk+1]=imag(eps[kk]-epsmes[kk]);
    
    return fun

def metal_multi(x, omega, epsmes):
    L = (len(x)-3)//4
    # Drude part
    eps_inf = x[0]; gam = x[1]; sig = x[2];
    eps = eps_inf*(ones(len(omega))+1j*zeros(len(omega)));
    eps = eps + 1j*sig/omega - 1j*sig/(omega + 1j*gam)
    # dielectric part
    for k in range(L):
        om_k = x[3+4*k] + 1j*x[4+4*k];
        sig_k = x[5+4*k] + 1j*x[6+4*k];
        eps = eps + 1j*sig_k/(omega-om_k) + 1j*conj(sig_k)/(omega+conj(om_k))

    fun = zeros(2*len(omega));
    fun[0::2] = real(eps - epsmes)
    fun[1::2] = imag(eps - epsmes)
    return fun

def dielec_multi(x, omega, epsmes):
    L = (len(x)-1)//4
    eps_inf = x[0];
    eps = eps_inf*(ones(len(omega))+1j*zeros(len(omega)));
    for k in range(L):
        om_k = x[1+4*k] + 1j*x[2+4*k];
        sig_k = x[3+4*k] + 1j*x[4+4*k];
        eps = eps + 1j*sig_k/(omega-om_k) + 1j*conj(sig_k)/(omega+conj(om_k))

    fun = zeros(2*len(omega));
    fun[0::2] = real(eps - epsmes)
    fun[1::2] = imag(eps - epsmes)
    return fun
