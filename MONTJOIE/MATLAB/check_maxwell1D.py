#!/usr/bin/python

from visuND import *
from pylab import *

Uprev = load1D('../Un0.dat');
Un = load1D('../Un1.dat');
Unext = load1D('../Un2.dat');

N = 500000
z = linspace(0, 1e-4, N);

dt = 0.1e-17
c = 299792458
gammaNL = 1e-32
cte = 1.0
#cte = 94096693978.1648
#gammaNL = 0
#cte = 0
mu = 4*pi*1e-7
epsilon = 1.0/(mu*c*c)
sigma = 10.0
eps_inf = 1.44

wk = (2*pi*c/6.5e-6)
wk2 = wk**2
alpha_k = 0.5
Tk = 1e-11

mixed_formulation = False

#wk2_bis = (2*pi*c/2.5e-6)**2
#alpha_k_bis = 0.325

print "N = ", N, len(Un)
Eprev = Uprev[0:N]
En = Un[0:N]
Enext = Unext[0:N]
Pprev = Uprev[N:2*N]
Pn = Un[N:2*N]
Pnext = Unext[N:2*N]
#Pprev_bis = Uprev[2*N:3*N]
#Pn_bis = Un[2*N:3*N]
#Pnext_bis = Unext[2*N:3*N]

dE_dz2 = calcule_acc(z, En);
#dE_dz2 = zeros(len(En))+1j*zeros(len(En))
#dz = z[1]
#dE_dz2[1:N-1] = (En[2:N] - 2.0*En[1:N-1] + En[0:N-2])/(dz*dz)
dE_dt2 = (Eprev - 2.0*En + Enext)/(dt*dt)
dE_dt = (Enext - Eprev)/(2.0*dt)
dP_dt = (Pnext - Pprev)/(2.0*dt)
dP_dt2 = (Pprev - 2.0*Pn + Pnext)/(dt*dt)
#dPbis_dt2 = (Pprev_bis - 2.0*Pn_bis + Pnext_bis)/(dt*dt)
#dP_dt2 = zeros(len(En))+1j*zeros(len(En))

Pnl_prev = abs(Eprev)**2*Eprev
Pnl = abs(En)**2*En
Pnl_next = abs(Enext)**2*Enext

P_NL = (Pnl_prev - 2.0*Pnl + Pnl_next)/(dt*dt)

if (mixed_formulation):
    Uprev = load1D('../Hn0.dat');
    Un = load1D('../Hn1.dat');
    Unext = load1D('../Hn2.dat');
    Hprev = Uprev[0:N]
    Hn = Un[0:N]
    Hnext = Unext[0:N]
    Qprev = Uprev[N:2*N]
    Qn = Un[N:2*N]
    Qnext = Unext[N:2*N]
    
    dH_dt = (Hnext - Hprev)/(2.0*dt)
    dQ_dt = (Qnext - Qprev)/(2.0*dt)
    dE_dz = calcule_vitesse(z, En);
    dH_dz = calcule_vitesse(z, Hn);
    P_NL = (Pnl_next - Pnl_prev)/(2.0*dt)

    rhsE = eps_inf/c*dE_dt + mu*sigma*c*En + 1.0/c*dP_dt - dH_dz + cte*cte*gammaNL/(c*epsilon)*P_NL
    print "erreur sur E ", norm(rhsE[1:]) / norm(dH_dz)
    rhsH = 1.0/c*dH_dt - dE_dz
    print "erreur sur H ", norm(rhsH[1:]) / norm(dE_dz)
    rhsP = 1.0/wk*dP_dt - Qn
    print "erreur sur P ", norm(rhsP[1:])/norm(Qn)
    rhsQ = 1.0/wk*dQ_dt + 1.0/(wk*Tk)*Qn + Pn - alpha_k*En
    print "erreur sur Q ", norm(rhsQ[1:])/norm(Pn)
else:
    rhsE = eps_inf/(c*c)*dE_dt2 + mu*sigma*dE_dt + 1.0/(c*c)*dP_dt2 - dE_dz2 + cte*cte*gammaNL/(c*c*epsilon)*P_NL;
    print "erreur sur E ", norm(rhsE[1:]) / norm(dE_dz2)
    if (Tk == 0):
        rhsP = 1.0/wk2*dP_dt2 + Pn - alpha_k*En
    else:
        rhsP = 1.0/wk2*dP_dt2 + 1.0/(wk2*Tk)*dP_dt + Pn - alpha_k*En
        
    print "erreur sur P ", norm(rhsP[1:])/norm(Pn)
    #rhsP = 1/wk2_bis*dPbis_dt2 + Pn_bis - alpha_k_bis*En
    #print "erreur sur Pbis ", norm(rhsP[1:])/norm(Pn)
