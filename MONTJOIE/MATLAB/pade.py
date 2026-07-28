from visuND import *
from numpy import *
from pylab import *

tau = 0.42466090014401
omega = 2.0*pi
f0 = lambda t : exp(-0.5*((t-4.0)/tau)**2)*cos(omega*t)

K = array(loadMat('../Kh.dat').todense())
invM = loadtxt('../invMh.dat')
A = dot(diag(invM), K)
L, V = eig(A)
invV = inv(V)
F = invM[0]*invV[:, 0]

num = find(imag(L)>1e-4)
Np = argsort(imag(L[num]))
n0 = num[Np]

w = loadtxt('../mat.dat')
C = loadtxt('../Denom.dat');
D = poly1d(array(C[::-1]))
C[1::2] = -C[1::2]; Nr = poly1d(array(C[::-1]))
R = lambda z : Nr(z) / D(z)
tabR = [lambda z : poly1d(array(w[::-1, 0]))(z) / D(z), lambda z : poly1d(array(w[::-1, 1]))(z) / D(z), lambda z : poly1d(array(w[::-1, 2]))(z) / D(z), lambda z : poly1d(array(w[::-1, 3]))(z) / D(z)]

ci = array([0.0694318442029738, 0.330009478207572, 0.669990521792428, 0.930568155797026])

gam = 0.136339
Ns = poly1d(loadtxt('../NumerS.dat')[::-1])
Ds = lambda z : (1.0-gam*z)**10
Rs = lambda z : Ns(z) / Ds(z)
ws = loadtxt('../A.dat')
tabRs = [lambda z : poly1d(array(ws[::-1, 0]))(z) / Ds(z), lambda z : poly1d(array(ws[::-1, 1]))(z) / Ds(z), lambda z : poly1d(array(ws[::-1, 2]))(z) / Ds(z), lambda z : poly1d(array(ws[::-1, 3]))(z) / Ds(z), lambda z : poly1d(array(ws[::-1, 4]))(z) / Ds(z)]

ci_s = array([0, 0.172673164646012, 0.5, 0.827326835353989, 1])
