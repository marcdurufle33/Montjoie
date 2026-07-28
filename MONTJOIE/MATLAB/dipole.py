from math import *
from pylab import *

class Dipole:
    def __init__(self):
        self.c0 = 299792458.0
        self.gamma = 1.125e13;
        self.Omega0 = 3.4586e15;
        self.Omega = self.Omega0 + 1j*self.gamma
        self.Lambda0 = 545e-9;
        self.z = 0;

    def __init__(self, z0):
        self.c0 = 1.0;
        self.gamma = 0.01;
        self.Omega0 = 2.0*pi;
        self.Omega = self.Omega0 + 1j*self.gamma
        self.Lambda0 = 1.0;
        self.z = z0;
    
    def Func(self, x):
        d_omega = x[0]+1j*x[1]
        k = (self.Omega + d_omega)/self.c0;
        z = self.z;
        coef = 1.5*self.gamma/(k*k*k)*(-k*k/z + 1j*k/(z*z) + 1.0 /(z*z*z))*exp(-1j*k*z);
        f = zeros([2]);
        f[0] = real(coef - d_omega);
        f[1] = imag(coef - d_omega);
        return f;
    
    def dFunc(self, x):
        h = 1e-6*abs(self.Omega);
        df = zeros([2, 2]);
        col0 = (self.Func(array([x[0]+h, x[1]])) - self.Func(array([x[0]-h, x[1]]))) / (2*h)
        col1 = (self.Func(array([x[0], x[1]+h])) - self.Func(array([x[0], x[1]-h]))) / (2*h)
        df[:, 0] = col0;
        df[:, 1] = col1;
        return df;


    def FindSol(self, x0):
        x = x0.copy()
        f = self.Func(x0)
        k = 0
        while ((norm(f) > 1.0) and (k < 10)):
            print "Residu = ", norm(f)
            df = self.dFunc(x)
            x = x - solve(df, f)
            f = self.Func(x)
            k = k+1
        
        return x

