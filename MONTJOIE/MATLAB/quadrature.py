#!/usr/bin/env python

from numpy import *
import math

def tgamma(n):

    prod = 1
    for i in range(n-1):
        prod = prod*(i+1)

    return prod

def GetJacobiPolynomial(n, alpha, beta):

    
    ab = zeros((n,2))
  
    ab[0,0] = (beta-alpha)/(alpha + beta + 2)
    ab[0,1] = tgamma(alpha+1) * tgamma(beta+1) / tgamma(alpha+beta+2)
    if (n > 1):
        ab[1,0] = double(pow(beta,2)-pow(alpha,2))/(( 2+alpha+beta)*((2+alpha+beta)+2))
        ab[1,1] = 4.0*(alpha+1)*(beta+1)/(pow(alpha+beta+2,2)*(alpha+beta+3))

    
    for i in range(n-2):
        nab = 2*(i+2)+alpha+beta
        ab[i+2,0] = double(pow(beta,2)-pow(alpha,2))/(nab*(nab+2))
        ab[i+2,1] = 4.0*(i+2+alpha)*(i+2+beta)*(i+2)*(i+2+alpha+beta)/(pow(nab,2)*(nab+1)*(nab-1))

    return ab

def SolveGauss(n, alpha, beta, ab):

    J = zeros((n, n))

    for i in range(n):
        J[i, i] = ab[i, 0]

    for i in range(1,n):
        J[i,i-1] = math.sqrt(ab[i,1])
        J[i-1,i] = J[i,i-1]
        
    [x, v] = linalg.linalg.eig(J)

    w = zeros(n)
    for i in range(n):
        x[i] = (x[i] + 1.0)/2
        w[i] = ab[0,1]*pow(v[0,i],2);

    return x, w

def ComputeGaussJacobi(r, alpha, beta):
    if (r == 0):
        return array([0.5]), array([1.0])
    n = r+1
    ab = GetJacobiPolynomial(n, alpha, beta);
    [x, w] = SolveGauss(n, alpha, beta, ab);
    return x, w
    #Sort(x, w);


def ComputeLobattoJacobi(r, alpha, beta):

    n = r+1

    ab = GetJacobiPolynomial(n, alpha, beta);

    ab[n-1,0] = double(alpha-beta)/(2*(n-2)+alpha+beta+2);
    ab[n-1,1] = 4.0*((n-2)+alpha+1)*((n-2)+beta+1)*((n-2)+alpha+beta+1)/((2*(n-2)+alpha+beta+1)*pow((2*(n-2)+alpha+beta+2),2))
  
    [x, w] = SolveGauss(n, alpha, beta, ab);

    points = zeros(n)
    poids = zeros(n)
    nb = 0
    for i in range(n):
        if (abs(x[i]) < 1e-12):
            points[0] = 0
            poids[0] = w[i]
        elif (abs(x[i] - 1.0) < 1e-12):
            points[r] = 1.0
            poids[r] = w[i]
        else:
            nb = nb+1
            points[nb] = x[i]
            poids[nb] = w[i]
    
        
    return points, poids

class LagrangeFunctions:
    def __init__(self, points):
        n = len(points)
        self.points = points.copy()
        self.cte_phi = zeros(n)
        for i in range(n):
            self.cte_phi[i] = 1.0
            for j in range(n):
                if (j != i):
                    self.cte_phi[i] *= points[i] - points[j];
	    
            self.cte_phi[i] = 1.0/self.cte_phi[i]

    def ComputeValuesPhiRef(self, x):
        denom = 0.0
        n = len(self.points)
        if (n == 1):
            return ones(n)
        
        phi = zeros(n)
        for i in range(n):
            if (abs(x-self.points[i]) <= 1e-15):
                phi = zeros(n)
                phi[i] = 1.0
                return phi
            else:
                term = self.cte_phi[i]/(x - self.points[i]);
                phi[i] = term;
                denom += term;
        
        denom = 1.0/denom;
        phi *= denom
        return phi

    def ComputeGradPhi(self):
        n = len(self.points)
        val_grad = zeros([n, n])
        if (n == 1):
            return val_grad

        for i in range(n):
            s = 0.0
            for j in range(n):
                if (j != i):
                    val_grad[j, i] = self.cte_phi[j]/(self.cte_phi[i]*(self.points[i]-self.points[j]));  
                    s += val_grad[j, i];

            val_grad[i, i] = -s;
            if (abs(val_grad[i, i]) <= 1e-15):
                val_grad[i, i] = 0.0
        
        return val_grad



