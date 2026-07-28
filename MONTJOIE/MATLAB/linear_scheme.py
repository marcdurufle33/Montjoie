from numpy import *
from quadrature import *

class PadeScheme:
    """ basic object that constructs arrays needed to implement Pade schemes
    Usage : p = PadeScheme(r)
            p = PadeScheme(r, True)
    
    r : order of the scheme (must be even)
    The second argument is optional, if true the weights are associated with the stable algorithm
    p : object containing the coefficients associated with Pade scheme
    p.Numer : numerator of R(z)
    p.Denom : denominator of R(z)
    p.ci : points where the source F must be evaluated
    p.omega_i : weights associated with points c_i
    p.real_root : real root (if present)
    p.complex_root : complex roots of denominator

    p.Advance is a method that can be used to compute X_n+1 from X by using the Pade scheme constructed
    """
    # constructor of the class PadeScheme
    def __init__(self, r, stable = True):
        #print "Ordre = ", r
        m = r//2
        n = m
        # coefficients of numerator and denominator of Pade approximant
        coefP = zeros(m+1)
        coefQ = zeros(m+1)
        for i in range(m+1):
            coefP[i] = float(tgamma(m+n+1-i) * tgamma(n+1)) / (tgamma(m+n+1)*tgamma(i+1)*tgamma(n+1-i))
            coefQ[i] = coefP[i]
            if (i%2 != 0):
                coefQ[i] = -coefP[i]

        self.stable_algo = stable
        self.Numer = poly1d(coefP[::-1])
        self.Denom = poly1d(coefQ[::-1])
        #print "Numerator = ", self.Numer
        #print "Denominator = ", self.Denom
        #print "Roots of denominator = ", self.Denom.r

        # Gauss-Legendre points
        points = ComputeGaussJacobi(m-1, 0, 0)[0]
        
        CoefC = zeros(2*m+1)
        for k in range(2*m+1):
            CoefC[k] = 1.0/(2.0**(k-1) * tgamma(k+1))

        # VanDerMonde matrix to find weights omega_i
        VDM = zeros([m, m])
        for i in range(len(points)):
            for j in range(len(points)):
                VDM[i, j] = (points[j]-0.5)**i

        self.ci = points
        #print "Points = ", self.ci
        
        invVDM = linalg.inv(VDM)
        self.omega_i = zeros([m, m])
        # loop over A^powL \Delta t^{powL+1}
        A = zeros([m, m])
        for powL in range(m):
            r = powL+1
            jmin = 2
            if (r%2 == 1):
                jmin = 1

            vec_f = zeros([2*m])
            for j in range(jmin, 2*m+1-r, 2):
                for i in range(0, r, 2):
                    vec_f[j-1] += coefP[i]*CoefC[r+j-i-1]

                for i in range(1, r, 2):
                    vec_f[j-1] -= coefP[i]*CoefC[r+j-i-1]

            if (stable):
                # for stable algorithm, we store the coefficients in a matrix A
                for j in range(m):
                    A[powL, j] = vec_f[j]
            else:                    
                rhs = zeros(m)
                for i in range(m):
                    rhs[i] = tgamma(i+1) * vec_f[i]
                    
                # weights omega_i for this power of A 
                self.omega_i[powL, :] = dot(invVDM, rhs)
                #print "Weights for r = ", r, " : ", self.omega_i[i,:]
                
        # roots of denominator
        self.complex_root = []
        for i in range(len(self.Denom.r)):
            if (imag(self.Denom.r[i]) == 0):
                self.real_root = real(self.Denom.r[i])
            else:
                z = self.Denom.r[i]
                if (imag(z) > 0):
                    self.complex_root.append(z)
                        
        if (stable):
            # polynomes du numerateur
            PolNumer = []
            PolNumer.append(poly1d([1.0]))
            Q = poly1d([1.0])
            z = poly1d([1.0, 0.0])
            for i in range(len(self.complex_root)-1, -1, -1):
                PolNumer.append(z*Q)
                Q = Q*poly1d([1.0/abs(self.complex_root[i])**2, 2.0*real(1.0/self.complex_root[i]), 1.0])
                PolNumer.append(Q)
                
            if (m%2 == 1):
                Q = Q*poly1d([1.0/self.real_root, 1.0])
                PolNumer.append(Q)
                
            # polynomials needed to expand the source are enumerated
            pol_source = []
            Denom = poly1d([1.0])
            num = len(PolNumer) - 2
            if (m%2 == 1):
                pol_source.append(PolNumer[num]); num -= 1
                
                # real root present, we multiply by (1 - z / real_root)
                Denom = poly1d([-1.0/self.real_root, 1.0])
                
            # loop on complex roots
            for i in range(len(self.complex_root)):
                pol_source.append(PolNumer[num]*Denom); num -= 1
                pol_source.append(PolNumer[num]*Denom); num -= 1
                
                Denom = Denom*poly1d([1.0/abs(self.complex_root[i])**2, -2.0*real(1.0/self.complex_root[i]), 1.0])
        
            # inverse of coefficients, to express x^i as a combination of pol_source
            coefBase = zeros([m, m])
            
            for j in range(m):
                for i in range(pol_source[j].coeffs.shape[0]):
                    coefBase[i, j] = pol_source[j][i]
            
            coefBase = linalg.inv(coefBase)
            polA = zeros([m, m])
            for powA in range(m):
                rhs = A[powA, :]
                
                for j in range(m):
                    for k in range(m):
                        polA[k, j] += coefBase[k, powA]*rhs[j]
            
            # computation of weights
            for powA in range(m):
                rhs = zeros(m)
                for i in range(m):
                    rhs[i] = tgamma(i+1) * polA[powA, i]
                
                self.omega_i[powA, :] = dot(invVDM, rhs)

    def SolveOperatorP2(self, a, F, sys, apply_mass):
        if (F.dtype == dtype('complex128')):
            Fr = real(F)
            Fi = imag(F)
            Yr = self.SolveOperatorP2(a, Fr, sys, apply_mass)
            Yi = self.SolveOperatorP2(a, Fi, sys, apply_mass)
            Y = Yr+1j*Yi
        else:
            N = len(F)
            if (apply_mass):
                Fc = sys.ApplyM(F) + 1j*zeros(N)
            else:
                Fc = F + 1j*zeros(N)
                
            v = sys.SolveComplex(a, Fc)
            b = a / (a - conj(a))
            Y = 2*real(b*v)
        
        return Y
        
    def Advance(self, tn, dt, Yold, sys):
        """ Method to compute X_n+1 from X_n with Pade scheme
        Usage :
            Xnext = pade.Advance(t, dt, Xn, sys)
        
        t : current time t_n
        dt : time step
        Xn : iterate X_n
        sys : object describing the linear ODE 
        the method returns X_n+1

        the object containing the linear ODE must contain ApplyM, SolveM, ApplyK
        (multiplication by M, M^{-1} and K)
        GetSource (computation of F(t))
        and SolveReal / SolveComplex (inversion by M + a K)
        """
        m = self.Numer.order
        EvalF = [0]*m
        for i in range(m):
            tcurrent = tn + self.ci[i]*dt
            EvalF[i] = sys.GetSource(tcurrent)
            
        Y = Yold.copy() 
        if (self.stable_algo):
            # stable algorithm
            num_source = 0
            if (m%2 == 1):
                Fn = zeros(len(Y))
                for i in range(m):
                    Fn += dt*self.omega_i[num_source, i]*EvalF[i]
                    
                num_source += 1
                Fn += sys.ApplyM(Y)
                Fn += dt/self.real_root*sys.ApplyK(Y)
                
                Y = sys.SolveReal(-dt/self.real_root, Fn)
                
            for z in self.complex_root:
                b = 2*dt*real(1.0 / z)
                a = dt*dt / abs(z)**2
                Fn = a*sys.ApplyK(Y)
                for i in range(m):
                    Fn += dt*dt*self.omega_i[num_source, i]*EvalF[i]
                    
                Fn = b*Y + sys.SolveM(Fn)
                Fn = sys.ApplyK(Fn) + sys.ApplyM(Y)
                
                for i in range(m):
                    Fn += dt*self.omega_i[num_source+1, i]*EvalF[i]
                    
                Y = self.SolveOperatorP2(-dt/z, Fn, sys, False)
                num_source += 2
        else:
            # Horner algorithm
            Fn = zeros(len(Y))
            if (Y.dtype == dtype('complex128')):
                Fn = zeros(len(Y)) + 1j*zeros(len(Y))
            
            # evaluation of phi_n + (P - Q) X_n
            last_coef = m
            if (m%2 == 0):
                last_coef = m-1
                for i in range(m):
                    Fn += dt*self.omega_i[m-1, i]*EvalF[i]
                    
                Fn = sys.SolveM(Fn)
                
            Fn += (self.Numer[last_coef] - self.Denom[last_coef])*Y
            for k in range(last_coef-1, -1, -1):
                Fn = dt*sys.ApplyK(Fn)
                for i in range(m):
                    Fn += dt*self.omega_i[k, i]*EvalF[i]
                    
                if (k > 0):
                    Fn = sys.SolveM(Fn)
                    
                if (k%2 == 1):
                    Fn += (self.Numer[k] - self.Denom[k])*Y

            # then we solve Q (X_{n+1} - X_n) = Fn
            if (m%2 == 1):
                Fn = sys.SolveReal(-dt/self.real_root, Fn)
                
            for k in range(len(self.complex_root)):
                apply_mass = True
                if ((k == 0) and (m%2 == 0)):
                    apply_mass = False

                Fn = self.SolveOperatorP2(-dt/self.complex_root[k], Fn, sys, apply_mass)

            Y += Fn          
        
        return Y

class LinearSdirkScheme:
    """  basic object that constructs arrays needed to implement Linear SDIRK schemes
     Usage : p = LinearSdirkScheme(r, extraS)

    r : order of the scheme (must be even)
    extraS : number of additional stages (0 by default)
    p : object containing the coefficients associated with the scheme
    p.Numer : numerator of R(z)
    p.Denom : denominator of R(z)
    p.gamma : unique pole of R(z)
    p.ci : points where the source F must be evaluated
    p.omega_i : weights associated with points c_i

    p.Advance is a method that can be used to compute X_n+1 from X by using the linear SDIRK scheme constructed
    """

    def Init(self, order, gamma, alpha):
        extraS = len(alpha)
        s = order-1
        self.gamma = gamma
        # denominator of R(z)
        self.Denom = poly1d([-gamma, 1.0])**(s+extraS)
        CoefN = zeros(s+2)
        # Taylor expansion of exponential
        for i in range(s+2):
            CoefN[i] = 1.0/tgamma(i+1)

        #print "Expo = ", CoefN
        # polynomial exp(z) * (1 - gamma z)^{s+extraS}
        self.Numer = poly1d(CoefN[::-1])*self.Denom

        # we truncate to remove higher order coefficients
        CoefN = zeros(s+1+extraS)
        for i in range(s+1+extraS):
            CoefN[i] = self.Numer[i]

        # we add parameters alpha
        for i in range(s+2, s+1+extraS):
            CoefN[i] += alpha[i-s-2]

        # Numerator of R(z)
        self.Numer = poly1d(CoefN[::-1])
        #print "Numerator = ", self.Numer
        #print "Denominator = ", self.Denom
        #print "Gamma = ", self.gamma

        # roots of numerator
        self.real_roots = []
        self.complex_roots = []
        for i in range(len(self.Numer.r)):
            if (imag(self.Numer.r[i]) == 0):
                self.real_roots.append(float(real(self.Numer.r[i])))
            else:
                z = self.Numer.r[i]
                if (imag(z) > 0):
                    self.complex_roots.append(z)
        
        #print "real roots = ", self.real_roots
        #print "complex roots = ", self.complex_roots
        nb_terms = self.Numer.order
        
        invFacto = zeros(2*(s+extraS+2))
        CoefC = zeros(2*(s+extraS+2))
        invFacto[0] = 1.0
        CoefC[0] = 1.0
        # Taylor expansion of exponential
        for k in range(1, 2*(s+extraS+2)):
            invFacto[k] = invFacto[k-1] / float(k)
            CoefC[k] = CoefC[k-1] / float(2*k)
        
        # coefficients alpha for inhomogeneous case
        nb_terms_add = nb_terms+1
        A = zeros([nb_terms_add, nb_terms_add])
        for r in range(1, nb_terms_add+1):
            for j in range(1, nb_terms_add+1):
                err = 0
                for i in range(0, min(r-1, self.Numer.order)+1):
                    k = r + j - i - 1
                    signK = 1.0
                    if (k%2 == 1):
                        signK = -signK
                        
                    err += (self.Denom[i] - signK*self.Numer[i])*CoefC[k]
                    
                A[r-1, j-1] = err
                
        orderN = self.Numer.order
        if (self.stable_algo):
            # for stable algorithm, we need to compute coefficient with another polynomial base 
            # (different from the canonical base)
            Pol = poly1d([-self.gamma, 1.0])
            
            z = poly1d([1.0, 0.0])
            Q = poly1d([1.0])
            PolNumer = []
            PolNumer.append(Q)
            for i in range(len(self.real_roots)-1, -1, -1):
                Q = Q*poly1d([-1.0/self.real_roots[i], 1.0])
                PolNumer.append(Q)
                
            for i in range(len(self.complex_roots)-1, -1, -1):
                PolNumer.append(z*Q)
                Q *= poly1d([1.0/abs(self.complex_roots[i])**2, -2.0*real(1.0/self.complex_roots[i]), 1.0])
                PolNumer.append(Q)
                
            # polynomials needed to expand the source are enumerated
            pol_source = []
            monoGamma = poly1d([1.0])
            num = len(PolNumer) - 2
            for i in range(len(self.complex_roots)):
                pol_source.append(PolNumer[num]*monoGamma); num -= 1
                pol_source.append(PolNumer[num]*monoGamma); num -= 1
                monoGamma = monoGamma*Pol*Pol
                
            for i in range(len(self.real_roots)):
                pol_source.append(PolNumer[num]*monoGamma); num -= 1
                monoGamma = monoGamma*Pol
                
            # inverse of coefficients, to express x^i as a combination of pol_source
            coefBase = zeros([orderN, orderN])
            
            for j in range(orderN):
                for i in range(pol_source[j].coeffs.shape[0]):
                    coefBase[i, j] = pol_source[j][i]
            
            coefBase = linalg.inv(coefBase)
            polA = zeros([orderN, nb_terms+1])
            for powA in range(nb_terms):
                rhs = A[powA, :]
                
                for j in range(nb_terms+1):
                    for k in range(orderN):
                        polA[k, j] += coefBase[k, powA]*rhs[j]

        # Gauss-Legendre points
        points = ComputeGaussJacobi(nb_terms, 0, 0)[0]
        self.ci = points
        #print "points", self.ci

        # VanDerMonde matrix to find weights omega_i
        VDM = zeros([len(points), len(points)])
        for i in range(len(points)):
            for j in range(len(points)):
                VDM[i, j] = (points[j]-0.5)**i

        invVDM = linalg.inv(VDM)
        self.omega_i = zeros([orderN, len(points)])
        for powA in range(orderN):

            rhs = zeros(len(points))
            if (self.stable_algo):
                for j in range(0, polA.shape[1]):
                    rhs[j] = polA[powA, j] / invFacto[j]
            else:
                for j in range(0, nb_terms+1):
                    rhs[j] = A[powA, j] / invFacto[j]
                
            # weights omega_i for this power of A 
            self.omega_i[powA, :] = dot(invVDM, rhs)
            #print "Weights for r = ", powA, " : ", self.omega_i[powA,:]

    # constructor of the class LinearSdirkScheme
    def __init__(self, r, extraS = 0, stable = True):
        alpha = zeros(extraS)
        self.stable_algo = stable
        if (r == 2):
            if (extraS == 0):
                gamma = 0.5
            else:
                print("Case not implemented")
        elif (r == 3):
            if (extraS == 0):
                gamma = 0.5 + 0.5/sqrt(3.0)
            else:
                print("Case not implemented")
        elif (r == 4):
            if (extraS == 0):
                gamma = 1.0 / sqrt(3.0) * cos(pi/18) + 0.5
            elif (extraS == 1):
                gamma = 0.394337567297407
            else:
                print("Case not implemented")
        elif (r == 6):
            if (extraS == 0):
                gamma = 0.47326839125829532445558852540261
            elif (extraS == 1):
                gamma = 0.284064638011799
            elif (extraS == 2):
                gamma = 0.204071
                alpha[0] = 1.9839430662e-4
            else:
                print("Case not implemented")
        elif (r==8):
            if (extraS == 1):
                gamma = 0.217049743094304
            elif (extraS == 2):
                gamma = 0.16689
                alpha[0] = 2.9259251764e-6
            elif (extraS == 3):
                gamma = 0.136339
                alpha[0] = 2.767416226e-6
                alpha[1] = -3.464398093e-6
            else:
                print("Case not implemented")
        elif (r==10):
            if (extraS == 2):
                gamma = 0.141940
                alpha[0] = 2.2982637210e-8
            elif (extraS == 3):
                gamma = 0.151706
                alpha[0] = 2.459114959e-8
                alpha[1] = -4.3140917546e-8
            else:
                print("Case not implemented")
        elif (r==12):
            if (extraS == 3):
                gamma = 0.132572
                alpha[0] = 1.644515143e-10
                alpha[1] = -2.89891484131e-10
            else:
                print("Case not implemented")
        else:
            print("Case not implemented")
        
        self.Init(r, gamma, alpha)

    def Advance(self, tn, dt, Yold, sys):
        """ Method to compute X_n+1 from X_n with Pade scheme
        Usage :
            Xnext = pade.Advance(t, dt, Xn, sys)
        
        t : current time t_n
        dt : time step
        Xn : iterate X_n
        sys : object describing the linear ODE 
        the method returns X_n+1

        the object containing the linear ODE must contain ApplyM, SolveM, ApplyK
        (multiplication by M, M^{-1} and K)
        GetSource (computation of F(t))
        and SolveReal / SolveComplex (inversion by M + a K)
        """
        EvalF = [0]*self.omega_i.shape[1]
        for i in range(self.omega_i.shape[1]):
            tcurrent = tn + self.ci[i]*dt
            EvalF[i] = sys.GetSource(tcurrent)
            
        Y = Yold.copy() 
        if (self.stable_algo):
            # stable algorithm
            num_source = 0
                
            for z in self.complex_roots:
                b = -2*dt*real(1.0 / z)
                a = dt*dt / abs(z)**2
                Fn = a*sys.ApplyK(Y)
                for i in range(self.omega_i.shape[1]):
                    Fn += dt*dt*self.omega_i[num_source, i]*EvalF[i]
                    
                Fn = b*Y + sys.SolveM(Fn)
                Fn = sys.ApplyK(Fn) + sys.ApplyM(Y)
                
                for i in range(self.omega_i.shape[1]):
                    Fn += dt*self.omega_i[num_source+1, i]*EvalF[i]
                    
                Y = sys.SolveReal(-dt*self.gamma, Fn)
                Fn = sys.ApplyM(Y)
                Y = sys.SolveReal(-dt*self.gamma, Fn)
                num_source += 2

            for z in self.real_roots:
                b = -dt / z
                Fn = b*sys.ApplyK(Y) + sys.ApplyM(Y)
                for i in range(self.omega_i.shape[1]):
                    Fn += dt*self.omega_i[num_source, i]*EvalF[i]

                Y = sys.SolveReal(-dt*self.gamma, Fn)                    
                num_source += 1

        else:
            # Horner algorithm
            last_coef = self.Numer.order
            Fn = self.Numer[last_coef]*Y
            for k in range(last_coef-1, -1, -1):
                Fn = dt*sys.SolveM(sys.ApplyK(Fn)) + self.Numer[k]*Y
                if ( k < self.omega_i.shape[0]):                    
                    Ku = zeros(len(Y))
                    for i in range(self.omega_i.shape[1]):
                        Ku += dt*self.omega_i[k, i]*EvalF[i]
                    
                    Fn += sys.SolveM(Ku)
                
            # then we solve Q X_{n+1} = Fn
            for k in range(self.Numer.order):
                Ku = sys.ApplyM(Fn)
                Fn = sys.SolveReal(-dt*self.gamma, Ku)
                
            Y = Fn          
        
        return Y

class LinearRkScheme:
    """  basic object that constructs arrays needed to implement Linear Runge-Kutta schemes
     Usage : p = LinearRkScheme(r, extraS)

    r : order of the scheme (must be even)
    extraS : number of additional stages (0 by default)
    p : object containing the coefficients associated with the scheme
    p.Numer : R(z)
    p.ci : points where the source F must be evaluated
    p.omega_i : weights associated with points c_i

    p.Advance is a method that can be used to compute X_n+1 from X
    """
    def Init(self, order, alpha):
        extraS = len(alpha)-order-1
        nb_terms = order
        nb_terms_add = order
        orderN = len(alpha)-1
        self.Numer = poly1d(alpha[::-1])
        # roots of numerator
        self.real_roots = []
        self.complex_roots = []
        for i in range(len(self.Numer.r)):
            if (imag(self.Numer.r[i]) == 0):
                self.real_roots.append(float(real(self.Numer.r[i])))
            else:
                z = self.Numer.r[i]
                if (imag(z) > 0):
                    self.complex_roots.append(z)

        # we store 1/k! et 1 / 2^k k!
        invFacto = zeros(2*nb_terms_add+1);
        CoefC = zeros(2*nb_terms_add+1)
        invFacto[0] = 1.0; CoefC[0] = 1.0
        for k in range(1, 2*nb_terms_add+1):
            invFacto[k] = invFacto[k-1] / float(k);
            CoefC[k] = CoefC[k-1] / float(2*k);

        if (not self.stable_algo):
            # using basis functions associated with Gauss points
            # to compute coefficients omega_i
            self.ci, weights = ComputeGaussJacobi(order-2, 0, 0)
            Q = LagrangeFunctions(self.ci)
            
	    # derivatives of these basis functions
            dPhi = Q.ComputeGradPhi();
            
            self.omega_i = zeros([orderN, len(self.ci)])
            self.omega_i[0,:] = weights
            
            ValPhi0 = Q.ComputeValuesPhiRef(0.0);
            
            Decomp = eye(len(self.ci))
            DerMat = zeros([orderN-1, len(self.ci)]);
            for l in range(1, orderN):
                qtilde = dot(Decomp, ValPhi0)
                DerMat[l-1, :] = qtilde
                for i in range(len(self.ci)):
                    dphi_tilde = qtilde[i];
                    for k in range(1, orderN):
                        if (l <= orderN-k):
                            self.omega_i[k, i] += dphi_tilde*alpha[k+l];
                
                Decomp = dot(Decomp, dPhi)
            
            return;
        
        # computation of coefficients alpha_j^r for the right hand side phi
        # phi = dt \sum_{r=1}^\infty (dt A)^{r-1} \sum_{j=1}^\infty \alpha_j^r dt^{j-1} F^{j-1}
        # where \alpha_j^r = \sum_{i=0}^{min(r-1, m)} (D_i - (-1)^k N_i) / 2^k k!
        # where k = r+j-i-1
        A = zeros([nb_terms_add, nb_terms_add]);
        for r in range(1, nb_terms_add+1):
            err = zeros(nb_terms_add+1-r)
            for j in range(1, nb_terms_add+2-r):
                # partie denominateur
                err[j-1] = CoefC[r+j-1];
                
                # numerateur
                for i in range(0, min(r-1, len(alpha)-1)+1):
                    k = r + j-i - 1;
                    signK = 1.0;
                    if (k%2 == 1):
                        signK = -signK;
                    
                    err[j-1] -= signK*alpha[i]*CoefC[k];
                
                A[r-1, j-1] = err[j-1];
        
        # calcul des ci : on prend les points de Gauss-Legendre
        points = ComputeGaussJacobi(order-2, 0, 0)[0]
        
        # Vandermonde Matrix for ci
        # les factorielles sont mises dans le second membre    
        VDM = zeros([len(points), len(points)])
        for i in range(len(points)):
            for j in range(len(points)):
                VDM[i, j] = (points[j] - 0.5)**i
        
        invVDM = linalg.inv(VDM)
        
        # series (1 - z / lambda_n) (1 - z / lambda_{n-1}) .. (1 - z/lambda_2)
        PolNumer = []
        PolNumer.append(poly1d([1.0]))
        Q = poly1d([1.0])
        z = poly1d([1.0, 0.0])
        for i in range(len(self.real_roots)-1, -1, -1):
            Q = Q*poly1d([-1.0/self.real_roots[i], 1.0])
            PolNumer.append(Q)
            
        for i in range(len(self.complex_roots)-1, -1, -1):
            PolNumer.append(z*Q)
            Q *= poly1d([1.0/abs(self.complex_roots[i])**2, -2.0*real(1.0/self.complex_roots[i]), 1.0])
            PolNumer.append(Q)

        # polynomials needed to expand the source are enumerated
        pol_source = []
        num = len(PolNumer) - 2
        for i in range(len(self.complex_roots)):
            pol_source.append(PolNumer[num]); num -= 1
            pol_source.append(PolNumer[num]); num -= 1

        for i in range(len(self.real_roots)):
            pol_source.append(PolNumer[num]); num -= 1
        
        # inverse of coefficients, to express x^i as a combination of pol_source
        coefBase = zeros([orderN, orderN]);
        for j in range(orderN):
            for i in range(pol_source[j].coeffs.shape[0]):
                coefBase[i, j] = pol_source[j][i]

        coefBase = linalg.inv(coefBase)
        polA = zeros([orderN, nb_terms-1])
        for powA in range(nb_terms):
            rhs = A[powA, :]
                
            for j in range(nb_terms-1):
                for k in range(orderN):
                    polA[k, j] += coefBase[k, powA]*rhs[j]

        # weights omega_i^r are computed
        self.ci = points
        self.omega_i = zeros([orderN, len(points)])
        for powA in range(orderN):

            rhs = zeros(len(points))
            for j in range(0, polA.shape[1]):
                rhs[j] = polA[powA, j] / invFacto[j]
            
            # weights omega_i for this power of A 
            self.omega_i[powA, :] = dot(invVDM, rhs)
    
    # constructor of the class LinearRkScheme
    def __init__(self, r, extraS = 0, stable = True):
        # for the first coefficients we put the exact values (1/k!)
        alpha = zeros([r+1+extraS])
        alpha[0] = 1.0
        for k in range(r):
            alpha[k+1] = alpha[k] / (k+1)

        self.stable_algo = stable
        if (r == 2):
            if (extraS == 1):
                alpha[3] = 1.451277982649155e-01
            elif (extraS == 2):
                alpha[3] = 1.665532314108146e-01
                alpha[4] = 2.327815361933148e-02
            elif (extraS == 3):
                alpha[3] = 1.618342913053687e-01;
                alpha[4] = 3.289792611743811e-02;
                alpha[5] = 2.839528016518102e-03;
            elif (extraS == 4):
                alpha[3] = 1.642981320398038e-01;
                alpha[4] = 3.657769285804588e-02;
                alpha[5] = 5.035250867609586e-03;
                alpha[6] = 3.001880509358407e-04;
            elif (extraS == 5):
                alpha[3] = 1.626462249413356e-01;
                alpha[4] = 3.762678272315501e-02;
                alpha[5] = 5.996644250417070e-03;
                alpha[6] = 5.826143210213330e-04;
                alpha[7] = 2.487327304531716e-05;
            elif (extraS == 6):
                alpha[3] = 1.627509585676844e-01;
                alpha[4] = 3.773348832445807e-02;
                alpha[5] = 6.387803046851333e-03;
                alpha[6] = 7.489561665296774e-04;
                alpha[7] = 5.356270766078865e-05;
                alpha[8] = 1.713109940102836e-06
            elif (extraS == 7):
                alpha[3] = 1.640094942014296e-01; 
                alpha[4] = 3.840429977823329e-02; 
                alpha[5] = 6.724597512047917e-03; 
                alpha[6] = 8.718626803227696e-04; 
                alpha[7] = 7.857554562878064e-05; 
                alpha[8] = 4.327975378833797e-06; 
                alpha[9] = 1.072985856243921e-07;
            elif (extraS == 8):
                alpha[3] = 1.649990588856614e-01; 
                alpha[4] = 3.927394350377206e-02; 
                alpha[5] = 7.055384479248899e-03; 
                alpha[6] = 9.695797812914759e-04; 
                alpha[7] = 9.943224646288322e-05; 
                alpha[8] = 7.129812259258231e-06; 
                alpha[9] = 3.148056880771953e-07; 
                alpha[10] = 6.324920988294407e-09;
            elif (extraS > 8):
                print("Case not implemented")
        elif (r == 4):
            if (extraS == 1):
                alpha[5] = 4.730163010446185e-03;
            elif (extraS == 2):
                alpha[5] = 6.541349497416528e-03;
                alpha[6] = 4.395282130923843e-04;
            elif (extraS == 3):
                alpha[5] = 7.241999849787970e-03;
                alpha[6] = 7.614940065988191e-04;
                alpha[7] = 3.521874589831831e-05;
            elif (extraS == 4):
                alpha[5] = 7.603292194142675e-03; 
                alpha[6] = 9.535828377031919e-04; 
                alpha[7] = 7.298469178025099e-05; 
                alpha[8] = 2.500124976522895e-06;
            elif (extraS == 5):
                alpha[5] = 7.817918289656257e-03; 
                alpha[6] = 1.075759999127459e-03; 
                alpha[7] = 1.026588721744709e-04; 
                alpha[8] = 6.038353896295552e-06; 
                alpha[9] = 1.628169027707504e-07;
            elif (extraS == 6):
                alpha[5] = 7.992535147077134e-03; 
                alpha[6] = 1.180030987873825e-03; 
                alpha[7] = 1.307878349087823e-04; 
                alpha[8] = 1.020785594818226e-05; 
                alpha[9] = 4.943966219870204e-07; 
                alpha[10] = 1.097077616437946e-08;
            elif (extraS == 7):
                alpha[5] = 9.619397138072583e-03; 
                alpha[6] = 3.970757223041604e-03; 
                alpha[7] = 1.979923031733034e-03; 
                alpha[8] = 6.726632799312973e-04; 
                alpha[9] = 1.385778310637994e-04; 
                alpha[10] = 1.585824201586086e-05; 
                alpha[11] = 7.742514686545619e-07;
            elif (extraS == 8):
                alpha[5] = 8.105487675563905e-03; 
                alpha[6] = 1.249316412377197e-03; 
                alpha[7] = 1.531845812394507e-04; 
                alpha[8] = 1.473468121845849e-05; 
                alpha[9] = 1.071860716775002e-06; 
                alpha[10] = 5.510748021396615e-08; 
                alpha[11] = 1.766727504578043e-09; 
                alpha[12] = 2.623218531216638e-11;
            elif (extraS > 8):
                print("Case not implemented")
        elif (r == 6):
            if (extraS == 1):
                alpha[7] = 2.070461615593214e-04;
            elif (extraS == 2):
                alpha[7] = 2.204061707466545e-04;
                alpha[8] = 1.942982735313673e-05;
            elif (extraS == 3):
                alpha[7] = 2.073919102492977e-04;
                alpha[8] = 2.499262304459253e-05;
                alpha[9] = 1.453234258464881e-06;
            elif (extraS == 4):
                alpha[7] = 2.358338644436141e-04; 
                alpha[8] = 4.056334413908446e-05; 
                alpha[9] = 4.775871882059528e-06; 
                alpha[10] = 2.442645091656458e-07;
            elif (extraS > 4):
                print("Case not implemented")
        elif (r == 8):
            if (extraS == 1):
                alpha[9] = 1.684112035592431e-06;
            elif (extraS == 2):
                alpha[9] = 2.288709306973234e-06;
                alpha[10] = 9.960040692054680e-08;
            elif (extraS == 3):
                alpha[9] = 2.528206540248994e-06
                alpha[10] = 1.724423811134767e-07
                alpha[11] = 5.449535772542617e-09
            elif (extraS == 4):
                alpha[9] = 2.638893313733145e-06
                alpha[10] = 2.150620166601062e-07
                alpha[11] = 1.123553506837818e-08
                alpha[12] = 2.690758844819519e-10
            elif (extraS == 5):
                alpha[9] = 2.703333893632985e-06
                alpha[10] = 2.435581983430564e-07
                alpha[11] = 1.631043038503232e-08
                alpha[12] = 6.905312067380033e-10
                alpha[13] = 1.342332862257654e-11
            elif (extraS == 6):
                alpha[9] = 2.711246141311401e-06 ; 
                alpha[10] = 2.500568374959440e-07 ; 
                alpha[11] = 1.817647917892119e-08 ; 
                alpha[12] = 9.481642471601341e-10 ; 
                alpha[13] = 3.089127728872379e-11 ; 
                alpha[14] = 4.655664953646905e-13 ;
            elif (extraS > 6):
                print("Case not implemented")
        elif (extraS > 0):
            print("Case not implemented")
        
        self.Init(r, alpha)

    def Advance(self, tn, dt, Yold, sys):
        """ Method to compute X_n+1 from X_n with Linear Runge-Kutta
        Usage :
            Xnext = rk.Advance(t, dt, Xn, sys)
        
        t : current time t_n
        dt : time step
        Xn : iterate X_n
        sys : object describing the linear ODE 
        the method returns X_n+1

        the object containing the linear ODE must contain ApplyM, SolveM, ApplyK
        (multiplication by M, M^{-1} and K)
        GetSource (computation of F(t))
        """
        m = self.omega_i.shape[1]
        EvalF = [0]*m
        for i in range(m):
            tcurrent = tn + self.ci[i]*dt
            EvalF[i] = sys.GetSource(tcurrent)
            
        Y = Yold.copy() 
        if (self.stable_algo):
            # stable algorithm
            num_source = 0
            for i in range(len(self.complex_roots)):
                b = -2.0*dt*real(1.0/self.complex_roots[i])
                a = dt*dt / abs(self.complex_roots[i])**2
                Fn = a*sys.ApplyK(Y)

                for i in range(m):
                    Fn += dt*dt*self.omega_i[num_source, i]*EvalF[i]

                Fn = sys.SolveM(Fn)
                Fn = sys.ApplyK(b*Y + Fn)

                for i in range(m):
                    Fn += dt*self.omega_i[num_source+1, i]*EvalF[i]

                Fn = sys.SolveM(Fn)
                Y = Y + Fn
                num_source += 2

            for i in range(len(self.real_roots)):
                b = -dt / self.real_roots[i]
                Fn = b*sys.ApplyK(Y)
                for i in range(m):
                    Fn += dt*self.omega_i[num_source, i]*EvalF[i]

                Fn = sys.SolveM(Fn)
                Y = Y + Fn
                num_source += 1
        else:
            # Horner algorithm
            # We evaluate U_{n+1} = Pol*U^n
            last_coef = self.Numer.order
    
            Fn = self.Numer[last_coef]*Y
            for k in range(last_coef-1, -1, -1):
                Fn = sys.ApplyK(Fn)
                if (k < self.omega_i.shape[0]):
                    for i in range(self.omega_i.shape[1]):
                        Fn += self.omega_i[k, i]*EvalF[i];
                
                Fn = sys.SolveM(Fn);
                Fn = dt*Fn + self.Numer[k]*Y	
            
            # updating U^n to the next value
            Y = Fn.copy()
        
        return Y

class ExampleOde:
    """
    Example of a class defining a linear ode that can be used
    by PadeScheme.Advance or LinearSdirkScheme.Advance 

    The considered ode is M dU/dt = K U(t) + h(t) F
    where M and K are dense matrices, F a vector and h a function    
    """
    
    # constructor of the class ExampleOde
    def __init__(self, M, K, F, h):
        self.M = M
        self.K = K
        self.F = F
        self.h = h
        
    def ApplyM(self, X):
        return dot(self.M, X)
        
    def SolveM(self, X):
        return linalg.solve(self.M, X)
        
    def ApplyK(self, X):
        return dot(self.K, X)
        
    def SolveReal(self, a, X):
        return linalg.solve(self.M + a*self.K, X)
        
    def SolveComplex(self, a, X):
        return linalg.solve(self.M + a*self.K, X)
        
    def GetSource(self, t):
        return self.h(t)*self.F

class ScalarOde:
    # constructor of the class ExampleOde
    def __init__(self, L, f):
        self.L = L
        self.f = f
        
    def ApplyM(self, X):
        return X
        
    def SolveM(self, X):
        return X
        
    def ApplyK(self, X):
        return self.L*X
        
    def SolveReal(self, a, X):
        return X/(1.0 + a*self.L)
        
    def SolveComplex(self, a, X):
        return X/(1.0 + a*self.L)
        
    def GetSource(self, t):
        return self.f(t)

def SolveOde(t0, tf, N, Y0, scheme, sys):
    """
    Returns the final solution y(tf) solving the following Cauchy problem :
     M dy/dt = K y + F(t)
     y(t0) = y0

    Usage : yfinal = SolveOde(t0, tf, N, Y0, pade, sys)
    
    t0 : initial time
    tf : final time
    N : number of iterations
    Y0 : initial condition
    scheme : which scheme to use (an instance of class PadeScheme or LinearSdirkScheme)
    sys : which ode to consider (e.g. an instance of ExampleOde)
    """
    Y = Y0.copy()
    t = t0
    dt = (tf - t0) / N
    for i in range(N):
        t = t0 + i*dt
        Y = scheme.Advance(t, dt, Y, sys)
        
    return Y

def AdvancePoly(xn, f, L, R, tabR, ci, tn, dt):
    xnext = R(L*dt)*xn
    for i in range(len(ci)):
        Fn = f(tn + ci[i]*dt)
        xnext = xnext + dt*tabR[i](L*dt)*Fn

    return xnext

def SolvePoly(x0, f, L, R, tabR, ci, t0, tf, N):
    Y = zeros(N+1) + 1j*zeros(N+1)
    Y[0] = x0
    t = t0
    dt = (tf - t0) / N
    for i in range(N):
        t = t0 + i*dt
        Y[i+1] = AdvancePoly(Y[i], f, L, R, tabR, ci, t, dt)
    
    return Y

from scipy.integrate import *

def SolveScalarOde(t0, tf, N, Y0, scheme, sys):
    """
    Returns the final solution y(tf) solving the following Cauchy problem :
     dy/dt = A y + F(t)
     y(t0) = y0

    Usage : yfinal = SolveOde(t0, tf, N, Y0, pade, sys)
    
    t0 : initial time
    tf : final time
    N : number of iterations
    Y0 : initial condition
    scheme : which scheme to use (an instance of class PadeScheme or LinearSdirkScheme)
    sys : which ode to consider (e.g. an instance of ExampleOde)
    """
    Y = Y0.copy()
    tabY = [Y0[0]]
    tabZ = [Y0[1]]
    t = [t0]
    dt = (tf - t0) / N
    for i in range(N):
        t.append(t0 + (i+1)*dt)
        Ynext = scheme.Advance(t[i], dt, Y, sys)
        tabY.append(Ynext[0])
        tabZ.append(Ynext[1])
        Y = Ynext.copy();
    
    return array(t), array(tabY), array(tabZ)

def SolveExact(x0, f, L, t0, tf, N):
    Y = zeros(N+1) + 1j*zeros(N+1)
    t = t0
    dt = (tf - t0) / N
    Y[0] = x0
    for i in range(N):
        t = t0 + (i+1)*dt
        Yr = quad(lambda s : real(exp(L*(t-s))*f(s)), t0, t)[0]
        Yi = quad(lambda s : imag(exp(L*(t-s))*f(s)), t0, t)[0]
        Y[i+1] = exp(L*t)*x0 + Yr + 1j*Yi

    return Y

