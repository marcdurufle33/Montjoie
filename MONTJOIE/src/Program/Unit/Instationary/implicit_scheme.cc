#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

//! function f(t, y) = y/(1+t^2)
class ScalarFunctionFirst : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp gamma;
  Real_wp coef_mass;
  VectReal_wp gamma_real, alpha_real;
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);
    f(0) = y(0)/(one + t*t);
  }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0) {}

  void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c) 
  {
    coef_mass = a; gamma = b;
  }
  
  // solving y - gamma f(t, y) = f
  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y)
  {
    // y (coef_mass - gamma/(1+t^2)) = f
    Real_wp one(1);
    y(0) = f(0) / (coef_mass - gamma/(one +t*t));
  }
  
  void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta, const VectReal_wp& gamma)
  {
    alpha_real = alpha; gamma_real = gamma;
  }

  void SolveOperatorReal(const Real_wp& t, const VectReal_wp& x, VectReal_wp& y, int num)
  {
    Real_wp one(1);
    y(0) = x(0) / (alpha_real(num) - gamma_real(num)/(one +t*t));    
  }

};

Real_wp err_order(0.03);

Real_wp GetRelativeError(const Real_wp& x, const Real_wp& y)
{
  return (abs(x-y)/abs(y));
}

Real_wp GetRelativeError(const VectReal_wp& x, const VectReal_wp& y)
{
  Real_wp err = 0, sum = 0;
  for (int i = 0; i < x.GetM(); i++)
    {
      err += square(x(i) - y(i));
      sum += y(i)*y(i);
    }
  
  err = sqrt(err); sum = sqrt(sum);
  return err/sum;
}

template<class TimeScheme>
void TestFirstOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  // first, we solve dy/dt = y/(1+t^2)
  Real_wp Tmax(3), zero(0);  
  VectReal_wp y0(1);
  y0(0) = Real_wp(1);
  ScalarFunctionFirst fct;
  RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);
  
  Real_wp y_dt = y0(0);
  Real_wp yref = exp(atan(Tmax));
  
  y0(0) = Real_wp(1);
  RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
  
  Real_wp y_dt_half = y0(0);
  
  Real_wp err = abs(y_dt - yref)/abs(yref);
  Real_wp err_half = abs(y_dt_half - yref)/abs(yref);
  Real_wp order = (log(err) - log(err_half))/log(2.0);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  //DISP(err); DISP(err_half); DISP(r); DISP(order);
  if (order < r*(1.0-err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }
}

class ScalarFunctionSecond  : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp alpha, bb;
  Real_wp coef_mass, coef_damping, coef_stiff;
  
  ScalarFunctionSecond() { alpha = Real_wp(1)/10; bb = Real_wp(1); }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);
    f(0) = -y(1) - alpha*y(0);
    f(1) = bb*y(0)/(one + t*t);
  }

  void EvaluateFunctionS(const Real_wp& t, const VectReal_wp& y, const VectReal_wp& yprime, VectReal_wp& f,
			 bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);    
    f(0) = -bb*y(0)/(one + t*t);    
        
    if (alpha != Real_wp(0))
      f(0) -= alpha*yprime(0);
  }
  
  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0) {}

  void ApplyOperatorSh(const Real_wp& a, const Real_wp& t, const VectReal_wp& u, 
                       const Real_wp& beta, VectReal_wp& f)
  {
    f(0) = beta*f(0) + a*alpha*u(0);
  }
  
  void ApplyOperatorKh(const Real_wp& a, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& f)
  {
    Real_wp one(1);
    f(0) = beta*f(0) - bb*a*u(0)/(one + t*t);
  }
  
  void Assemble(VectReal_wp& u) const {}
  
  void SolveOperatorDh(VectReal_wp& Dh) {}
  
  void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c) 
  {
    coef_mass = a; coef_damping = b; coef_stiff = c;
  }
  
  void ApplyOperatorDh(const Real_wp& a, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& f)
  {
    f(0) = beta*f(0) + a*u(0);
  }
  
  void ApplyOperatorDhMinusdtSh(const Real_wp& a, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& f)
  {
    f(0) = beta*f(0) + a*(u(0) - alpha*coef_damping*u(0)); 
  }
  
  void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& f) {}
  
  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& u, VectReal_wp& f)
  {
    Real_wp one(1);
    f(0) = u(0) / (coef_mass + coef_damping*alpha + coef_stiff*bb/(one+t*t)); 
  }
  
  
};

template<class TimeScheme>
void TestSecondOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  Real_wp zero(0), Tmax(1);
  RungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(4);
      
  {
    ScalarFunctionSecond fct;
  
    // trying scalar problem
    VectReal_wp yinit(2);
    Real_wp yinit_0, yinit_p0;
    yinit_0 = Real_wp(2)/5; yinit_p0 = -Real_wp(4)/5;
    yinit_0 = Real_wp(1); yinit_p0 = -fct.alpha;
    
    yinit(0) = yinit_0; yinit(1) = -fct.alpha*yinit_0 - yinit_p0;
    VectReal_wp y_rk0(yinit);

    RunTimeScheme(zero, Tmax, dt/2, y_rk0, fct, rk);
    
    VectReal_wp y0(1), y0_prime(1);
    y0(0) = yinit_0; y0_prime(0) = yinit_p0;
    
    RunSecondOrderScheme(zero, Tmax, dt, y0, y0_prime, fct, scheme);
    
    VectReal_wp y_rk1(yinit);
    
    RunTimeScheme(zero, Tmax, dt/4, y_rk1, fct, rk);
    
    VectReal_wp y1(1), y1_prime(1);
    y1(0) = yinit_0; y1_prime(0) = yinit_p0;
    
    RunSecondOrderScheme(zero, Tmax, dt/2, y1, y1_prime, fct, scheme);
    
    Real_wp err = GetRelativeError(y0(0), y_rk0(0));
    Real_wp err_half = GetRelativeError(y1(0), y_rk1(0));
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    if (order < r*(1.0-err_order))
      {
        DISP(err); DISP(err_half); DISP(r); DISP(order);
        cout << "Scheme incorrect " << endl;
        abort();
      }  
  }
}


//! function f(t, y) = y + t, and M = 2
class LinearScalarFunctionFirst : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp K, M, invM;
  Real_wp gamma;
  Real_wp coef_mass;
  VectComplex_wp coef_facto;
  VectReal_wp coef_facto_real;
  Matrix<Real_wp> mat_lu; IVect pivot;
  
  LinearScalarFunctionFirst()
  {
    K = 0.5; M = 1.0; invM = 1.0;
    // K = 1.0; M = 2.0; invM = 0.5;
  }

  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    if (invert_mass)
      f(0) = invM*K*y(0);
    else
      f(0) = K*y(0);
    
    if (source)
      AddPrimitiveTimeSource(Real_wp(1), t, 0, f);
  }

  // f = beta f + alpha Kh y
  void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    f(0) = beta*f(0) + alpha*K*y(0);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    f(0) = beta*f(0) + M*alpha*y(0);
  }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0) {}

  void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c) 
  {
    coef_mass = M*a; gamma = b*K;
  }

  void SolveOperatorDh(Vector<Real_wp>& f)
  {
    f(0) *= invM;
  }

  void Assemble(Vector<Real_wp>& Vh) const{}
  
  void SolveMassMatrix(Vector<Real_wp>& X){}

  void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& x, 
		       const Real_wp& beta, Vector<Real_wp>& y){ y = beta*y + alpha*x;}

  //
  Real_wp GetFinalTimeSource() const {return Real_wp(100);}
  
  // solving y - gamma f(t, y) = f
  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y)
  {
    // y (coef_mass - gamma) = f
    y(0) = f(0) / (coef_mass - gamma);
  }


  void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				const VectComplex_wp& gamma)    
  {
    coef_facto = Real_wp(1) / (M*alpha - gamma*K);
  }
  
  void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			    int num_system)
  {
    Y(0) = X(0)*coef_facto(num_system);
  }

  void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
			     const VectReal_wp& gamma)    
  {
    coef_facto_real = Real_wp(1) / (M*alpha - gamma*K);
  }
  
  void SolveOperatorReal(const Real_wp& t, const VectReal_wp& X, VectReal_wp& Y,
			 int num_system)
  {
    Y(0) = X(0)*coef_facto_real(num_system);
  }

  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			      int nb_deriv, Vector<Real_wp>& b_src)
  {
    // we specify the source here (t)
    //b_src(0) += alpha*t;
    //b_src(0) += alpha*0.5*t;
    b_src(0) += alpha*exp(-t/4);
  }

  void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& Agauss, const Real_wp& dt, const VectReal_wp& Y,
				    const VectReal_wp& ki)
  {
    mat_lu = Agauss; mat_lu *= -K*dt;
    for (int i = 0; i < Agauss.GetM(); i++)
      mat_lu(i, i) += M;

    GetLU(mat_lu, pivot);
  }

  void SolveOperatorDhPlusAijKh(const VectReal_wp& tn, const Matrix<Real_wp>& A, const Real_wp& dt,
				const VectReal_wp& Y, VectReal_wp& ki)
  {
    VectReal_wp b(1);
    for (int i = 0; i < tn.GetM(); i++)
      {
	ki(i) = K*Y(0); b(0) = 0;
	AddPrimitiveTimeSource(1.0, tn(i), 0, b);
	ki(i) += b(0);
      }
    
    SolveLU(mat_lu, pivot, ki);
  }
  
};


//! function f(t, y) = M^-1 K y
class LinearVectorialFunctionFirst : public VirtualOdeSystem<Real_wp>
{
private:
  Matrix<Real_wp> M, K, invM;
  Vector<Real_wp> F;
  Matrix<Real_wp> facto_real;
  Vector<Matrix<Complex_wp> > facto_cplx;
  Vector<Matrix<Real_wp> > facto_real_vec;
  IVect pivot;
  
public :
  int GetM() { return 5; }
  
  LinearVectorialFunctionFirst()
  {
    srand(0);
    
    M.Reallocate(5, 5);
    K.Reallocate(5, 5);
    M.FillRand(); M *= 1e-9;
    //M.SetIdentity();
    K.FillRand(); K *= 1e-10;

    F.Reallocate(5);
    F.FillRand(); F *= 1e-9;
    
    //DISP(M); DISP(K);
    //M.WriteText("M.dat"); K.WriteText("K.dat");
    //F.WriteText("F.dat");
    invM = M;
    GetInverse(invM);
    
  }

  Real_wp GetSourceTime(const Real_wp& t)
  {
    return exp(-t/4);
  }


  VectReal_wp GetExactSolution(const Real_wp& Tmax, const VectReal_wp& yinit, int n)
  {
    Matrix<Real_wp> A;
    A.Reallocate(5, 5);
    Mlt(invM, K, A);  
    Matrix<Real_wp> B;
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(n);
    VectReal_wp yfinal(5), Ftilde(5), x(5);
    yfinal.Zero();
    Mlt(invM, F, Ftilde);
    for (int k = 0; k <= gauss.GetOrder(); k++)
      {
	Real_wp s = gauss.Points(k)*Tmax;
	B = A;
	B *= -s;
	GetExponential(B);
	Mlt(B, Ftilde, x);
	Real_wp poids = Tmax*gauss.Weights(k)*GetSourceTime(s);
	Add(poids, x, yfinal);
      }
    
    A *= Tmax;
    GetExponential(A);

    x = yinit + yfinal;
    Mlt(A, x, yfinal);
    return yfinal;
  }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    VectReal_wp Ky(5);
    Mlt(K, y, Ky);
    if(source)
      AddPrimitiveTimeSource(Real_wp(1), t, 0, Ky);
    
    if (invert_mass)
      Mlt(invM, Ky, f);
    else
      f = Ky;
  }

  // f = beta f + alpha Kh y
  void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    MltAdd(alpha, K, y, beta, f);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    MltAdd(alpha, M, y, beta, f);
  }

  void SolveOperatorDh(Vector<Real_wp>& f)
  {
    VectReal_wp y(f);
    Mlt(invM, y, f);
  }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0) {}

  void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c) 
  {
    facto_real = M;
    facto_real *= a;
    Add(-b, K, facto_real);
    GetInverse(facto_real);
  }

  Real_wp GetFinalTimeSource() const {return Real_wp(100);}

  void Assemble(Vector<Real_wp>& Vh) const{}

  void SolveMassMatrix(Vector<Real_wp>& X){ VectReal_wp Y(X); Mlt(invM, Y, X); }
  
  void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t, 
		       const Vector<Real_wp>& x, 
		       const Real_wp& beta, 
		       Vector<Real_wp>& y)
  { 
    MltAdd(alpha, M, x, beta, y); 
  }

  // solving y - gamma f(t, y) = f
  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y)
  {
    Mlt(facto_real, f, y);
  }


  void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				const VectComplex_wp& gamma)    
  {
    facto_cplx.Reallocate(alpha.GetM());
    for (int k = 0; k < alpha.GetM(); k++)
      {
	facto_cplx(k).Reallocate(M.GetM(), M.GetM());
	for (int i = 0; i < M.GetM(); i++)
	  for (int j = 0; j < M.GetM(); j++)
	    facto_cplx(k)(i, j) = Complex_wp(M(i, j), 0);
	
	facto_cplx(k) *= alpha(k);
	Add(-beta(k), K, facto_cplx(k));
	GetInverse(facto_cplx(k));    
      }
  }
  
  void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			    int num_system)
  {
    Mlt(facto_cplx(num_system), X, Y);
  }

  void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
				const VectReal_wp& gamma)    
  {
    facto_real_vec.Reallocate(alpha.GetM());
    for (int k = 0; k < alpha.GetM(); k++)
      {
	facto_real_vec(k).Reallocate(M.GetM(), M.GetM());
	for (int i = 0; i < M.GetM(); i++)
	  for (int j = 0; j < M.GetM(); j++)
	    facto_real_vec(k)(i, j) = M(i, j);
	
	facto_real_vec(k) *= alpha(k);
	Add(-beta(k), K, facto_real_vec(k));
	GetInverse(facto_real_vec(k));    
      }
  }
  
  void SolveOperatorReal(const Real_wp& t, const VectReal_wp& X, VectReal_wp& Y,
			 int num_system)
  {
    Mlt(facto_real_vec(num_system), X, Y);
  }

  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			      int nb_deriv, Vector<Real_wp>& b_src)
  {
    Real_wp pulse = GetSourceTime(t);
    b_src += alpha*pulse*F;
  }

  void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& Agauss, const Real_wp& dt, const VectReal_wp& Y,
				    const VectReal_wp& ki)
  {
    int N = K.GetM()*Agauss.GetM();
    facto_real.Reallocate(N, N);
    for (int i = 0; i < Agauss.GetM(); i++)
      for (int j = 0; j < Agauss.GetN(); j++)
	{
	  Real_wp coef = -Agauss(i, j)*dt;
	  int offset_row = i*K.GetM();
	  int offset_col = j*K.GetM();
	  for (int k = 0; k < K.GetM(); k++)
	    for (int l = 0; l < K.GetM(); l++)
	      facto_real(offset_row + k, offset_col + l) = coef*K(k, l);
	}
    
    for (int i = 0; i < Agauss.GetM(); i++)
      for (int k = 0; k < K.GetM(); k++)
	for (int l = 0; l < K.GetM(); l++)
	  facto_real(k + i*K.GetM(), l + i*K.GetM()) += M(k, l);
    
    GetLU(facto_real, pivot);
  }

  
  void SolveOperatorDhPlusAijKh(const VectReal_wp& tn, const Matrix<Real_wp>& A, const Real_wp& dt,
				const VectReal_wp& Y, VectReal_wp& ki_)
  {
    VectReal_wp ky(K.GetM()), kt(K.GetM()); ky.Zero();
    Mlt(K, Y, ky);

    VectReal_wp ki(ki_.GetM());
    for (int i = 0; i < tn.GetM(); i++)
      {
	Real_wp pulse = GetSourceTime(tn(i));
	kt = ky + pulse*F;
	for (int j = 0; j < K.GetM(); j++)
	  ki(j + i*K.GetM()) = kt(j);
      }
    
    SolveLU(facto_real, pivot, ki);

    for (int i = 0; i < tn.GetM(); i++)
      for (int j = 0; j < K.GetM(); j++)
	ki_(i + j*tn.GetM()) = ki(j + i*K.GetM());
  }
  
};

template<class TimeScheme>
void TestLinearFirstOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  // first, we solve 2 dy/dt = y + t  (scalar equation)
  Real_wp Tmax(3), zero(0);  

  {
    VectReal_wp y0(1);
    y0(0) = Real_wp(1);
    LinearScalarFunctionFirst fct;
    RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);
    
    Real_wp y_dt = y0(0);
    //Real_wp yref = (3.0*exp(0.5*Tmax) - (Tmax+2.0));
    Real_wp yref = Real_wp(7)/3*exp(0.5*Tmax) - Real_wp(4)/3*exp(-Tmax/4);
    //Real_wp yref = exp(0.5*Tmax);
    y0(0) = Real_wp(1);
    RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
    
    Real_wp y_dt_half = y0(0);
    
    Real_wp err = abs(y_dt - yref)/abs(yref);
    Real_wp err_half = abs(y_dt_half - yref)/abs(yref);
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    // DISP(err); DISP(err_half); DISP(r); DISP(order);
    if (order < r*(1.0-err_order))
      {
	DISP(err); DISP(err_half); DISP(r); DISP(order);
	cout << "Scheme incorrect " << endl;
	abort();
      }
  }
  
  // then, we solve M dy/dt = Ay + F (vectorial equation)
  {
    LinearVectorialFunctionFirst fct;
    VectReal_wp yinit(fct.GetM());
    yinit.FillRand(); yinit *= 1e-9;
    VectReal_wp y0 = yinit;
    RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);

    VectReal_wp y_dt = y0;
    VectReal_wp yref = fct.GetExactSolution(Tmax, yinit, 100);
    
    y0 = yinit;
    RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
    
    VectReal_wp y_dt_half = y0;
  
    Real_wp err = Norm2(y_dt - yref)/Norm2(yref);
    Real_wp err_half = Norm2(y_dt_half - yref)/Norm2(yref);
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    //DISP(err); DISP(err_half); DISP(r); DISP(order);
    if (order < r*(1.0-err_order))
      {
	DISP(err); DISP(err_half); DISP(r); DISP(order);
	cout << "Scheme incorrect " << endl;
	abort();
      }
  }
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  Real_wp threshold = 1e4*epsilon_machine;
  DISP(threshold);
  
  bool double_prec = false;
  if (threshold > 1e-16)
    double_prec = true;

  {
    cout << "Testing Linear Sdirk Scheme ..." << endl;
    LinearSdirkScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/2;
    if (double_prec)
      dt = 0.2;

    scheme.SetOrder(8, 3, LinearSdirkScheme_Iterator<Real_wp>::STABLE_ALTERNATIVE);
    //scheme.SetOrder(3, 0, LinearSdirkScheme_Iterator<Real_wp>::POLYNOMIAL);
    TestLinearFirstOrderScheme(scheme, 8, dt);
  }
  
  return FinalizeMontjoie();
  
  /* {
    cout << "Testing teta-scheme ..." << endl;
    TetaScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10000;
    if (double_prec)
      dt = 0.01;
    
    TestSecondOrderScheme(scheme, 2, dt);
    } */

  {
    cout << "Testing Pade Scheme ..." << endl;
    PadeScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10;
    if (double_prec)
      dt = 0.01;

    scheme.SetOrder(2);
    TestLinearFirstOrderScheme(scheme, 2, dt);

    scheme.SetOrder(4);
    TestLinearFirstOrderScheme(scheme, 4, dt);

    scheme.SetOrder(6);
    TestLinearFirstOrderScheme(scheme, 6, dt);
    
    scheme.SetOrder(8);
    TestLinearFirstOrderScheme(scheme, 8, dt);
  }
  
  {
    cout << "Testing Gauss-Runge-Kutta Scheme ..." << endl;
    GaussRungeKutta_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10;
    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(2);
    TestLinearFirstOrderScheme(scheme, 2, dt);

    scheme.SetOrder(4);
    TestLinearFirstOrderScheme(scheme, 4, dt);

    scheme.SetOrder(6);
    TestLinearFirstOrderScheme(scheme, 6, dt);
    
    scheme.SetOrder(8);
    TestLinearFirstOrderScheme(scheme, 8, dt);
  }  
  
  {
    cout << "Testing SDIRK schemes ..." << endl;
    SdirkScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10000;
    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(2);
    TestFirstOrderScheme(scheme, 2, dt);

    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt);

    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt);

    scheme.SetOrder(5);
    TestFirstOrderScheme(scheme, 5, dt);
  }

  {
    cout << "Testing DIRK schemes ..." << endl;
    DirkScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/1000;
    if (double_prec)
      dt = 0.01;
    
    // to be corrected
    scheme.SetOrder(4);
    TestLinearFirstOrderScheme(scheme, 2, dt);
  }

  {
    cout << "Testing Implicit Adams schemes ... " << endl;
    AdamsImplicit_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10000;
    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(1);
    TestFirstOrderScheme(scheme, 1, dt);    

    scheme.SetOrder(2);
    TestFirstOrderScheme(scheme, 2, dt);    

    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt);    

    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt);    
    
    if (double_prec)
      dt = 0.02;
    
    scheme.SetOrder(5);
    TestFirstOrderScheme(scheme, 5, dt);    

    scheme.SetOrder(6);
    TestFirstOrderScheme(scheme, 6, dt);    
  }

  {
    cout << "Testing Milne-Simpson schemes ..." << endl;
    MilneSimpson_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10000;
    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(1);
    TestFirstOrderScheme(scheme, 1, dt);

    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt);

    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt);

    scheme.SetOrder(5);
    TestFirstOrderScheme(scheme, 5, dt);

    scheme.SetOrder(6);
    TestFirstOrderScheme(scheme, 6, dt);
  }
  
  {
    cout << "Testing BDF schemes ..." << endl;
    BackwardDifferentiation_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10000;
    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(1);
    TestFirstOrderScheme(scheme, 1, dt);

    scheme.SetOrder(2);
    TestFirstOrderScheme(scheme, 2, dt);
    
    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt);
    
    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt);
    
    scheme.SetOrder(5);
    TestFirstOrderScheme(scheme, 5, dt);
    
    scheme.SetOrder(6);
    TestFirstOrderScheme(scheme, 6, dt);    
  }
  
  cout << "All tests passed successfully" << endl;

  return 0;
}
