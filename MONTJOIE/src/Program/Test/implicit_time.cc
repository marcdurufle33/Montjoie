#define MONTJOIE_WITH_TWO_DIM

#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

bool inhg_case = true;

//! function f(t, y) = y + t, and M = 2
class ScalarFunctionFirst : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp K, M, invM;
  Real_wp gamma;
  Real_wp coef_mass;
  VectComplex_wp coef_facto;
  
  ScalarFunctionFirst()
  {
    K = 0.5; M = 1.0; invM = 1.0;
    // K = 1.0; M = 2.0; invM = 0.5;
  }

  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    f(0) = K*y(0);
    if (source)
      AddPrimitiveTimeSource(Real_wp(1), t, 0, f);

    if (invert_mass)
      f(0) *= invM;
  }

  void EvaluateDerivativeFunction(const Real_wp& t, int n,
                                  const VectReal_wp& y, VectReal_wp& df,
                                  bool invert_mass = true, bool source = true)
  {
    df(0) = K*y(0);
    if (source)
      df(0) += pow(-Real_wp(1)/4, n)*exp(-t/4);
    
    if (invert_mass)
      df(0) *= invM;
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
    DISP(coef_facto);
  }
  
  void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			    int num_system)
  {
    Y(0) = X(0)*coef_facto(num_system);
  }

  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			      int nb_deriv, Vector<Real_wp>& b_src)
  {
    // we specify the source here (t)
    //b_src(0) += alpha*t;
    //b_src(0) += alpha*0.5*t;
    
    if (inhg_case)
      b_src(0) += alpha*pow(Real_wp(-1)/4, nb_deriv)*exp(-t/4);
  }

};


//! function f(t, y) = M^-1 K y
class VectorialFunctionFirst : public VirtualOdeSystem<Real_wp>
{
private:
  Matrix<Real_wp> M, K, invM;
  Vector<Real_wp> F;
  Matrix<Real_wp> facto_real;
  Vector<Matrix<Complex_wp> > facto_cplx;
    
public :
  int GetM() { return 5; }
  
  VectorialFunctionFirst()
  {
    M.Reallocate(5, 5);
    K.Reallocate(5, 5);
    M.FillRand(); M *= 1e-9;
    //M.SetIdentity();
    K.FillRand(); K *= 1e-10;

    F.Reallocate(5);
    F.FillRand(); F *= 1e-9;
    
    //DISP(M); DISP(K);
    M.WriteText("M.dat"); K.WriteText("K.dat");
    F.WriteText("F.dat");
    invM = M;
    GetInverse(invM);
    
  }

  Real_wp GetSourceTime(const Real_wp& t, int n)
  {
    return exp(-t/4)*pow(-Real_wp(1)/4, n);
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
    if (inhg_case)
      {
        Mlt(invM, F, Ftilde);
        for (int k = 0; k <= gauss.GetOrder(); k++)
          {
            Real_wp s = gauss.Points(k)*Tmax;
            B = A;
            B *= -s;
            GetExponential(B);
            Mlt(B, Ftilde, x);
            Real_wp poids = Tmax*gauss.Weights(k)*GetSourceTime(s, 0);
            Add(poids, x, yfinal);
          }
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

  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			      int nb_deriv, Vector<Real_wp>& b_src)
  {
    if (inhg_case)
      {
        Real_wp pulse = GetSourceTime(t, nb_deriv);
        b_src += alpha*pulse*F;
      }
    //DISP(t); DISP(pulse);
  }
  
};

Real_wp err_order(0.03);

void RunScalarScheme(const Real_wp& t0, const Real_wp& tf, const Real_wp& dt, VectReal_wp& y)
{
  Real_wp A = 0.5, Q(1);
  int nb_iter = toInteger(round(tf - t0) / dt);
  DISP(nb_iter); DISP(nb_iter*dt-tf);
  
  // Fourth-order
  /* Q = 1.0 - 0.5*dt*A + dt*dt/12*A*A;
  for (int n = 0; n < nb_iter; n++)
    {
      Real_wp t = n*dt+0.5*dt;
      Real_wp F = 0.5*t, dF = 0.5, F2 = 0;
      Real_wp phi = dt*(F + dt*dt/24*F2) - dt*dt*dt/12*A*dF;
      Real_wp b = dt*A*y(0) + phi;
      y(0) += b/Q;
    }
    return; */
  
  // Sixth-order
  Q = 1.0 - 0.5*dt*A + dt*dt/10*A*A - pow(dt*A, 3)/120;
  for (int n = 0; n < nb_iter; n++)
    {
      Real_wp t = n*dt+0.5*dt;
      //Real_wp F = 0.5*t, dF = 0.5, F2 = 0, F3 = 0, F4 = 0;
      Real_wp F = exp(-t/4), dF = -F/4, F2 = -dF/4, F3 = -F2/4, F4 = -F3/4;
      Real_wp phi = dt*(F + dt*dt/24*(F2 + dt*dt/80*F4)) - dt*dt*dt/12*A*(dF + dt*dt/40*F3) + dt*dt*dt/60*A*A*(F+dt*dt/8*F2); 
      //DISP(phi);
      Real_wp b = dt*A*(y(0) + dt*dt/60*A*A*y(0)) + phi; //DISP(b);
      y(0) += b/Q; //DISP(y);
      //int test_input; cout << "waiting" << endl; cin >> test_input;
    }
}

template<class TimeScheme>
void TestFirstOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  // first, we solve 2 dy/dt = y + t  (scalar equation)
  Real_wp Tmax(3), zero(0);  

  {
    VectReal_wp y0(1);
    y0(0) = Real_wp(1);
    ScalarFunctionFirst fct;
    RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);
    //RunScalarScheme(zero, Tmax, dt, y0);
    
    Real_wp y_dt = y0(0);
    //Real_wp yref = (3.0*exp(0.5*Tmax) - (Tmax+2.0));
    Real_wp yref(0);
    if (inhg_case)
      yref = Real_wp(7)/3*exp(0.5*Tmax) - Real_wp(4)/3*exp(-Tmax/4);
    else
      yref = exp(0.5*Tmax); // cas sans la source
    
    DISP(yref); DISP(y_dt);
    y0(0) = Real_wp(1);
    RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
    //RunScalarScheme(zero, Tmax, dt/2, y0);
    
    Real_wp y_dt_half = y0(0); DISP(y_dt_half);
    
    Real_wp err = abs(y_dt - yref)/abs(yref);
    Real_wp err_half = abs(y_dt_half - yref)/abs(yref);
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    DISP(err); DISP(err_half); DISP(r); DISP(order);
    if (order < r*(1.0-err_order))
      {
	DISP(err); DISP(err_half); DISP(r); DISP(order);
	cout << "Scheme incorrect " << endl;
	abort();
      }
  }
  
  // then, we solve M dy/dt = Ay + F (vectorial equation)
  {
    VectorialFunctionFirst fct;
    VectReal_wp yinit(fct.GetM());
    yinit.FillRand(); yinit *= 1e-9; yinit.WriteText("Y0.dat");
    VectReal_wp y0 = yinit; DISP(yinit);
    RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);
    
    VectReal_wp y_dt = y0;
    VectReal_wp yref = fct.GetExactSolution(Tmax, yinit, 100);
    y_dt.WriteText("Y.dat");
    
    y0 = yinit;
    RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
    
    VectReal_wp y_dt_half = y0;
  
    Real_wp err = Norm2(y_dt - yref)/Norm2(yref);
    Real_wp err_half = Norm2(y_dt_half - yref)/Norm2(yref);
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    DISP(err); DISP(err_half); DISP(r); DISP(order);
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

  if (argc < 4)
    {
      cout <<"Provide the parameter N, order and extra stages." << endl;
      cout << "Usage: ./test.x N order extraS" << endl;
    }

  Real_wp threshold = 1e4*epsilon_machine;
  DISP(threshold);
   
  int N = atoi(argv[1]);
  int r = atoi(argv[2]);
  int s = atoi(argv[3]);

  //PadeScheme_Iterator<Real_wp> scheme;
  //scheme.SetOrder(r);

  // cas series de Taylor
  //TaylorSeries_Iterator<Real_wp> scheme;
  //scheme.SetOrder(r);
  
  // cas Runge-Kutta lineaire
  LinearRungeKutta_Iterator<Real_wp> scheme;
  scheme.SetOrder(r, s, false);
  
  //LinearSdirkScheme_Iterator<Real_wp> scheme;
  //scheme.SetOrder(r, s, true);

   //LinearSdirkScheme_Iterator<Real_wp> scheme;
   //scheme.SetOrder(r, s, scheme.STABLE_DERIVATIVE);

  //LinearSdirkScheme_Iterator<Real_wp> scheme;
  //scheme.SetOrder(r, s, scheme.SDIRK);

  //LinearSdirkScheme_Iterator<Real_wp> scheme;
  //scheme.SetOrder(r, s, scheme.STABLE_WEIGHTS);
  
  Real_wp dt = Real_wp(3)/N;
  TestFirstOrderScheme(scheme, r, dt);
}

