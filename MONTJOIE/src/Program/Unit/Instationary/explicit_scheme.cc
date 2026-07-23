#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

//! function f(t, y) = y/(1+t^2)
class ScalarFunctionFirst : public VirtualOdeSystem<Real_wp>
{
public :
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);
    f(0) = y(0)/(one + t*t);
  }

  void EvaluateDerivativeFunction(const Real_wp& t, int n, const VectReal_wp& y, VectReal_wp& f,
				  bool invert_mass = true, bool source = true)
  {
    Real_wp one(1); f(0) = Real_wp(0);
    switch (n)
      {
      case 0 : f(0) = y(0)/(one + t*t); break;
      case 1 : f(0) = -Real_wp(2)*t*y(0)/square(one + t*t); break;
      case 2 : f(0) = (Real_wp(-2) + Real_wp(6)*t*t)*y(0) 
          / pow( one + t*t, Real_wp(3)); break;
      case 3 : f(0) = Real_wp(24)*t*(one -t*t)*y(0) / pow(one + t*t, Real_wp(4)); break;
      }    
  }
  
  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}
};

//! function f(t, y) = A y
class VectorialFunctionFirst  : public VirtualOdeSystem<Real_wp>
{
public :
  Matrix<Real_wp> A;
  
  VectorialFunctionFirst(Matrix<Real_wp>& A_) { A = A_; }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Mlt(A, y, f);
  }

  void EvaluateDerivativeFunction(const Real_wp& t, int n, const VectReal_wp& y, VectReal_wp& f,
				  bool invert_mass = true, bool source = true)
  {
    Mlt(A, y, f);
  }
  
  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}
};

#ifdef MONTJOIE_WITH_MPFR
Real_wp err_order(0.03, MONTJOIE_DEFAULT_PRECISION);
#else
Real_wp err_order(0.03);
#endif

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

//! function f(t, y) for the n-body system
class AutonomousFunctionFirst  : public VirtualOdeSystem<Real_wp>
{
public :
  VectReal_wp masse;
  Real_wp G;
  Vector<ofstream> trajectoire;
  
  AutonomousFunctionFirst()
  {
    int nb_body = 3;
    G = 1.0;
    masse.Reallocate(nb_body);
    masse(0) = 0.2;
    masse(1) = 1.0;
    masse(2) = 0.1;
    
    /* int nb_body = 6;
    G = 2.95912208286e-4;
    masse.Reallocate(nb_body);
    masse(0) = 0.000954786104043;
    masse(1) = 0.000285583733151;
    masse(2) = 0.0000437273164546;
    masse(3) = 0.0000517759138449;
    masse(4) = 1.0/1.3e8;
    masse(5) = 1.00000597682; 
    */
      
    /*trajectoire.Reallocate(nb_body);
    for (int i = 0; i < nb_body; i++)
      {
        string name = "P" + to_str(i) + ".txt";
        trajectoire(i).precision(15);
        trajectoire(i).open(name.data());        
        }*/
  }

  int GetM() { return 6*masse.GetM(); }
    
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& pos_vit, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    int N = masse.GetM();
    Vector<R3> pos(N);
    
    // imposing dpos_i/dt = vi, and retrieving positions
    int offset = 0;
    for (int i = 0; i < N; i++)
      {
        // position of body i
        pos(i).Init(pos_vit(offset), pos_vit(offset+1), pos_vit(offset+2));
        
        // dpos_i/dt = vi
        f(offset) = pos_vit(offset+3);
        f(offset+1) = pos_vit(offset+4);
        f(offset+2) = pos_vit(offset+5);
        
        offset += 6;
      }
    
    // imposing dvi/dt = -G \sum_j m_j (xi - xj) / || xi - xj||^3
    R3 diff_x, acc_i; Real_wp r2, r3;
    offset = 0;
    for (int i = 0; i < N; i++)
      {
        acc_i.Fill(0);
        for (int j = 0; j < N; j++)
          if (j != i)
            {
              diff_x = pos(i) - pos(j);
              r2 = diff_x(0)*diff_x(0) + diff_x(1)*diff_x(1) + diff_x(2)*diff_x(2);
              r3 = r2*sqrt(r2);
              Add(-masse(j)/r3, diff_x, acc_i);
            }
        
        //if (i == N-1)
        //acc_i.Fill(0);
        
        f(offset+3) = G*acc_i(0);
        f(offset+4) = G*acc_i(1);
        f(offset+5) = G*acc_i(2);

        offset += 6;
      }
  }

  void EvaluateFunctionS(const Real_wp& t, const VectReal_wp& pos, const VectReal_wp& vit, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    int N = masse.GetM();
    Vector<R3> x(N);
    
    int offset = 0;
    for (int i = 0; i < N; i++)
      {
        // position of body i
        x(i).Init(pos(offset), pos(offset+1), pos(offset+2));
        offset += 3;
      }
    
    // f = -G \sum_j m_j (xi - xj) / || xi - xj||^3
    R3 diff_x, acc_i; Real_wp r2, r3;
    offset = 0;
    for (int i = 0; i < N; i++)
      {
        acc_i.Fill(0);
        for (int j = 0; j < N; j++)
          if (j != i)
            {
              diff_x = x(i) - x(j);
              r2 = diff_x(0)*diff_x(0) + diff_x(1)*diff_x(1) + diff_x(2)*diff_x(2);
              r3 = r2*sqrt(r2);
              Add(-masse(j)/r3, diff_x, acc_i);
            }
        
        //if (i == N-1)
        //acc_i.Fill(0);
        
        f(offset) = G*acc_i(0);
        f(offset+1) = G*acc_i(1);
        f(offset+2) = G*acc_i(2);

        offset += 3;
      }
  }

  void Extract(const VectReal_wp& y, VectReal_wp& y0, VectReal_wp& y0_prime)
  {
    int N = y.GetM()/6;
    int offset = 0;
    y0.Reallocate(3*N); y0_prime.Reallocate(3*N);
    for (int i = 0; i < N; i++)
      {
        y0(3*i) = y(offset);
        y0(3*i+1) = y(offset+1);
        y0(3*i+2) = y(offset+2);

        y0_prime(3*i) = y(offset+3);
        y0_prime(3*i+1) = y(offset+4);
        y0_prime(3*i+2) = y(offset+5);
        
        offset += 6;
      }
  }

  void Regroup(const VectReal_wp& y0, const VectReal_wp& y0_prime, VectReal_wp& y, VectReal_wp& yref)
  {
    int N = y0.GetM()/3;
    int offset = 0;
    y.Reallocate(6*N);
    y.Fill(0);
    for (int i = 0; i < N; i++)
      {
        y(offset) = y0(3*i);
        y(offset+1) = y0(3*i+1);
        y(offset+2) = y0(3*i+2);
        
	if (y0_prime.GetM() > 0)
	  {
	    y(offset+3) = y0_prime(3*i);
	    y(offset+4) = y0_prime(3*i+1);
	    y(offset+5) = y0_prime(3*i+2);
	  }
	else
	  {
	    yref(offset+3) = 0;
	    yref(offset+4) = 0;
	    yref(offset+5) = 0;
	  }
	
        offset += 6;
      }
  }

  void EvaluateDerivativeFunction(const Real_wp& t, int n, const VectReal_wp& y, VectReal_wp& f,
				  bool invert_mass = true, bool source = true)
  {
    f.Fill(0);
  }


  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void GiveIterate(int n, const Real_wp& t, VectReal_wp& Y) 
  {
    /*int N = masse.GetM();
    int offset = 0;
    for (int i = 0; i < N; i++)
      {
        trajectoire(i) << Y(offset) << " " << Y(offset+1) << " " << Y(offset+2) << endl;
        offset += 6;
        } */
  }
  
  void GiveFinalIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
    //for (int i = 0; i < masse.GetM(); i++)
    //trajectoire(i).close();
  }
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}

};

class LinearFunctionFirst  : public VirtualOdeSystem<Real_wp>
{
public :
  bool source_taken;
  
  LinearFunctionFirst()
  {
    source_taken = true;
  }
  
  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& f) 
  {
    Real_wp one(1), coef_x = Real_wp(3)/10, coef_y = Real_wp(35)/52;
    Real_wp gaussian = alpha*exp(-Real_wp(36)*square(t-one));
    f(0) += coef_x*gaussian;
    f(1) += coef_y*gaussian;    
  }
  
  void SolveMassMatrix(VectReal_wp&) {}
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
                        bool invert_mass = true, bool add_source = true)
  {
    Real_wp one(1), coef_x = Real_wp(3)/10, coef_y = Real_wp(35)/52;
    f(0) = y(1);
    f(1) = -y(0);
    if (add_source && source_taken)
      {
        Real_wp gaussian = exp(-Real_wp(36)*square(t-one));
        f(0) += coef_x*gaussian;
        f(1) += coef_y*gaussian;
      }
  }

  void EvaluateDerivativeFunction(const Real_wp& t, int n, const VectReal_wp& y, VectReal_wp& f,
				  bool invert_mass = true, bool source = true)
  {
    Real_wp one(1), coef_x = Real_wp(3)/10, coef_y = Real_wp(35)/52;
    f(0) = y(1);
    f(1) = -y(0);
    
    if (!source_taken)
      return;
    
    Real_wp gaussian = exp(-Real_wp(36)*square(t-one));
    Real_wp pulse(gaussian);
    
    switch (n)
      {
      case 1 : pulse = -Real_wp(72)*(t-one)*gaussian; break;
      case 2 : pulse = Real_wp(72)*(Real_wp(72)*square(t-one) - one)*gaussian; break;
      case 3 : pulse = (-Real_wp(5184)*(Real_wp(72)*square(t-one) - one)*(t-one) 
                        + Real_wp(72)*(Real_wp(144)*(t-one)))*gaussian; break;
      }

    f(0) += coef_x*pulse;
    f(1) += coef_y*pulse;
  }
  
  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  // returns the frequency associated with the source h
  Real_wp GetFrequency() const { return Real_wp(1); }
  
  // computes pulse = d^n h/dt (t), n is always equal to 0
  void SourceOnlyTime(const Real_wp& t, int n, Real_wp& pulse)
  {
    Real_wp one(1);
    pulse = exp(-Real_wp(36)*square(t-one));
    //pulse = Real_wp(0);
  }
  
  // sets Dirichlet condition u_i = f(t)
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0)
  {
  }
  
  // returns the final time after which the source h(t) is vanishing
  Real_wp GetFinalTimeSource() const { return Real_wp(2); }
  
  // fills a vector with B
  void FillSource(VectReal_wp& b) const
  {
    Real_wp coef_x = Real_wp(3)/10, coef_y = Real_wp(35)/52;
    b(0) = coef_x; b(1) = coef_y;
  }
  
  // inversion by the mass matrix : replaces b by M^{-1} b
  // in pratice the system considered is M dy/dt = h(t) B + A y
  // M being the mass matrix, therefore we need to compute M^-1 B
  void SolveOperatorDh(VectReal_wp& b) { }
  
};


template<class TimeScheme>
void TestLinearScheme(TimeScheme& scheme, int r, Real_wp dt,
                      const Real_wp& threshold, bool check_order = true, bool add_source = true)
{
  // testing linear scheme of the type
  // dy/dt = h(t) B + A y
  // with A a matrix and B a vector
  Real_wp Tmax(5), zero(0);
  VectReal_wp y0(2);
  LinearFunctionFirst fct;
  fct.source_taken = add_source;
  
  y0(0) = Real_wp(1);
  y0(1) = Real_wp(1)/2;
  VectReal_wp yinit(y0);
  
  RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);

  if (check_order)
    {
      VectReal_wp y_0(y0);
      
      y0 = yinit;
      RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
      VectReal_wp y1(y0);
      
      y0 = yinit;
      RunTimeScheme(zero, Tmax, dt/4, y0, fct, scheme);
      VectReal_wp y2(y0);
      
      Real_wp err = GetRelativeError(y_0, y1);
      Real_wp err_half = GetRelativeError(y1, y2);
      
      DISP(err);
      Real_wp order = (log(err) - log(err_half))/log(2.0);
      cout << "Order expected : " << r << ", Order measured : " << order << endl;
      if (order < r*(1.0- err_order))
        {
          DISP(err); DISP(err_half); DISP(r); DISP(order);
          cout << "Scheme incorrect " << endl;
          abort();
        }
    }
  else
    {
      // testing with a Runge-Kutta scheme
      RungeKutta_Iterator<Real_wp> rk;
      rk.SetOrder(6);
      
      VectReal_wp y_rk(yinit);
      RunTimeScheme(zero, Tmax, dt/100, y_rk, fct, rk);
      
      if (GetRelativeError(y0, y_rk) > threshold)
        {          
          DISP(y0); DISP(y_rk);
          DISP(GetRelativeError(y0, y_rk));

          cout << "scheme incorrect" << endl;
          abort();
        }
    }

}

template<class TimeScheme>
void TestFirstOrderScheme(TimeScheme& scheme, int r, Real_wp dt, Real_wp dt_bis)
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
  if (order < r*(1.0-err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }
  
  // then, a linear system dy/dt = A y
  Matrix<Real_wp> A(4, 4);
  A.FillRand(); Mlt(Real_wp(1)/RAND_MAX, A);
  for (int i = 0; i < A.GetM(); i++)
    {
      for (int j = i+1; j < A.GetM(); j++)
        A(i, j) = -A(j, i);
      
      A(i, i) = 0.0;
    }
    
  VectorialFunctionFirst fctA(A);
  y0.Reallocate(A.GetM());
  y0.FillRand(); Mlt(Real_wp(1)/RAND_MAX, y0);
  VectReal_wp yinit(y0);

  VectReal_wp Lreal, Limag; Matrix<Real_wp> P; Vector<int> pivot;
  Matrix<Complex_wp> Ptrue(A.GetM(), A.GetM()), Plu;
  GetEigenvaluesEigenvectors(A, Lreal, Limag, P);  
  int num = 0;
  while (num < A.GetM())
    {
      if (Limag(num) == 0)
        {
          for (int j = 0; j < A.GetM(); j++)
            Ptrue(j, num) = P(j, num);
          
          num++;
        }
      else
        {
          for (int j = 0; j < A.GetM(); j++)
            {
              Ptrue(j, num) = Complex_wp(P(j, num), P(j, num+1));
              Ptrue(j, num+1) = Complex_wp(P(j, num), -P(j, num+1));
            }
          
          num += 2;
        }
    }
  
  VectComplex_wp ytmp(A.GetM()), ysol(A.GetM());
  Plu = Ptrue; GetLU(Plu, pivot);
  
  // exact solution
  VectReal_wp y_exact(y0);  
  for (int i = 0; i < A.GetM(); i++)
    ytmp(i) = Complex_wp(y0(i), zero);
  
  SolveLU(Plu, pivot, ytmp);
  
  for (int i = 0; i < ytmp.GetM(); i++)
    ytmp(i) = exp(Complex_wp(Lreal(i), Limag(i))*Tmax)*ytmp(i);
  
  ysol.Fill(Complex_wp(zero, zero));
  Mlt(Ptrue, ytmp, ysol);
  for (int i = 0; i < A.GetM(); i++)
    y_exact(i) = real(ysol(i));
  
  // numerical solution for dt
  RunTimeScheme(zero, Tmax, dt, y0, fctA, scheme);
  
  err = GetRelativeError(y0, y_exact);

  // numerical solution for dt/2
  y0 = yinit;
  RunTimeScheme(zero, Tmax, dt/2, y0, fctA, scheme);
  
  err_half = GetRelativeError(y0, y_exact);

  order = (log(err) - log(err_half))/log(2.0);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  if (order < r*(1.0-err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }
  
  // then, an autonomous system dy/dt = f(y)
  AutonomousFunctionFirst fct_auto;
  int N = fct_auto.GetM();
  
  y0.Reallocate(N); y0.Fill(0);
  /* y0(0) = -3.5023653; y0(1) = -3.8169847; y0(2) = -1.5507963;
     y0(3) = 0.00565429; y0(4) = -0.00412490; y0(5) = -0.00190589;
     y0(6) = 9.0755314; y0(7) = -3.0458353; y0(8) = -1.6483708;
     y0(9) = 0.00168318; y0(10) = 0.00483525; y0(11) = 0.00192462;
     y0(12) = 8.3101420; y0(13) = -16.2901086; y0(14) = -7.2521278;
     y0(15) = 0.00354178; y0(16) = 0.00137102; y0(17) = 0.00055029;
     y0(18) = 11.4707666; y0(19) = -25.7294829; y0(20) = -10.8169456;
     y0(21) = 0.00288930; y0(22) = 0.00114527; y0(23) = 0.00039677;
     y0(24) = -15.5387357; y0(25) = -25.2225594; y0(26) = -3.1902382;
     y0(27) = 0.00276725; y0(28) = -0.00170702; y0(29) = -0.00136504; */
  
  y0(0) = 1.0; y0(1) = 0.0; y0(2) = 0.0;
  y0(3) = 0.0; y0(4) = 1.1; y0(5) = 0.04;
  y0(6) = 0.0; y0(7) = 0.0; y0(8) = 0.0;
  y0(9) = 0.0; y0(10) = 0.0; y0(11) = 0.0;
  y0(12) = 0.5; y0(13) = 0.8; y0(14) = 0.2;
  y0(15) = -0.04; y0(16) = 0.95; y0(17) = 0.2;

  Tmax = Real_wp(1);
  dt = dt_bis;
  yinit = y0;
  RunTimeScheme(zero, Tmax, dt, y0, fct_auto, scheme);
  VectReal_wp y_0(y0);
  
  y0 = yinit;
  RunTimeScheme(zero, Tmax, dt/2, y0, fct_auto, scheme);
  VectReal_wp y1(y0);

  y0 = yinit;
  RunTimeScheme(zero, Tmax, dt/4, y0, fct_auto, scheme);
  VectReal_wp y2(y0);
  
  err = GetRelativeError(y_0, y1);
  err_half = GetRelativeError(y1, y2);
  
  order = (log(err) - log(err_half))/log(2.0);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  if (order < r*(1.0- err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }
}

class ScalarFunctionSecond  : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp alpha;
  
  ScalarFunctionSecond() { alpha = Real_wp(0); }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);
    f(0) = -y(1) - alpha*y(0);
    f(1) = y(0)/(one + t*t);
  }

  void EvaluateFunctionS(const Real_wp& t, const VectReal_wp& y, const VectReal_wp& yprime, VectReal_wp& f,
			 bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);    
    f(0) = -y(0)/(one + t*t);    
        
    if (alpha != Real_wp(0))
      f(0) -= alpha*yprime(0);
  }
  
  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}
};


template<class TimeScheme>
void TestSecondOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  Real_wp zero(0), Tmax(1);
  RungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(min(r, 8));
      
  {
    ScalarFunctionSecond fct;
  
    // trying scalar problem
    VectReal_wp yinit(2);
    yinit(0) = Real_wp(2)/5; yinit(1) = -Real_wp(4)/5;
    
    VectReal_wp y_rk0(yinit);

    if (r > 8)
      RunTimeScheme(zero, Tmax, dt/20, y_rk0, fct, rk);
    else
      RunTimeScheme(zero, Tmax, dt/2, y_rk0, fct, rk);
    
    VectReal_wp y0(1), y0_prime(1);
    y0(0) = yinit(0); y0_prime(0) = -yinit(1);
    
    RunSecondOrderScheme(zero, Tmax, dt, y0, y0_prime, fct, scheme);
    
    VectReal_wp y_rk1(yinit);
    
    if (r > 8)
      RunTimeScheme(zero, Tmax, dt/50, y_rk1, fct, rk);
    else
      RunTimeScheme(zero, Tmax, dt/4, y_rk1, fct, rk);
    
    VectReal_wp y1(1), y1_prime(1);
    y1(0) = yinit(0); y1_prime(0) = -yinit(1);
    
    RunSecondOrderScheme(zero, Tmax, dt/2, y1, y1_prime, fct, scheme);
    
    VectReal_wp y0_(2), y1_(2);
    y0_.Fill(0); y1_.Fill(0);
    y0_(0) = y0(0); 
    if (y0_prime.GetM() > 0)
      y0_(1) = -y0_prime(0);
    else
      y_rk0(1) = 0.0;
    
    y1_(0) = y1(0);
    if (y1_prime.GetM() > 0)
      y1_(1) = -y1_prime(0);
    else
      y_rk1(1) = 0.0;
    
    Real_wp err = GetRelativeError(y0_, y_rk0);
    Real_wp err_half = GetRelativeError(y1_, y_rk1);
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    DISP(err); DISP(err_half);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    if (order < r*(1.0-err_order))
      {
        DISP(err); DISP(err_half); DISP(r); DISP(order);
        cout << "Scheme incorrect " << endl;
        abort();
      }  
  }
  
  {
    // then, an autonomous system d^2 y/dt^2 = f(y)
    AutonomousFunctionFirst fct_auto;
    int N = fct_auto.GetM();
    
    VectReal_wp yinit;
    yinit.Reallocate(N); yinit.Fill(0);
    yinit(0) = 1.0; yinit(1) = 0.0; yinit(2) = 0.0;
    yinit(3) = 0.0; yinit(4) = 1.1; yinit(5) = 0.04;
    yinit(6) = 0.0; yinit(7) = 0.0; yinit(8) = 0.0;
    yinit(9) = 0.0; yinit(10) = 0.0; yinit(11) = 0.0;
    yinit(12) = 0.5; yinit(13) = 0.8; yinit(14) = 0.2;
    yinit(15) = -0.04; yinit(16) = 0.95; yinit(17) = 0.2;
    
    VectReal_wp y_rk0(yinit);
    if (r > 8)
      RunTimeScheme(zero, Tmax, dt/20, y_rk0, fct_auto, rk);
    else
      RunTimeScheme(zero, Tmax, dt/2, y_rk0, fct_auto, rk);
    
    VectReal_wp y_rk1(yinit);
    if (r > 8)
      RunTimeScheme(zero, Tmax, dt/40, y_rk1, fct_auto, rk);
    else
      RunTimeScheme(zero, Tmax, dt/4, y_rk1, fct_auto, rk);
    
    VectReal_wp y0, y0_prime;
    fct_auto.Extract(yinit, y0, y0_prime);

    RunSecondOrderScheme(zero, Tmax, dt, y0, y0_prime, fct_auto, scheme);

    VectReal_wp y1, y1_prime;
    fct_auto.Extract(yinit, y1, y1_prime);

    RunSecondOrderScheme(zero, Tmax, dt/2, y1, y1_prime, fct_auto, scheme);
    
    VectReal_wp y0_, y1_;
    fct_auto.Regroup(y0, y0_prime, y0_, y_rk0);
    fct_auto.Regroup(y1, y1_prime, y1_, y_rk1);
  
    Real_wp err = GetRelativeError(y0_, y_rk0);
    Real_wp err_half = GetRelativeError(y1_, y_rk1);
    
    Real_wp order = (log(err) - log(err_half))/log(2.0);
    cout << "Order expected : " << r << ", Order measured : " << order << endl;
    if (order < r*(1.0- err_order))
      {
        DISP(err); DISP(err_half); DISP(r); DISP(order);
        cout << "Scheme incorrect " << endl;
        abort();
      }
  }
  
}

void GenerateSdpMatrix(const Real_wp& coef, Matrix<Real_wp>& A)
{
  Real_wp zero(0), one(1);
  int N = A.GetM();
  Matrix<Real_wp> B(N, N);
  B.FillRand(); Mlt(coef, B);
  A.Fill(zero);
  MltAdd(one, SeldonNoTrans, B, SeldonTrans, B, zero, A);
}

class LinearFunctionSecond  : public VirtualOdeSystem<Real_wp>
{
public :
  Matrix<Real_wp> Dh, Kh, Sh, Dh_lu, DhMinusdtSh, DhPlusdtSh;
  IVect pivot_Dh, pivot_DhPlusdtSh;
  VectReal_wp F;
  
  LinearFunctionSecond()
  {
    Real_wp one(1);
    int N = 4;
    F.Reallocate(N); Dh.Reallocate(N, N); Kh.Reallocate(N, N); Sh.Reallocate(N, N);
    Real_wp coef = one/RAND_MAX;
    F.FillRand(); Mlt(coef, F);
    
    GenerateSdpMatrix(coef, Dh);
    for (int i = 0; i < Dh.GetM(); i++)
      Dh(i, i) += one;
    
    GenerateSdpMatrix(coef, Kh);
    GenerateSdpMatrix(coef, Sh);
    Mlt(-one, Kh);
    
    Dh_lu = Dh;
    GetLU(Dh_lu, pivot_Dh);
  }
  
  int GetM() { return Kh.GetM(); }
  
  void Init(const Real_wp& dt)
  {
    DhMinusdtSh = Dh;
    Add(-dt/2, Sh, DhMinusdtSh); 
    
    DhPlusdtSh = Dh;
    Add(dt/2, Sh, DhPlusdtSh); 
    
    GetLU(DhPlusdtSh, pivot_DhPlusdtSh);
  }
  
  void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& X)
  {
    Real_wp one(1), beta(50);
    Real_wp gaussian = exp(-beta*square(t-one));    
    Real_wp pulse = gaussian;
    if (n == 2)
      pulse = -2*beta*(t-one)*gaussian;
    else if (n == 3)
      pulse = -2*beta*gaussian + square(2*beta*(t-one))*gaussian;
    else if (n > 3)
      {
        abort();
      }

    pulse *= alpha;
    for (int i = 0; i < F.GetM(); i++)
      X(i) += pulse*F(i);
  }
  
  void Extract(const VectReal_wp& Y, VectReal_wp& U, VectReal_wp& V)
  {
    int N = Kh.GetM();
    U.Reallocate(N); V.Reallocate(N);
    for (int i = 0; i < N; i++)
      {
        U(i) = Y(i);
        V(i) = Y(N+i);
      }
  }
  
  void Regroup(const VectReal_wp& U, const VectReal_wp& V, VectReal_wp& Y)
  {
    int N = Kh.GetM();
    Y.Reallocate(2*N);
    for (int i = 0; i < N; i++)
      {
        Y(i) = U(i);
        Y(N+i) = V(i);
      }
  }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& Y, VectReal_wp& ProdY,
			bool invert_mass = true, bool source = true)
  {
    Real_wp zero(0), one(1);
    int N = Kh.GetM();
    VectReal_wp U, V, ProdU(N), ProdV(N);
    ProdU.Fill(zero); ProdV.Fill(zero);
    Extract(Y, U, V);

    // Dh dV/dt - Kh U + Sh V = F
    Mlt(Kh, U, ProdV);
    MltAdd(-one, Sh, V, one, ProdV);
    if (source)
      AddScalarTimeSource(one, t, 1, ProdV);
    
    SolveLU(Dh_lu, pivot_Dh, ProdV);
    
    // dU/dt = V
    Copy(V, ProdU);
    
    Regroup(ProdU, ProdV, ProdY);
  }
  
  void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, Sh, u, beta, v);
  }

  void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, Kh, u, beta, v);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, Dh, u, beta, v);
  }

  void Assemble(VectReal_wp& u) const {}
  
  void SolveOperatorDh(VectReal_wp& u) { SolveLU(Dh_lu, pivot_Dh, u); }
  
  void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, DhMinusdtSh, u, beta, v);
  }
  
  void SolveOperatorDhPlusdtSh(VectReal_wp& u) { SolveLU(DhPlusdtSh, pivot_DhPlusdtSh, u); }

  //void SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& u, VectReal_wp& v)
  // {
  //  VectReal_wp x(u);
  //  SolveLU(DhPlusdtSh, pivot_DhPlusdtSh, x);
  //  Add(alpha, x, v);
  // }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}

};

template<class TimeScheme>
void TestLinearSecond(TimeScheme& scheme, int r, Real_wp dt)
{
  RungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(4);
  Real_wp one(1), zero(0), Tmax(2);
  Real_wp coef = one/RAND_MAX;
  Tmax = Real_wp(1);
  
  {
    LinearFunctionSecond fct;
    
    VectReal_wp yinit(2*fct.GetM());
    yinit.FillRand(); Mlt(coef, yinit);
    
    VectReal_wp y_rk0(yinit);
        
    RunTimeScheme(zero, Tmax, dt/2, y_rk0, fct, rk);
    
    VectReal_wp y0, y0_prime;
    fct.Extract(yinit, y0, y0_prime);
    fct.Init(dt);
    
    RunSecondOrderScheme(zero, Tmax, dt, y0, y0_prime, fct, scheme);
    
    VectReal_wp y_rk1(yinit);
    
    RunTimeScheme(zero, Tmax, dt/4, y_rk1, fct, rk);
    
    VectReal_wp y1, y1_prime;
    fct.Extract(yinit, y1, y1_prime);
    fct.Init(dt/2);

    RunSecondOrderScheme(zero, Tmax, dt/2, y1, y1_prime, fct, scheme);

    VectReal_wp u_rk0, u_rk1, v_rk0, v_rk1;
    fct.Extract(y_rk0, u_rk0, v_rk0);
    fct.Extract(y_rk1, u_rk1, v_rk1);
    
    Real_wp err = GetRelativeError(y0, u_rk0);
    Real_wp err_half = GetRelativeError(y1, u_rk1);
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

class LinearFunctionSplit  : public VirtualOdeSystem<Real_wp>
{
public :
  Matrix<Real_wp> Dh, Bh, RhScal, RhVec, ShScal, ShVec;
  Matrix<Real_wp> Dh_lu, DhMinusdtSh, DhPlusdtSh;
  Matrix<Real_wp> Bh_lu, BhMinusdtSh, BhPlusdtSh;
  IVect pivot_Dh, pivot_DhPlusdtSh, pivot_Bh, pivot_BhPlusdtSh;
  VectReal_wp Fscal, Fvec;
  
  LinearFunctionSplit()
  {
    Real_wp one(1);
    int Nvec = 3, Nscal = 2;
    Fscal.Reallocate(Nscal); Fvec.Reallocate(Nvec);
    Dh.Reallocate(Nscal, Nscal); ShScal.Reallocate(Nscal, Nscal);
    Bh.Reallocate(Nvec, Nvec); ShVec.Reallocate(Nvec, Nvec);
    
    Real_wp coef = one/RAND_MAX;
    Fscal.FillRand(); Mlt(coef, Fscal);
    Fvec.FillRand(); Mlt(coef, Fvec);
    
    GenerateSdpMatrix(coef, Dh);
    for (int i = 0; i < Dh.GetM(); i++)
      Dh(i, i) += one;
    
    GenerateSdpMatrix(coef, Bh);
    for (int i = 0; i < Bh.GetM(); i++)
      Bh(i, i) += one;

    GenerateSdpMatrix(coef, ShScal);
    GenerateSdpMatrix(coef, ShVec);
        
    Dh_lu = Dh;
    GetLU(Dh_lu, pivot_Dh);
    Bh_lu = Bh;
    GetLU(Bh_lu, pivot_Bh);

    RhScal.Reallocate(Nscal, Nvec); RhScal.FillRand(); Mlt(coef, RhScal);
    RhVec.Reallocate(Nvec, Nscal);
    Transpose(RhScal, RhVec); Mlt(-one, RhVec);
  }
  
  int GetM() { return Dh.GetM()+Bh.GetM(); }
  int GetNscal() { return Dh.GetM(); }
  
  void Init(const Real_wp& dt)
  {
    DhMinusdtSh = Dh;
    Add(-dt/2, ShScal, DhMinusdtSh); 
    
    DhPlusdtSh = Dh;
    Add(dt/2, ShScal, DhPlusdtSh); 
    
    GetLU(DhPlusdtSh, pivot_DhPlusdtSh);

    BhMinusdtSh = Bh;
    Add(-dt/2, ShVec, BhMinusdtSh); 
    
    BhPlusdtSh = Bh;
    Add(dt/2, ShVec, BhPlusdtSh); 
    
    GetLU(BhPlusdtSh, pivot_BhPlusdtSh);
  }
  
  void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& X)
  {
    Real_wp one(1), beta(50);
    Real_wp gaussian = exp(-beta*square(t-one));    
    Real_wp pulse = gaussian;
    if (n == 1)
      pulse = -2*beta*(t-one)*gaussian;
    else if (n == 2)
      pulse = -2*beta*gaussian + square(2*beta*(t-one))*gaussian;
    else if (n > 2)
      {
        abort();
      }

    pulse *= alpha;
    for (int i = 0; i < Fscal.GetM(); i++)
      X(i) += pulse*Fscal(i);
  }

  void AddVectorialTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& X)
  {
    Real_wp one(1), beta(50);
    Real_wp gaussian = exp(-beta*square(t-one));    
    Real_wp pulse = gaussian;
    if (n == 1)
      pulse = -2*beta*(t-one)*gaussian;
    else if (n == 2)
      pulse = -2*beta*gaussian + square(2*beta*(t-one))*gaussian;
    else if (n > 2)
      {
        abort();
      }

    pulse *= alpha;
    for (int i = 0; i < Fvec.GetM(); i++)
      X(i) += pulse*Fvec(i);
  }

  void Extract(const VectReal_wp& Y, VectReal_wp& U, VectReal_wp& V)
  {
    int Nscal = Dh.GetM(), Nvec = Bh.GetM();
    U.Reallocate(Nscal); V.Reallocate(Nvec);
    for (int i = 0; i < Nscal; i++)
      U(i) = Y(i);

    for (int i = 0; i < Nvec; i++)
      V(i) = Y(Nscal+i);
  }
  
  void Regroup(const VectReal_wp& U, const VectReal_wp& V, VectReal_wp& Y)
  {
    int Nscal = Dh.GetM(), Nvec = Bh.GetM();
    Y.Reallocate(Nscal+Nvec);
    for (int i = 0; i < Nscal; i++)
        Y(i) = U(i);

    for (int i = 0; i < Nvec; i++)
      Y(Nscal+i) = V(i);
  }

  void EvaluateFunction(const Real_wp& t, const VectReal_wp& Y, VectReal_wp& ProdY,
			bool invert_mass = true, bool source = true)
  {
    Real_wp zero(0), one(1);
    int Nscal = Dh.GetM(), Nvec = Bh.GetM();
    VectReal_wp U, V, ProdU(Nscal), ProdV(Nvec);
    ProdU.Fill(zero); ProdV.Fill(zero);
    Extract(Y, U, V);

    // Dh dU/dt - RhS V + Sh U = Fs    
    Mlt(RhScal, V, ProdU);
    MltAdd(-one, ShScal, U, one, ProdU);
    if (source)
      AddScalarTimeSource(one, t, 0, ProdU);
    
    SolveLU(Dh_lu, pivot_Dh, ProdU);
    
    // Bh dV/dt - RhV U + ShV V = Fv
    Mlt(RhVec, U, ProdV);
    MltAdd(-one, ShVec, V, one, ProdV);
    if (source)
      AddVectorialTimeSource(one, t, 0, ProdV);
    
    SolveLU(Bh_lu, pivot_Bh, ProdV);
    
    Regroup(ProdU, ProdV, ProdY);
  }


  void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, ShScal, u, beta, v);
  }

  void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, ShVec, u, beta, v);
  }

  void ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                             const Real_wp& beta, VectReal_wp& v, bool extrapol = true)
  {
    MltAdd(alpha, RhScal, u, beta, v);
  }

  void ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v, bool extrapol = true)
  {
    MltAdd(alpha, RhVec, u, beta, v);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, Dh, u, beta, v);
  }

  void Assemble(VectReal_wp& u) const {}
  
  void SolveOperatorDh(VectReal_wp& u) { SolveLU(Dh_lu, pivot_Dh, u); }
  
  void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, DhMinusdtSh, u, beta, v);
  }

  void SolveOperatorBh(VectReal_wp& u) { SolveLU(Bh_lu, pivot_Bh, u); }
  
  void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v)
  {
    MltAdd(alpha, BhMinusdtSh, u, beta, v);
  }
  
  void SolveOperatorDhPlusdtSh(VectReal_wp& u) { SolveLU(DhPlusdtSh, pivot_DhPlusdtSh, u); }

  void SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& u, VectReal_wp& v)
  {
    VectReal_wp x(u);
    SolveLU(DhPlusdtSh, pivot_DhPlusdtSh, x);
    Add(alpha, x, v);
  }

  void SolveOperatorBhPlusdtSh(VectReal_wp& u) { SolveLU(BhPlusdtSh, pivot_BhPlusdtSh, u); }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}


};

template<class TimeScheme>
void TestLinearFirstSplit(TimeScheme& scheme, int r, Real_wp dt)
{
  RungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(4);
  Real_wp one(1), zero(0), Tmax(2);
  Real_wp coef = one/RAND_MAX;
  Tmax = Real_wp(2);
  
  {
    LinearFunctionSplit fct;
    
    VectReal_wp yinit(fct.GetM());
    yinit.FillRand(); Mlt(coef, yinit);
    //yinit.Fill(zero);
    
    VectReal_wp y_rk0(yinit);
        
    RunTimeScheme(zero, Tmax, dt/2, y_rk0, fct, rk);
    
    VectReal_wp y0, y0_prime;
    fct.Extract(yinit, y0, y0_prime);
    fct.Init(dt);
    
    RunFirstOrderScheme(zero, Tmax, dt, y0, y0_prime, fct, scheme);
    
    VectReal_wp y_rk1(yinit);
    
    RunTimeScheme(zero, Tmax, dt/4, y_rk1, fct, rk);
    
    VectReal_wp y1, y1_prime;
    fct.Extract(yinit, y1, y1_prime);
    fct.Init(dt/2);

    RunFirstOrderScheme(zero, Tmax, dt/2, y1, y1_prime, fct, scheme);
    
    VectReal_wp u_rk0, u_rk1, v_rk0, v_rk1;
    fct.Extract(y_rk0, u_rk0, v_rk0);
    fct.Extract(y_rk1, u_rk1, v_rk1);
    
    Real_wp err = GetRelativeError(y0, u_rk0);
    Real_wp err_half = GetRelativeError(y1, u_rk1);
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

class ExponentialMatrix
{
public :
  Matrix<Complex_wp> Plu, Ptrue;
  Vector<Complex_wp> Lambda;
  IVect pivot;
  
  void Init(const Matrix<Real_wp>& A_)
  {
    Matrix<Real_wp> A(A_);
    VectReal_wp Lreal, Limag; Matrix<Real_wp> P;
    Ptrue.Reallocate(A.GetM(), A.GetM());
    GetEigenvaluesEigenvectors(A, Lreal, Limag, P);  
    int num = 0;
    while (num < A.GetM())
      {
        if (Limag(num) == 0)
          {
            for (int j = 0; j < A.GetM(); j++)
              Ptrue(j, num) = P(j, num);
            
            num++;
          }
        else
          {
            for (int j = 0; j < A.GetM(); j++)
              {
                Ptrue(j, num) = Complex_wp(P(j, num), P(j, num+1));
                Ptrue(j, num+1) = Complex_wp(P(j, num), -P(j, num+1));
              }
            
            num += 2;
          }
      }
    
    Lambda.Reallocate(A.GetM());
    for (int i = 0; i < A.GetM(); i++)
      Lambda(i) = Complex_wp(Lreal(i), Limag(i));
      
    Plu = Ptrue; GetLU(Plu, pivot);
  }
  
  void Evaluate(const Real_wp& Tmax, const VectReal_wp& y0, VectReal_wp& y_exact)
  {
    int N = Ptrue.GetM();
    VectComplex_wp ytmp(N), ysol(N);
    Real_wp zero(0);
    
    // exact solution
    for (int i = 0; i < N; i++)
      ytmp(i) = Complex_wp(y0(i), zero);
    
    SolveLU(Plu, pivot, ytmp);
    
    for (int i = 0; i < ytmp.GetM(); i++)
      ytmp(i) = exp(Lambda(i)*Tmax)*ytmp(i);
  
    ysol.Fill(Complex_wp(zero, zero));
    Mlt(Ptrue, ytmp, ysol);
    for (int i = 0; i < N; i++)
      y_exact(i) = real(ysol(i));
  }
  
};

class SumLinearOperator  : public VirtualOdeSystem<Real_wp>
{
public :
  ExponentialMatrix expA, expB;
  
  SumLinearOperator(const Matrix<Real_wp>& A, const Matrix<Real_wp>& B)
  {
    expA.Init(A); expB.Init(B);
  }
  
  void SolveSplitOperator(const Real_wp& t, const Real_wp& dt, 
                          const VectReal_wp& Un, VectReal_wp& Un_next, int num)
  {
    if (num == 0)
      {
        expA.Evaluate(dt, Un, Un_next);
      }
    else
      {
        expB.Evaluate(dt, Un, Un_next);
      }
  }  

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0) {}
  
};


template<class Scheme>
void TestSplittingScheme(Scheme& scheme, int r, const Real_wp& dt, const Real_wp& dt_bis)
{
  Matrix<Real_wp> A(4, 4), B(4, 4), C(4, 4);
  A.FillRand(); Mlt(Real_wp(1)/RAND_MAX, A);
  B.FillRand(); Mlt(Real_wp(1)/RAND_MAX, B);

  for (int i = 0; i < A.GetM(); i++)
    {
      for (int j = i+1; j < A.GetM(); j++)
        {
          A(i, j) = -A(j, i);
          B(i, j) = -B(j, i);
        }
      
      A(i, i) = 0.0;
      B(i, i) = 0.0;
    }
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      C(i, j) = A(i, j) + B(i, j);
  
  Real_wp Tmax(2), zero(0);
  ExponentialMatrix fct_ex;
  fct_ex.Init(C);

  SumLinearOperator fctA(A, B);
  VectReal_wp y0, y_exact(A.GetM());
  y0.Reallocate(A.GetM());
  y0.FillRand(); Mlt(Real_wp(1)/RAND_MAX, y0);

  fct_ex.Evaluate(Tmax, y0, y_exact);

  VectReal_wp yinit(y0);  
  // numerical solution for dt
  RunTimeScheme(zero, Tmax, dt, y0, fctA, scheme);
  
  Real_wp err = GetRelativeError(y0, y_exact);

  // numerical solution for dt/2
  y0 = yinit;
  RunTimeScheme(zero, Tmax, dt/2, y0, fctA, scheme);
  
  Real_wp err_half = GetRelativeError(y0, y_exact);

  Real_wp order = (log(err) - log(err_half))/log(2.0);
  //DISP(err); DISP(err_half);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  if (order < r*(1.0-err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
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
  
  /*  {
    cout << endl << "Testing multistep schemes" << endl;
    SymmetricMultistep_Iterator<Real_wp> scheme;
    cout << "Rentrez N tel que dt = 1/N" << endl;
    int N; cin >> N;
    Real_wp dt = Real_wp(1)/N;
    
    VectReal_wp theta;
    scheme.SetOrder(2, theta);
    TestSecondOrderScheme(scheme, 2, dt);

    theta.Reallocate(1); theta(0) = 0.2;
    scheme.SetOrder(4, theta);
    TestSecondOrderScheme(scheme, 4, dt/20);

    theta.Reallocate(2); theta(0) = 0.2; theta(1) = 0.7;
    scheme.SetOrder(6, theta);
    TestSecondOrderScheme(scheme, 6, dt/20);

    theta.Reallocate(3); theta(0) = 0.2; theta(1) = 0.7; theta(2) = 1.2;
    scheme.SetOrder(8, theta);
    TestSecondOrderScheme(scheme, 8, dt/20);
  }
  */
  
  {
    cout << "Testing explicit Runge-Kutta Nystrom schemes ..." << endl;
    cout << "Rentrez N tel que dt = 1/N" << endl;
    
    Real_wp dt = Real_wp(1)/100;
    RungeKuttaNystrom_Iterator<Real_wp> scheme;

    scheme.SetOrder(2);
    TestSecondOrderScheme(scheme, 2, dt);

    scheme.SetOrder(3);
    TestSecondOrderScheme(scheme, 3, dt);

    scheme.SetOrder(4);
    TestSecondOrderScheme(scheme, 4, dt);
    
    scheme.SetOrder(5);
    TestSecondOrderScheme(scheme, 5, dt);
    
    dt = Real_wp(1) / 20;
    
    scheme.SetOrder(6);
    TestSecondOrderScheme(scheme, 6, dt);

    dt = Real_wp(1) / 10;
    
    scheme.SetOrder(7);
    TestSecondOrderScheme(scheme, 7, dt);
    
    if (!double_prec)
      {
	dt = Real_wp(1) / 50;
	scheme.SetOrder(8);
	TestSecondOrderScheme(scheme, 8, dt);
	
	scheme.SetOrder(10);
	TestSecondOrderScheme(scheme, 10, dt);
	
	dt = Real_wp(1) / 10;
	
	scheme.SetOrder(12);
	TestSecondOrderScheme(scheme, 11, dt);
      }
  }


  {
    cout << endl << "Testing multistep schemes" << endl;
    SymmetricMultistep_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10;
    
    VectReal_wp theta;
    scheme.SetOrder(2, theta);
    TestSecondOrderScheme(scheme, 2, dt);

    theta.Reallocate(1); theta(0) = pi_wp/2;
    scheme.SetOrder(4, theta);
    TestSecondOrderScheme(scheme, 4, dt/20);

    if (epsilon_machine < 1e-24)
      {
	theta.Reallocate(2);
	theta(0) = pi_wp/2; theta(1) = 2*pi_wp/3;
	scheme.SetOrder(6, theta);
	TestSecondOrderScheme(scheme, 6, dt/20);
	
	theta.Reallocate(3);
	theta(0) = pi_wp/6; theta(1) = pi_wp/2; theta(2) = 2*pi_wp/3;
	scheme.SetOrder(8, theta);
	TestSecondOrderScheme(scheme, 8, dt/20);
      }
  }
  
  {
    cout << "Testing explicit Runge-Kutta Nystrom schemes ..." << endl;
    Real_wp dt = Real_wp(1)/100;
    RungeKuttaNystrom_Iterator<Real_wp> scheme;
    
    if (double_prec)
      dt = 0.01;

    scheme.SetOrder(2);
    TestSecondOrderScheme(scheme, 2, dt);

    scheme.SetOrder(3);
    TestSecondOrderScheme(scheme, 3, dt);

    scheme.SetOrder(4);
    TestSecondOrderScheme(scheme, 4, dt);
    
    scheme.SetOrder(5);
    TestSecondOrderScheme(scheme, 5, dt);
    
    dt = Real_wp(1)/40;
    if (double_prec)
      dt = 0.05;

    scheme.SetOrder(6);
    TestSecondOrderScheme(scheme, 6, dt);

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.05;

    scheme.SetOrder(2, true);
    TestSecondOrderScheme(scheme, 2, dt);
    
    scheme.SetOrder(4, true);
    TestSecondOrderScheme(scheme, 4, dt);

    dt = Real_wp(1)/40;
    if (double_prec)
      dt = 0.05;
    
    scheme.SetOrder(6, true);
    TestSecondOrderScheme(scheme, 6, dt);

    if (double_prec)
      dt = 0.5;

    if (!double_prec)
      {
	dt = Real_wp(1)/20;
        scheme.SetOrder(8, true);
        TestSecondOrderScheme(scheme, 8, dt);

        scheme.SetOrder(12);
        TestSecondOrderScheme(scheme, 12, dt);

	dt = Real_wp(1)/50;
	scheme.SetOrder(10);
	TestSecondOrderScheme(scheme, 10, dt);

	dt = Real_wp(1)/8;	
	scheme.SetOrder(12, false, scheme.SHARP_QURESHI_GRAZIER);
	TestSecondOrderScheme(scheme, 12, dt);
      }
  }
    
  if (!double_prec)
  {
    cout << "Testing optimal modified equation schemes ..." << endl;
    
    OptimalModifiedEquation_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;
    
    // testing stable algorithm
    scheme.SetOrder(2, 0, false, 10);
    TestLinearScheme(scheme, 2, dt, threshold, true);

    scheme.SetOrder(4, 6, false, 10);
    TestLinearScheme(scheme, 4, dt, threshold, true);

    dt = Real_wp(1)/10;
    scheme.SetOrder(6, 8, false, 50);
    TestLinearScheme(scheme, 6, dt, threshold, true, false);

    scheme.SetOrder(8, 8, false, 1000);
    TestLinearScheme(scheme, 8, dt, threshold, true, false);

    scheme.SetOrder(10, 8, false, 20000);
    TestLinearScheme(scheme, 10, dt, threshold, true, false);

    // testing unstable algorithm
    dt = Real_wp(1)/100;
    scheme.SetOrder(6, 0, true, 10);
    TestLinearScheme(scheme, 6, dt, threshold, true, false);

    scheme.SetOrder(8, 0, true, 500);
    TestLinearScheme(scheme, 8, dt, threshold, true, false);

    scheme.SetOrder(4, 8, true, 100);
    TestLinearScheme(scheme, 4, dt, threshold, true, false);

    scheme.SetOrder(6, 12, true, 1000);
    TestLinearScheme(scheme, 6, dt, threshold, true, false);

    dt = Real_wp(1)/4;
    
    scheme.SetOrder(8, 8, true, 10000);
    TestLinearScheme(scheme, 8, dt, threshold, true, false);
    
    scheme.SetOrder(12, 12, true, 100000);
    TestLinearScheme(scheme, 10, dt, threshold, true, false);
  }

  {
    cout << "Testing splitting schemes" << endl;
    SplitScheme_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/10;
    
    if (!double_prec)
      dt = Real_wp(1)/50;
    
    scheme.SetOrder(1);
    TestSplittingScheme(scheme, 1, dt, dt);

    scheme.SetOrder(2);
    TestSplittingScheme(scheme, 2, dt, dt);

    scheme.SetOrder(4);
    TestSplittingScheme(scheme, 4, dt, dt);

    if (!double_prec)
      dt = Real_wp(1)/20;

    scheme.SetOrder(6);
    TestSplittingScheme(scheme, 6, dt, dt);

    scheme.SetOrder(8);
    TestSplittingScheme(scheme, 8, dt, dt);
  }

  {
    cout << endl << "Testing modified equation schemes" << endl;
    ModifiedEquationIterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;
    
    scheme.SetOrder(2);
    TestLinearSecond(scheme, 2, dt);

    scheme.SetOrder(4);
    TestLinearSecond(scheme, 4, dt);
  }
  
  {
    cout << endl << "Testing modified equation schemes on first-order formulation" << endl;
    ModifiedEquationSystemIterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;

    scheme.SetOrder(2);
    TestLinearFirstSplit(scheme, 2, dt);
    
    scheme.SetOrder(4);
    TestLinearFirstSplit(scheme, 4, dt);
  }
  
  {
    cout << "Testing explicit Runge-Kutta schemes ..." << endl;
    
    // explicit Runge-Kutta schemes
    RungeKutta_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/1000;
    if (double_prec)
      dt = 0.002;

    scheme.SetOrder(1);    
    TestFirstOrderScheme(scheme, 1, dt, dt);

    dt = Real_wp(1)/200;
    if (double_prec)
      dt = 0.005;

    scheme.SetOrder(2);    
    TestFirstOrderScheme(scheme, 2, dt, dt);

    if (double_prec)
      dt = 0.01;

    scheme.SetOrder(3);    
    TestFirstOrderScheme(scheme, 3, dt, dt);

    scheme.SetOrder(4);    
    TestFirstOrderScheme(scheme, 4, dt, dt);

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.02;

    scheme.SetOrder(5);    
    TestFirstOrderScheme(scheme, 5, dt, dt);

    if (double_prec)
      dt = 0.04;

    scheme.SetOrder(6);    
    TestFirstOrderScheme(scheme, 6, dt, dt);

    if (double_prec)
      dt = 0.06;

    //scheme.SetOrder(7);    
    //TestFirstOrderScheme(scheme, 7, dt, max(Real_wp(1)/10000, dt/4));
    
    dt = Real_wp(1)/20;
    if (double_prec)
      dt = 0.08;

    scheme.SetOrder(8);    
    TestFirstOrderScheme(scheme, 8, dt, max(Real_wp(1)/10000, dt/2));

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.02;

    scheme.SetOrder(5, scheme.FEHLBERG);    
    TestFirstOrderScheme(scheme, 5, dt, dt/2);

    if (double_prec)
      dt = 0.03;

    //scheme.SetOrder(5, scheme.DOPRI);    
    //TestFirstOrderScheme(scheme, 5, dt, dt/2);

    dt = Real_wp(1)/20;
    if (double_prec)
      dt = 0.08;

    scheme.SetOrder(8, scheme.FEHLBERG);    
    TestFirstOrderScheme(scheme, 8, dt, dt);

    scheme.SetOrder(8, scheme.DOPRI);    
    TestFirstOrderScheme(scheme, 8, dt, dt);    
  }  

  if (double_prec)
    {
      cout << "Testing Tal-ezer scheme ..." << endl;
      
      Talezer_Iterator<Real_wp> scheme;
      Real_wp dt(1);
      
      scheme.SetSpectralRadius(Real_wp(1));
      scheme.SetOrderIntegration(20);
      scheme.SetOrder(15, dt);
      
      TestLinearScheme(scheme, 1, dt, threshold, false);
    }

  {
    cout << "Testing Adams-Bashforth schemes ..." << endl;
    
    AdamsBashforth_Moulton_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/500;
    if (double_prec)
      dt = 0.001;
    
    scheme.SetOrder(1);
    TestFirstOrderScheme(scheme, 1, dt, max(Real_wp(1)/10000, dt/8));

    dt = Real_wp(1)/200;
    if (double_prec)
      dt = 0.002;
    
    scheme.SetOrder(2);
    TestFirstOrderScheme(scheme, 2, dt, dt);

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.01;

    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt, dt);
    
    dt = Real_wp(1)/50;
    
    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt, dt/4);

    scheme.SetOrder(5);
    TestFirstOrderScheme(scheme, 5, dt, dt/4);

    if (double_prec)
      dt = 0.02;

    //scheme.SetOrder(6);
    //TestFirstOrderScheme(scheme, 6, dt, dt);    

    if (!double_prec)
      {
        //scheme.SetOrder(7);
        //TestFirstOrderScheme(scheme, 7, dt, dt);    
      }
  }
  
  {
    cout << "Testing Adams-Bashforth-Moulton schemes ..." << endl;
    
    AdamsBashforth_Moulton_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/1000;

    if (double_prec)
      dt = 0.001;
    
    scheme.SetOrder(1, scheme.BASHFORTH_MOULTON);
    TestFirstOrderScheme(scheme, 1, dt, dt);

    dt = Real_wp(1)/200;
    if (double_prec)
      dt = 0.002;

    scheme.SetOrder(2, scheme.BASHFORTH_MOULTON);
    TestFirstOrderScheme(scheme, 2, dt, max(Real_wp(1)/10000, dt/2));

    if (double_prec)
      dt = 0.005;

    scheme.SetOrder(3, scheme.BASHFORTH_MOULTON);
    TestFirstOrderScheme(scheme, 3, dt, dt);

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.005;

    scheme.SetOrder(4, scheme.BASHFORTH_MOULTON);
    TestFirstOrderScheme(scheme, 4, dt, dt);

    if (double_prec)
      dt = 1.0/60;

    //scheme.SetOrder(5, scheme.BASHFORTH_MOULTON);
    //TestFirstOrderScheme(scheme, 5, dt, dt/2);

    if (!double_prec)
      {
        //scheme.SetOrder(6, scheme.BASHFORTH_MOULTON);
        //TestFirstOrderScheme(scheme, 6, dt, dt);    

        //scheme.SetOrder(7, scheme.BASHFORTH_MOULTON);
        //TestFirstOrderScheme(scheme, 7, dt, dt);    
      }
  }
  

  {
    cout << "Testing Taylor schemes ..." << endl;
    
    TaylorSeries_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;
    
    scheme.SetOrder(1);
    TestLinearScheme(scheme, 1, dt, threshold, true);

    scheme.SetOrder(2);
    TestLinearScheme(scheme, 2, dt, threshold, true);

    scheme.SetOrder(3);
    TestLinearScheme(scheme, 3, dt, threshold, true);

    scheme.SetOrder(4);
    TestLinearScheme(scheme, 4, dt, threshold, true);
  }
  
  {
    cout << "Testing low-storage explicit Runge-Kutta schemes ..." << endl;
    
    // low-storage Runge-Kutta schemes
    LowStorageRK_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;

    scheme.SetOrder(1);
    TestFirstOrderScheme(scheme, 1, dt, dt);

    scheme.SetOrder(2);
    TestFirstOrderScheme(scheme, 2, dt, dt);

    if (double_prec)
      dt = 0.001;
    
    scheme.SetOrder(3);
    TestFirstOrderScheme(scheme, 3, dt, dt);

    if (double_prec)
      dt = 0.01;
    
    scheme.SetOrder(3, 4);
    TestFirstOrderScheme(scheme, 3, dt, dt);

    dt = Real_wp(1)/50;
    if (double_prec)
      dt = 0.02;
    
    scheme.SetOrder(4);
    TestFirstOrderScheme(scheme, 4, dt, dt);

    scheme.SetOrder(4, 4);
    TestFirstOrderScheme(scheme, 4, dt, dt);

    scheme.SetOrder(4, 5);
    TestFirstOrderScheme(scheme, 4, dt, dt/2);
  }
  
  {
    cout << "Testing Butcher schemes ..." << endl;
    
    MultiStepButcher_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/200;
    
    scheme.SetOrder(1);
    TestLinearScheme(scheme, 1, dt, threshold, true);

    scheme.SetOrder(2);
    TestLinearScheme(scheme, 2, dt, threshold, true);

    dt = Real_wp(1)/100;
    if (double_prec)
      dt = 0.02;

    scheme.SetOrder(3);
    TestLinearScheme(scheme, 3, dt, threshold, true);

    scheme.SetOrder(4);
    TestLinearScheme(scheme, 4, dt, threshold, true);

    scheme.SetOrder(1, scheme.PECE);
    TestLinearScheme(scheme, 1, dt, threshold, true);

    scheme.SetOrder(2, scheme.PECE);
    TestLinearScheme(scheme, 2, dt, threshold, true);

    scheme.SetOrder(3, scheme.PECE);
    TestLinearScheme(scheme, 3, dt, threshold, true);

    scheme.SetOrder(4, scheme.PECE);
    TestLinearScheme(scheme, 4, dt, threshold, true);

  }
  
  {
    cout << "Testing Nystrom schemes ..." << endl;
    
    Nystrom_Iterator<Real_wp> scheme;
    Real_wp dt = Real_wp(1)/100;
    
    scheme.SetOrder(1, scheme.NYSTROM);
    TestFirstOrderScheme(scheme, 1, dt, dt);

    scheme.SetOrder(2, scheme.NYSTROM);
    TestFirstOrderScheme(scheme, 2, dt, dt);

    scheme.SetOrder(3, scheme.NYSTROM);
    TestFirstOrderScheme(scheme, 3, dt, dt);

    if (double_prec)
      dt = 0.01;

    //scheme.SetOrder(4, scheme.NYSTROM);
    //TestFirstOrderScheme(scheme, 4, dt, max(Real_wp(1)/10000, dt/1000));

    //scheme.SetOrder(5, scheme.NYSTROM);
    //TestFirstOrderScheme(scheme, 5, dt, dt);

    //scheme.SetOrder(6, scheme.NYSTROM);
    //TestFirstOrderScheme(scheme, 6, dt, dt);    

    //scheme.SetOrder(7, scheme.NYSTROM);
    //TestFirstOrderScheme(scheme, 7, dt, dt);    
  }
    
  cout << "All tests passed successfully" << endl;

  return 0;
}
