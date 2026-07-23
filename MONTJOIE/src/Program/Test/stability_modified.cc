#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

// computes roots of a polynomial P, real parts of roots are contained in R
// imaginary roots are contained in Rimag
Real_wp GetRoots(const UnivariatePolynomial<Real_wp>& P, VectReal_wp& R, VectReal_wp& Rimag)
{
  int first_coef = P.GetOrder();
  for (int i = 0; i < P.GetM(); i++)
    if (abs(P(i)) > 1e3*epsilon_machine)
      {
	first_coef = i;
	break;
      }
  
  int last_coef = 0;
  Real_wp one(1), zero(0), coef;
  for (int i = P.GetM()-1; i >= 0; i--)
    if (P(i) != zero)
      {
	last_coef = i;
	break;
      }
  
  // companion matrix
  int N = last_coef - first_coef;
  Matrix<Real_wp> A, V;
  A.Reallocate(N, N);
  A.Fill(zero);
  for (int i = 0; i < N-1; i++)
    A(i+1, i) = one;
  
  coef = one/P(last_coef);
  for (int i = first_coef; i < last_coef; i++)
    A(i-first_coef, N-1) = -coef*P(i);
  
  // roots are equal to eigenvalues of the companion matrix
  GetEigenvaluesEigenvectors(A, R, Rimag, V);
  
  return P(first_coef);
}

// returns CFL associated with polynomial P
// this polynomial is used in the first-order formulation
// (U^{n+1} - U{n-1})/2 = i P(i x)
// the CFL number is the maximal x such that |P(x)| <= 1
Real_wp GetCFL(UnivariatePolynomial<Real_wp>& P)
{
  VectReal_wp roots, root_imag;
  Real_wp epsilon = epsilon_machine*1e10, one(1);

  // computing minimal x such that P = 1 + epsilon
  P(0) -= one + epsilon;
  GetRoots(P, roots, root_imag);

  Real_wp cfl = 1e10;
  for (int i = 0; i < roots.GetM(); i++)
    if ((root_imag(i) == Real_wp(0)) && (roots(i) > epsilon_machine) && (roots(i) < cfl))
      cfl = roots(i);

  // taking minimum of previous x and minimal x such that P = -1 - epsilon
  P(0) += 2.0*(one+epsilon);
  GetRoots(P, roots, root_imag);
  for (int i = 0; i < roots.GetM(); i++)
    if ((root_imag(i) == Real_wp(0)) && (roots(i) > epsilon_machine) && (roots(i) < cfl))
      cfl = roots(i);  
  
  P(0) = 0;
  return cfl;
}

Real_wp GetTauCFL(const UnivariatePolynomial<Real_wp>& P, const VectReal_wp& tau)
{
  Real_wp x0 = tau(tau.GetM()-1);
  Real_wp dx = Real_wp(1)/10;
  if (tau.GetM() >= 2)
    dx = (x0 - tau(tau.GetM()-2)) / 100;
  
  Real_wp target = -P.Evaluate(x0);
  Real_wp x1 = x0 +dx, f;
  while (true)
    {
      f = P.Evaluate(x1);
      if (target < 0)
        {
          if ( f < target)
            break;          
        }
      else
        {
          if (f > target)
            break;
        }
      
      x0 = x1;
      x1 += dx;
    }
  
  //DISP(x0); DISP(x1); DISP(target); DISP(f);
  // dichotomy
  Real_wp c;
  while (abs(x1-x0) > max(1e-32, 10*epsilon_machine))
    {
      c = (x0+x1)/2;
      f = P.Evaluate(c);
      if (target < 0)
        {
          if (f < target)
            x1 = c;
          else
            x0 = c;
        }
      else
        {
          if (f > target)
            x1 = c;
          else
            x0 = c;
        }
    }
  
  //DISP(c); DISP(f);
  return c;
}

//! class to compute approximation of sinus with optimal stability
/*!
  Optimal stability is obtained by computing tangent points tau_i
  such that P(tau_i) = -1 or 1
  and P(tau_i)' = 0
  this last condition provides a non-linear system which can be solved
  with Newton's method
  
  P is searched in the form
  P = Esin + x^{2r+1} R
  where R is a polynomial with even orders
  and Esin the Taylor expansion of sinus :
  Esin = x - x^3/6 + x^5/120 ... + (-1)^{r+1} x^{2r-1} / (2r-1)!
*/
template<class T>
class ApproxSinus
{
public :
  // Taylor expansion of sinus, derivative and second derivative
  UnivariatePolynomial<T> ESin, dESin, d2ESin;
  // order : order of Taylor expansion, n : number of tangent points searched
  int n, order;
  // right hand side (containing alternatively -1 and +1 depending on the tangent point)
  Vector<T> v_rhs;
  
  void Init(int m, int k)
  {    
    order = m;
    n = k;
    ESin.SetOrder(2*m-1);
    UnivariatePolynomial<T> x, monome;
    T one(1);
    x.SetOrder(1);
    x(1) = one;
    monome = x;
    ESin = x;
    T coef = one;
    for (int i = 1; i < m; i++)
      {
	coef *= -one/((2*i)*(2*i+1));
	monome *= x*x;
	ESin += coef*monome;
      }
    
    DerivatePolynomial(ESin, dESin);
    DerivatePolynomial(dESin, d2ESin);
    
    // filling v_rhs with -1 and +1
    // epsilon is introduce to ensure that the polynomial
    // is strictly inside interval [-1, 1]
    v_rhs.Reallocate(n);
    int pos = order%2;
    Real_wp epsilon = 0;
    for (int i = 0; i < n; i++)
      {
	if (pos == 0)
	  v_rhs(i) = -one + epsilon;
	else
	  v_rhs(i) = one - epsilon;
	
	pos = 1-pos;
      }
  }
  
  int GetNbUnknowns()
  {
    return n;
  }
  
  // initial solution for tangent points
  // tangents points should be close to roots of Chebychev polynomials
  void FindInitGuess(Vector<T>& tau)
  {
    tau.Reallocate(n);
    /* int order_tcheb = (2*(1+n)-1); DISP(order_tcheb);
    Real_wp dtheta = pi_wp/order_tcheb;
    int j = 0; tau.Fill(0);
    if (order%2 == 0)
      for (int i = order_tcheb/2-1; i >= order_tcheb/2-n; i--)
        tau(j++) = order_tcheb*cos(i*dtheta);
    else
      for (int i = order_tcheb/2; i >= order_tcheb/2-n+1; i--)
        tau(j++) = order_tcheb*cos(i*dtheta);
    */
    if (order%2 == 0)
      for (int i = 0; i < n; i++)
        tau(i) = pi_wp/2 + (i+1)*pi_wp;
    else
      for (int i = 0; i < n; i++)
        tau(i) = pi_wp/2 + i*pi_wp;
  }
  
  // computes VanderMonde matrix
  void GetVandermondeMatrix(const Vector<T>& tau, Matrix<T>& Vdm)
  {
    Vdm.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Vdm(i, j) = pow(tau(i), 2*(order+j)+1);
  }
  
  // computes differential of Vandermonde matrix with respect to tau
  void GetDerivativeVdm(const Vector<T>& tau, Matrix<T>& dV)
  {
    dV.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	dV(i, j) = (2*(order+j)+1)*pow(tau(i), 2*(order+j));
  }
  
  // computes polynomial expansion of R from tangent points tau
  void GetDecompositionR(const Vector<T>& tau, Vector<T>& R)
  {
    // evaluating Vandermonde matrix
    Matrix<T> Vdm;
    GetVandermondeMatrix(tau, Vdm);
    // factorizing it
    IVect pivot;
    GetLU(Vdm, pivot);
    R.Reallocate(n);
    // then solving v_i = Esin(tau_i) + tau_i^{2r+1} \sum R_k  tau_i^{2 k}
    // unknowns are coefficients R_k
    for (int i = 0; i < n; i++)
      R(i) = v_rhs(i) - ESin.Evaluate(tau(i));
    
    SolveLU(Vdm, pivot, R);
  }
  
  // evaluates non-linear system to solve
  // The non-linear system is obtained with imposing P'(tau_i) = 0
  // => Esin'(tau_i) + tau_i^{2r} \sum R_k  2 (order+1+k) tau_i^{2 k} = 0
  void EvaluateFunction(const Vector<T>& tau, Vector<T>& f)
  {
    Vector<T> R;
    GetDecompositionR(tau, R);
    f.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
	f(i) = dESin.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  f(i) += R(j)*(2*(order+j)+1)*pow(tau(i), 2*(order+j));
      }
  }
  
  // numerical jacobian matrix
  void EvaluateJacobianNum(const Vector<T>& tau, Matrix<T>& df_num)
  {
    df_num.Reallocate(n, n);
    Vector<T> x(tau), fp, fm; Real_wp h = 1e-6;
    for (int j = 0; j < n; j++)
      {
	x(j) += h;
	EvaluateFunction(x, fp);
	x(j) -= 2*h;
	EvaluateFunction(x, fm);
	x(j) += h;
	
	for (int i = 0; i < n; i++)
	  df_num(i, j) = (fp(i) - fm(i))/(2*h);
      }
  }
    
  // evaluates jacobian matrix
  void EvaluateJacobian(const Vector<T>& tau, Matrix<T>& df)
  {
    Vector<T> f;
    EvaluateJacobian(tau, f, df);
  }
  
  // evaluates analytical jacobian
  void EvaluateJacobian(const Vector<T>& tau, Vector<T>& f, Matrix<T>& df)
  {
    Vector<T> R, dR;
    Matrix<T> Vdm;
    GetVandermondeMatrix(tau, Vdm);
    GetInverse(Vdm);
    f.Reallocate(n);
    R.Reallocate(n);
    for (int i = 0; i < n; i++)
      f(i) = v_rhs(i) - ESin.Evaluate(tau(i));
    
    Mlt(Vdm, f, R);
    
    for (int i = 0; i < n; i++)
      {
	f(i) = dESin.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  f(i) += R(j)*(2*(order+j)+1)*pow(tau(i), 2*(order+j));
      }
    
    df.Reallocate(n, n);
    df.Fill(0);
    Matrix<T> dV;
    GetDerivativeVdm(tau, dV);
    
    T psi_seconde;
    for (int i = 0; i < n; i++)
      {
	psi_seconde = d2ESin.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  psi_seconde += R(j)*(2*(order+j)+1)*(2*(order+j))*pow(tau(i), 2*(order+j)-1);	
	
	df(i, i) = psi_seconde;
	for (int j = 0; j < n; j++)
	  {
	    for (int l = 0; l < n; l++)
	      df(i, j) -= dV(i, l)*Vdm(l, j)*f(j);
	  }
      }    
  }
  
};

template<class T>
class RecurrenceSinus
{
public :
  UnivariatePolynomial<Real_wp> Ptarget;
  
  void EvaluateFunction(const VectReal_wp& x, VectReal_wp& f)
  {
    // expressing relation recurrence
    UnivariatePolynomial<Real_wp> Pn_prev, Pn, Pn_next, monome;
    monome.SetOrder(1);
    monome(0) = Real_wp(0); monome(1) = Real_wp(1);
    Pn_prev.SetOrder(0);
    Pn_prev(0) = -Real_wp(1);
    Pn = x(0)*monome;
    for (int k = 2; k < 2*x.GetM(); k++)
      {
        Pn_next = -x(k/2)*monome*Pn + Pn_prev;
        Pn_prev = Pn;
        Pn = Pn_next;
      }
    
    for (int i = 1; i < 2*x.GetM(); i += 2)
      {
        if (((i-1)/2)%2 == 0)
          f(i/2) = Pn(i) - Ptarget(i);
        else
          f(i/2) = Ptarget(i) + Pn(i);
      }
  }
  
  void EvaluateJacobian(const VectReal_wp& x, Matrix<Real_wp>& df)
  {
    int n = x.GetM();
    VectReal_wp y(x), fm(n), fp(n);
    Real_wp h = pow(epsilon_machine, Real_wp(1)/3);
    for (int j = 0; j < x.GetM(); j++)
      {
        y(j) += h;
        EvaluateFunction(y, fp);
        
        y(j) -= Real_wp(2)*h;
        EvaluateFunction(y, fm);
        
        for (int i = 0; i < x.GetM(); i++)
          df(i, j) = (fp(i) - fm(i)) / (Real_wp(2)*h);
        
        y(j) += h;
      }
    //DISP(df);
  }
  
  void FindRecurrence(const UnivariatePolynomial<Real_wp>& P,
                      ofstream& file_out)
  {
    Ptarget = P;
    int n = (P.GetOrder()+1)/2;
    VectReal_wp x(n);
    x(0) = Real_wp(1)/(2*n-1);
    for (int i = 1; i < x.GetM(); i++)
      x(i) = Real_wp(2)*x(0) - Real_wp(rand())*1e-15;
    
    // algo de Newton
    /* bool test_loop = true;
    Real_wp epsilon = 1e6*epsilon_machine;
    Real_wp residu = 1e30, residu_prec;
    VectReal_wp f(n); Matrix<Real_wp> df(n, n); IVect pivot(n);
    this->EvaluateFunction(x, f);    
    residu = Norm2(f);
    cout << "residu = " << residu << endl;
    while (test_loop)
      {
        this->EvaluateJacobian(x, df);
        GetLU(df, pivot);
        SolveLU(df, pivot, f);
        for (int i = 0; i < x.GetM(); i++)
          x(i) -= f(i);
        
        this->EvaluateFunction(x, f);
        residu_prec = residu;
        residu = Norm2(f);
        if ((residu_prec < residu) || (residu < epsilon))
          test_loop = false;	  

        cout << "residu = " << residu << endl;
      }
    */
    
    // utilisation de MinPack
    for (int p = 0; p < 1000; p++)
      {
        VectReal_wp f(n), scale_eqn; Matrix<Real_wp> fjac(n, n);
        
        //x(0) = Real_wp(2)/(2*n-1) * Real_wp(rand())/RAND_MAX;
        //for (int i = 1; i < x.GetM(); i++)
        //x(i) = Real_wp(3)*x(0)* Real_wp(rand())/RAND_MAX;
        
        for (int i = 0; i < x.GetM(); i++)
          x(i) = Real_wp(2)*(-2.0*Real_wp(rand())/RAND_MAX + Real_wp(1))/n;
        
        x(0) = abs(x(0));
        /*x(0) = 1.4203497533155446970375474554575e-2;
        x(1) = 1.9868024500204107268732669640581e-1;
        x(2) = -9.4820308141555142132786719160676e-2;
        x(3) = 7.8738297696164570500229093814705e-2;
        x(4) = -2.9429868117380800229324848282508e-2;
        x(5) = 4.9286591051519149195452304933472e-1;
        x(6) = -7.1998319813054873440506341355291e-2;
        x(7) = 4.7525811468298390208922031787061e-1;
        x(8) = -1.0387063829304748583885021822945e-1;
        x(9) = 4.037306893550181743979349504751e-2;
        x(10) = 0.0;
        this->EvaluateFunction(x, f);
        DISP(Norm2(f));
        */
        IVect Control(10); VectReal_wp RControl(10);
        Control(0) = 1000; Control(3) = 1; Control(4) = 0;
        RControl(0) = 1e10*epsilon_machine; RControl(1) = 1.0;
        int err = SolveMinpack(*this, x, f, fjac, scale_eqn, Control, RControl);
        
        this->EvaluateFunction(x, f);
        if (Norm2(f) < 1e-15)
          {
            DISP(x); DISP(Norm2(f));
            for (int i = 0; i < x.GetM(); i++)
              file_out << "alpha(" << i << ") = " << x(i) << ";\n";
            
            break;
          }
      }
  }
  
};

//! class to compute approximation of cosinus with optimal stability
/*!
  Optimal stability is obtained by computing tangent points tau_i
  such that P(tau_i) = -1 or 1
  and P(tau_i)' = 0
  this last condition provides a non-linear system which can be solved
  with Newton's method
  
  P is searched in the form
  P = Ecos + x^{2r} R
  where R is a polynomial with even orders
  and Ecos the Taylor expansion of cosinus :
  Ecos = 1 - x^2/2 + x^4/24 ... + (-1)^{r+1} x^{2r} / (2r)!
*/
template<class T>
class ApproxCosinus
{
public :
  // Taylor expansion of cosinus, derivative and second derivative
  UnivariatePolynomial<T> ECos, dECos, d2ECos;
  // order : order of Taylor expansion, n : number of tangent points searched
  int n, order;
  // right hand side (containing alternatively -1 and +1 depending on the tangent point)
  Vector<T> v_rhs;
  
  void Init(int m, int k)
  {    
    order = m;
    n = k;
    ECos.SetOrder(2*m-1);
    UnivariatePolynomial<T> x, monome;
    T one(1), zero(0);
    x.SetOrder(1);
    x(1) = one;
    monome.SetOrder(0); monome(0) = one;
    ECos.Fill(zero);
    ECos(0) = one;
    T coef = one;
    for (int i = 1; i <= m; i++)
      {
	coef *= -one/((2*i-1)*(2*i));
	monome *= x*x;
	ECos += coef*monome;
      }
    
    DerivatePolynomial(ECos, dECos);
    DerivatePolynomial(dECos, d2ECos);
    
    // filling v_rhs with -1 and +1
    // epsilon is introduce to ensure that the polynomial
    // is strictly inside interval [-1, 1]
    v_rhs.Reallocate(n);
    int pos = order%2;
    Real_wp epsilon = 0;
    for (int i = 0; i < n; i++)
      {
	if (pos == 0)
	  v_rhs(i) = one - epsilon;
	else
	  v_rhs(i) = -one + epsilon;
	
	pos = 1-pos;
      }
  }
  
  int GetNbUnknowns()
  {
    return n;
  }
  
  // initial solution for tangent points
  // tangents points should be close to k pi
  void FindInitGuess(Vector<T>& tau)
  {
    tau.Reallocate(n);
    // for odd orders, the instability occurs at the first tangent point pi
    // for even orders, the instability occurs at the second tangent point 2 pi
    if (order%2 == 1)
      for (int i = 0; i < n; i++)
	tau(i) = pi_wp + pi_wp*i;
    else
      for (int i = 0; i < n; i++)
	tau(i) = pi_wp + pi_wp*(i+1);
  }
  
  // computes Vandermonde matrix
  void GetVandermondeMatrix(const Vector<T>& tau, Matrix<T>& Vdm)
  {
    Vdm.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	Vdm(i, j) = pow(tau(i), 2*(order+j)+2);
  }
  
  // computes differential of Vandermonde matrix with respect to tau
  void GetDerivativeVdm(const Vector<T>& tau, Matrix<T>& dV)
  {
    dV.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	dV(i, j) = (2*(order+j)+2)*pow(tau(i), 2*(order+j)+1);
  }
  
  // computes polynomial expansion of R from tangent points tau
  void GetDecompositionR(const Vector<T>& tau, Vector<T>& R)
  {
    // evaluating Vandermonde matrix
    Matrix<T> Vdm;
    GetVandermondeMatrix(tau, Vdm);
    // factorizing it
    IVect pivot;
    GetLU(Vdm, pivot);
    R.Reallocate(n);
    // then solving v_i = ECos(tau_i) + tau_i^{2r+2} \sum R_k  tau_i^{2 k}
    // unknowns are coefficients R_k
    for (int i = 0; i < n; i++)
      R(i) = v_rhs(i) - ECos.Evaluate(tau(i));
    
    SolveLU(Vdm, pivot, R);
  }
  
  // evaluates non-linear system to solve
  // The non-linear system is obtained with imposing P'(tau_i) = 0
  // => ECos'(tau_i) + tau_i^{2r+1} \sum R_k  2 (order+k)+2) tau_i^{2 k} = 0
  void EvaluateFunction(const Vector<T>& tau, Vector<T>& f)
  {
    Vector<T> R;
    GetDecompositionR(tau, R);
    f.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
	f(i) = dECos.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  f(i) += R(j)*(2*(order+j)+2)*pow(tau(i), 2*(order+j)+1);
      }
  }
  
  // numerical jacobian matrix
  void EvaluateJacobianNum(const Vector<T>& tau, Matrix<T>& df_num)
  {
    df_num.Reallocate(n, n);
    Vector<T> x(tau), fp, fm; Real_wp h = 1e-6;
    for (int j = 0; j < n; j++)
      {
	x(j) += h;
	EvaluateFunction(x, fp);
	x(j) -= 2*h;
	EvaluateFunction(x, fm);
	x(j) += h;
	
	for (int i = 0; i < n; i++)
	  df_num(i, j) = (fp(i) - fm(i))/(2*h);
      }
  }
    
  // evaluates jacobian matrix
  void EvaluateJacobian(const Vector<T>& tau, Matrix<T>& df)
  {
    Vector<T> f;
    EvaluateJacobian(tau, f, df);
  }
  
  // evaluates analytical jacobian
  void EvaluateJacobian(const Vector<T>& tau, Vector<T>& f, Matrix<T>& df)
  {
    Vector<T> R, dR;
    Matrix<T> Vdm;
    GetVandermondeMatrix(tau, Vdm);
    GetInverse(Vdm);
    f.Reallocate(n);
    R.Reallocate(n);
    for (int i = 0; i < n; i++)
      f(i) = v_rhs(i) - ECos.Evaluate(tau(i));
    
    Mlt(Vdm, f, R);
    
    for (int i = 0; i < n; i++)
      {
	f(i) = dECos.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  f(i) += R(j)*(2*(order+j)+2)*pow(tau(i), 2*(order+j)+1);
      }
    
    df.Reallocate(n, n);
    df.Fill(0);
    Matrix<T> dV;
    GetDerivativeVdm(tau, dV);
    
    T psi_seconde;
    for (int i = 0; i < n; i++)
      {
	psi_seconde = d2ECos.Evaluate(tau(i));
	for (int j = 0; j < n; j++)
	  psi_seconde += R(j)*(2*(order+j)+2)*(2*(order+j)+1)*pow(tau(i), 2*(order+j));	
	
	df(i, i) = psi_seconde;
	for (int j = 0; j < n; j++)
	  {
	    for (int l = 0; l < n; l++)
	      df(i, j) -= dV(i, l)*Vdm(l, j)*f(j);
	  }
      }    
  }
  
};

void PolyFit(const VectReal_wp& x, const VectReal_wp& y, int n,
	     UnivariatePolynomial<Real_wp>& P)
{
  // methode QR
  int m = x.GetM();
  Matrix<Real_wp> Vdm(m, n+1);
  Vdm.Fill(Real_wp(1));
  for (int i = 0; i < m; i++)
    for (int j = 1; j <= n; j++)
      Vdm(i, j) = pow(x(i), j);
  
  Vector<Real_wp> tau;
  GetQR(Vdm, tau);
  
  VectReal_wp coef(y);
  SolveQR(Vdm, tau, coef);
  
  P.SetOrder(n);
  for (int j = 0; j <= n; j++)
    P(j) = coef(j);
}

void FindLawEfficiency(int order, int m_max, const string& file_name)
{
  ifstream file_in(file_name.data());
  
  if (!file_in.is_open())
    {
      cout << "Unable to open file " << file_name << endl;
      abort();
    }
  
  string ligne; int r, i = -1; Real_wp eff;
  bool section_ok = false;
  VectReal_wp step_x, eval;
  while (!file_in.eof())
    {
      getline(file_in, ligne);
      
      if (file_in.fail())
	break;
      
      if (ligne.size() > 4)
	{
	  if (ligne.substr(0, 3) == string("m ="))
	    {
	      int m = to_num<int>(ligne.substr(4, 1));
	      if (m == order)
		{
		  section_ok = true;
		  i++;
		}
	      else
		section_ok = false;
	    }
	  
	  if (section_ok)
	    {
	      if (ligne.substr(0, 3) == string("r ="))
		{
		  r = to_num<int>(ligne.substr(4, ligne.size()-4));
		  step_x.PushBack(Real_wp(1) / (2*(r + order)-1));	      
		}
	      
	      if (ligne.substr(0, 12) == string("Efficacite ="))
		{
		  eff = to_num<Real_wp>(ligne.substr(13, ligne.size()-13));
		  eval.PushBack(eff);
		}
	    }
	}
    }
  
  cout.precision(32);
  DISP(step_x); DISP(eval);
  
  file_in.close();

  UnivariatePolynomial<Real_wp> P;
  PolyFit(step_x, eval, 25, P);
  DISP(P);
}

template<class T>
T EvaluateEvalSin(const VectReal_wp& alpha, const T& x)
{
  T eval = 1.0;
  for (int k = alpha.GetM()-1; k >= 1; k--)
    eval = -T(alpha(k))*x*x*eval + T(1);
  
  eval *= x;
  return eval;
}

template<class T>
T EvaluateRecSin(const VectReal_wp& alpha, const T& x)
{
  T Q0 = 1.0;
  T Q1 = -T(alpha(0))*x, Qnext;
  for (int k = 1; k < alpha.GetM(); k++)
    {
      Qnext = T(alpha(k))*x*Q1 + Q0;
      Q0 = Q1; Q1 = Qnext;
      
      Qnext = -T(alpha(k))*x*Q1 + Q0;
      Q0 = Q1; Q1 = Qnext;
    }
  
  return -Q1;
}

template<class T>
T EvaluateRecSin4(int r, const T& x)
{
  T xsi = 0.5*(T(2)*r+1)/sqrt(T(r*(r+1)));
  T coef_xsi = xsi/(T(2)*r+1);
  
  T Q0 = 1, Q1 = coef_xsi*x, Qnext;
  for (int i = 1; i <= r; i++)
    {
      Qnext = -2*coef_xsi*x*Q1 + Q0;
      Q0 = Q1; Q1 = Qnext;

      Qnext = 2*coef_xsi*x*Q1 + Q0;
      Q0 = Q1; Q1 = Qnext;      
    }
  
  return Q1/xsi;
}

void CheckDoubleCFL(const VectReal_wp& tau, const Real_wp& cfl, const VectReal_wp& alpha)
{
  return;
  bool recurrence = true;
  VectReal_wp beta;
  /*
  DISP(sin(Real_wp(0.4)));
  int r = alpha.GetM()-5;

  DISP(alpha); DISP(r);
  beta.Reallocate(r+5);
	    switch (r)
	      {
	      case 1 :
		beta(0) = to_num<Real_wp>("5.5462054282148278034367069010781e-1");
		beta(1) = to_num<Real_wp>("4.5363382341602965850743180707892e-1");
		beta(2) = to_num<Real_wp>("1.4178440193067762573900598300572e-1");
		beta(3) = to_num<Real_wp>("-2.0979055655164342167501261941812e-1");
		beta(4) = to_num<Real_wp>("1.5840540082698138121279345227151e-1");
		beta(5) = to_num<Real_wp>("-9.8653612443528024127889313045831e-2");
		break;
	      case 2 :
		beta(0) = to_num<Real_wp>("4.3760135822253414895771360927587e-1");
		beta(1) = to_num<Real_wp>("-8.3631254886620617208720550189317e-2");
		beta(2) = to_num<Real_wp>("2.9877478091430204142717369460918e-1");
		beta(3) = to_num<Real_wp>("4.0487172367966396967827560433562e-1");
		beta(4) = to_num<Real_wp>("-9.8145353288666396520247973806684e-2");
		beta(5) = to_num<Real_wp>("1.5675429068734770812302608945921e-1");
		beta(6) = to_num<Real_wp>("-1.1622554532856085445722047368389e-1");
		break;
	      case 3 :
		beta(0) = to_num<Real_wp>("5.9557460646229699802429897608598e-2");
		beta(1) = to_num<Real_wp>("-4.9986760394208784616160888188979e-2");
		beta(2) = to_num<Real_wp>("5.729341901557947077785286038705e-1");
		beta(3) = to_num<Real_wp>("-5.6777384181713962416624236877106e-2");
		beta(4) = to_num<Real_wp>("1.582972225538703635436514677692e-1");
		beta(5) = to_num<Real_wp>("-1.4678186817581440509558824062881e-1");
		beta(6) = to_num<Real_wp>("5.9915306017002097575561189318235e-1");
		beta(7) = to_num<Real_wp>("-1.3639592077417859475184849673575e-1");
		break;
	      case 4 :
		beta(0) = to_num<Real_wp>("6.6833858919114633214837360668518e-2");
		beta(1) = to_num<Real_wp>("-5.1969249648370063292483204907096e-2");
		beta(2) = to_num<Real_wp>("3.0601916678616665991902243165809e-1");
		beta(3) = to_num<Real_wp>("2.4713675471060738872583744147318e-1");
		beta(4) = to_num<Real_wp>("-6.3333549389545541232300213496116e-2");
		beta(5) = to_num<Real_wp>("2.7465944538321828515677345332194e-1");
		beta(6) = to_num<Real_wp>("2.6800103014397551773078051573676e-1");
		beta(7) = to_num<Real_wp>("7.1967340546913629768925666219499e-2");
		beta(8) = to_num<Real_wp>("-1.1931479745208050999139345067477e-1");
		break;
	      case 6 :
		beta(0) = to_num<Real_wp>("4.2598135471478686970806284089629e-2");
		beta(1) = to_num<Real_wp>("-3.7823767249716001418989228033099e-2");
		beta(2) = to_num<Real_wp>("4.2053173977281746605937521886161e-1");
		beta(3) = to_num<Real_wp>("-4.4043175596230162458537358850622e-2");
		beta(4) = to_num<Real_wp>("1.5024336202668584129606480298488e-1");
		beta(5) = to_num<Real_wp>("2.3123279143084998432008913172275e-1");
		beta(6) = to_num<Real_wp>("3.2396528584918277121379384850532e-1");
		beta(7) = to_num<Real_wp>("-1.4694794337104915004404981693123e-1");
		beta(8) = to_num<Real_wp>("1.0665713758152413987001839383556e-1");
		beta(9) = to_num<Real_wp>("-8.4033837253175126290524236890858e-2");
		beta(10) = to_num<Real_wp>("3.7620271337631550481952960706065e-2");
		break;
	      case 8 :
		beta(0) = to_num<Real_wp>("3.1950635567016969381315851945034e-2");
		beta(1) = to_num<Real_wp>("-2.9845572536789141857221170812021e-2");
		beta(2) = to_num<Real_wp>("4.976486688903574382840742977027e-1");
		beta(3) = to_num<Real_wp>("7.3429268558378318772773274232592e-2");
		beta(4) = to_num<Real_wp>("1.5992240140708422186927481210951e-1");
		beta(5) = to_num<Real_wp>("-4.773782106274187797644124991281e-2");
		beta(6) = to_num<Real_wp>("3.6893988180036987766849460487039e-1");
		beta(7) = to_num<Real_wp>("-4.3486807806047725694972757807044e-2");
		beta(8) = to_num<Real_wp>("1.145269578790083225186833290496e-1");
		beta(9) = to_num<Real_wp>("-2.0859728364568566301845372526954e-1");
		beta(10) = to_num<Real_wp>("6.6082874370738628504100490392772e-2");
		beta(11) = to_num<Real_wp>("6.1948439986443482492760627742939e-2");
		beta(12) = to_num<Real_wp>("-4.478164340813285094438838424412e-2");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}				
	      }
            */
  /* VectReal_wp beta(r+1);
  beta(0) = Real_wp(1)/(2*r+1);
  for (int i = 1; i <= r; i++)
  beta(i) = Real_wp(2)/(2*r+1); */
  
  int r = tau.GetM()+1; DISP(r);
  VectReal_wp step_x(10001), eval_y(10001), eval_y_opt(10001);  
  Vector<double> eval_y_double(10001);
  for (int i = 0; i < step_x.GetM(); i++)
    {
      step_x(i) = Real_wp(i)/(step_x.GetM()-1)*cfl;
      eval_y(i) = EvaluateRecSin4(r, step_x(i));
      eval_y_opt(i) = EvaluateEvalSin(alpha, step_x(i));
      eval_y_double(i) = EvaluateEvalSin(alpha, double(step_x(i)));
    }
  
  step_x.WriteText("x.dat");
  eval_y.WriteText("eval_sin4.dat");
  eval_y_opt.WriteText("eval_sin_opt.dat");
  eval_y_double.WriteText("eval_sin_double.dat");

  return;
  
  DISP(EvaluateRecSin(beta, Real_wp(0.4)));
  double y, cfl_eval = double(cfl);
  double err(0); Real_wp yref;
  for (int i = 0; i < tau.GetM(); i++)
    {
      if (recurrence)
        {
          y = EvaluateRecSin(beta, double(tau(i)));
          yref = EvaluateRecSin(beta, tau(i));
        }
      else
        {
          y = EvaluateEvalSin(alpha, double(tau(i)));
          yref = EvaluateEvalSin(alpha, tau(i));
        }
        
      DISP(tau(i)); DISP(y); DISP(1.0-abs(y));
      DISP(Real_wp(1)-abs(yref));
      if (1.0 - abs(y) < -2e-14)
        {
          if (double(tau(i)) < cfl_eval)
            {
              cfl_eval = double(tau(i));
              err = 1.0-abs(y);
            }
        }
    }
  
  if (recurrence)
    y = EvaluateRecSin(beta, double(cfl));
  else
    y = EvaluateEvalSin(alpha, double(cfl));
  
  DISP(cfl); DISP(y); DISP(1.0-abs(y));
  DISP(cfl_eval); DISP(err);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  srand(time(NULL));
  
  if (argc < 3)
    {
      cout << "Entrez l'ordre et le nombre de termes additionnels" << endl;
      return -1;
    }
  
  int K = atoi(argv[1]);
  int m_max = atoi(argv[2]);
  
  //FindLawEfficiency(K, m_max, "optimal_coef.dat");
  //return 0;
  
  // calcul de la CFL pour l'equation modifiee
  if (false)
    {
    ofstream file_out("cfl.dat");     
    file_out.precision(32);
    for (int N = 0; N < K; N++)
      {
	
	// CFL formulation ordre 1 avec une approximation de Taylor du sinus
	ApproxSinus<Real_wp> eqn;
	eqn.Init(N+1, 0);
	Real_wp cfl = GetCFL(eqn.ESin);
	
	cout << "efficacite = " << cfl / (2*N+1) << endl;
	file_out << "m = " << N << " CFL = " << cfl << endl;
        file_out << "Efficacite = " << cfl / (2*N+1) << endl;
	
	// CFL formulation ordre 2 avec une approximation de Taylor du cosinus
	//ApproxCosinus<Real_wp> eqn;
	//eqn.Init(N+1, 0);
	//Real_wp cfl = GetCFL(eqn.ECos);
	
	//DISP(N); DISP(cfl);
	//cout << "efficacite = " << cfl / (2*(N+1)) << endl;
	//file_out << "m = " << N << " CFL = " << cfl << endl;
        //file_out << "Efficacite = " << cfl / (2*(N+1)) << endl;
      }
     
    file_out.close();
    return 0;
  }
  
  // calcul du polynome optimal
  ofstream file_out("cfl.dat");
  file_out.precision(16);
  string nom = "reccurence" + to_str(K) + ".dat";
  ofstream file_alpha_out(nom.data());
  file_alpha_out.precision(64);
  for (int r = m_max; r <= m_max; r++)
    {
      // cas du sinus
      ApproxSinus<Real_wp> eqn;
      eqn.Init(K, r);

      // cas du cosinus
      //ApproxCosinus<Real_wp> eqn;
      //eqn.Init(K, r);
      
      // algo de Newton
      Vector<Real_wp> tau, f, R;
      Matrix<Real_wp> df; IVect pivot;
      eqn.FindInitGuess(tau);
      bool test_loop = true;
      Real_wp epsilon = 1e6*epsilon_machine;
      Real_wp residu, residu_prec;
      eqn.EvaluateJacobian(tau, f, df);
      residu = Norm2(f);
      while (test_loop)
	{
          GetLU(df, pivot);
	  SolveLU(df, pivot, f);
	  for (int i = 0; i < tau.GetM(); i++)
	    tau(i) -= f(i);
	  
          residu_prec = residu;
          eqn.EvaluateJacobian(tau, f, df);
	  residu = Norm2(f);
          cout << "Residu = " << residu << endl;
	  if ((residu_prec < residu) || (residu < epsilon))
            test_loop = false;
        }
      
      eqn.GetDecompositionR(tau, R);
      UnivariatePolynomial<Real_wp> Popt;
      
      // cas du sinus
      Popt = eqn.ESin;
      Popt.ResizeOrder(2*(K+r)-1);
      for (int i = 0; i < r; i++)
        Popt(2*(K+i)+1) = R(i);
      
      /* for (int i = 0; i < tau.GetM(); i++)
        {
          DISP(i); DISP(tau(i));
          DISP(Popt.Evaluate(tau(i)));
          } */

      // cas du cosinus
      //Popt = eqn.ECos;
      //Popt.ResizeOrder(2*(K+r));
      //for (int i = 0; i < r; i++)
      //Popt(2*(K+i)+2) = R(i);
      
      VectReal_wp alpha(K+r);
      alpha.Fill(0);
      
      file_alpha_out << "m = " << K << endl;
      file_alpha_out << "r = " << r << endl;
      //RecurrenceSinus<Real_wp> rec;
      //rec.FindRecurrence(Popt, file_alpha_out);      
      
      for (int i = 1; i < K+r; i++)
        {
          alpha(i) = abs(Popt(2*i+1)/Popt(2*i-1)); 
          file_alpha_out << "alpha(" << i-1 << ") = to_num<Real_wp>(\"" << abs(Popt(2*i+1)/Popt(2*i-1)) << "\");\n"; 
        }
      
      DISP(Popt);
      //Real_wp cfl = GetCFL(Popt);
      Real_wp cfl = GetTauCFL(Popt, tau);
      DISP(cfl);
      cout << "m = " << r << ", Efficacite = " << cfl / (2*(r+K)-1) << endl;
      file_out << "K = " << K << " r = " << r << " CFL = " << cfl << endl;
      file_out << R << endl;
      
      CheckDoubleCFL(tau, cfl, alpha);

      file_alpha_out << "CFL = " << cfl << endl;
      file_alpha_out << "Efficacite = " << cfl / (2*(r+K)-1) << endl;
      file_alpha_out << endl;
    }
  
  file_out.close();
  
  return FinalizeMontjoie();
}
