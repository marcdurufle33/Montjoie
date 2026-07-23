#include "Share/MontjoieCommon.hxx"
#include "Solver/NonLinearEquations.hxx"
#include "Solver/NonLinearEquations.cxx"

using namespace Montjoie;

// #define USE_LOCAL_MAXIMA

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

// returns the amplification factor for a symmetric multistep scheme
// defined through coefficients alpha, beta
Real_wp GetAmplificationFactor(const VectReal_wp& alpha, const VectReal_wp& beta,
                               const Real_wp& z)
{
  UnivariatePolynomial<Real_wp> P;
  P.SetOrder(alpha.GetM()-1);
  
  for (int i = 0; i < alpha.GetM(); i++)
    P(i) = alpha(i) - beta(i)*z ;
  
  //DISP(P);
  VectReal_wp R, Rimag;
  GetRoots(P, R, Rimag);
  
  //DISP(R); DISP(Rimag);
  Real_wp G = 1.0;
  for (int i = 0; i < R.GetM(); i++)
    G = max(G, Seldon::abs(Complex_wp(R(i), Rimag(i))));
  
  return G;
}

// returns the CFL of the symmetric multistep scheme
Real_wp GetCFL_Scheme(const VectReal_wp& alpha, const VectReal_wp& beta, bool display = false)
{
  Real_wp one(1), zero(0);
  Real_wp dz = one/100000;
  Real_wp z0 = -one/100000;
  //Real_wp target = one + 1e3*epsilon_machine;
  Real_wp target = one + 1e-12;
  Real_wp G0 = GetAmplificationFactor(alpha, beta, z0);
  Real_wp dz_target = 1e-14; int nb_iter_max(50);

#ifdef MONTJOIE_WITH_MPFR
  target = one + 1e-15;
  dz_target = 1e4*epsilon_machine;
  nb_iter_max = 200;
#endif
  
  if (G0 > target)
    {
      // instable algorithm
      return zero;
    }
  
  // searching first point with |G| > 1
  Real_wp z1 = z0, G1 = G0;
  Real_wp Gprev = G0, zprev = z0;
  while (G1 < target)
    {      
      if (z1 <= -one/1000)
        dz = one/10000;
      else if (z1 <= -one/100)
        dz = one/500;
      else if (z1 <= -one/50)
        dz = one/200;
      else if (z1 <= -one/20)
        dz = one/100;
      else if (z1 < -one/10)
        dz = one/50;
      else if (z1 < -one)
        dz = one/5;
      else if (z1 < -Real_wp(10))
        dz = one;
      
      zprev = z0; z0 = z1; Gprev = G0; G0 = G1; z1 -= dz;
      G1 = GetAmplificationFactor(alpha, beta, z1);
      
#ifdef USE_LOCAL_MAXIMA
      if ((G0 > G1) && (G0 > Gprev) && (z0 < -one))
        {
          if (display)
            cout << "local maxima at z = " << z0 << endl;
          //DISP(G0); DISP(G1); DISP(Gprev);
          // the local maximum is found by dichotomy
          Real_wp z_min = z1, z_max = zprev, z_med = z0, Gmed = G0;
          int nb_iter = 0;
          while ((abs(z_max-z_min) > dz_target) && (nb_iter < nb_iter_max))
            {
              Real_wp za = 0.5*(z_min+z_med), zb = 0.5*(z_max+z_med);
              Real_wp Ga = GetAmplificationFactor(alpha, beta, za);
              Real_wp Gb = GetAmplificationFactor(alpha, beta, zb);
              if (Gmed > max(Ga, Gb))
                {
                  z_min = za;
                  z_max = zb;
                }
              else if (Ga > max(Gmed, Gb))
                {
                  z_max = z_med;
                  z_med = za;
                  Gmed = Ga;
                }
              else
                {
                  z_min = z_med;
                  z_med = zb;
                  Gmed = Gb;
                }
              
              nb_iter++;
            }
          
          if (display)
            cout << "value at z = " << z_med << " : " << Gmed << endl;
          
          if (Gmed >= target)
            {
              G1 = Gmed;
              z0 = zprev;
              z1 = z_med;
            }
          else if (Gmed > one)
            {
              //return z_med;
            }
                    
          //DISP(z_med); DISP(Gmed);
        }
#endif
      
      // maximal searched value
      if (z1 < -Real_wp(100))
	return z1;
    }
  
  if (display)
    cout << "Value at z = " << z1 << " : " << G1 << endl;

  // then dichotomy method to find solution of |G| = 1
  Real_wp z, G; int nb_iter = 0;
  while ((abs(z1-z0) > dz_target) && (nb_iter < nb_iter_max))
    {
      z = 0.5*(z0 + z1);
      G = GetAmplificationFactor(alpha, beta, z);
      if (G > target)
	z1 = z;
      else
	z0 = z;
      
      nb_iter++;
    }
  
  z = 0.5*(z0 + z1);
  return z;
}


template<class T>
class SchemeMinimizationCFL : public VirtualMinimizedFunction<T>
{
  int r;
 
public : 
  SchemeMinimizationCFL(int order)
  {
    r = order;
    this->n = (order-2)/2;
    cout << "Order of the time scheme " << r << endl;
    cout << "Number of free parameters " << this->n << endl;
  }

  virtual void FindInitGuess(Vector<T>& x)
  {
    // random initial guess
    x.Reallocate(this->n);
    for (int i = 0; i < this->n; i++)
      x(i) = 2.0*pi_wp*T(rand())/RAND_MAX;
  }

  // finds coefficients alpha and beta from coefficients theta
  bool FindCoefAlphaBeta(const Vector<T>& theta, const T& teta_impl,
			 Vector<T>& alpha, Vector<T>& beta)
  {
    T one(1);
    // on forme le polynome rho
    UnivariatePolynomial<T> rho, zeta, pol;
    zeta.SetOrder(1);
    zeta(0) = 0; zeta(1) = one;
    pol.SetOrder(2);
    rho = (zeta-1)*(zeta-1);
    for (int i = 0; i < theta.GetM(); i++)
      {
        pol(0) = one; pol(1) = -2.0*cos(theta(i)); pol(2) = one;
        rho *= pol;
      }
    
    beta.Reallocate(r+1);
    alpha.Reallocate(r+1);
    beta.Fill(0);
    
    // cas du schema explicite
    int N = beta.GetM()-2;
    // on explicite le systeme lineaire que doit resoudre beta
    Matrix<T> sys(N, N);
    Vector<T> row(N+2), rhs(N);
    row.Fill(one);
    for (int s = 2; s < N+2; s++)
      {
	// row(i) contient i^{s-2}
	// on considere ici que alpha_0 = alpha_N+1 = 1
	// terme alpha_{N+1} (N+1)^s
        rhs(s-2) = row(N+1)*(N+1)*(N+1);
	
	// boucle pour i allant de 1 a N
        for (int i = 1; i < N+1; i++)
          {
	    // terme beta_i i^{s-2} 
            sys(s-2, i-1) = row(i);
	    
	    // terme alpha_i i^s
            rhs(s-2) += rho(i)*row(i)*i*i;
          }
        
	// on divise par s(s-1) car on a au final l'equation :
	// \sum i^{s-2} beta_i = \sum i^s / [s(s-1)] \alpha_i
        rhs(s-2) /= s*(s-1);
        
	// dans le cas implicite, beta_0 = beta_N+1 = teta_impl
	// si s=2, on retranche beta_0+beta_N+1 = 2 teta_impl
	// si s=1, on retranche beta_N+1 (N+1)^s-2
	if (s == 2) 
	  rhs(0) -= 2 * teta_impl;
	else
	  rhs(s-2) -= row(N+1) * teta_impl;
	
	// on multiplie row(i) par i pour le s suivant
        for (int i = 1; i < N+2; i++)
          row(i) *= i;        
      }
    
    // on resout le systeme lineaire
    IVect pivot(N);
    GetLU(sys, pivot);
    
    SolveLU(sys, pivot, rhs);
    
    // on remplit alpha et beta
    alpha.Reallocate(N+2);
    alpha(0) = one; alpha(N+1) = one;
    for (int i = 1; i < N+1; i++)
      alpha(i) = rho(i);
    
    beta(0) = teta_impl; beta(N+1) = teta_impl;
    for (int i = 1; i < N+1; i++)
      beta(i) = rhs(i-1);
    
    return true;
  }
  
  // calcul de la CFL en fonction des parametres libres
  virtual void EvaluateFunction(const Vector<T>& x, T& feval)
  {
    Vector<T> alpha, beta;
    // cas explicite
    bool success = FindCoefAlphaBeta(x, Real_wp(0), alpha, beta);
    if (success)
      feval = GetCFL_Scheme(alpha, beta);
    else
      feval = T(0);
  }
  
  // calcul de la CFL et de son gradient en fonction des parametres
  virtual void EvaluateFunctionGradient(const Vector<T>& x, T& feval, Vector<T>& fjac)
  {
    // cas explicite
    VectReal_wp alpha, beta;
    bool success = FindCoefAlphaBeta(x, Real_wp(0), alpha, beta);
    if (!success)
      {
        feval = T(0);
        fjac.Fill(0);
        return;
      }
    
    feval = GetCFL_Scheme(alpha, beta);
    // derivee numerique
    VectReal_wp xp(x), xm(x);
    Real_wp h = 1e-6, fm, fp;
    for (int i = 0; i < this->n; i++)
      {
	xp(i) = x(i) + h;
	xm(i) = x(i) - h;

	FindCoefAlphaBeta(xp, Real_wp(0), alpha, beta);
	fp = GetCFL_Scheme(alpha, beta);

	FindCoefAlphaBeta(xm, Real_wp(0), alpha, beta);
	fm = GetCFL_Scheme(alpha, beta);
	
	fjac(i) = (fp - fm)/(2.0*h);
	
	xp(i) = x(i); xm(i) = x(i);
      }
  }

  void WriteCoefficients(const VectReal_wp& alpha, const VectReal_wp& beta)
  {
    ofstream file_out("coef.dat"); file_out.precision(16);
    for (int i = 0; i < alpha.GetM(); i++)
      file_out << "alpha(" << i << ") = " << alpha(i) <<"; ";
    
    file_out << endl;
    for (int i = 0; i < beta.GetM(); i++)
      file_out << "beta(" << i << ") = " << beta(i) <<"; ";  
    
    file_out << endl; file_out.close();
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  time_t n = time(NULL); DISP(n);
  srand(n);
  //srand(1394483403);
  
  // class implementing function to minimize
  Vector<Real_wp> xmin; 
  int order = 4;
  SchemeMinimizationCFL<Real_wp> fct(order);
  fct.SetGslAlgorithm(fct.SIMPLEX);

  VectReal_wp alpha, beta;

  VectReal_wp theta(fct.GetM());
  
  //theta(0) = 2.0*pi_wp/5;   theta(1) = 4.0*pi_wp/5;
  //theta(0) = 0.1;
  //theta(0) = 0.112142565901302;
  
  
  //fct.FindCoefAlphaBeta(theta, Real_wp(0.10050251256281401), alpha, beta);
  //DISP(alpha); DISP(beta);  
  //DISP(GetCFL_Scheme(alpha, beta));
  
  //exit(0);

  // calcul du facteur d'amplification lorque qu'on fait varier z
  /*
  VectReal_wp ParamZ;
  Linspace(Real_wp(0), Real_wp(10), 10000, ParamZ);
  
  VectReal_wp Gfactor(ParamZ.GetM());
  for (int i = 0; i < ParamZ.GetM(); i++)
    Gfactor(i) = GetAmplificationFactor(alpha, beta, -ParamZ(i));

  Gfactor.WriteText("factor.dat");
  exit(0); */


  // calcul de la cfl pour deux parametres
  Vector<Real_wp> Cfl;
  VectReal_wp PhiParam, ThetaParam;
  Linspace(Real_wp(1)/20, pi_wp, 200, PhiParam);
  Linspace(Real_wp(0), Real_wp(1), 200, ThetaParam);
  Cfl.Reallocate(PhiParam.GetM());
  Cfl.Fill(0);  

  int j = rank_proc;
  {
    for (int i = 0; i < ThetaParam.GetM(); i++)
      {
        theta(0) = PhiParam(j);
        fct.FindCoefAlphaBeta(theta, ThetaParam(i), alpha, beta);
        //Cfl(i, j) = GetCFL_Scheme(alpha, beta);
        Cfl(i) = GetCFL_Scheme(alpha, beta);
        DISP(i); DISP(j); DISP(Cfl(i));
      }
  }
  
  Vector<int64_t> xtmp;
  VectReal_wp AllCfl;
  if (rank_proc == 0)
    AllCfl.Reallocate(ThetaParam.GetM()*Cfl.GetM());
  
  MpiGather(MPI_COMM_WORLD, Cfl, xtmp, AllCfl, Cfl.GetM(), 0);
  
  if (rank_proc == 0)
    {
      ThetaParam.WriteText("theta_t.dat");
      PhiParam.WriteText("phi_t.dat");
      AllCfl.WriteText("cfl_t.dat");
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // cas avec deux parametres
  /* VectReal_wp ThetaParam;
  Matrix<Real_wp> Cfl;
  Linspace(Real_wp(0), pi_wp, 500, ThetaParam);

  Cfl.Reallocate(ThetaParam.GetM(), ThetaParam.GetM());
  Cfl.Fill(0);
  for (int i = 0; i < ThetaParam.GetM(); i++)
    {
      DISP(i);
      for (int j = 0; j < ThetaParam.GetM(); j++)
        {
          theta(0) = ThetaParam(i);
          theta(1) = ThetaParam(j);
          fct.FindCoefAlphaBeta(theta, alpha, beta);
          Cfl(i, j) = GetCFL_Scheme(alpha, beta);
        }  
    }
  
  ThetaParam.WriteText("theta_c.dat");
  Cfl.WriteText("cfl_c.dat");
  */
#ifndef MONTJOIE_WITH_MPFR
  Real_wp fmin(0), feval; 
  for (int k = 0; k < 10; k++)
    {
      // launching a minimization at each step
      // the initial guess is set as random in FindInitGuess
      Vector<Real_wp> xsol;
      
      feval = MinimizeParametersGsl(fct, xsol);
      
      if (feval < fmin)
	{
	  // displays the new minimum reached
	  xmin = xsol;
	  fmin = feval;
	  DISP(xmin); DISP(sqrt(abs(fmin)));
	  fct.FindCoefAlphaBeta(xmin, Real_wp(0), alpha, beta);
          DISP(GetCFL_Scheme(alpha, beta));
          DISP(alpha); DISP(beta);
	}
    }

  DISP(xmin);
  DISP(sqrt(abs(fmin)));
#endif
    
  return FinalizeMontjoie();
}
