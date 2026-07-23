#include "Share/MontjoieCommon.hxx"
#include "Solver/NonLinearEquations.hxx"
#include "Solver/NonLinearEquations.cxx"

using namespace Montjoie;

int rank_proc(0), nb_proc(1);

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


// returns the amplification factor for a symmetric multistep scheme
// defined through coefficients alpha, beta
void GetRoots(const VectReal_wp& alpha, const VectReal_wp& beta,
              const Real_wp& z, VectComplex_wp& all_roots)
{
  UnivariatePolynomial<Real_wp> P;
  P.SetOrder(alpha.GetM()-1);
  
  for (int i = 0; i < alpha.GetM(); i++)
    P(i) = alpha(i) - beta(i)*z ;
  
  //DISP(P);
  VectReal_wp R, Rimag;
  GetRoots(P, R, Rimag);

  all_roots.Reallocate(R.GetM());
  for (int i = 0; i < R.GetM(); i++)
    all_roots(i) = Complex_wp(R(i), Rimag(i));
}

void FindPermutationCloseRoots(const VectComplex_wp& prev_roots, const VectComplex_wp& new_roots,
                               IVect& permut)
{
  int N = prev_roots.GetM();
  Vector<bool> root_selected(N);
  root_selected.Fill(false);
  permut.Reallocate(N);
  for (int j = 0; j < N; j++)
    {
      Real_wp dist_min = 1e300;
      int pos = -1;
      for (int k = 0; k < N; k++)
        if (!root_selected(k))
          {
            Real_wp dist = abs(prev_roots(j) - new_roots(k));
            if (dist < dist_min)
              {
                pos = k;
                dist_min = dist;
              }
          }

      permut(j) = pos;
      root_selected(pos) = true;
    }
}

Real_wp GetDifferenceAngle(const Complex_wp& z1, const Complex_wp& z2)
{
  Real_wp angle1 = arg(z1);
  Real_wp angle2 = arg(z2);
  if (angle2 < angle1-pi_wp)
    angle2 += 2*pi_wp;
  
  if (angle2 > angle1 + pi_wp)
    angle2 -= 2*pi_wp;
  
  return angle2 - angle1;
}


Real_wp GetIntersectionZ(const Complex_wp& pt1, const Complex_wp& pt2,
                         const Real_wp& v1, const Real_wp& v2)
{
  if (v1 == v2)
    return Real_wp(1);
  
  Real_wp angle1 = arg(pt1);
  Real_wp angle2 = arg(pt2);
  Real_wp zmin = 1e300;
  Real_wp zi = (angle2 - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = zi;
  
  zi = (angle2 + 2*pi_wp - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = min(zmin, zi);
  
  zi = (angle2 - 2*pi_wp - angle1) / (v1 - v2);
  if (zi > 0)
    zmin = min(zmin, zi);          
  
  return zmin;
}


// returns the CFL of the symmetric multistep scheme
Real_wp GetCFL_Scheme(const VectReal_wp& alpha, const VectReal_wp& beta, bool display = false)
{
  Real_wp one(1), zero(0);
  Real_wp dz_min(1e-6), dz_max(1);
  Real_wp z0 = -dz_min;
  if (epsilon_machine < 1e-30)
    {
      z0 = -Real_wp(1e-12);
      dz_min = Real_wp(1e-8);
    }
  
  Real_wp target = one + 1e-12;
  VectComplex_wp all_roots, new_roots;
  IVect permut;
  GetRoots(alpha, beta, z0, all_roots);
  Real_wp G0 = all_roots.GetNormInf(), dz_target = 1e-12;
  int nb_iter_max(50);
  
  if (epsilon_machine < 1e-30)
    target = one + 1e-15;  
  
  if (G0 > target)
    {
      // instable algorithm
      return zero;
    }
  
  // searching first point with |G| > 1
  Real_wp z1 = z0, G1 = G0;
  while (G1 < target)
    {      
      z0 = z1; G0 = G1;
      
      // roots at z-dz are evaluated to estimate the velocity angle for each root
      GetRoots(alpha, beta, z0-dz_min, new_roots);     
      FindPermutationCloseRoots(all_roots, new_roots, permut);
      VectReal_wp angle_velocity(all_roots.GetM());
      for (int i = 0; i < all_roots.GetM(); i++)
        {
          Real_wp diff_angle = GetDifferenceAngle(all_roots(i), new_roots(permut(i)));
          angle_velocity(i) = diff_angle / dz_min;
        }

      // setting the new dz by computing the angle where eventual intersections may occur
      Real_wp dist_min = 1e300;
      for (int i = 0; i < all_roots.GetM(); i++)
        for (int j = 0; j < all_roots.GetM(); j++)
          if (i != j)
            {
              Real_wp dist = GetIntersectionZ(all_roots(i), all_roots(j),
                                              angle_velocity(i), angle_velocity(j));
              
              dist_min = min(dist_min, dist);
            }
      
      Real_wp dz = dist_min/2;
      if (dz < dz_min)
        dz = dz_min;
      else if (dz > dz_max)
        dz = dz_max;
      
      // evaluating the roots for the next point and amplification factor
      z1 -= dz;      
      GetRoots(alpha, beta, z1, all_roots);
      G1 = all_roots.GetNormInf();

      // maximal searched value
      if (z1 < -Real_wp(100))
	return -Real_wp(100);
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
    if (rank_proc == 0)
      {
        cout << "Order of the time scheme " << r << endl;
        cout << "Number of free parameters " << this->n << endl;
      }
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
    
    // symmetric coefficients : 
    // alpha = (alpha_N, alpha_N-1, ..., alpha_0, alpha_1, ..., alpha_N)
    // beta = (beta_N, beta_N-1, ..., beta_0, beta_1, ..., beta_N)
    // taille du systeme : on resout que la moitie des inconnues beta_0, beta_N-1
    int N = r/2;

    // on explicite le systeme lineaire que doit resoudre beta_0, ..., beta_N-1
    Matrix<T> sys(N, N);
    Vector<T> powI(N+1), rhs(N);
    sys.Zero();
    powI.Zero();
    rhs.Zero();
    
    // premiere equation : \alpha_0 + 2 \sum_{i=1}^N alpha_i = 0
    // deja verifiee par 1 est racine double de rho

    // deuxieme equation : 2 \sum_{i=1}^N i^2/2 alpha_i = \beta_0 + 2 \sum_{i=1}^N \beta_i
    sys(0, 0) = one;  
    
    // le terme 2 \beta_N est reverse dans le second membre car beta_N = teta_implicite
    rhs(0) = -2*teta_impl;
    
    // terme N^2 alpha_N
    rhs(0) += N*N*rho(0); powI(N) = N*N;
    for (int i = 1; i < N; i++)
      {
        // terme 2 beta_i
        sys(0, i) = T(2);

        // terme i^2 alpha_i
        rhs(0) += i*i*rho(N-i);

        // on initialise powI(i) = i^2 pour la boucle suivante
        powI(i) = i*i;
      }
    
    // equation pour obtenir un ordre s
    // 1/[s(s-1)] \sum_{i=1}^N i^s alpha_i = \sum_{i=1}^N i^{s-2} \beta_i 
    int num = 1;
    for (int s = 4; s <= r; s += 2)
      {
        // coef contient 1 / (s(s-1))
        T coef = one / (s*(s-1));
        
	// powI(i) contient i^{s-2}
        
        // terme \beta_N N^{s-2} est reverse dans le second membre car beta_N = teta_implicite
        rhs(num) = -powI(N)*teta_impl;
        
        // terme N^s alpha_N / (s (s-1))
        rhs(num) += powI(N) * N * N * rho(0) * coef;
        for (int i = 1; i < N; i++)
          {
	    // terme beta_i i^{s-2} 
            sys(num, i) = powI(i);
	    
	    // terme alpha_i i^s / [s(s-1)]
            rhs(num) += rho(N-i) * powI(i) * i * i * coef;
          }
	
	// on multiplie row(i) par i*i pour le s+2 suivant
        for (int i = 1; i <= N; i++)
          powI(i) *= i*i;        
        
        num++;
      }
    
    // on resout le systeme lineaire
    IVect pivot(N);
    GetLU(sys, pivot);
    
    SolveLU(sys, pivot, rhs);
    
    // on remplit alpha et beta
    for (int i = 0; i < alpha.GetM(); i++)
      alpha(i) = rho(i);    
    
    beta(0) = teta_impl; beta(2*N) = teta_impl;
    beta(N) = rhs(0);
    for (int i = 1; i < N; i++)
      {
        beta(N+i) = rhs(i);
        beta(N-i) = rhs(i);
      }
    
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
  cout.precision(16);

#ifdef SELDON_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#endif

  /* if (argc != 2)
    {
      cout << "Entrez le parametre i " << endl;
      abort();
    }

    int num_angle = atoi(argv[1]); */
  if (rank_proc == 0)
    {
      cout << "Epsilon machine = " << epsilon_machine << endl;
      //cout << "Angle = " << num_angle << endl;
    }
  
  time_t n = time(NULL); //DISP(n);
  srand(n);
  //srand(1394483403);
  
  // class implementing function to minimize
  Vector<Real_wp> xmin; 
  int order = 8;
  SchemeMinimizationCFL<Real_wp> fct(order);
  fct.SetGslAlgorithm(fct.SIMPLEX);
  
  VectReal_wp alpha, beta;
  VectReal_wp theta(fct.GetM());
  theta.Zero();
  
  // exemple de calcul de la CFL en un point particulier
  //theta(0) = pi_wp/2;

  //theta(0) = pi_wp/2;
  //theta(1) = 2*pi_wp/3;

  theta(0) = pi_wp/6;
  theta(1) = pi_wp/2;
  theta(2) = 2*pi_wp/3;
  
  fct.FindCoefAlphaBeta(theta, Real_wp(0), alpha, beta);
  DISP(alpha); DISP(beta);  
  DISP(GetCFL_Scheme(alpha, beta));
  
  // calcul du facteur d'amplification lorque qu'on fait varier z
  
  VectReal_wp ParamZ;
  Linspace(Real_wp(0.0), Real_wp(4), 10001, ParamZ);
  
  VectReal_wp Gfactor(ParamZ.GetM());
  Matrix<Real_wp> Groot_real(ParamZ.GetM(), order);
  Matrix<Real_wp> Groot_imag(ParamZ.GetM(), order);
  VectComplex_wp all_roots, prev_roots(order); IVect permut(order);
  permut.Fill();
  for (int i = 0; i < ParamZ.GetM(); i++)
    {
      if (i%100 == 0)
        DISP(i);

      GetRoots(alpha, beta, -ParamZ(i), all_roots);
      Gfactor(i) = all_roots.GetNormInf();
      
      if (i > 0)
        FindPermutationCloseRoots(prev_roots, all_roots, permut);
      for (int j = 0; j < order; j++)
        {
          Groot_real(i, j) = real(all_roots(permut(j)));
          Groot_imag(i, j) = imag(all_roots(permut(j)));
          prev_roots(j) = Complex_wp(Groot_real(i, j), Groot_imag(i, j));
        }
    }
  
  Gfactor.WriteText("factor.dat");
  Groot_real.WriteText("real_part.dat");
  Groot_imag.WriteText("imag_part.dat");

  exit(0);

  // calcul de la cfl pour un parametre
  /* Vector<Real_wp> Cfl;
  VectReal_wp PhiParam;
  Linspace(Real_wp(0.583364)*pi_wp, Real_wp(0.583368)*pi_wp, 1001, PhiParam);
  Cfl.Reallocate(PhiParam.GetM());
  Cfl.Fill(0);  
  
  for (int i = 0; i < PhiParam.GetM(); i++)
    {
      theta(0) = PhiParam(i);
      theta(1) = PhiParam(i);
      fct.FindCoefAlphaBeta(theta, Real_wp(0), alpha, beta);
      Cfl(i) = GetCFL_Scheme(alpha, beta);
      DISP(i); DISP(Cfl(i));
    }
  
  Cfl.WriteText("cfl.dat");  
  return FinalizeMontjoie(); */
  
  // calcul de la cfl pour deux parametres
  /* Vector<double> Cfl;
  VectReal_wp PhiParam, ThetaParam;
  Linspace(Real_wp(0), pi_wp, 201, PhiParam);
  Linspace(Real_wp(0), pi_wp, 201, ThetaParam);
  Cfl.Reallocate(ThetaParam.GetM());
  Cfl.Zero();  
  
  int j = rank_proc;
  //for (int j = 0; j < PhiParam.GetM(); j++)
  {
    for (int i = 0; i < ThetaParam.GetM(); i++)
      {
        theta(0) = PhiParam(j);
        theta(1) = ThetaParam(i);
        theta(2) = num_angle*pi_wp/180; //DISP(theta);
        fct.FindCoefAlphaBeta(theta, Real_wp(0), alpha, beta);
        //Cfl(i, j) = GetCFL_Scheme(alpha, beta);
        Cfl(i) = toDouble(GetCFL_Scheme(alpha, beta));
        //DISP(i); DISP(j); DISP(Cfl(i));
      }
  }
  
  Vector<int64_t> xtmp;
  Vector<double> AllCfl;
  if (rank_proc == 0)
    AllCfl.Reallocate(ThetaParam.GetM()*PhiParam.GetM());
  
  MpiGather(MPI_COMM_WORLD, Cfl, xtmp, AllCfl, Cfl.GetM(), 0);
  
  if (rank_proc == 0)
    {
      //ThetaParam.WriteText("theta_mp.dat");
      //PhiParam.WriteText("phi_mp.dat");
      AllCfl.WriteText("cfl_" + to_str(num_angle) + ".dat");
    }
  
    MPI_Barrier(MPI_COMM_WORLD); */
  
#ifndef MONTJOIE_WITH_MPFR
  // boucle de minimization
  /* Real_wp fmin(0), feval; 
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
  DISP(sqrt(abs(fmin))); */
#endif
    
  return FinalizeMontjoie();
}
