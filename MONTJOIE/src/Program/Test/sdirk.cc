#include "Solver/MontjoieSolver.hxx"
#include "Output/MontjoieOutput.hxx"
#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

bool find_cfl = false;

void CheckScheme(const VectReal_wp& param, int s);

// binomial coefficients C_n^k
// we fill the row for a given n
void GetCnk(int n, Vector<int>& row)
{
  row.Reallocate(n+1);
  row.Fill(1);

  if (n <= 1)
    return;

  for (int i = 2; i <= n; i++)
    {
      for (int j = i-1; j >= 1; j--)
        row(j) = row(j-1) + row(j);
    }
}

// we compute Pmod = |P(i \sqrt{x})|^2
void GetSquareModulus(const UnivariatePolynomial<Real_wp>& P,
		      UnivariatePolynomial<Real_wp>& Pmod)
{
  UnivariatePolynomial<Real_wp> Pr, Pi, z;
  z.SetOrder(1); z(0) = Real_wp(0); z(1) = Real_wp(1);
  
  // computing Pr(i \sqrt{z}) (to avoid high degrees)
  Pr.SetOrder(P.GetOrder()/2);  
  for (int i = 0; i <= P.GetOrder(); i += 2)
    {
      Pr(i/2) = P(i);
      if (i%4 == 2)
	Pr(i/2) = -P(i);
    }

  // computing \sqrt{z} Pi(i \sqrt{z}) (to avoid high degrees)
  Pi.SetOrder((P.GetOrder()-1)/2);  
  for (int i = 1; i <= P.GetOrder(); i += 2)
    {
      Pi((i-1)/2) = P(i);
      if (i%4 == 3)
	Pi((i-1)/2) = -P(i);
    }

  // modulus |P(i sqrt(z))|^2
  Pmod = Pr*Pr + z*Pi*Pi;
  
}

// returns the CFL associated with the stability function R = P / Q
Real_wp GetCFL(const UnivariatePolynomial<Real_wp>& P,
	       const UnivariatePolynomial<Real_wp>& Q,
               bool display = false)
{
  // finding the values z for which |P|^2 = |Q|^2
  UnivariatePolynomial<Real_wp> Pmod, Qmod, Pdiff;

  GetSquareModulus(P, Pmod);
  GetSquareModulus(Q, Qmod);
  
  if (display) { DISP(Pmod); DISP(Qmod); }
  bool mod_different = false;
  if (Pmod.GetM() == Qmod.GetM())
    {
      for (int i = 0; i < Pmod.GetM(); i++)
	{
	  if (abs(Pmod(i) - Qmod(i)) > 10*sqrt(epsilon_machine))
	    {
	      if (!mod_different)
		{
		  // here we manage the bad slope at the origin
		  if (Pmod(i) > Qmod(i))
		    return Real_wp(0);
		}
	      
	      mod_different = true;
	    }
	}
    }
  else
    mod_different = true;

  if (display) { DISP(mod_different); }
  // if the modulus of P and Q are equal, the CFL is infinite
  if (!mod_different)
    return Real_wp(1e300);

  // otherwise we find the roots of |P|^2 - |Q|^2
  UnivariatePolynomial<Real_wp> R;
  VectReal_wp Li, Lr;
  R = Pmod - Qmod;
  
  for (int i = 0; i <= R.GetOrder(); i++)
    if (abs(R(i)) <= 1e4*epsilon_machine)
      R(i) = Real_wp(0);
  
  //DISP(epsilon_machine);
  SolvePolynomialEquation(R, Lr, Li);

  Sort(Lr, Li);
  if (display) { DISP(R); DISP(Lr); DISP(Li); }
  
  Real_wp cfl(1e300);
  int num_root = 0;
  while (num_root < Lr.GetM())
    {
      // we select only real roots different from 0
      if ((abs(Li(num_root)) <= epsilon_machine) && (Lr(num_root) >= Real_wp(1e-12)))
	{
	  int multiplicity = 1;
	  while ((num_root < Lr.GetM()-1) && (abs(Lr(num_root+1) - Lr(num_root)) <= 1e4*epsilon_machine))
	    {
	      multiplicity++;
	      num_root++;
	    }
          
          // for odd multiplicity, the CFL is lower than this root
	  if (multiplicity%2 == 1)
	    cfl = min(cfl, sqrt(Lr(num_root)));
	}

      num_root++;
    }

  return cfl;
}

// on calcule la fonction de stabilite R = Numer/Denom
// en fonction de s (la methode est d'ordre s+1) et des parametres
// param(0) = gamma
// param(1:) = les constantes alpha, beta, etc
void GetPolynomial_Sdirk(const VectReal_wp& param, int s,
                         UnivariatePolynomial<Real_wp>& Numer,
                         UnivariatePolynomial<Real_wp>& Denom,
                         UnivariatePolynomial<Real_wp>& Pol)
{
  // premiere parametre : gamma
  Real_wp Gamma = param(0);
  int p = param.GetM() + s;
  
  // Pol = 1 - gamma z
  Pol.SetOrder(1); Pol(0) = Real_wp(1); Pol(1) = -Gamma;
  
  // Numer : developpement limite de l'exponentielle
  Numer.SetOrder(s+1);
  Real_wp invFacto(1);
  Numer(0) = Real_wp(1);
  for (int i = 0; i <= s; i++)
    {      
      Numer(i+1) = invFacto;
      invFacto *= Real_wp(1)/(i+2);
    }
  
  // Denom = (1-gamma z)^p
  Denom = Pow(Pol, p);
  
  // Pol = (1 - gamma z)^p exp(z) en enlevant les termes de degre superieur
  Pol = Numer*Denom;

  // on ajoute les parametres additionels
  for (int i = s+1; i < p; i++)
    Pol(i+1) += param(i-s);
  
  Numer = Pol;
  Numer.ResizeOrder(p);
}

void ComputeErrorSource(const UnivariatePolynomial<Real_wp>& NumerR,
                        const UnivariatePolynomial<Real_wp>& DenomR,
                        int order, int nb_coef, VectReal_wp& coef)
{
  int nb_terms_add = order+1;
  // calcul de 1/k! et 1 / 2^k k!
  VectReal_wp invFacto(2*nb_terms_add+1), CoefC(2*nb_terms_add+1);
  invFacto(0) = Real_wp(1); CoefC(0) = Real_wp(1);
  for (int k = 1; k <= 2*nb_terms_add; k++)
    {
      invFacto(k) = invFacto(k-1) / Real_wp(k);
      CoefC(k) = CoefC(k-1) / Real_wp(2*k);
    }
  
  // terms with X^{0} (to see the error constants for the homogeneous case
  VectReal_wp err(nb_terms_add);
  err.Zero();
  /* 
  for (int r = 0; r < nb_terms_add; r++)
    {
      Real_wp sum(0);
      for (int i = 0; i <= min(r, NumerR.GetOrder()); i++)
        {
          Real_wp signK(1);
          if ((r-i)%2 == 1)
            signK = -signK;
          
          sum += (DenomR(i) - signK*NumerR(i))*CoefC(r-i);
        }
      
      err(r) = sum;
    }
  
    cout << "Errors for the homogeneous case" << endl << err << endl; */
  
  // computation of coefficients alpha_j^r for the right hand side phi
  // phi = dt \sum_{r=1}^\infty (dt A)^{r-1} \sum_{j=1}^\infty \alpha_j^r dt^{j-1} F^{j-1}
  // where \alpha_j^r = \sum_{i=0}^{min(r-1, m)} (D_i - (-1)^k N_i) / 2^k k!
  // where k = r+j-i-1
  //cout << endl << "Coefficients alpha for inhomogeneous case" << endl;
  VectReal_wp xtmp;
  for (int r = 1; r <= nb_terms_add; r++)
    {
      err.Zero();
      for (int j = 1; j <= nb_terms_add+1-r; j++)
        {
          for (int i = 0; i <= min(r-1, NumerR.GetOrder()); i++)
            {
              int k = r + j-i - 1;
              Real_wp signK(1);
              if (k%2 == 1)
                signK = -signK;
              
              err(j-1) += (DenomR(i) - signK*NumerR(i))*CoefC(k);
            }

          xtmp.PushBack(err(j-1));
        }
      
      //cout << "Terms with (dt A)^" << r-1 << endl << err << endl;
    }

  coef.Reallocate(nb_coef);
  for (int i = 0; i < nb_coef; i++)
    coef(i) = xtmp(xtmp.GetM()-1-i);
}


// on calcule la CFL d'un schema de type SDIRK
// en fonction de s (la methode est d'ordre s+1) et des parametres
// param(0) = gamma
// param(1:) = les constantes alpha, beta, etc
Real_wp GetCFL_Sdirk(const VectReal_wp& param, int s, bool print_pol = false)
{
  UnivariatePolynomial<Real_wp> Denom, Numer, Pol;
  GetPolynomial_Sdirk(param, s, Numer, Denom, Pol);
  if (print_pol)
    {
      DISP(Numer); DISP(Denom);
    }
  
  Real_wp cfl = GetCFL(Numer, Denom);
  return cfl;
}


Real_wp GetError_Sdirk(const VectReal_wp& xb, int s, Real_wp& cfl)
{
  UnivariatePolynomial<Real_wp> NumerR, DenomR, Pol;
  UnivariatePolynomial<Real_wp> exp_taylor;
  
  GetPolynomial_Sdirk(xb, s, NumerR, DenomR, Pol);
  //DISP(NumerR); DISP(DenomR); DISP(Pol);
  cfl = GetCFL(NumerR, DenomR);
  //CheckScheme(xb, s);
  //int test_input; cout << "waiting" << endl; cin >> test_input;
  
  if (xb.GetM() >= 2)
    Pol(s+2) = -xb(1);
  
  Real_wp invFacto(1);
  for(int j = 0; j <= s; j++)
    invFacto *= Real_wp(1)/(j+2);
  
  //DISP(invFacto);
  Real_wp err = abs(invFacto + Pol(s+2));
  //DISP(err); DISP(cfl);

  //VectReal_wp coef;
  //ComputeErrorSource(NumerR, DenomR, s, xb.GetM(), coef);
  //err = Norm1(coef);
  
  //int test_input; cout << "Waiting" << endl; cout << endl; cin >> test_input;
  return err;
    
  // val_ref: developpement limite de l'exponentielle
  //exp_taylor.SetOrder(s + 2);
  //exp_taylor(0) = invFacto;
  //for (int i = 0; i <= s+1; i++)
  //{      
  //exp_taylor(i+1) = invFacto;
  //invFacto *= Real_wp(1)/(i+2);
  //}
  //DISP(exp_taylor);
  
  // We evaluate error in the imaginary axis to be minimize
  /*VectReal_wp z;
    Linspace(Real_wp(0), Real_wp(50), 10000, z);
    
    for (int i = 0; i < 1; i++)
    {
    Complex_wp zi = Iwp*z(i);
    Complex_wp val_ref = exp_taylor(s+1)*zi + exp_taylor(s);
    
    // Horner algorithm
    int order = NumerR.GetOrder();
    Complex_wp num = NumerR(order);
    Complex_wp denom = DenomR(order);
    for (int j = s; j >= 0; j--)
    {
    num = num*zi + NumerR(j);
    denom = denom*zi + DenomR(j);
    val_ref = val_ref*zi + exp_taylor(j);
    }
    err = abs(num/denom - val_ref);
    feval += toDouble(err);
    }*/
}


// classe de base pour definir une sequence N-D de parametres
class SequenceParameter
{
public:
  // nombre de parametres (dimension)
  virtual int GetNbParam() const = 0;
  // nombre de points totaux a balayer
  virtual int GetNbPoints() const = 0;
  // on est sur le point k a balayer
  virtual void SetParameters(int k, VectReal_wp& x) = 0;
  
};


//! class for minimizing error with with free parameters
/*!
  This class implements the computation of the free parameters gamma, alpha
  ... which minimize the global error e^z - R(z)
*/
template<class T>
class SdirkErrorMinimization : public VirtualMinimizedFunction<double>
{
  // pour changer s, appeler SetOrder
  // order = s+1
  int s;
  int type_init;
  enum {RANDOM, USER};
  Vector<double> params_init;
  Vector<double> coef_scale;

public : 
  //! default constructor
 SdirkErrorMinimization()
  {
    s = 3;
    type_init = USER;
    this->n = 1;
    //this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX;
    this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX2_RAND;
    //this->type_algo = VirtualMinimizedFunction<T>::CG;
    //this->type_algo = VirtualMinimizedFunction<T>::CG_FR;
    //this->type_algo = VirtualMinimizedFunction<T>::STEEPEST_DESCENT;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS2;

    this->step_size = 0.001;
  }
  
  void SetInitialGuess(const Vector<double>& params, const Vector<T>& coefs)
  {
    type_init = USER;
    params_init = params;
    coef_scale.Reallocate(coefs.GetM());
    for (int i = 0; i < coef_scale.GetM(); i++)
      {
        coef_scale(i) = toDouble(coefs(i));
        params_init(i) /= coef_scale(i);
      }
  }

  //! initial guess for free parameters
  void FindInitGuess(Vector<double>& params)
  {
    if (type_init == USER)
      {
	params = params_init;
	return;
      }

    params.Reallocate(this->n);
    //x(0) = 2.0*T(rand())/RAND_MAX; params(1) = 2.0*T(rand())/RAND_MAX;
    for (int i = 0; i < this->n; i++)
      params(i) = double(rand())/RAND_MAX;
    //params(i) = 3.0*T(rand())/RAND_MAX - 1.5;
    
  }

  //! initializes the number of free parameters
  void SetOrder(int s, int n)
  {
    this->s = s;
    this->n = n; //number of free parameters
  }

  //! computes feval = first_non_zero_coeff(e^z - R(z)) from free parameters contained in params
  void EvaluateFunction(const Vector<double>& param, double& feval)
  {
    int nbParam = param.GetM();
    Vector<T> xb(nbParam);
    for (int i = 0; i < nbParam; i++)
      xb(i) = param(i)*coef_scale(i);

    Real_wp cfl;
    Real_wp err = GetError_Sdirk(xb, this->s, cfl);
    if (find_cfl)
      {
        feval = -toDouble(cfl);
      }
    else
      {
        feval = toDouble(err);
        
        // to penalize non A-stable schemes, we add 1e300 if the scheme is not A-stable
        if(cfl < Real_wp(1e100))
          feval = toDouble(1e300);	  
      }

    //DISP(param); DISP(find_cfl); DISP(err); DISP(feval);
    //cout << "waiting" << endl; int test_input; cin >> test_input;
  }

  //! computes the gradient of EvaluateFunction
  void EvaluateFunctionGradient(const Vector<double>& param, 
				double& feval_d, Vector<double>& fjac)
  {
  }

};

// function to compute optimal parameters regardind cfl and approximation error
void ComputeCFLOptimal(SdirkErrorMinimization<Real_wp>& fct,
                       SequenceParameter& seq, int nb_selection, int s,
                       const VectReal_wp& coef_scale, VectReal_wp& xmin_opt)
{
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int nb_proc(1), rank_proc(0);
#endif

  // points are distributed to the different processors
  int N = seq.GetNbPoints(), nb_pts_proc, offset_proc;
  GetParallelDistributionPoints(nb_proc, rank_proc, N, nb_pts_proc, offset_proc);
  
  Vector<int64_t> xtmp;
  
  // main loop over points to find best candidates
  Vector<double> Cfl(nb_selection); Cfl.Fill(1e300);
  if (find_cfl)
    Cfl.Zero();
  
  Vector<int> NumParam(nb_selection); NumParam.Fill(-1);
  VectReal_wp xmin(fct.GetM()); xmin.Zero();
  int old_percent = 0, percent = 0;
  //DISP(nb_selection); DISP(N);
  for (int k = 0; k < N; k++)
    {
      if ((k >= offset_proc) && (k < offset_proc+nb_pts_proc))
        {
          if (rank_proc == 0)
            {
              old_percent = percent;
              percent = toInteger(round(100*(k-offset_proc) / Real_wp(nb_pts_proc)));
              if (percent != old_percent)
                {
                  cout << percent << "%  ";
                }
            }

	  seq.SetParameters(k, xmin);
          
          // we check that the parameters are all different
         
          Real_wp cfl;
	  Real_wp err = GetError_Sdirk(xmin, s, cfl);
          //DISP(xmin); DISP(err); DISP(cfl);
	  // a quel rang se situe cette cfl ?
	  int rank = nb_selection;
          if (find_cfl)
            {
              for (int p = nb_selection-1; p >= 0; p--)
                if (cfl > Cfl(p))
                  rank = p;              
            }
          else if (cfl >= Real_wp(1e299))
            for (int p = nb_selection-1; p >= 0; p--)
              if (err < Cfl(p))
                rank = p;
          
	  if (rank < nb_selection)
	    {
	      // on insere ce nouveau optimum
	      for (int p = nb_selection-2; p >= rank; p--)
		{
		  NumParam(p+1) = NumParam(p);
		  Cfl(p+1) = Cfl(p);
		}
	      
	      NumParam(rank) = k;
	      if (find_cfl)
                Cfl(rank) = toDouble(cfl);
              else
                Cfl(rank) = toDouble(err);
	    }
	}
    }
  
  DISP(Cfl); DISP(NumParam);
  cout << "Starting optimization" << endl;
  
  // then optimization is performed for each candidate  
  Vector<double> x_opt(xmin.GetM());
#ifdef MONTJOIE_WITH_GSL
  double fmin(-1), feval; 
  for (int k = 0; k < nb_selection; k++)
    if (NumParam(k) != -1)
      {      
        Vector<double> xsol(xmin.GetM()); xsol.Zero();
        
        VectReal_wp xsol_(xmin.GetM()); xsol_.Zero();
        seq.SetParameters(NumParam(k), xsol_);
        
        for (int k = 0; k < xsol.GetM(); k++)
          xsol(k) = toDouble(xsol_(k));
        
        fct.SetInitialGuess(xsol, coef_scale);
        
        MinimizeParametersGsl(fct, xsol);
        
        fct.EvaluateFunction(xsol, feval);
        //DISP(xsol); DISP(feval);
        //if(fmin<0)
        //fmin = feval;
        
        //DISP(fmin);
        //DISP(xsol);
        
        if (feval < fmin || fmin < 0)
          {
            x_opt = xsol;
            fmin = feval;
            //DISP(x_opt);
            // DISP(fmin);
          }
      }
#else
  if (NumParam(0) != -1)
    seq.SetParameters(NumParam(0), xmin);
  
  x_opt.Reallocate(xmin.GetM());
  for (int k = 0; k < xmin.GetM(); k++)
    x_opt(k) = toDouble(xmin(k));
  //DISP(x_opt);
#endif  

  for (int k = 0; k < xmin.GetM(); k++)
    xmin(k) = x_opt(k)*coef_scale(k);

  // parameters and cfl are gathered
  Vector<double> AllCfl;
  Vector<double> AllParam;
  double cfl_opt(0); Real_wp cfl;
  
  Real_wp cfl_sdirk = GetError_Sdirk(xmin, s, cfl);
  cfl_opt = toDouble(cfl_sdirk);
  if (find_cfl)
    cfl_opt = toDouble(cfl);
  
  if (rank_proc == 0)
    {
      AllCfl.Reallocate(nb_proc);
      AllParam.Reallocate(nb_proc*xmin.GetM());
    }
  
#ifdef SELDON_WITH_MPI
  MpiGather(MPI_COMM_WORLD, &cfl_opt, xtmp, AllCfl.GetData(), 1, 0);
  MpiGather(MPI_COMM_WORLD, x_opt, xtmp, AllParam, xmin.GetM(), 0);
#else
  cout << "Recompile with MPI" << endl;
  abort();
#endif
  
  if (rank_proc == 0)
    {
      Vector<Vector<double> > all_optParam(nb_proc);
      int num = 0;
      for (int i = 0; i < nb_proc; i++)
        {
          all_optParam(i).Reallocate(xmin.GetM());
          for (int j = 0; j < all_optParam(i).GetM(); j++)
            all_optParam(i)(j) = AllParam(num++);
        }
      
      Sort(AllCfl, all_optParam);
      
      for (int k = nb_proc-1; k >= 0; k--)
        {
	  cout << "For processor " << k << endl;
          if (find_cfl)
            cout << "CFL = " << AllCfl(k) << "\t" ;// endl;
          else
            cout << "Error = " << AllCfl(k) << "\t" ;// endl;
          cout << "For parameters = " << all_optParam(k) << "\t";// endl;
	  double feval;
	  fct.EvaluateFunction(all_optParam(k), feval);
	  cout << "error = " << feval << endl;
        }

      int k0 = 0;
      if (find_cfl)
        k0 = all_optParam.GetM() - 1;

      xmin_opt.Reallocate(all_optParam(k0).GetM());
      for (int i = 0; i < all_optParam(k0).GetM(); i++)
        xmin_opt(i) = all_optParam(k0)(i);
    }
  
#ifdef SELDON_WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}


// on calcule la CFL pour les parametres de la sequence
// s : la methode est d'ordre s+1
// seq : definition de la sequence de parametres a tester
// CflP : CFL obtenue pour toutes les valeurs de la sequence
void ComputeCFLParam(SequenceParameter& seq, Vector<double>& CflP, Vector<double>& ErrP, int s,
                     const VectReal_wp& coef_scale)
{
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int nb_proc(1), rank_proc(0);
#endif

  // points are distributed to the different processors
  int N = seq.GetNbPoints(), nb_pts_proc, offset_proc;
  GetParallelDistributionPoints(nb_proc, rank_proc, N, nb_pts_proc, offset_proc);

  srand(0);
  VarRandomGenerator var;
  IVect permut(N); permut.Fill();
  var.GenerateRandomPermutation(N, permut);
  
  // in the root processor, we retrieve the number of points
  // computed on each proc and offset
  Vector<int64_t> xtmp;
  IVect NbPointsPerProc, OffsetPointsProc;
  // nmax is the maximum number of points (among all proc)
  int nmax = N;
  NbPointsPerProc.Reallocate(nb_proc);
  OffsetPointsProc.Reallocate(nb_proc);
  if (nb_proc > 1)
    {
#ifdef SELDON_WITH_MPI
      MPI_Gather(&nb_pts_proc, 1, MPI_INTEGER,
                 NbPointsPerProc.GetData(), 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      MPI_Gather(&offset_proc, 1, MPI_INTEGER,
                 OffsetPointsProc.GetData(), 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      nmax = NbPointsPerProc.GetNormInf();
      
      MPI_Bcast(&nmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#else
      cout << "Compile with MPI" << endl;
      abort();
#endif
    }
  else
    {
      NbPointsPerProc(0) = N;
      OffsetPointsProc(0) = 0;
    }

  // main loop over points
  Vector<double> Cfl(nmax), Err(nmax); Cfl.Zero(); Err.Zero();
  VectReal_wp param(seq.GetNbParam()); param.Zero();
  int old_percent = 0, percent = 0;
  for (int k = 0; k < N; k++)
    {
      if ((k >= offset_proc) && (k < offset_proc+nb_pts_proc))
        {
          if (rank_proc == 0)
            {
              old_percent = percent;
              percent = toInteger(round(100*(k-offset_proc) / Real_wp(nb_pts_proc)));
              if (percent != old_percent)
                {
                  cout << percent << "%  "; cout.flush();
                }
            }

          seq.SetParameters(permut(k), param);          
          
	  Real_wp cfl;
          Real_wp err = GetError_Sdirk(param, s, cfl);
	  Cfl(k-offset_proc) = toDouble(cfl);
          Err(k-offset_proc) = toDouble(err);
	  //just to write an example of params, cfl, and poly obtained
          if (k == 67498)
            {
              //DISP(k); DISP(param); DISP(cfl);
              //GetCFL_Sdirk(param, s, true);
            }
        }
    }

  if (rank_proc == 0)
    cout << endl;
  
  Vector<double> AllCfl, AllErr;
  if (rank_proc == 0)
    {
      AllCfl.Reallocate(nmax*nb_proc);
      AllErr.Reallocate(nmax*nb_proc);
    }

#ifdef SELDON_WITH_MPI
  MpiGather(MPI_COMM_WORLD, Cfl, xtmp, AllCfl, Cfl.GetM(), 0);
  MpiGather(MPI_COMM_WORLD, Err, xtmp, AllErr, Cfl.GetM(), 0);
#else
  AllCfl = Cfl;
  AllErr = Err;
#endif
  
  if (rank_proc == 0)
    {
      // we reorder with the good ordering (zeros are removed)
      CflP.Reallocate(N); CflP.Zero();
      ErrP.Reallocate(N); ErrP.Fill(1e300);
      Real_wp err_min = 1e300; int num_min = -1;
      for (int p = 0; p < nb_proc; p++)
        for (int k = 0; k < NbPointsPerProc(p); k++)
          {
            int num = permut(OffsetPointsProc(p) + k);
            CflP(num) = AllCfl(nmax*p+k);
            ErrP(num) = AllErr(nmax*p+k);
            if ( (CflP(num) > 1e100) && (ErrP(num) < err_min))
              {
                num_min = num;
                err_min = ErrP(num);
              }
          }            
      
      if (num_min != -1)
        {
          seq.SetParameters(num_min, param);
          for (int i = 0; i < param.GetM(); i++)
            param(i) /= coef_scale(i);
          
          cout << "Minimal errror " << ErrP(num_min) << " reached for param = " << param << endl;
        }
    }
  
#ifdef SELDON_WITH_MPI  
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

// cas avec 4 parametres
class SequenceGrid4D : public SequenceParameter
{
protected:
  int nx, ny, nz, nt;
  VectReal_wp param_x, param_y, param_z, param_t;
  
public :
  SequenceGrid4D(const VectReal_wp& x, const VectReal_wp& y, const VectReal_wp& z, const VectReal_wp& t)
    : nx(x.GetM()), ny(y.GetM()), nz(z.GetM()), nt(t.GetM()), param_x(x), param_y(y), param_z(z), param_t(t)
  {
  }

  int GetNbParam() const { return 4; }
  
  int GetNbPoints() const { return nx*ny*nz*nt; }

  void SetParameters(int num, VectReal_wp& x)
  {
    int nzt = nz*nt, nyzt = ny*nzt;
    int i = num/nyzt, j = (num-i*nyzt)/nzt, k = (num-i*nyzt-j*nzt)/nt, l = num%nt;
    x(0) = param_x(i); x(1) = param_y(j); x(2) = param_z(k); x(3) = param_t(l);
  }

};

// cas avec 3 parametres
class SequenceGrid3D : public SequenceParameter
{
protected:
  int nx, ny, nz;
  VectReal_wp param_x, param_y, param_z;
  
public :
  SequenceGrid3D(const VectReal_wp& x, const VectReal_wp& y, const VectReal_wp& z)
    : nx(x.GetM()), ny(y.GetM()), nz(z.GetM()), param_x(x), param_y(y), param_z(z)
  {
  }
  
  int GetNbParam() const { return 3; }
  
  int GetNbPoints() const { return nx*ny*nz; }

  void SetParameters(int num, VectReal_wp& x)
  {
    int nyz = ny*nz;
    int i = num/nyz, j = (num-i*nyz)/nz, k = num%nz;
    x(0) = param_x(i); x(1) = param_y(j); x(2) = param_z(k);
  }

};


// cas avec 2 parametres
class SequenceGrid2D : public SequenceParameter
{
protected:
  int nx, ny;
  VectReal_wp param_x, param_y;
  
public :
  SequenceGrid2D(const VectReal_wp& x, const VectReal_wp& y)
    : nx(x.GetM()), ny(y.GetM()), param_x(x), param_y(y)
  {
  }

  int GetNbParam() const { return 2; }
  
  int GetNbPoints() const { return nx*ny; }

  void SetParameters(int k, VectReal_wp& x)
  {
    int i = k/ny, j = k%ny;
    x(0) = param_x(i); x(1) = param_y(j);
  }

};


// cas avec 1 parametres
class SequenceGrid1D : public SequenceParameter
{
protected:
  int nx;
  VectReal_wp param_x;
  
public :
  SequenceGrid1D(const VectReal_wp& x)
    : nx(x.GetM()), param_x(x)
  {
  }

  int GetNbParam() const { return 1; }
  
  int GetNbPoints() const { return nx; }

  void SetParameters(int i, VectReal_wp& x)
  {
    x(0) = param_x(i);
  }
  
};

void CheckScheme(const VectReal_wp& param, int s)
{
  UnivariatePolynomial<Real_wp> NumerR;
  UnivariatePolynomial<Real_wp> DenomR, Pol;
  GetPolynomial_Sdirk(param, s, NumerR, DenomR, Pol);
  DISP(param); DISP(NumerR); DISP(DenomR); DISP(Pol);

  VectReal_wp z;
  ofstream file_out("exp_app.dat"); file_out.precision(15);
  Linspace(Real_wp(0), Real_wp(100), 100000, z);
  Real_wp err1, z1, err2, z2, Gmax(0);
  for (int i = 0; i < z.GetM(); i++)
    {
      Complex_wp zi = Iwp*z(i);
      Complex_wp val_ref = exp(zi);
      
      // Horner algorithm
      int order = NumerR.GetOrder();
      Complex_wp num = NumerR(order);
      Complex_wp denom = DenomR(order);
      for (int j = order-1; j >= 0; j--)
        {
          num = num*zi + NumerR(j);
          denom = denom*zi + DenomR(j);
        }
      
      Gmax = max(Gmax, abs(num / denom));
      file_out << z(i) << " " <<  realpart(num / denom) << " " << imagpart(num/denom) << '\n';

      if (i == 1)
        {
          cout << "Erreur pour " << zi << " = " << abs(num/denom - val_ref) << endl;
          err1 = abs(num/denom - val_ref); z1 = z(i);
        }

      if (i == 2)
        {
          cout << "Erreur pour " << zi << " = " << abs(num/denom - val_ref) << endl;
          err2 = abs(num/denom - val_ref); z2 = z(i);
          Real_wp rnum = abs((log(err1) - log(err2)) / (log(z2) - log(z1)));
          cout << "Measured order = " << rnum-1 << endl;
          cout << "Constante d'erreur = " << err2 / pow(z2, toInteger(round(rnum))) << endl;
        }
    }

  cout << "Maximal value of G = " << Gmax << endl;
  Real_wp cfl = GetCFL(NumerR, DenomR, true);
  cout << "CFL = " << cfl << endl;

  file_out.close();
}

// computes coefficient of the source
void GetCoefficientSource(const UnivariatePolynomial<Real_wp>& P,
                          const UnivariatePolynomial<Real_wp>& Q,
                          int order, VectReal_wp& ci, Vector<VectReal_wp>& wi)
{
  DISP(order);
  VectReal_wp b; Matrix<Real_wp> A;
  A.Reallocate(order+1, order+1);
  b.Reallocate(order+1);
  b.Zero();
  A.Zero();
  
  // calcul de 1/k!
  VectReal_wp invFacto(order+1);
  invFacto(0) = Real_wp(1);
  for (int k = 1; k <= order; k++)
    invFacto(k) = invFacto(k-1) / Real_wp(k);
  
  DISP(invFacto);
  // contribution du numerateur :  -P X^n
  for (int i = 0; i <= P.GetOrder(); i++)
    if (i <= order)
      b(i) = -P(i);
  
  // contribution du denominateur : Q X^{n+1}
  for (int i = 0; i <= Q.GetOrder(); i++)
    {
      for (int k = 0; k <= order; k++)
        if ((i+k) <= order)
          {
            b(i+k) += Q(i)*invFacto(k);
            for (int j = 1; j <= k; j++)
              A(i+k, j) += Q(i)*invFacto(k);
          }
    }
  
  DISP(A);
  DISP(b);

  // calcul des ci : on prend les points de Gauss-Lobatto
  VectReal_wp weights;
  ComputeGaussLobatto(ci, weights, order-1);
  
  DISP(ci);

  // Matrice de Vandermonde pour les ci
  // les factorielles sont mises dans le second membre
  Matrix<Real_wp> VDM(ci.GetM(), ci.GetM());  
  for (int i = 0; i < ci.GetM(); i++)
    for (int j = 0; j < ci.GetM(); j++)
      VDM(i, j) = pow(ci(j), i);
  
  Vector<int> pivot; VectReal_wp rhs(ci.GetM());
  GetLU(VDM, pivot);

  wi.Reallocate(min(P.GetOrder(), order));
  for (int powA = 0; powA <= order-1; powA++)
    {
      rhs.Zero();
      for (int j = 1; j <= order-powA; j++)
        rhs(j-1) = A(powA+j, j) / invFacto(j-1);

      // forme alternative
      VectReal_wp rhs_vec(order-powA);
      rhs_vec.Zero();
      for (int j = 1; j <= order-powA; j++)
	{
	  int r = powA + 1;
	  for (int i = 0; i <= powA; i++)
	    rhs_vec(j-1) += Q(i) * invFacto(r + j - i - 1);

	  rhs_vec(j-1) /= invFacto(j-1);
	}
      
      DISP(powA); DISP(rhs); DISP(rhs_vec);
      SolveLU(VDM, pivot, rhs);
      if (powA <= wi.GetM()-1)
        wi(powA) = rhs;
      
      cout << "Pour powA = " << powA << endl;
      cout << "wi = " << rhs << endl;
      if (powA <= wi.GetM()-1)
	cout << "wi(powA) = " << wi(powA) << endl;
      cout << endl << endl;
    }
}
  

class LinearSdirkScheme : public VirtualTimeScheme<Real_wp>
{
public:  
  Real_wp gamma, dt;
  UnivariatePolynomial<Real_wp> Numer;
  Vector<VectReal_wp> wi; VectReal_wp ci;

  VectReal_wp Un, Fn, KhUn;
  Vector<VectReal_wp> EvalF;
  
  void Init(const UnivariatePolynomial<Real_wp>& P, const UnivariatePolynomial<Real_wp>& Q,
            const Real_wp& gam_, int order)
  {
    gamma = gam_;    
    Numer = P;
    
    GetCoefficientSource(P, Q, order+1, ci, wi);
  }
  
  void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                           VectReal_wp& u0, VirtualOdeSystem<Real_wp>& sys)
  {
    Un = u0;
    Fn = Un; KhUn = Un;
    dt = deltat;
    sys.FactorizeOperatorDhPlusGammaKh(Real_wp(1), gamma*deltat, gamma*deltat);

    EvalF.Reallocate(ci.GetM());
    for (int i = 0; i < ci.GetM(); i++)
      EvalF(i) = Fn;
  }

  void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<Real_wp>& sys)
  {
    for(int i = 0; i < ci.GetM(); i++)
      {
        Real_wp tcurrent = t + ci(i)*dt;
        // compute y = y + alpha f(t)	
        EvalF(i).Zero();
        sys.AddPrimitiveTimeSource(Real_wp(1), tcurrent, 0, EvalF(i));
        //sys.Assemble(EvalF(i));
      }
    
    //DISP(Un);
    int last_coef = Numer.GetOrder();
    Fn = Numer(last_coef)*Un;
    for (int k = last_coef-1; k >= 0; k--)
      {
        sys.EvaluateFunction(t, Fn, KhUn, false, false);
        Fn = dt*KhUn + Numer(k)*Un;
        if (k < wi.GetM())
	  for(int i = 0; i < wi(k).GetM(); i++)
	    Fn += dt*wi(k)(i)*EvalF(i);
      }
    
    for (int k = 0; k < Numer.GetOrder(); k++)
      {
        sys.SolveOperatorDhPlusGammaKh(t, Fn, KhUn);
        Fn = KhUn;
      }
    
    Un = Fn;
    //DISP(Un);
    //int test_input; cout << " waiting " <<endl; cin >> test_input;
  }

  void Clear(){}
  void ClearFirst(Vector<Real_wp>& Ufinal) { Ufinal = Un;}

  void ChangeTimeStep(const Real_wp& deltat) {}
    
  Vector<Real_wp>& GetIterate() { return Un; }
  const Vector<Real_wp>& GetIterate() const { return Un; }
    
  int GetNumberOfIterates() const { return 1;}
  Vector<Real_wp>& GetIterate(int k) { return Un; }
  
  
};

//! function f(t, y) = y + exp(t)/(1+t^2)
class ScalarFunctionFirst : public VirtualOdeSystem<Real_wp>
{
public :
  Real_wp gamma;
  Real_wp coef_mass;
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Real_wp one(1);
    f(0) = y(0);
    if (source)
      f(0) += exp(t) / (Real_wp(1) + t*t);
  }
  
  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& F)
  {
    F(0) += alpha*exp(t) / (Real_wp(1) + t*t);
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
    // y (coef_mass - gamma) = f
    Real_wp one(1);
    y(0) = f(0) / (coef_mass - gamma);
  }
  
};

template<class TimeScheme>
void TestFirstOrderScheme(TimeScheme& scheme, int r, Real_wp dt)
{
  // first, we solve dy/dt = y
  Real_wp Tmax(3), zero(0);  
  VectReal_wp y0(1);
  y0(0) = Real_wp(1);
  ScalarFunctionFirst fct;
  RunTimeScheme(zero, Tmax, dt, y0, fct, scheme);
  
  Real_wp y_dt = y0(0);
  Real_wp yref = (Real_wp(1) + atan(Tmax))*exp(Tmax);
  DISP(yref); DISP(y_dt);
  
  y0(0) = Real_wp(1);
  RunTimeScheme(zero, Tmax, dt/2, y0, fct, scheme);
  
  Real_wp y_dt_half = y0(0);
  
  Real_wp err = abs(y_dt - yref)/abs(yref);
  Real_wp err_half = abs(y_dt_half - yref)/abs(yref);
  Real_wp order = (log(err) - log(err_half))/log(2.0);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  //DISP(err); DISP(err_half); DISP(r); DISP(order);
  Real_wp err_order(0.03);
  if (order < r*(1.0-err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }
}


// case where no additional stages are added
// Gamma solves an equation
void FindSdirkMinimal(int order)
{
  UnivariatePolynomial<Real_wp> PolGamma, Numer, Denom;  
  PolGamma.SetOrder(order);
  Numer.SetOrder(order+1);
  Denom.SetOrder(order);

  // on genere les coefs C_n^k directement
  IVect row;
  GetCnk(order, row);

  Real_wp invFacto(1);
  Numer(0) = Real_wp(1);
  for (int i = 0; i <= order; i++)
    {      
      // Numer : developpement limite de l'exponentielle
      Numer(i+1) = abs(invFacto);
      
      // Denom : (1- gamma z)^order
      Denom(i) = row(i);
      if (i%2 == 1)
        Denom(i) = -row(i);

      // on remplit l'equation que doit satisfaire gamma
      PolGamma(order-i) = row(i) * invFacto;
      invFacto *= -Real_wp(1)/(i+2);
    }
  
  DISP(Numer); DISP(Denom); DISP(PolGamma);
  VectReal_wp R, Rimag;
  SolvePolynomialEquation(PolGamma, R, Rimag);
  
  //DISP(R); DISP(Rimag);
  
  // boucle sur les racines
  Real_wp gamma_only(0);
  for (int p = 0; p < R.GetM(); p++)
    {
      Real_wp gam = R(p);
      
      // vrai denominateur avec gamma
      UnivariatePolynomial<Real_wp> DenomR(Denom);
      Real_wp powG (1);
      for (int i = 0; i <= order; i++)
        {
          DenomR(i) *= powG;
          powG *= gam;
        }
      
      // vrai numerateur en enlevant les zeros
      // dus aux constantes alpha, beta, etc
      UnivariatePolynomial<Real_wp> NumerR;
      NumerR = Numer*DenomR;
      
      Real_wp asymptote = NumerR(order) / DenomR(order);
      if (abs(asymptote) <= Real_wp(1))
        {
          cout << "Gamma = " << gam << endl;
          cout << "Asymptote = " <<  asymptote << endl;
          cout << "Constantes d'erreur = " << endl;
	  Real_wp invFacto(1);
	  for (int i = 0; i <= order; i++)
	    invFacto *= Real_wp(1)/(i+1);
	  
          for (int i = order+1; i < NumerR.GetM(); i++)
            {
	      if (i == (order+1))
		cout << NumerR(i) << "  ";
	      else
		cout << NumerR(i)+invFacto << "  ";

	      invFacto *= Real_wp(1)/(i+1);
	    }
	  
          cout << endl;
          
          DISP(NumerR);
          NumerR.ResizeOrder(order);
          cout << "Numerateur de R(z) = " << NumerR << endl;
	  for (int j = order; j >= 0; j--)
	    cout << NumerR(j) << ", ";

	  cout << endl;
          cout << "Denominateur de R(z) = " << DenomR << endl;
          
	  LinearSdirkScheme scheme;
          scheme.Init(NumerR, DenomR, gam, order);
          
          TestFirstOrderScheme(scheme, order+1, Real_wp(1)/100);
          
          Real_wp cfl = GetCFL(NumerR, DenomR);
          DISP(cfl);
          if (cfl > 1e100)
            gamma_only = gam;
          
	  VectReal_wp Rn, Rn_imag;
	  SolvePolynomialEquation(NumerR, Rn, Rn_imag);
	  VectComplex_wp Rcplx(Rn.GetM());
	  for (int i = 0; i < Rn.GetM(); i++)
            Rcplx(i) = Complex_wp(Rn(i), Rn_imag(i));
          
          cout << "Racines du numerateur " << Rcplx << endl;
	  int test_input; cout << "waiting for the next root..." << endl; cin >> test_input;
	}
    }

  cout << endl << endl;
  cout.precision(32);
  if (gamma_only != Real_wp(0))
    cout << "Only acceptable value of gamma = " << gamma_only << endl;
  else
    cout << "No acceptable value of gamma = " << endl;
}

void FillParameter(const R2& pt, const Real_wp& gam, const VectReal_wp& coef_scale, VectReal_wp& param)
{
  param(0) = gam; param(1) = pt(0)*coef_scale(1); param(2) = pt(1)*coef_scale(2); 
}

bool ProjectPointStabilityRegion(int s, int n, const Real_wp& gam, R2& pt,
				 const VectReal_wp& coef_scale, const R2& normale)
{
  VectReal_wp param(n);
  FillParameter(pt, gam, coef_scale, param);
  
  Real_wp cfl = GetCFL_Sdirk(param, s);
  if (cfl < 1e100)
    {
      pt -= Real_wp(1e-12)*normale;
      FillParameter(pt, gam, coef_scale, param);
      cfl = GetCFL_Sdirk(param, s);
    }
  
  if (cfl < 1e100)
    {
      cout << "Impossible " << endl;
      cout << "the point " << pt << " is not in the stability domain" << endl;
      abort();
    }

  Real_wp coef(1); bool test_loop = true;
  R2 ptA = pt, ptB = pt + normale;
  while (test_loop)
    {
      ptB = pt + coef*normale;
      FillParameter(ptB, gam, coef_scale, param);  
      cfl = GetCFL_Sdirk(param, s);
      if (cfl > 1e100)
	coef *= Real_wp(2);
      else
	test_loop = false;
      
      if (coef > Real_wp(16))
	{
	  cout << "Point not on the boundary ? " << pt << endl;
	  DISP(pt); DISP(normale); DISP(coef);
	  return false;
	}
    }
  
  // bisection method between ptA and ptB
  R2 ptC;
  while (ptA.Distance(ptB) > 1e-12)
    {
      ptC = Real_wp(0.5)*(ptA + ptB);
      FillParameter(ptC, gam, coef_scale, param);  
      cfl = GetCFL_Sdirk(param, s);
      if (cfl < 1e100)
	ptB = ptC;
      else
	ptA = ptC;
    }

  pt = ptA;
  return true;
}

void RefineBoundary(const Mesh<Dimension2>& mesh, int n, int s, const Real_wp& gam,
                    const VectReal_wp& coef_scale, Mesh<Dimension2>& mesh_dest)
{
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int nb_proc(1), rank_proc(0);
#endif

  // edges are distributed to the different processors
  int N = mesh.GetNbBoundaryRef(), nb_pts_proc, offset_proc;
  GetParallelDistributionPoints(nb_proc, rank_proc, N, nb_pts_proc, offset_proc);

  // main loop over points
  Vector<bool> VertexTreated(mesh.GetNbVertices());
  VertexTreated.Fill(false);
  for (int num = 0; num < 2; num++)
    {
      for (int i = offset_proc; i < offset_proc+nb_pts_proc; i++)
        {
          int n0 = mesh.BoundaryRef(i).numVertex(0);
          int n1 = mesh.BoundaryRef(i).numVertex(1);
          R2 pt0 = mesh.Vertex(n0), pt1 = mesh.Vertex(n1);
          R2 normale(pt1(1)-pt0(1), pt0(0)-pt1(0));
          R2 milieu = Real_wp(0.5)*(pt0 + pt1);
          
          int num_elem = mesh.BoundaryRef(i).numElement(0);
          R2 center;
          int nb_vert = mesh.Element(num_elem).GetNbVertices();
          for (int j = 0; j < nb_vert; j++)
            center += mesh.Vertex(mesh.Element(num_elem).numVertex(j));
          
          Mlt(Real_wp(1)/nb_vert, center);
          R2 diff = milieu - center;
          if (DotProd(diff, normale) < 0)
            normale = -normale;
          
          if ((num == 1) || (nb_vert == 4))
            {
              if (!VertexTreated(n0))
                {
                  VertexTreated(n0) = true;
                  bool pt_found = ProjectPointStabilityRegion(s, n, gam, pt0, coef_scale, normale);
                  if (!pt_found)
                    VertexTreated(n0) = false;
                  
                  mesh_dest.Vertex(n0) = pt0;                  
                }
              
              if (!VertexTreated(n1))
                {
                  VertexTreated(n1) = true;
                  bool pt_found = ProjectPointStabilityRegion(s, n, gam, pt1, coef_scale, normale);
                  if (!pt_found)
                    VertexTreated(n1) = false;
                  
                  mesh_dest.Vertex(n1) = pt1;
                }
            }
        }

#ifdef SELDON_WITH_MPI
      if (nb_proc > 1)
        {
          MPI_Status status; Vector<int64_t> xtmp;
          if (rank_proc == 0)
            {
              for (int i = 1; i < nb_proc; i++)
                {
                  int nb_pts_treated = 0; Vector<int> list_points;
                  VectReal_wp new_points;
                  MPI_Recv(&nb_pts_treated, 1, MPI_INTEGER, i, 23, MPI_COMM_WORLD, &status);
                  if (nb_pts_treated > 0)
                    {
                      list_points.Reallocate(nb_pts_treated);
                      new_points.Reallocate(2*nb_pts_treated);
                      MPI_Recv(list_points.GetData(), list_points.GetM(), MPI_INTEGER, i, 24, MPI_COMM_WORLD, &status);
                      MpiRecv(MPI_COMM_WORLD, new_points, xtmp, 2*nb_pts_treated, i, 25, status);
                    }
                  
                  for (int i = 0; i < nb_pts_treated; i++)
                    {
                      int nv = list_points(i);
                      if (!VertexTreated(nv))
                        {
                          VertexTreated(nv) = true;
                          mesh_dest.Vertex(nv)(0) = new_points(2*i);
                          mesh_dest.Vertex(nv)(1) = new_points(2*i+1);
                        }
                    }
                }          
            }
          else
            {
              int nb_pts_treated = 0;
              for (int i = 0; i < mesh.GetNbVertices(); i++)
                if (VertexTreated(i))
                  nb_pts_treated++;
              
              Vector<int> list_points(nb_pts_treated);
              VectReal_wp new_points(2*nb_pts_treated);
              nb_pts_treated = 0;
              for (int i = 0; i < mesh.GetNbVertices(); i++)
                if (VertexTreated(i))
                  {
                    list_points(nb_pts_treated) = i;
                    new_points(2*nb_pts_treated) = mesh_dest.Vertex(i)(0);
                    new_points(2*nb_pts_treated+1) = mesh_dest.Vertex(i)(1);
                    nb_pts_treated++;
                  }
              
              MPI_Ssend(&nb_pts_treated, 1, MPI_INTEGER, 0, 23, MPI_COMM_WORLD);
              if (nb_pts_treated > 0)
                {
                  MPI_Ssend(list_points.GetData(), nb_pts_treated, MPI_INTEGER, 0, 24, MPI_COMM_WORLD);
                  MpiSsend(MPI_COMM_WORLD, new_points, xtmp, 2*nb_pts_treated, 0, 25);
                }
            }
          
          if (num == 0)
            MpiBcast(MPI_COMM_WORLD, VertexTreated, xtmp, VertexTreated.GetM(), 0);
        }
#endif

    }

  if (rank_proc == 0)
    {
      mesh_dest.ForceCoherenceMesh(true);
      mesh_dest.Write("contour.mesh");
    }
        
  //int test_input; cout << "waiting" << endl; cin >> test_input;
}
		    
int main(int argc, char** argv)
{
  if (argc < 3)
    {
      cout << "Provide the paramter s, the number of additional stages" << endl;
      cout << "Usage : ./test.x 10 n" << endl;
      abort();
    }
  
  InitMontjoie(argc, argv);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  int s = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  if (n == 0)
    {
      // cas sans etape additionnelle
      FindSdirkMinimal(s);
      return FinalizeMontjoie();
    }

  // pour les constantes alpha, beta, on divise par 1/k! la valeur max du parametre
  Real_wp invFacto(1);
  for (int i = 0; i <= s; i++)
    invFacto *= Real_wp(1)/(i+2);
  
  VectReal_wp coef_scale(n);
  coef_scale.Fill(1.0);
  for (int i = 1; i < n; i++)
    {
      coef_scale(i) = invFacto;
      invFacto *= Real_wp(1)/(s+2+i);
    }
  
  if (argc == 5)
    {
      if (rank_proc == 0)
        cout << "on raffine le contour du domaine de stabilite" << endl;
      
      string name_mesh(argv[3]);
      Mesh<Dimension2> mesh, mesh_dest;
      mesh.Read(name_mesh);

      mesh_dest.Read(name_mesh);
      mesh_dest.ClearElements();
      
      if (rank_proc == 0)
        cout << "Mesh read successfully" << endl;
      
      Real_wp gam = atof(argv[4]);
      
      RefineBoundary(mesh, n, s, gam, coef_scale, mesh_dest);

      return FinalizeMontjoie();
    }

  // cas ou on teste directement le schema SDIRK trouve
  /* LinearSdirkScheme_Iterator<Real_wp> sdirk;
     sdirk.SetOrder(s+1, n);
  
     LinearSdirkScheme scheme;
     scheme.Init(sdirk.GetNumeratorStabilityFunction(),
     sdirk.GetDenominatorStabilityFunction(),
     sdirk.GetGammaCoefficient(), s);
  
     // TestFirstOrderScheme(scheme, s+1, Real_wp(1)/100);
     DISP(GetCFL(sdirk.GetNumeratorStabilityFunction(),
     sdirk.GetDenominatorStabilityFunction()));
  
     VectReal_wp coef;
     ComputeErrorSource(sdirk.GetNumeratorStabilityFunction(),
     sdirk.GetDenominatorStabilityFunction(), s, 3, coef);
  
     DISP(coef);
  
     return FinalizeMontjoie(); */

  
  // cas ou on teste un jeu de parametres
  /* VectReal_wp param(3);  

  // choix optimaux s=5, l=2
  // ./test.x 5 2 0.204071 0.204072 0.999907304 0.999907307
  // param(0) = 0.204071; param(1) = 1.9839430662e-4;
  
  // choix optimaux s=7, l=2
  // ./test.x 7 2 0.16688964 0.16689033 1.061759726 1.061759730
  // param(0) = 0.16689; param(1) = 2.9259251764e-6;
  
  // choix optimaux s=9, l=2
  // ./test.x 9 2 0.141938 0.141944 0.91738 0.91742
  // param(0) = 0.141940; param(1) = 2.2982637210e-8;
  
  // choix optimaux s=7, l=3
  //param(0) = 0.136339; param(1) = 2.766997072e-6; param(2) = -3.464371455249e-6;

  // choix optimaux s=9, l=3
  //param(0) = 0.151706; param(1) = Real_wp(0.9816111)/39916800; param(2) =  Real_wp(-20.664613068541)/479001600;
  
  // choix optimaux s=11, l=3
  param(0) = 0.132571; param(1) = 1.64451085e-10; param(2) = -2.89888832092e-10;
  
  CheckScheme(param, s); 
  
  return FinalizeMontjoie();
  
  */
  
  // cas ou on balaie pour trouver les meilleurs parametres
  if (argc != 7+2*n)
    {
      cout << "Provide the range of all parameters" << endl;
      cout << "Usage : ./test.x 10 n gamma alpha ..." << endl;
      abort();
    }


  VectReal_wp param_max(n), param_min(n);
  for (int p = 0; p < n; p++)
    {
      param_min(p) = to_num<Real_wp>(string(argv[3+2*p]));
      param_max(p) = to_num<Real_wp>(string(argv[3+2*p+1]));
    }

  int off_param = 3+2*n;
  int N = atoi(argv[off_param]), Nx = atoi(argv[off_param+1]), Ny = atoi(argv[off_param+2]), Nz = atoi(argv[off_param+3]);
  if (rank_proc == 0)
    { DISP(N); DISP(Nx); DISP(Ny); DISP(Nz); }

  invFacto = Real_wp(1);
  for (int i = 0; i <= s; i++)
    invFacto *= Real_wp(1)/(i+2);

  for (int i = 1; i < n; i++)
    {
      if (rank_proc == 0)
	{ DISP(i); DISP(Real_wp(1)/invFacto); }

      param_min(i) *= invFacto;
      param_max(i) *= invFacto;
      invFacto *= Real_wp(1)/(s+2+i);
    }

  if (rank_proc == 0)
    {
      cout << "Parametres max = " << param_max << endl;
      cout << "Parametres min = " << param_min << endl;
    }

  
  SdirkErrorMinimization<Real_wp> fct;
  fct.SetOrder(s, n);

  Vector<double> CflP, ErrP;
  if (n == 1)
    {
      VectReal_wp Gamma, xmin; 
      Linspace(param_min(0), param_max(0), N, Gamma);
      
      SequenceGrid1D grid(Gamma);
      //ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      ComputeCFLOptimal(fct, grid, 4, s, coef_scale, xmin);
      //if (rank_proc == 0)
      //CheckScheme(xmin, s);
    }
  else if (n == 2)
    {      
      // calcul de la CFL pour une sequence de parametres
      VectReal_wp Gamma, Alpha, xmin; 
      Linspace(param_min(0), param_max(0), N, Gamma);
      Linspace(param_min(1), param_max(1), Nx, Alpha);
      
      SequenceGrid2D grid(Gamma, Alpha);
      //ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      ComputeCFLOptimal(fct, grid, 4, s, coef_scale, xmin);
      //if (rank_proc == 0)
      //CheckScheme(xmin, s);
    }
  else if (n == 3)
    {
      // calcul de la CFL pour une sequence de parametres
      VectReal_wp Gamma, Alpha, Beta, xmin; 
      Linspace(param_min(0), param_max(0), N, Gamma);
      Linspace(param_min(1), param_max(1), Nx, Alpha);
      Linspace(param_min(2), param_max(2), Ny, Beta);
      
      SequenceGrid3D grid(Gamma, Alpha, Beta);
      ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      //ComputeCFLOptimal(fct,grid, 4, s, coef_scale, xmin);
      //if (rank_proc == 0)
      //CheckScheme(xmin, s);
    }
  else if (n == 4)
    {
      // calcul de la CFL pour une sequence de parametres
      VectReal_wp Gamma, Alpha, Beta, Eta, xmin; 
      Linspace(param_min(0), param_max(0), N, Gamma);
      Linspace(param_min(1), param_max(1), Nx, Alpha);
      Linspace(param_min(2), param_max(2), Ny, Beta);
      Linspace(param_min(3), param_max(3), Nz, Eta);
      
      SequenceGrid4D grid(Gamma, Alpha, Beta, Eta);
      //ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      //ComputeCFLOptimal(fct,grid, 4, s, coef_scale, xmin);
      if (rank_proc == 0)
        CheckScheme(xmin, s);
    }
  else
    {
      cout << "Case not treated" << endl;
      abort();
    }
  
  if (rank_proc == 0)
    {
      DISP(CflP.GetNormInf());
      CflP.WriteText("CflR.dat");
      ErrP.WriteText("ErrR.dat");
      
      // pour un choix de parametre
      // on verifie le schema
      //VectReal_wp param(2);
      //param(0) = 0.3525; param(1) = 0.010625;
      //CheckScheme(param, s);
    }
  
  return FinalizeMontjoie();
}
