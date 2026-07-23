#include "Solver/MontjoieSolver.hxx"
#include "Output/MontjoieOutput.hxx"
#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

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

Real_wp GetGmax(const UnivariatePolynomial<Real_wp>& P,
                const UnivariatePolynomial<Real_wp>& Q)
{
  if (P.GetOrder() != Q.GetOrder())
    {
      cout << "Case not handled" << endl;
      abort();
    }

  if (abs(P(P.GetOrder())) > abs(Q(Q.GetOrder())))
    return Real_wp(1e30);
  
  // on cherche le max de |P|/|Q|
  UnivariatePolynomial<Real_wp> Pmod, Qmod, dPmod, dQmod, diff;
  
  GetSquareModulus(P, Pmod);
  GetSquareModulus(Q, Qmod);
  
  DerivatePolynomial(Pmod, dPmod);
  DerivatePolynomial(Qmod, dQmod);

  // points critiques de  |P|^2/|Q|^2
  diff = Pmod*dQmod - dPmod*Qmod;

  for (int i = 0; i <= diff.GetOrder(); i++)
    if (abs(diff(i)) <= 1e4*epsilon_machine)
      diff(i) = Real_wp(0);

  VectReal_wp Li, Lr;
  SolvePolynomialEquation(diff, Lr, Li);

  Real_wp G = 1.0;
  for (int i = 0; i < Lr.GetM(); i++)
    if ((abs(Li(i)) <= 1e4*epsilon_machine) && (Lr(i) >= 0))
      G = max(G, sqrt(abs(Pmod.Evaluate(Lr(i)) / Qmod.Evaluate(Lr(i)))));

  return G;
}

// on calcule la fonction de stabilite R = Numer/Denom
// en fonction de s (la methode est d'ordre s+1) et des parametres
void GetPolynomial_Sdirk(const Real_wp& xsi, const Real_wp& Gamma,
                         const VectReal_wp& param, int order,
                         UnivariatePolynomial<Real_wp>& Numer,
                         UnivariatePolynomial<Real_wp>& Denom,
                         UnivariatePolynomial<Real_wp>& Pol)
{
  int p = param.GetM() + order;
  
  // Pol = 1 - gamma z
  Pol.SetOrder(1); Pol(0) = Real_wp(1); Pol(1) = -Gamma;
  
  // Numer : developpement limite de l'exponentielle
  Numer.SetOrder(order);
  Numer(0) = Real_wp(1);
  for (int i = 1; i <= order; i++)
    Numer(i) = Numer(i-1) / Real_wp(i) * xsi;

  //DISP(xsi); DISP(Numer); DISP(Pol); DISP(p);
  // Denom = (1-gamma z)^p
  Denom = Pow(Pol, p); //DISP(Gamma); DISP(p); DISP(Pol); DISP(Denom);
  
  // Pol = (1 - gamma z)^p exp(z) en enlevant les termes de degre superieur
  Pol = Numer*Denom;

  //DISP(Pol);
  // on ajoute les parametres additionels
  for (int i = order; i < p; i++)
    Pol(i+1) += param(i-order);
  
  Numer = Pol;
  Numer.ResizeOrder(p);
  //DISP(Numer);
}

R2 GetBounds_Sdirk(const Real_wp& xsi, const Real_wp& Gamma,
                   int order, int n)
{
  int p = n + order;
  
  // Pol = 1 - gamma z
  UnivariatePolynomial<Real_wp> Pol, Numer, Denom;
  Pol.SetOrder(1); Pol(0) = Real_wp(1); Pol(1) = -Gamma;
  
  // Numer : developpement limite de l'exponentielle
  Numer.SetOrder(order);
  Numer(0) = Real_wp(1);
  for (int i = 1; i <= order; i++)
    Numer(i) = Numer(i-1) / Real_wp(i) * xsi;

  //DISP(xsi); DISP(Numer); DISP(Pol); DISP(p);
  // Denom = (1-gamma z)^p
  Denom = Pow(Pol, p); //DISP(Gamma); DISP(p); DISP(Pol); DISP(Denom);
  
  // Pol = (1 - gamma z)^p exp(z) en enlevant les termes de degre superieur
  Pol = Numer*Denom;

  Real_wp coef_num = Pol(p);
  Real_wp coef_denom = Denom(p);
  DISP(coef_num); DISP(coef_denom);
  
  Real_wp xmin = -coef_denom - coef_num;
  Real_wp xmax = coef_denom - coef_num;

  return R2(xmin, xmax);
  //DISP(Numer);
}


// on calcule la CFL d'un schema de type SDIRK
// en fonction de s (la methode est d'ordre s+1) et des parametres
// param(0) = gamma
// param(1:) = les constantes alpha, beta, etc
Real_wp GetCFL_Sdirk(const Real_wp& xsi, const Real_wp& Gamma,
                     const VectReal_wp& param, int order, bool print_pol = false)
{
  UnivariatePolynomial<Real_wp> Denom, Numer, Pol;
  GetPolynomial_Sdirk(xsi, Gamma, param, order, Numer, Denom, Pol);
  if (print_pol)
    {
      DISP(Numer); DISP(Denom);
    }
  
  Real_wp cfl = GetCFL(Numer, Denom);
  return cfl;
}


Real_wp GetError_Sdirk(const Real_wp& xsi, const Real_wp& Gamma,
                       const VectReal_wp& param, int order, Real_wp& cfl)
{
  UnivariatePolynomial<Real_wp> NumerR, DenomR, Pol;
  UnivariatePolynomial<Real_wp> exp_taylor;
  
  GetPolynomial_Sdirk(xsi, Gamma, param, order, NumerR, DenomR, Pol);
  cfl = GetCFL(NumerR, DenomR);
  
  Real_wp err = GetGmax(NumerR, DenomR);
  if (err <= Real_wp(1)+1e4*epsilon_machine)
    {
        Real_wp invFacto = xsi;
        for(int j = 0; j < order; j++)
          invFacto *= xsi/(j+2);

        err = abs(invFacto - param(0));
    }
  
  return err;
}


// classe de base pour definir une sequence N-D de parametres
class SequenceParameter
{
public:
  Real_wp Gamma, xsi;
  
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

  Real_wp Gamma, xsi;
  
  //! default constructor
 SdirkErrorMinimization()
  {
    s = 3; Gamma = Real_wp(0); xsi = Real_wp(0);
    type_init = USER;
    this->n = 1;
    //this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX;
    this->type_algo = VirtualMinimizedFunction<T>::SIMPLEX2_RAND;
    //this->type_algo = VirtualMinimizedFunction<T>::CG;
    //this->type_algo = VirtualMinimizedFunction<T>::CG_FR;
    //this->type_algo = VirtualMinimizedFunction<T>::STEEPEST_DESCENT;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS;
    //this->type_algo = VirtualMinimizedFunction<T>::BFGS2;

    this->step_size = 0.1;
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
    for (int i = 0; i < this->n; i++)
      params(i) = double(rand())/RAND_MAX;
  }

  //! initializes the number of free parameters
  void SetOrder(int s, int n, const Real_wp& G, const Real_wp& x)
  {
    this->s = s;
    this->n = n; //number of free parameters
    this->Gamma = G;
    this->xsi = x;
    
  }

  //! computes feval = first_non_zero_coeff(e^z - R(z)) from free parameters contained in params
  void EvaluateFunction(const Vector<double>& param, double& feval)
  {
    int nbParam = param.GetM();
    Vector<T> xb(nbParam);
    for (int i = 0; i < nbParam; i++)
      xb(i) = param(i)*coef_scale(i);

    Real_wp cfl;
    Real_wp err = GetError_Sdirk(xsi, Gamma, xb, this->s, cfl);
    feval = toDouble(err);
    DISP(xb); DISP(feval);
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
	  Real_wp err = GetError_Sdirk(fct.xsi, fct.Gamma, xmin, s, cfl);
          //DISP(xmin); DISP(err); DISP(cfl);
	  // a quel rang se situe cette cfl ?
	  int rank = nb_selection;
          //if (cfl >= Real_wp(1e299))
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
              Cfl(rank) = toDouble(err);
	    }
	}
    }
  
  DISP(Cfl); DISP(NumParam); DISP(coef_scale);
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
        
        MinimizeParametersGsl(fct, xsol, 1e-10, 1000);
        
        fct.EvaluateFunction(xsol, feval);
        DISP(xsol); DISP(feval);
        //if(fmin<0)
        //fmin = feval;
        
        DISP(fmin);
        DISP(xsol);
        
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
  
  Real_wp cfl_sdirk = GetError_Sdirk(fct.xsi, fct.Gamma, xmin, s, cfl);
  cfl_opt = toDouble(cfl_sdirk);
  
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
          cout << "Error = " << AllCfl(k) << "\t" ;// endl;
          cout << "For parameters = " << all_optParam(k) << "\t";// endl;
	  double feval;
	  fct.EvaluateFunction(all_optParam(k), feval);
	  cout << "error = " << feval << endl;
        }

      int k0 = 0;
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
          Real_wp err = GetError_Sdirk(seq.xsi, seq.Gamma, param, s, cfl);
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


int main(int argc, char** argv)
{
  if (argc < 4)
    {
      cout << "Provide the parameter s, the number of additional stages and the quadrature point" << endl;
      cout << "Usage : ./test.x s n num_point" << endl;
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
  int num_point = atoi(argv[3]);
  
  LinearSdirkScheme_Iterator<Real_wp> scheme;
  scheme.SetOrder(s+1, n, LinearSdirkScheme_Iterator<Real_wp>::STABLE_ALTERNATIVE);

  Real_wp xsi = (1.0-scheme.GetQuadraturePoint(num_point));
  DISP(xsi);
  
  /*VectReal_wp param(n); param.Zero();
  Real_wp cfl;
  Real_wp err = GetError_Sdirk(Real_wp(0.5), scheme.GetGammaCoefficient(), param, s, cfl);
  DISP(err); DISP(cfl);*/
    
  // cas ou on balaie pour trouver les meilleurs parametres
  if (argc != 8+2*n)
    {
      cout << "Provide the range of all parameters" << endl;
      cout << "Usage : ./test.x s n num_point alpha ..." << endl;
      abort();
    }


  VectReal_wp param_max(n), param_min(n), coef_scale(n);
  coef_scale.Fill(Real_wp(1));
  for (int p = 0; p < n; p++)
    {
      param_min(p) = to_num<Real_wp>(string(argv[4+2*p]));
      param_max(p) = to_num<Real_wp>(string(argv[4+2*p+1]));
    }
  
  int off_param = 4+2*n;
  int N = atoi(argv[off_param]), Nx = atoi(argv[off_param+1]), Ny = atoi(argv[off_param+2]), Nz = atoi(argv[off_param+3]);
  if (rank_proc == 0)
    { DISP(N); DISP(Nx); DISP(Ny); DISP(Nz); }
  
  Real_wp invFacto = xsi;
  for (int i = 0; i < s; i++)
    invFacto *= xsi/(i+2);
  
  DISP(invFacto);
  for (int i = 0; i < n; i++)
    {
      if (rank_proc == 0)
	{ DISP(i); DISP(Real_wp(1)/invFacto); }

      coef_scale(i) = invFacto;
      param_min(i) *= invFacto;
      param_max(i) *= invFacto;
      invFacto *= xsi/(s+2+i);
    }

  if (rank_proc == 0)
    {
      cout << "Parametres max = " << param_max << endl;
      cout << "Parametres min = " << param_min << endl;
    }
  
  
  SdirkErrorMinimization<Real_wp> fct;
  fct.SetOrder(s, n, scheme.GetGammaCoefficient(), xsi);

  R2 bounds = GetBounds_Sdirk(xsi, scheme.GetGammaCoefficient(), s, n);
  DISP(bounds);
  bounds *= Real_wp(1) / coef_scale(n-1);
  DISP(bounds);
  
  Vector<double> CflP, ErrP;
  if (n == 1)
    {
      VectReal_wp Gamma, xmin; 
      Linspace(param_min(0), param_max(0), N, Gamma);
      
      SequenceGrid1D grid(Gamma);
      grid.xsi = xsi;
      grid.Gamma = scheme.GetGammaCoefficient();
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
      grid.xsi = xsi;
      grid.Gamma = scheme.GetGammaCoefficient();
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
      grid.xsi = xsi;
      grid.Gamma = scheme.GetGammaCoefficient();
      //ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      ComputeCFLOptimal(fct,grid, 4, s, coef_scale, xmin);
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
      grid.xsi = xsi;
      grid.Gamma = scheme.GetGammaCoefficient();
      ComputeCFLParam(grid, CflP, ErrP, s, coef_scale);
      //ComputeCFLOptimal(fct,grid, 4, s, coef_scale, xmin);
    }
  else
    {
      cout << "Case not treated" << endl;
      abort();
    }
  
  if (rank_proc == 0)
    {
      DISP(CflP.GetNormInf());
      //CflP.WriteText("CflR.dat");
      //ErrP.WriteText("ErrR.dat");
    }
  
  return FinalizeMontjoie();
}
