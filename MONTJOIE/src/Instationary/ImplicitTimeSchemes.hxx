#ifndef MONTJOIE_FILE_IMPLICIT_TIME_SCHEMES_HXX

namespace Montjoie
{
  
  //! Base class for teta-scheme
  /*!
    This scheme can be used to solve the following second-order equation
    D_h d^2 U/dt^2 - K_h U + S_h dU/dt = F(t)
    The operator D_h is linear (i.e. it can be represented as a matrix)
    the operators K_h and S_h can be non-linear
    
    The implemented method is a theta-scheme :
    D_h (U^{n+1} - 2 U^n + U^{n-1})/dt^2 + Sh (U^{n+1} - U^{n-1}) / (2 dt)
    - Kh (\theta U^{n-1} + (1-2 \theta) U^n + \theta U^{n+1}) \, = \, F^n    
    
    Usually theta is taken equal to 1/4
   */
  template<class T>
  class TetaScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    Real_wp dt, dt2; //!< time step
    Real_wp teta; //!< teta-parameter
    Vector<T> U1, Ah_u1, Ah_u0, Prod_Uh, U0;  
    bool fourth_order_init;

  public :
    TetaScheme_Iterator();
    
    void SetTheta(const Real_wp& theta);
    Real_wp GetTheta() const;
    void SetFourthOrderInitialScheme();
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                              Vector<T>& du0_dt, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearSecond(Vector<T>&, Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };


  //! Base class for SDIRK (Singly Diagonally Implicit Runge-Kutta)
  template<class T>
  class SdirkScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    //! time step
    Real_wp dt;
    //! value of diagonal
    Real_wp gamma;
    //! value of Runge-Kutta coefficients A 
    Matrix<Real_wp> A;
    //! value of b and c
    Vector<Real_wp> b, c;
    //! Esdirk ?
    bool explicit_first_stage;
    //! alternative type (a -> 0, b -> 1) (e.g. S33a/S33b)
    int alternative_type;
    //! intermediary unknowns
    Vector<Vector<T> > Zn;
    //! main unknown
    Vector<T> Y;
    //! evaluation of functional
    Vector<T> Fn;
    //! combination of Zn
    Vector<T> Un_tmp;
    
  public :
    //! default constructor
    SdirkScheme_Iterator();
    
    void SetOrder(int r, bool esdirk_scheme = false, int type = 0);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                             Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();    
    const Vector<T>& GetIterate() const;
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
  //! base class listing order equations to satisfy for SDIRK schemes
  class SdirkNonLinearEquation : public VirtualLeastSquaresFunction<Real_wp>
  {
  public :
    Real_wp gamma; int r, s;
    UnivariatePolynomial<Real_wp> Numer;
    
    void Init(const UnivariatePolynomial<Real_wp>& P, const UnivariatePolynomial<Real_wp>& Q,
              const Real_wp& gam_, int order, Matrix<Real_wp>& A, VectReal_wp& b, VectReal_wp& c);
    
    void FillAB(const VectReal_wp& x, Matrix<Real_wp>& A, VectReal_wp& b, VectReal_wp& c);

    void GetStability(const Matrix<Real_wp>& A, const VectReal_wp& b,
                      UnivariatePolynomial<Real_wp>& N);

    void FindInitGuess(VectReal_wp& x);
    void EvaluateF(const VectReal_wp& x, VectReal_wp& f);

    void EvaluateJacobian(const VectReal_wp& x, VectReal_wp& f, Matrix<Real_wp, General, ColMajor>& df);
  
  };
  
  
  //! Base class for for LINEAR_SDIRK (Linear Singly Diagonally Implicit Runge-Kutta)
  template<class T>
  class LinearSdirkScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    //! time step
    Real_wp dt;
    // Value of diagonal DenomR = (1-gamma z)^{s-1 + extraS)
    Real_wp gamma;
    //! order of the scheme
    int order;
    
    //Stability function R = NumerR/DenomR of linear Sdirk
    // NumerR = DenomR*(1 + z + z^2/2! + ... + z^{s+1}/(s+1)!) + O(z^{s+extraS})
    UnivariatePolynomial<Real_wp> NumerR;
    //DenomR = (1-gamma z)^{s+extraS)
    UnivariatePolynomial<Real_wp> DenomR;

    //! Gauss lobatto points
    VectReal_wp ci;
    //! vect wi for each powL
    Matrix<Real_wp> mat_wi, polA;

    //! vectors used to compute U^n+1 from U^n
    //! main unknown
    Vector<T> Y;
    //! evaluation of functional F_n = DenomR*Y_(n+1) -  NumerR*Y_n
    Vector<T> Fn;
    Vector<Vector<T> > EvalF;  
    //! intermediary vectors
    Vector<T> KhUn;
    
    VectReal_wp real_roots; VectComplex_wp complex_roots;
    bool stable_algo, alternative_source;    
    bool sdirk_algo, lobatto_points;
    
    //! value of Runge-Kutta coefficients A 
    Matrix<Real_wp> A_rk;
    //! value of b and c
    Vector<Real_wp> b_rk, c_rk;
    
  protected:
    // case where no additional stages are added
    // We assume we know optimal gamma with DenomR = (1-gamma)^{s-1}
    // order is the scheme order
    void FindSdirkMinimal(int order_, Real_wp gamma_);
    

    // on calcule la fonction de stabilite R = Numer/Denom
    // en fonction de s (la methode est d'ordre s+1) et des parametres
    // param(0) = gamma
    // param(1:) = les constantes alpha, beta, etc
    void GetPolynomial_Sdirk(const VectReal_wp& param, int s,
			     UnivariatePolynomial<Real_wp>& Numer,
			     UnivariatePolynomial<Real_wp>& Denom,
			     UnivariatePolynomial<Real_wp>& Pol);

    void AddParamAlternativeSource(int i, UnivariatePolynomial<Real_wp>& Q);
    
    //Compute ci and mat_wi use to evaluate source function for Linear Sdirk Scheme
    void ComputeCoefForLinearSdirkSourceFunction();

  public :
    // type of algorithms
    enum {POLYNOMIAL, STABLE_WEIGHTS, SDIRK, POLY_ALTERNATIVE, STABLE_ALTERNATIVE};
    
    //! default constructor
    LinearSdirkScheme_Iterator();
    LinearSdirkScheme_Iterator(int s, int extraS, int stable_);

    inline const Matrix<Real_wp>& GetWeights() const { return mat_wi; }

    inline const Real_wp& GetGammaCoefficient() const { return gamma; }
    inline const Real_wp& GetQuadraturePoint(int i) const { return ci(i); }
    
    inline const UnivariatePolynomial<Real_wp>& GetNumeratorStabilityFunction() const { return NumerR; }
    inline const UnivariatePolynomial<Real_wp>& GetDenominatorStabilityFunction() const { return DenomR; }
    
    void SetOrder(int s, int extraS=0, int stable_ = 0);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                             Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void WriteCoefficients(string file_name);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();    
    const Vector<T>& GetIterate() const;
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };


  //! Base class for DIRK (Diagonally Implicit Runge-Kutta)
  template<class T>
  class DirkScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    //! time step
    Real_wp dt;
    //! values of diagonal
    //VectReal_wp gamma;
    //! value of Runge-Kutta coefficients A 
    Matrix<Real_wp> A;
    //! value of b and c
    Vector<Real_wp> b, c;
   //! intermediary unknowns
    Vector<Vector<T> > Zn;
    //! main unknown
    Vector<T> Y;
    //! evaluation of functional
    Vector<T> Fn;
    //! combination of Zn
    Vector<T> Un_tmp;
    
  public :
    //! default constructor
    DirkScheme_Iterator();
    
    void SetOrder(int);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                             Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();    
    const Vector<T>& GetIterate() const;
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
    //! Base class for the Pade scheme
  template<class T>
  class PadeScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    //! time step
    Real_wp dt;
    unsigned int order;
    //! Gauss lobatto points
    VectReal_wp ci;
    //! vect wi for each powL
    Matrix<Real_wp> mat_wi;
    
    //! real root if order = 4m+2
    Real_wp real_root;
    
    //! complex roots for order >= 4
    VectComplex_wp complex_root;

    //! coefficients in P - Q
    VectReal_wp coefPmQ;
    
    //! vectors used to compute U^n+1 from U^n
    
    //! main unknown
    Vector<T> Y;
    //! evaluation of functional F_n = Q_r*Y_(n+1) -  P_s*Y_n
    Vector<T> Fn;
    Vector<Vector<T> > EvalF;  
    //! intermediary vectors
    Vector<T> Un_tmp, KhUn;

    VectReal_wp Freal;
    //VectReal_wp Fimag;
    VectComplex_wp Fn_cplx, KhUn_cplx;

    // should we use stable algorithm ?
    bool stable_algo;

    // Compute n!
    Real_wp fac(int n)
    {
      Real_wp one(1);
      if(n == 0 || n == 1)
	return one;
      return n*fac(n-1);
    }
    
    void FindPadeCoeff(int order);

    void SolveOperatorP2(const Real_wp&, const Vector<T>&, const Complex_wp&,
			 int k, Vector<T>&, VirtualOdeSystem<T>& sys, bool apply_mass);
    
    //Compute ci and mat_wi use to evaluate source function for pade Scheme
    void ComputeCoefForPadeSourceFunction();

  public :
    //! default constructor
    PadeScheme_Iterator();
    PadeScheme_Iterator(int, bool s);
    
    const Matrix<Real_wp>& GetWeights() const { return mat_wi; }
    
    void SetOrder(int, bool s = false);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                             Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();    
    const Vector<T>& GetIterate() const;
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
  //! Base class for Gauss-Runge-Kutta schemes
  template<class T>
  class GaussRungeKutta_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    Real_wp dt;
    int order;
    Vector<T> Y, ki;
    Matrix<Real_wp> AgaussRK;
    Vector<Real_wp> BgaussRK, CgaussRK;

  public :
    //! default constructor
    GaussRungeKutta_Iterator();
    
    void SetOrder(int);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                             VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };
  
  
  //! Implementation of implicit Adams schemes
  /*!
    This scheme can solve a general first-order differential system
    dy/dt = f(t, y)
    
    Adams-Moulton (type_algo = MOULTON) schemes are written as
    (named also implicit Adams)
    
    y^{n+1} = y^n + dt \sum_{i=0}^order beta_i f(t^{n+1-i}, y^{n+1-i})
   */
  template<class T>
  class AdamsImplicit_Iterator : public AdamsBashforth_Moulton_Iterator<T>
  {
  protected :
    
  public :
    AdamsImplicit_Iterator();
    
    void SetOrder(int order);

    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
  };
  

  //! Implementation of implicit MilneSimpson schemes
  /*!
    This scheme can solve a general first-order differential system
    dy/dt = f(t, y)
    
    MilneSimpson (type_algo = MILNE) schemes are written as
    
    y^{n+1} = y^{n-1} + dt \sum_{i=0}^order beta_i f(t^{n+1-i}, y^{n+1-i})
   */
  template<class T>
  class MilneSimpson_Iterator : public Nystrom_Iterator<T>
  {
  protected :
    
  public :    
    MilneSimpson_Iterator();
        
    void SetOrder(int order);

    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
  };
  
  
  //! Backward differentiation formulas
  /*!
    This scheme can solve a general first-order differential system
    dy/dt = f(t, y)
    
    They are written in the following form :
    \sum \beta_k y_{n+1-k} = f(t^n+dt, y^{n+1})
   */
  template<class T>
  class BackwardDifferentiation_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    Real_wp dt;
    VectReal_wp beta;
    int order;
    Vector<Vector<T> > Yn_history;
    Vector<T> Y, Y_star, F_star, Fy;
    
  public : 

    void SetOrder(int);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                             VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
  //! Implicit Symmetric multistep scheme
  template<class T>
  class ImplicitSymmetricMultistep_Iterator : public VirtualTimeScheme<T>
  {
  protected :
    Real_wp t_begin_ener, t_end_ener, deltat_ener;
    WriteOnTheGoWithTinyBuffer<Real_wp,4> output_ener;
    int size_buffer_ener;
    string file_output_ener;
    VectReal_wp alpha, beta;
    Vector<Vector<T> > yk, fk;
    Vector<T> Y_star, F_star;
    Real_wp dt, energy;
    RungeKuttaNystrom_Iterator<T> rk_nystrom;    

  public :
    ImplicitSymmetricMultistep_Iterator();
    
    void GetAlphaBeta(VectReal_wp& A, VectReal_wp& B);
    void FindCoefAlphaBeta(const VectReal_wp& theta, const Real_wp& phi);
    
    void SetOrder(int order, const VectReal_wp& theta, const Real_wp& phi);
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                              Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearSecond(Vector<T>&, Vector<T>&);
   
    void SetParametersOutputEnergy(const Vector<string>& parameters);
    void SetEnergyFile(const Vector<string>& parameters);

    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;    
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    Real_wp ComputeEnergy(const Real_wp& t, VirtualOdeSystem<T>& sys, TinyVector<Real_wp,4>& energy);

    void WriteEnergy(const Real_wp& t, VirtualOdeSystem<T>& sys);
    
    size_t GetMemorySize() const;
    
  };
  
}

#define MONTJOIE_FILE_IMPLICIT_TIME_SCHEMES_HXX
#endif

