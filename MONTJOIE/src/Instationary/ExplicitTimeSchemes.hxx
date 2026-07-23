#ifndef MONTJOIE_FILE_EXPLICIT_TIME_SCHEMES_HXX

namespace Montjoie
{

  template<class T>
  void RunTimeScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
		     Vector<T>& Y0, VirtualOdeSystem<T>& sys, VirtualTimeScheme<T>& scheme,
                     int num_iter0 = 0);
  
  //! base class for explicit Runge-Kutta schemes
  /*!
    We consider explicit Runge-Kutta schemes for the resolution of the system :
    dy/dt = f(t, y)
    
   These schemes can be written in the following form
   k_i = f(t^n + c_i dt, y^n + dt \sum_j a_{i, j} k_j)
   y^{n+1} = y^n + dt \sum_i b_i k_i
   
   The coefficients c_i, b_i and a_{i, j} are stored in internal arrays C, B, A
   and are automatically filled when SetOrder is called

   See the comments of RunTimeScheme to see how to use these schemes
   */
  template<class T>
  class RungeKutta_Iterator : public VirtualTimeScheme<T>, public NonLinearEquations_Newton<Real_wp>
  {
  protected :
    VectReal_wp B; //!< coefficients B 
    VectReal_wp C; //!< coefficients C
    Matrix<Real_wp> A; //!< coefficient matrix A
    VectReal_wp Bembed; //!< embedded coefficients B
    Real_wp dt; //!< time step
    Vector<Vector<T> > ki; //!< intermediary vectors
    Vector<T> Y, Y_tilde; //!< main iterate and  intermediary vector
    int type_algo, order_rk;
    
  public :
    enum {BUTCHER, FEHLBERG, DOPRI};
    
    RungeKutta_Iterator();
    
    void SetOrder(int order, int type = 0);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    Real_wp GetCFL(bool real_p) const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    void GetCoefficients(Matrix<Real_wp>& a, VectReal_wp& b, VectReal_wp& c) const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  protected : 
    // variables used to refine coefficients in multiple precision
    Matrix<Real_wp> mat_jacobian_lu; IVect pivot_jacobian;
    
  public :
    // optimisation of RK scheme 
    void ComputeScheme(const VectReal_wp& X, VectReal_wp& F);
    Real_wp GetNorm2Vector(const VectReal_wp& scheme);

    void ComputeAndFactoriseDiff(const VectReal_wp& X, const VectReal_wp& scale);
    void SolveDifferential(const VectReal_wp& F, VectReal_wp& X);
    
    void EvaluateFunction(const VectReal_wp&, VectReal_wp&);
    void EvaluateJacobian(const VectReal_wp&, Matrix<Real_wp>&);
    
  };


  template<class T>
  class LinearRungeKutta_Iterator : public VirtualTimeScheme<T>
  {
  protected:
    //! Polynomial stability function
    UnivariatePolynomial<Real_wp> Pol;

    //! reals roots of the polynomial
    VectReal_wp real_roots;
    //! complex roots of the polynomial
    VectComplex_wp complex_roots;

    bool stable_algo, stable_source;
    Real_wp dt;
    Vector<T> Y, Fn;
    Vector<T> KhUn;
    int order, extraS;
    
    //! Quadrature points (Gauss Lobatto or Gauss Legendre)
    VectReal_wp ci;

    //! vect wi for each powL
    Matrix<Real_wp> mat_wi, polA, DerMat;
    
    Vector<Vector< T> >   EvalF;
    
    //Compute ci and mat_wi use to evaluate source function for Linear Runge Kutta Scheme
    void ComputeCoefForLinearRungeKuttaSourceFunction();

  public:
    LinearRungeKutta_Iterator();
    
    void SetOrder(int order, int extraS = 0, bool s = true,
                  bool s_source = false);
    
    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    const VectReal_wp& GetPoints();
    const Matrix<Real_wp>& GetWeights();
    const Matrix<Real_wp>& GetDerivativeMatrix();
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);
    Real_wp GetCFL(bool real_p) const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };
  

  //! base class for Talezer schemes
  /*!
    This scheme can only be applied to the following linear system of equations
    dy/dt = f(t, y)
    with f(t, y) = h(t) B + A y
    with a matrix A, a vector B and h a scalar function
    
    The scheme is obtained by exploiting the Jacobi-Anger expansion of exp(i k R cos teta)
    involving Bessel functions. We denote by R an upper bound of the eigenvalues of A, the scheme
    is then equal to :

    y^{n+1} = \sum_k epsilon_k J_k(dt R) Q_k(A/R) y^n + \sum_k b_k Q_k(A/R) B
    
    with epsilon_k = 1 if k is equal to 0, 2 otherwise
    J_k is the bessel function of order k
    Q_k are "imaginary" Tchebycheff polynomials determined by the recurrence relation
    Q_{k+1}(X) = 2 X Q_k + Q_{k-1}  (here X = A/R)
    and Q_0 = I, Q_1 = X
    
    the coefficients b_k are obtained as the integrals :
    b_k = \int_{0}^{dt} epsilon_k J_k(tau R) h((n+1)dt - tau) dtau
    these coefficients are computed with a Gauss-Legendre integration
    
    Usually, in this scheme, instead of decreasing dt, it is better to increase the
    order (k goes from 0 to order). By increasing slightly k, we can obtain quickly
    a machine precision solution because of spectral convergence.
    
    In the system of equations, you need to implement the following functions :
    
    // returns the frequency associated with the source h
    Real_wp GetFrequency();
    
    // computes pulse = d^n h/dt (t), n is always equal to 0
    void SourceOnlyTime(const Real_wp& t, int n, Real_wp& pulse);
    
    // sets Dirichlet condition u_i = f(t)
    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y);
    
    // computes result = A y
    // this function is always called with add_source = false (only A y is needed)
    void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& result,
                          bool invert_mass = true, bool add_source = true);
    
    // returns the final time after which the source h(t) is vanishing
    Real_wp GetFinalTimeSource();
    
    // fills a vector with B
    void FillSource(VectReal_wp& b);
    
    // inversion by the mass matrix : replaces b by M^{-1} b
    // in pratice the system considered is M dy/dt = h(t) B + A y
    // M being the mass matrix, therefore we need to compute M^-1 B
    void SolveOperatorDh(VectReal_wp& b);
     
    
    In order to use the scheme, you can write :

    Talezer_Iterator<Real_wp> scheme;
    
    // setting the value R (upper bound of spectral radius of A)
    // this function must be called before SetOrder
    scheme.SetSpectralRadius(4.5);
    
    // providing the order to use for quadrature formula involved in coefficients b_k
    scheme.SetOrderIntegration(12);

    // providing the order (here the order is the number of elements 
    // to take into account in the sum, k goes from 0 to order-1)
    scheme.SetOrder(20, dt);
    
    // then you can call RunTimeScheme if you want
    RunTimeScheme(t0, tf, dt, u0, sys, scheme);
   */
  template<class T>
  class Talezer_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    VectReal_wp coef; //!< coefficients
    Real_wp Rmax; //!< parameter R (max eigenvalue of the operator, usually)
    Vector<Vector<Real_wp> > coef_bessel; //!< evaluation of bessel functions 
    int order_integration_bessel; //!< number of bessel functions
    Vector<T> coef_bk; //!< coefficients
    Globatto<Real_wp> lob; //!< quadrature formula for source
    Vector<T> Y, By, Qn, Qnp1, Bsource, Bn, Bnp1; //!< intermediary vectors
    Real_wp dt; //!< time step
    
  public :

    Talezer_Iterator();
    
    void SetOrder(int order, const Real_wp& deltat);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;    
    void SetSpectralRadius(const Real_wp& R);
    void SetOrderIntegration(int r);
    
    void ComputeCoefficientsSource(const Real_wp& t0, const Real_wp& deltat, VirtualOdeSystem<T>& sys);

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

  
  //! base class for Adams-Bashforth-Moulton schemes
  /*!
    This scheme can solve a general first-order differential system
    dy/dt = f(t, y)
    
    Adams-Bashforth (type_algo = BASHFORTH) schemes are written as :
    
    y^{n+1} = y^n + dt \sum_{i=0}^order beta_star_i f(t^{n-i}, y^{n-i})
    
    
    Adams-Moulton (type_algo = MOULTON) schemes are written as
    (named also implicit Adams)
    
    y^{n+1} = y^n + dt \sum_{i=0}^order beta_i f(t^{n+1-i}, y^{n+1-i})
    
    However, this scheme is implicit, therefore actually implemented
    in ImplicitTimeSchemes.cxx (class AdamsImplicit_Iterator)
    
    Adams-Bashforth-Moulon (type_algo = BASHFORTH_MOULTON)
    is a predictor-corrector scheme, that can be written as :
    
    y^{predict} = y^n + dt \sum_{i=0}^order beta_star_i f(t^{n-i}, y^{n-i})
    y^{n+1} = y^n + dt beta_0 f(t^{n+1}, y^{predict})
                 + dt \sum_{i=1}^order beta_i f(t^{n+1-i}, y^{n+1-i})
    
    this last-scheme is explicit
   */
  template<class T>
  class AdamsBashforth_Moulton_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    VectReal_wp beta_star; //!< coefficients
    VectReal_wp beta; //!< coefficients
    Vector<Vector<T> > fi;
    Vector<T> Y_star, F_star, Y;
    Real_wp dt;
    int type_algo;
    RungeKutta_Iterator<T> rk_scheme;
    
  public :
    // different algorithms
    enum {BASHFORTH, MOULTON, BASHFORTH_MOULTON};
    
    AdamsBashforth_Moulton_Iterator();
    
    void SetOrder(int order, int type = 0);

    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
  //! base class for Taylor schemes for a linear operator
  /*!
    This scheme can only be applied to the following linear system of equations
    dy/dt = f(t, y)
    with f(t, y) = B(t) + A y
    with a matrix A, a vector B(t), only B depends on t
    
    we derive the equation n times to obtain all the derivatives of y
    dy/dt = B(t^n) + A y
    d^2 y/dt^2 = dB/dt(t^n) + A dy/dt
    d^3 y/dt^3 = d^2 B/dt^2 + A d^2 y/dt^2
    ...
    d^k y/dt^k = d^k B/dt^k + A d^{k-1} y/dt^{k-1}
    
    The solution y^(n+1) is obtained with the Taylor expansion :
    y^{n+1} = y^n + \sum_{k=1}^K  d^k y/dt^k  dt^k/k!

    In the system of equations, you need to implement the following function :
    
    // should compute result = d^k B/dt^k + A u
    void EvaluateDerivativeFunction(const Real_wp& t, int k,
                                 const VectReal_wp& u, VectReal_wp& result);
    
    This scheme is equivalent to Runge-Kutta schemes when B = 0, and when the
    order r is lower or equal to 4. When r >= 5, the schemes are different because Runge-Kutta
    schemes comprise more than r+1 stages.
   */
  template<class T>
  class TaylorSeries_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    int order; //!< order of the time scheme
    Vector<Real_wp> dtk_div_factork; //!< coefficients dt^k / k!
    Vector<T> Y, Y_tilde, Fy; //!< intermediary vectors
    Real_wp dt; //!< time step
    
  public :

    TaylorSeries_Iterator();
    TaylorSeries_Iterator(int order_);
    
    void SetOrder(int r);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();    
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };

  
  //! base class for low-storage Runge-Kutta schemes
  /*!
    low-storage Runge-Kutta schemes using 2N storage for the resolution of
    dy/dt = f(t, y)
    However, since Montjoie does not implement an additive version of the evaluation of f :
    u <- u + alpha f(t0, y0)
    but a non-additive version :
    u <- f(t0, y0)
    this class stores 3 vectors instead of 2 vectors for the optimal schemes
    
    The scheme can be written as the following algorithm :
    U = y^n
    Do i = 1, s
      dU = alpha_i dU + dt f(t + c_i dt, U)
      U = U + beta_i dU
    End
    y^n+1 = U

    The coefficients alpha, beta, and c are derived in the papers of Carpenter/Kennedy
   */
  template<class T>
  class LowStorageRK_Iterator  : public VirtualTimeScheme<T>
#ifdef MONTJOIE_WITH_MULTIPLE
                               , public NonLinearEquations_Newton<Real_wp>
#endif
  {
  protected :
    VectReal_wp coef_alpha; //!< coefficients
    VectReal_wp coef_beta; //!< coefficients
    VectReal_wp coef_C; //!< coefficients
    Vector<T> Y, Y_tilde, Fy; //!< intermediary vectors
    Real_wp dt; //!< time step
    
  public :
    LowStorageRK_Iterator();
    
    void SetOrder(int r, int nb_steps = 0);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    Real_wp GetCFL(bool real_p) const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
#ifdef MONTJOIE_WITH_MULTIPLE
  protected : 
    // variables used to refine coefficients in multiple precision
    VectReal_wp param_scheme;
    Matrix<Real_wp> mat_jacobian_lu; IVect pivot_jacobian;
    
  public :
    template<class Vector1, class Vector2>
    void FillGlobalVector(const Vector1& A, const Vector1& B, Vector2& X);

    template<class Vector1, class Vector2>
    void FillVectorAB(const Vector1& X, Vector2& A, Vector2& B);

    template<class Vector1, class Vector2>
    void GetCoefficientsC(const Vector1& A, const Vector1& B, Vector2& C);

    void GetCoefficientsAB(const VectReal_wp& A, const VectReal_wp& B,
                           Matrix<Real_wp>& a, VectReal_wp& b);
    
    void ComputeScheme(const VectReal_wp& X, VectReal_wp& F);
    Real_wp GetNorm2Vector(const VectReal_wp& scheme);

    void ComputeAndFactoriseDiff(const VectReal_wp& X, const VectReal_wp& scale);
    void SolveDifferential(const VectReal_wp& F, VectReal_wp& X);
    
    void EvaluateFunction(const VectReal_wp&, VectReal_wp&);
    void EvaluateJacobian(const VectReal_wp&, Matrix<Real_wp>&);
#endif
    
  };

  
  //! base class for multistep butcher's schemes
  /*!
    This scheme can only be applied to the following linear system of equations
    dy/dt = f(t, y)
    with f(t, y) = B(t) + A y
    with a matrix A, a vector B(t), only B depends on t
    
    we derive the equation n times to obtain all the derivatives of y
    dy/dt = B(t^n) + A y^n
    d^2 y/dt^2 = dB/dt(t^n) + A dy/dt
    d^3 y/dt^3 = d^2 B/dt^2 + A d^2 y/dt^2
    ...
    d^k y/dt^k = d^k B/dt^k + A d^{k-1} y/dt^{k-1}
    
    The scheme is written as :
    y^{n+1} = y^n + \sum_{k=0}^r dt^k alpha_k d^k y/dt^k(t^n) 
                  + \sum_{k=0}^r dt^k beta_k d^k y/dt^k(t^{n-1})
    
    A correction step can be performed :
    y^{n+1} = y^n + \sum_{k=0}^r dt^k gamma_k d^k y/dt^k(t^{n+1})
                  + \sum_{k=0}^r dt^k delta_k d^k y/dt^k(t^{n})
    where the values d^k y/dt^k(t^{n+1}) are replaced by the y^{n+1} evaluated
    during the first step (also called predictive step)
    
    This scheme is described in the Butcher's book (p 119)

    In the system of equations, you need to implement the following function :
    
    // should compute result = d^k B/dt^k + A u
    void EvaluateDerivativeFunction(const Real_wp& t, int k,
                                    const VectReal_wp& u, VectReal_wp& result);
   */
  template<class T>
  class MultiStepButcher_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    int order; //!< order of the time scheme
    VectReal_wp alpha; //!< coefficients
    VectReal_wp beta; //!< coefficients
    VectReal_wp gamma; //!< coefficients
    VectReal_wp delta; //!< coefficients
    int type_algo;    
    Vector<Vector<T> > fnm1, fnm2;
    Vector<T> Y_star, Y;
    Real_wp dt;
    
  public :
    enum {NO_CORRECTION, PECE};
    
    MultiStepButcher_Iterator();
    
    void SetOrder(int r, int type = 0);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };


  //! base class for operator-splitting schemes
  /*!
    This scheme can be used to solve the following class of equations :
    du/dt = A u + B u
  */
  template<class T>
  class SplitScheme_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    int order;
    Real_wp dt;
    Vector<T> Un, Un_next;
    Vector<Real_wp> coef_ci, coef_di;
    
  public :
    SplitScheme_Iterator();
    
    void SetOrder(int r);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };
  
  
  //! base class for Nystrom schemes
  /*!
    This class can solve differential equations of the form
    dy/dt = f(t, y)
    
    Nystrom schemes can be written as
    y^{n+1} = y^{n-1} + dt \sum beta_i f(t^{n-i}, y^{n-i})
    
    For a second-order scheme, we have beta_0 = 2, i.e. the classical leap-frog scheme
   */
  template<class T>
  class Nystrom_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    VectReal_wp beta_star; //!< coefficients
    VectReal_wp beta; //!< coefficients
    Vector<Vector<T> > fi;
    Vector<T> Y_star, F_star, Y, Ynm1;
    Real_wp dt;
    int type_algo;
    RungeKutta_Iterator<T> rk_scheme;
    
  public :
    // different algorithms
    enum {NYSTROM, MILNE, MILNE_PREDICTOR};
    
    Nystrom_Iterator();
    
    void SetOrder(int order, int type = 0);

    void SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			     Vector<T>& u0, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearFirst(Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };


  //! base class for Runge-Kutta-Nystrom schemes
  /*!
    These schemes solve the following system :
    d^2 y/dt^2 = f(t, y, dy/dt)
    
    A special case appears when f does not depend on dy/dt :
    d^2 y/dt^2 = f(t, y)
    
    By default, it will be considered that f does not depend on dy/dt    
   */
  template<class T>
  class RungeKuttaNystrom_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    Vector<Vector<T> > ki;
    Vector<T> Y, Yprime, Ytilde, Yprime_tilde;
    bool f_dependent_from_yprime; int order_rk;
    Matrix<Real_wp> A, Abar;
    Vector<Real_wp> B, Bbar, C;
    Real_wp dt;
    int type_pair;

  public :
    enum {OPTIMAL_PAIR, DORMAND_PRINCE, SHARP_QURESHI_GRAZIER};
    
    RungeKuttaNystrom_Iterator();
    
    void SetOrder(int order, bool dependent = false, int type = OPTIMAL_PAIR);
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                              Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys);
    
    void GetCoefficients(Matrix<Real_wp>& coefA, Matrix<Real_wp>& coefAbar,
			 Vector<Real_wp>& coefB, Vector<Real_wp>& coefBbar,
			 Vector<Real_wp>& coefC) const;
    
    void Clear();
    void ClearSecond(Vector<T>&, Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);    
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;    

    Vector<T>& GetDeriveIterate();
    const Vector<T>& GetDeriveIterate() const;    

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };
  
  
  //! Symmetric multistep scheme
  template<class T>
  class SymmetricMultistep_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    VectReal_wp alpha, beta;
    Vector<Vector<T> > yk, fk;
    Vector<T> Y_star, F_star;
    Real_wp dt;
    RungeKuttaNystrom_Iterator<T> rk_nystrom;
    
  public :
    SymmetricMultistep_Iterator();
    
    void GetAlphaBeta(VectReal_wp& A, VectReal_wp& B);
    void FindCoefAlphaBeta(const VectReal_wp& theta);
    
    void SetOrder(int order, const VectReal_wp& theta);
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                              Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys);
    
    void Clear();
    void ClearSecond(Vector<T>&, Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;    
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };
  
}

#define MONTJOIE_FILE_EXPLICIT_TIME_SCHEMES_HXX
#endif
