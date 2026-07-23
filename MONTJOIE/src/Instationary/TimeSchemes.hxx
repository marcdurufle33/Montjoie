#ifndef MONTJOIE_FILE_TIME_SCHEMES_HXX

namespace Montjoie
{

  //! available time schemes
  class TimeSchemeEnum
  {
  public :
    // LEAP_FROG : (u^(n+1) - u^(n-1))/(2 dt) = g(u^n)    
    // LEAP_FROG_EXACT : (u^(n+1) - u^(n-1))/2 = -i sin(i dt g(u^n))    
    // LEAP_FROG_ORDER2 : (u^(n+1) - 2 u^n + u^(n-1)) / dt^2 = g(u^n)    
    // LEAP_FROG_SYSTEM : (u^(n+1) - u^n)/dt = f(v^(n+1/2))
    //                    (v^(n+3/2) - v^n)/dt = g(u^(n+1))
    
    // OPTIMAL_LEAP_FROG : (u^(n+1) - 2 u^n + u^n+1) = Tn(dt^2 g(u^n))
    //            with Tn an optimal polynomial for the best CFL number
    // OPTIMAL_LEAP_FROG_PML : (u^(n+1) - u^(n-1))/2 = -i Tn(i dt g(u^n))
    //            with Tn an optimal polynomialfor the best CFL number 
    enum { LEAP_FROG, LEAP_FROG_EXACT, LEAP_FROG_ORDER2, LEAP_FROG_SYSTEM,
           OPTIMAL_LEAP_FROG, OPTIMAL_LEAP_FROG_PML,
           SPLIT_EXACT, SYMPLECTIC, NYSTROM, RUNGE_KUTTA_NYSTROM, 
           ADAMS_BASHFORTH, MULTISTEP_BUTCHER, TALEZER, TAYLOR_SERIES,
           RUNGE_KUTTA_LOW_STORAGE, RUNGE_KUTTA, LINEAR_RUNGE_KUTTA, RUNGE_KUTTA_EXACT,
           GAUSS_RUNGE_KUTTA, SDIRK, LINEAR_SDIRK, DIRK, CRANK_NICOLSON, ESDIRK, THETA_SCHEME,
           PREDICTOR, ADAMS_IMPLICIT, MILNE_SIMPSON, BDF, SYMMETRIC_MULTISTEP,
           LOCAL_IMPERIALE_SCHEME, LOCAL_PIPERNO_SCHEME, PADE_SCHEME,
	   IMPLICIT_SYMMETRIC_MULTISTEP, LOCAL_IMPLICIT_RK, LOCAL_IMEX_RK};  
    
  };


  //! base class to overload if you want to specify a new set of ordinary differential
  //! equations to be solved by time schemes implemented in Montjoie
  template<class T>
  class VirtualOdeSystem
  {
  public:
    // when you call RunTimeScheme, this method is called such that you 
    // can know the number of iterations that will be performed
    virtual void GiveNumberIterations(const Real_wp& dt, int nb_max_iter);
    
    // when you call RunTimeScheme, this method is called to give you the solution U^n
    virtual void GiveIterate(int n, const Real_wp&, Vector<T>&);
    
    // and the final iterate
    virtual void GiveFinalIterate(int n, const Real_wp&, Vector<T>&);
    virtual void GiveVectorialIterate(int n, const Real_wp&, Vector<T>&);


    // evaluation of main function f in system y' = f(t, y)
    virtual void EvaluateFunction(const Real_wp& t, const Vector<T>& y, Vector<T>& f,
				  bool invert_mass = true, bool source = true);

    // evaluation of main function f in system y'' = f(t, y, y')
    virtual void EvaluateFunctionS(const Real_wp& t, const Vector<T>& y, const Vector<T>& yp, Vector<T>& f,
				   bool invert_mass = true, bool source = true);
    
    virtual void SetDirichletCondition(const Real_wp& t, int n, Vector<T>&, Real_wp alpha = 1.0);
    virtual void SetDirichletConditionSource(const Real_wp& t, int n, Vector<T>&, Real_wp alpha = 1.0);
    
    // derivatives of f
    virtual void EvaluateDerivativeFunction(const Real_wp& t, int n,
					    const Vector<T>& y, Vector<T>& df,
					    bool invert_mass = true, bool source = true);
    
    // mass matrix for M y' = f(t, y)
    virtual void SolveMassMatrix(Vector<T>& X);    
    virtual void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
				 const Vector<T>& x, const Real_wp& beta, Vector<T>& y);
    
    // for Tal-ezer schemes  : source in separated variables (time and space)
    virtual Real_wp GetFinalTimeSource() const;
    virtual void FillSource(Vector<T>&) const;
    
    virtual void SourceOnlyTime(const Real_wp& t0, int n, T& pulse);

    // operators for scalar equation D_h d^2 U/dt^2 + K_h U = F
    virtual void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
				 const Real_wp& beta, Vector<T>& Yh);

    virtual void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
				 const Real_wp& beta, Vector<T>& Prod_Uh);
    
    virtual void SolveOperatorDh(Vector<T>& Y);

    virtual void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
					  const Real_wp& beta, Vector<T>& Prod_Uh);

    virtual void SolveOperatorDhPlusdtSh(Vector<T>& Y);    
    virtual void SolveOperatorDhPlusdtSh(const Real_wp&, const Vector<T>& X, Vector<T>& Y);

    virtual void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
				 const Real_wp& beta, Vector<T>& Prod_Uh);
    
    // operators for split first-order formulation
    virtual void ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
					  const Real_wp& beta, Vector<T>& Vh, bool extrapol = true);

    virtual void ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
				       const Real_wp& beta, Vector<T>& Vh, bool extrapol = true);
        
    virtual void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
					  const Real_wp& beta, Vector<T>& Prod_Uh);
    
    virtual void SolveOperatorBh(Vector<T>& Y);

    virtual void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
					  const Real_wp& beta, Vector<T>& Vh);
    
    virtual void SolveOperatorBhPlusdtSh(Vector<T>& Y);
    
    // sources
    virtual void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
					int nb_deriv, Vector<T>& b_src);
    
    virtual void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, Vector<T>& Y);
    virtual void AddVectorialTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, Vector<T>& Y);
    
    // for split-schemes
    virtual void SolveSplitOperator(const Real_wp& t, const Real_wp& dt,
				    const Vector<T>& Un, Vector<T>& Un_next, int num);
    
    // for implicit-schemes (only one system to solve)
    virtual void FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta,
						const Real_wp& gamma);    
    
    virtual void SolveOperatorDhPlusGammaKh(const Real_wp& t, const Vector<T>& X, Vector<T>& Y);

    // DIRK case (with several real coefficients)
    virtual void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
				       const VectReal_wp& gamma);    

    virtual void SolveOperatorReal(const Real_wp& t, const Vector<T>& X, Vector<T>& Y,
				   int num_system);
    
    // Gauss-Runge-Kutta case
    virtual void FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& AgaussRK, const Real_wp& dt,
					      const Vector<T>& Y, const Vector<T>& ki);
    
    virtual void SolveOperatorDhPlusAijKh(const VectReal_wp& t, const Matrix<Real_wp>& A,
					  const Real_wp& dt, const Vector<T>& Y, Vector<T>& ki);
    
    
    // For the Pade schemes (complex roots)
    virtual void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
					  const VectComplex_wp& gamma);

    virtual void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y, int);
    
  };

  
  //! base class for any time scheme
  template<class T>
  class VirtualTimeScheme
  {
  public:
    virtual ~VirtualTimeScheme();
    
    // initializes the scheme with the initial condition U0 and time step dt
    virtual void SetInitialCondition(const Real_wp& t, const Real_wp& dt,
				     Vector<T>& U0, VirtualOdeSystem<T>& sys);

    virtual void SetInitialConditionS(const Real_wp& t, const Real_wp& dt,
                                      Vector<T>& U0, Vector<T>& dU0, VirtualOdeSystem<T>& sys);
    
    virtual void SetParametersOutputEnergy(const Vector<string>& parameters);
    virtual void SetEnergyFile(const Vector<string>& parameters);
    
    // advances the scheme, ie computes U^n+1 from previous iterates
    virtual void Advance(const Real_wp& t, int n, VirtualOdeSystem<T>& sys) = 0;
    
    // changes the time step
    virtual void ChangeTimeStep(const Real_wp&) = 0;
    
    // returns the iterate U^n
    virtual Vector<T>& GetIterate() = 0;
    virtual const Vector<T>& GetIterate() const = 0;

    virtual Vector<T>& GetVectorialIterate();
    virtual int GetNumberOfVectorialIterates() const;
    virtual Vector<T>& GetVectorialIterate(int k);

    // returns the number of stored vectors needed to advance to the next time step
    virtual int GetNumberOfIterates() const = 0;
    // returns the stored vector k needed to advance to the next time step
    virtual Vector<T>& GetIterate(int k) = 0;

    // releases memory used by vectors
    virtual void Clear() = 0;
    // releases memory and returns in Un the current iterate
    virtual void ClearFirst(Vector<T>& Un);
    virtual void ClearSecond(Vector<T>& Un, Vector<T>& dUn);
    
    virtual Real_wp GetCFL(bool real_p) const;
    virtual UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    virtual size_t GetMemorySize() const;
    virtual void GetMemoryUsed(map<string, size_t>& var) const;
    
  };
  
  
  //! Class for solving evolution systems with any time-scheme
  template<class T>
  class All_TimeScheme : public InputDataProblem_Base
  {
  protected :
    //! type of time scheme
    int type_time_scheme;
    
    //! order of time scheme
    int order_time_scheme;

    //! Pointer to the selected time scheme
    VirtualTimeScheme<T>* time_scheme;
    
    //! number of iterates used by time scheme
    int nb_iterates;
    
  public :
    All_TimeScheme();
    ~All_TimeScheme();
    
    VirtualTimeScheme<T>& GetTimeScheme() const;

    size_t GetMemorySize() const;
    void GetMemoryUsed(map<string, size_t>& var) const;
    
    int GetTimeSchemeType() const;
    int GetOrder() const;
    
    void SetInputData(const string& keyword, const Vector<string>& parameters);

    bool FirstOrderScheme() const;
    bool StaggeredScheme() const;
    bool IsImplicitScheme() const;
    bool LocalTimeStepping() const;
    void SetTimeScheme(int r, const Vector<string>& type);
        
    void CheckPresenceScheme() const;
    void ChangeTimeStep(const Real_wp& dt);

    Real_wp GetCflScheme(bool real_p) const;
    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void SetInitialCondition(const Real_wp& t0, const Real_wp& dt, Vector<T>& U0,
                             VirtualOdeSystem<T>& sys);
    
    void SetInitialCondition(const Real_wp& t0, const Real_wp& dt, Vector<T>& U0,
                             Vector<T>& V0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& z, int n, VirtualOdeSystem<T>& sys);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;
    
    Vector<T>& GetVectorialIterate();
    
    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    void Clear();
    void Clear(Vector<T>& U);
    void Clear(Vector<T>& U, Vector<T>& V);
    
  };
  
}

#define MONTJOIE_FILE_TIME_SCHEMES_HXX
#endif

