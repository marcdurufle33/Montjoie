#ifndef MONTJOIE_FILE_MODIFIED_EQUATION_HXX

namespace Montjoie
{
  //! modified equation scheme for second-order wave equation
  /*!
    This scheme can be applied to the following linear system :
    D_h d^2 U/dt^2 - K_h U + S_h dU/dt = F(t)
    
    The second-order time scheme can be written as :
    D_h (U^{n+1} - 2 U^n + U^{n-1})/dt^2 + S_h (U^{n+1} - U^n)/(2 dt) - K_h U^n = F(t^n)
    The fourth-order scheme is obtained via the "modified equation" technique :
    
    - The term (U^{n+1} - 2 U^n + U^{n-1})/dt^2 is replaced by 
    (U^{n+1} - 2 U^n + U^{n-1})/dt^2 + dt^2/12 d^4 U/dt^4  
    which is a fourth-order approximation of d^2 U/dt^2
    
    - The term S_h (U^{n+1} - U^n)/(2 dt) is replaced by
    S_h (U^{n+1} - U^n)/(2 dt) + dt^2/6 d^3 U/dt^3
    
    - The values of d^3 U/dt^3, and d^4 U/dt^4 are obtained by derivating the equation :
    D_h d^3 U/dt^3 = K_h dU/dt - S_h d^2 U/dt^2 + d^2/dt^2 ( F(t) )
    D_h d^4 U/dt^4 = K_h d^2 U/dt^2 - S_h d^3 U/dt^3 + d^3/dt^3 ( F(t) )
    the value of dU/dt being obtained with an uncentered approximation with U^n, U^{n-1}, U^{n-2}
    
    This scheme is compatible with the main function RunSecondOrderScheme,
    and can be used in the following way :
    
    ModifiedEquationIterator<Real_wp> scheme;
    
    // only order 2 and 4 are implemented
    scheme.SetOrder(4);
    
    // you need to prepare the resolution of operators D_h + dt/2 S_h (for second order),
    // and operator D_h as well for fourth order
    SystemEquation sys;    
    Real_wp dt = 0.1;
    sys.InitComputation(dt);

    // then the scheme can be used with RunSecondOrderScheme for example
    RunSecondOrderScheme(t0, tf, dt, U0, dU0_dt, sys, scheme);
    
    // on exit, only U0 contains U(tf), dU0_dt is empty
    // in SystemEquation, you need to implement the following functions
    // (example given with real numbers):

    // computes U = U + alpha d^{n-1} F/dt^{n-1}
    void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& U);
    
    // computes v = beta v + alpha S_h u
    void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                          const Real_wp& beta, VectReal_wp& v);
    
    // computes v = beta v + alpha K_h u
    void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                         const Real_wp& beta, VectReal_wp& v);

    // computes v = beta v + alpha D_h u
    void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                       const Real_wp& beta, VectReal_wp& v);

    // assembles vector u (useful for parallel computation)
    void Assemble(VectReal_wp& u);
  
    // replaces u by D_h^{-1} u 
    void SolveOperatorDh(VectReal_wp& u);

    // computes v = beta v + (D_h - dt/2 S_h) u  
    void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& u,
                                const Real_wp& beta, VectReal_wp& v);

    // replaces u by (D_h + dt/2 S_h)^{-1} u
    void SolveOperatorDhPlusdtSh(VectReal_wp& u);
  
    // sets Dirichlet condition u_i = d^n g_i/dt^n
    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& u, Real_wp alpha = 1);
   */
  template<class T>
  class ModifiedEquationIterator  : public VirtualTimeScheme<T>
  {
  protected :
    int order;
    Vector<T> U0, U1, Ah_u1, U2, Prod_Uh;
    Vector<T> dU0_dt, dU_dt2, dU_dt3, dU_dt4;
    Real_wp dt, dt2;
    
  public :
    ModifiedEquationIterator();
    
    void SetOrder(int);
    int GetOrder() const;
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                              Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    void Clear();
    void ClearSecond(Vector<T>&, Vector<T>&);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    Real_wp GetCFL(bool real_p) const;
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    size_t GetMemorySize() const;
    
  };


  template<class T>
  void RunSecondOrderScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
			    Vector<T>& U0, Vector<T>& dU0_dt, VirtualOdeSystem<T>& sys,
			    VirtualTimeScheme<T>& scheme);
  

  //! modified equation scheme for first-order wave equation
  /*!
    This scheme is intended to solve a first-order split system of the following form :
    D_h dU/dt - R_h^S V + S_h U = F^S(t)
    B_h dV/dt - R_h^V U + S_h^V V = F^V(t)
    where B_h, D_h, R_h^V, R_h^V, S_h, S_h^V are matrices independent of time, 
    and F^S(t), F^V(t) are vectors depending on time

    This system is solved with a leap-frog method for second-order method
    D_h (U^(n+1) - U^n)/dt - R_h^S V^{n+1/2} + S_h (U^{n+1} + U^n)/2 = F^S(t^{n+1/2})
    B_h (V^{n+1/2} - V^{n-1/2})/dt - R_h^V U^n + S_h^V (V^{n+1/2} + V^{n-1/2})/2 = F^V(t^n)
    
    The fourth-order scheme is obtained with a modified equation approach
    (U^(n+1) - U^n)/dt is replaced by (U^(n+1) - U^n)/dt + dt^2/24 d^3 U/dt^3
    (U^(n+1) + U^n)/dt is replaced by (U^(n+1) + U^n)/2 + dt^2/8 d^2 U/dt^2
    and similarly for terms in V
    
    The derivatives dU/dt, dV/dt, d^2 V/dt^2, d^2 U/dt^2, d^3 V/dt^3 are evaluated by
    using the equation :
    
    dU/dt = D_h^{-1} ( R_h^S V - S_h U + F^S)
    dV/dt = B_h^{-1} ( R_h^S U - S_h^V V + F^V)
    d^2 U/dt^2 = D_h^{-1} ( R_h^S dV/dt - S_h dU/dt + d F^S / dt)
    d^2 V/dt^2 = B_h^{-1} ( R_h^S dU/dt - S_h^V dV/dt + d F^V / dt)
    ...
    
    The values of U at time t^{n+1/2} and V at time t^{n+1} are obtained 
    with a second-order uncentered approximation :
    V^{n+1} = 3/2 V^{n+1/2} - 1/2 V^{n-1/2}
    U^{n+1/2} = 3/2 U^{n} - 1/2 U^{n-1}

    This scheme is compatible with the main function RunFirstOrderScheme,
    and can be used in the following way :
    
    ModifiedEquationSystemIterator<Real_wp> scheme;
    
    // only order 2 and 4 are implemented
    scheme.SetOrder(4);
    
    // you need to prepare the resolution of operators D_h + dt/2 S_h (for second order),
    // and operator D_h as well for fourth order
    SystemEquation sys;    
    Real_wp dt = 0.1;
    sys.InitComputation(dt);

    // then the scheme can be used with RunSecondOrderScheme for example
    RunFirstOrderScheme(t0, tf, dt, U0, V0, sys, scheme);
    
    // on exit U0 contains the solution at time t = tf
   */
  template<class T>
  class ModifiedEquationSystemIterator  : public VirtualTimeScheme<T>
  {
  protected :
  public :
    int order; Real_wp dt, dt2;
    Vector<T> Un, Vn, Prod_Uh, Prod_Vh;
    Vector<T> dU_dt, dU_dt2, dU_dt3, dV_dt, dV_dt2, dV_dt3,
      Un_m1, Vn_m05, Vn_p1, Un_demi;

  public :
    ModifiedEquationSystemIterator();
    
    void SetOrder(int);
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                              Vector<T>& v0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    void Clear();
    void ClearSecond(Vector<T>& U, Vector<T>& V);

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    Real_wp GetCFL(bool real_p) const;
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;
    Vector<T>& GetVectorialIterate();

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    int GetNumberOfVectorialIterates() const;
    Vector<T>& GetVectorialIterate(int k);

    size_t GetMemorySize() const;
    
  };

  template<class T>
  void RunFirstOrderScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
			   Vector<T>& U0, Vector<T>& V0, VirtualOdeSystem<T>& sys, VirtualTimeScheme<T>& scheme);

  
  template<class T>
  class OptimalModifiedEquationSecond_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    Vector<T> Un;
    
  public :
    void SetOrder(int) {}
    
    void SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                              Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys) {}
    
    void Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys) {}
    
    void ChangeTimeStep(const Real_wp&) {}
    
    void Clear() {}
    void ClearSecond(Vector<T>&, Vector<T>&) {}

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const {UnivariatePolynomial<Real_wp> P; return P; }
    
    Vector<T>& GetIterate() { return Un; }
    const Vector<T>& GetIterate() const { return Un; }
    
    int GetNumberOfIterates() const { return 2; }
    Vector<T>& GetIterate(int k) { return Un; }

    size_t GetMemorySize() const { return 0; }
        
  };  


  //! Optimal modified equation on leap-frog scheme
  /*!
    The following schemes solves equation :
    dU / dt = A U + F(t) 
    where A is skew-symmetric, F is the source term
    The used second-order scheme is the leap frog scheme :
    (U^{n+1} - U^{n-1}) / (2 dt) = A U^n + F(t^n)
   */
  template<class T>
  class OptimalModifiedEquation_Iterator  : public VirtualTimeScheme<T>
  {
  protected :
    int order, nb_additional_terms, nb_rk_step_first, order_rk_first;
    Real_wp coef_xsi, xsi, alpha1, alpha2;
    Real_wp dt;
    bool fourth_order;
    Vector<T> Unm1, Y, Q0, Q1, Qkp1;
    bool unstable_algo;
    Vector<Real_wp> alpha;
    
  public :
    OptimalModifiedEquation_Iterator();
    
    void SetOrder(int r, int m = 0, bool horner = false, int nb_step = 10, int order_rk = 4);
    
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
  
}    

#define MONTJOIE_FILE_MODIFIED_EQUATION_HXX
#endif
