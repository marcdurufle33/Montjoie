#ifndef MONTJOIE_FILE_MODIFIED_EQUATION_CXX

namespace Montjoie
{
  
  //! running a second-order time scheme to solve d^2 u /dt = f(u, du/dt)
  /*!
    \param[in] t_begin inital time t0
    \param[in] t_end final time tf
    \param[in] deltat time step
    \param[inout] U0 solution u at t = 0, on output final solution u(tf)
    \param[inout] dU0_dt solution du/dt at t = 0, on output final solution du/dt(tf)
    \param[inout] sys object containing the definition of f
    \param[inout] scheme time scheme to use
   */
  template<class T>
  void RunSecondOrderScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
			    Vector<T>& U0, Vector<T>& dU0_dt, VirtualOdeSystem<T>& sys,
			    VirtualTimeScheme<T>& scheme)
  {
    Real_wp dt, t; 
    int nb_max_iter;
    if (deltat != Real_wp(0))
      nb_max_iter = toInteger(ceil(abs(t_end - t_begin)/deltat ) );
    else
      return;
    
    // if t_end <= t_begin, we exit 
    if (nb_max_iter < 1)
      return;
    
    dt = deltat; t = t_begin;
    sys.GiveNumberIterations(dt, nb_max_iter);

    // first we initialize the scheme
    scheme.SetInitialConditionS(t, deltat, U0, dU0_dt, sys);
       
    // iterations in time
    // main loop
    for (int n_time = 0; n_time < nb_max_iter; n_time++)
      {
	t = t_begin + n_time*dt;
        // the iterate is provided 
	sys.GiveIterate(n_time, t, scheme.GetIterate());
	
        scheme.Advance(t, n_time, sys);
      }
    
    // end of the scheme
    t = t_begin + nb_max_iter*dt;
    sys.GiveFinalIterate(nb_max_iter, t, scheme.GetIterate());
    scheme.ClearSecond(U0, dU0_dt);
  }


  /****************************
   * ModifiedEquationIterator *
   ****************************/


  //! default constructor
  template<class T>
  ModifiedEquationIterator<T>::ModifiedEquationIterator()
  {
    dt = 0; dt2 = 0;
    order = 2;
  }
  
  
  //! sets the order
  template<class T>
  void ModifiedEquationIterator<T>::SetOrder(int r)
  {
    order = r;
    if ((order != 2) && (order != 4))
      {
        cout << "The order " << r << " is not imlemented " << endl;
        abort();
      }
  }
  

  template<class T>
  int ModifiedEquationIterator<T>::GetOrder() const
  {
    return order;
  }
  
  
  //! initialisationof time schemes with the initial conditions u(0) and du/dt(0)
  template<class T>
  void ModifiedEquationIterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                         Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys)
  {
    T zero, one;  dt = deltat; dt2 = dt*dt;
    SetComplexZero(zero); SetComplexOne(one);
    int N = u0.GetM();
    if (order == 2 )
      {
	U0.Reallocate(N);
	U1.Reallocate(N); Ah_u1.Reallocate(N); Ah_u1.Fill(zero);
	Copy(u0, U0); Copy(U0, U1); u0.Clear();
	
        sys.SetDirichletCondition(t, 0, U0);
        
	// first step, we use a one-step algorithm to compute U^1 from initial condition U^0
	// (we assume that potential sources are null near t = initial_time)
	// basic second order initialization :
	// U^1 = U^0 + dt dU/dt(0) + dt^2/2 d^2 U/dt^2(0)
	// computing D_h d^2 U /dt^2 = Kh U - Sh dU/dt
	sys.ApplyOperatorSh(Real_wp(1), t, du0_dt, Real_wp(0), U1);
	sys.ApplyOperatorKh(Real_wp(1), t, U0, -Real_wp(1), U1);
	sys.SolveOperatorDh(U1);
	// then forming U1
	Mlt(T(0.5*dt), U1); Add(one, du0_dt, U1);
	Mlt(dt, U1); Add(one, U0, U1);

        // non-homogeneous Dirichlet condition
        sys.SetDirichletCondition(t+dt, 0, U1);
	
	Ah_u1.Fill(zero); dU0_dt.Clear();

      }
    else if (order == 4)
      {
	U0.Reallocate(N); dU0_dt.Reallocate(N);
	Copy(u0, U0); Copy(du0_dt, dU0_dt); u0.Clear(); du0_dt.Clear();
	U1.Reallocate(N); U1.Fill(0); U2.Reallocate(N); U2.Fill(0);
	Copy(U0, U1);
	dU_dt2.Reallocate(N); dU_dt2.Fill(0);
	dU_dt3.Reallocate(N); dU_dt3.Fill(0);
	dU_dt4.Reallocate(N); dU_dt4.Fill(0);
	
        sys.SetDirichletCondition(t-dt, 0, U0);
            
	// first step, we use a one-step algorithm to compute U^0 and U^2 from initial condition U^1
	// (we suppose that sources are null near t = initial_time)
	// basic fourth order initialization :
	// U^2 = U^1 + dt dU/dt(0) + dt^2/2 d^2 U/dt^2(0) 
	//      + dt^3/6 d^3 U/dt^3(0) + dt^4/24 d^4 U/dt^4
	// U^0 = U^1 - dt dU/dt(0) + dt^2/2 d^2 U/dt^2(0) 
	//      - dt^3/6 d^3 U/dt^3(0) + dt^4/24 d^4 U/dt^4
	// computing d^2 U /dt^2 = Kh U - Sh dU/dt
	sys.ApplyOperatorSh(Real_wp(1), t, dU0_dt, Real_wp(0), U2);
	sys.ApplyOperatorKh(Real_wp(1), t, U0, -Real_wp(1), U2);
	Copy(U2, dU_dt2); sys.SolveOperatorDh(dU_dt2);
	sys.SetDirichletCondition(t, 2, dU_dt2);
        
	sys.ApplyOperatorSh(Real_wp(1), t, dU_dt2, Real_wp(0), U2);
	sys.ApplyOperatorKh(Real_wp(1), t, dU0_dt, -Real_wp(1), U2);
	Copy(U2, dU_dt3); sys.SolveOperatorDh(dU_dt3);
	sys.SetDirichletCondition(t, 3, dU_dt3);
	
	sys.ApplyOperatorSh(Real_wp(1), t, dU_dt3, Real_wp(0), U2);
	sys.ApplyOperatorKh(Real_wp(1), t, dU_dt2, -Real_wp(1), U2);
	Copy(U2, dU_dt4); sys.SolveOperatorDh(dU_dt4);
	sys.SetDirichletCondition(t, 4, dU_dt4);
	
	// then forming U0 and U2
	Mlt(one/24, dU_dt4); Mlt(one/6, dU_dt3); Mlt(one/2, dU_dt2);
	Copy(dU_dt4, U2); Mlt(dt, U2); Add(one, dU_dt3, U2); Mlt(dt, U2);
	Add(one, dU_dt2, U2); Mlt(dt, U2); Add(one, dU0_dt, U2);
	Mlt(dt, U2); Add(one, U1, U2);
	
	Copy(dU_dt4, U0); Mlt(dt, U0); Add(-one, dU_dt3, U0); Mlt(dt, U0);
	Add(one, dU_dt2, U0); Mlt(dt, U0); Add(-one, dU0_dt, U0);
	Mlt(dt, U0); Add(one, U1, U0);
	
	sys.SetDirichletCondition(t, 0, U1);
	sys.SetDirichletCondition(t+dt, 0, U2);
	
	Prod_Uh.Reallocate(N); Prod_Uh.Fill(0);
      }
  }
	
  
  //! vectors needed to advance the scheme are cleared
  template<class T>
  void ModifiedEquationIterator<T>::Clear()
  {
    U0.Clear(); U1.Clear(); Ah_u1.Clear();
    U2.Clear(); Prod_Uh.Clear(); dU0_dt.Clear(); dU_dt2.Clear(); dU_dt3.Clear(); dU_dt4.Clear();
  }
  
  
  //! vectors needed to advance the scheme are cleared
  template<class T>
  void ModifiedEquationIterator<T>::ClearSecond(Vector<T>& Ufinal, Vector<T>& dUfinal_dt)
  {
    Ah_u1.Clear(); U2.Clear();
    Prod_Uh.Clear(); dU0_dt.Clear(); dU_dt2.Clear(); dU_dt3.Clear(); dU_dt4.Clear();
    
    dUfinal_dt.Clear();
    if (order == 2)
      Ufinal = U0;
    else
      Ufinal = U1;

    U0.Clear(); U1.Clear();
  }
  
  
  //! returns stability function of the scheme
  template<class T>
  UnivariatePolynomial<Real_wp> ModifiedEquationIterator<T>::GetStabilityFunction() const
  {
    UnivariatePolynomial<Real_wp> Pol;
    Pol.SetOrder(0);
    Pol(0) = GetCFL(true);
    return Pol;
  }


  //! returns the CFL of the scheme
  template<class T>
  Real_wp ModifiedEquationIterator<T>::GetCFL(bool real_p) const
  {
    if (order == 2)
      return 2.0;
    else if (order == 4)
      return 2.0*sqrt(Real_wp(3));
    
    return Real_wp(0);
  }
  

  //! returns u(t^n)
  template<class T>
  Vector<T>& ModifiedEquationIterator<T>::GetIterate()
  {
    if (order == 2)
      return U0;
    else
      return U1;
  }


  //! returns u(t^n)
  template<class T>
  const Vector<T>& ModifiedEquationIterator<T>::GetIterate() const
  {
    if (order == 2)
      return U0;
    else
      return U1;
  }
  
  
  //! changes time step during the simulation
  template<class T>
  void ModifiedEquationIterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int ModifiedEquationIterator<T>::GetNumberOfIterates() const
  {
    if (order == 4)
      return 3;
    
    return 2;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& ModifiedEquationIterator<T>::GetIterate(int k)
  {
    switch (k)
      {
      case 0 : return U0;
      case 1 : return U1;        
      case 2 : return U2;
      }
    
    return U0;
  }

  
  //! computes u^{n+1} from u^{n-1} and u^n, and overwrites old vectors
  template<class T>
  void ModifiedEquationIterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (order == 2)
      {
	///////////////////////////
	// ORDER TWO ON EQUATION //
	///////////////////////////
	
	// we are solving the following equation :
	// D_h d^2 U /dt^2  -  Kh U + Sh dU/dt = F
	
	// scheme of the type
	// Dh ( U^{n+1} - 2 U^n + U^{n-1} ) / dt^2 - K_h U^n 
        // + \Sigma_h (U^{n+1} - U^{n-1}) / (2 dt) = F^n
	
	// now we compute Ah_u1 = F^n + K_h U^n
	//U0.WriteText("Un_prev.dat");	U1.WriteText("Un.dat");
	//glob_chrono.Start(VirtualTimer::STIFFNESS);
        sys.ApplyOperatorKh(Real_wp(1), t+dt, U1, Real_wp(0), Ah_u1);
        //glob_chrono.Stop(VirtualTimer::STIFFNESS);

	// we add the source term Fh(n)
	//glob_chrono.Start(VirtualTimer::OUTPUT);
	sys.AddScalarTimeSource(Real_wp(1), t+dt, 1, Ah_u1);
	//glob_chrono.Stop(VirtualTimer::OUTPUT);
	
	// we compute 2 Dh U^n + dt^2 (F^n + K_h U^n) -  (Dh - dt Sh/2) U^{n-1}
	//glob_chrono.Start(VirtualTimer::MASS);
	sys.ApplyOperatorDhMinusdtSh(-Real_wp(1), t+dt, U0, dt2, Ah_u1);
	sys.ApplyOperatorDh(Real_wp(2), t+dt, U1, Real_wp(1), Ah_u1);
	//glob_chrono.Stop(VirtualTimer::MASS);

	// now solving (Dh + dt/2 Sh) U^{n+1} = Ah_u1
	//glob_chrono.Start(VirtualTimer::FLUX);
	Copy(U1, U0);
        Copy(Ah_u1, U1);
	sys.SetDirichletConditionSource(t+2*dt, 0, U1);
	//glob_chrono.Stop(VirtualTimer::FLUX);

	//glob_chrono.Start(VirtualTimer::SOLVE);
	sys.SolveOperatorDhPlusdtSh(U1); //U1.WriteText("Un_next.dat");
        //glob_chrono.Stop(VirtualTimer::SOLVE);

        /*	if (n_time%100 == 0)
	  {
	    glob_chrono.DisplayTime(VirtualTimer::STIFFNESS, " perform the product with stiffness matrix");
	    glob_chrono.DisplayTime(VirtualTimer::OUTPUT, " add the source terms");
	    glob_chrono.DisplayTime(VirtualTimer::MASS, " perform the product with mass matrix");
	    glob_chrono.DisplayTime(VirtualTimer::SOLVE, " solve the mass matrix");
	    glob_chrono.DisplayTime(VirtualTimer::FLUX, " copy vectors and Dirichlet");
	  }
        */
      }
    else if (order == 4)
      {

	////////////////////////////////////////////////
	// MODIFIED EQUATION ON SECOND ORDER EQUATION //
	////////////////////////////////////////////////
		
	// we solve  d^2 U/dt2 + sigma_h du/dt - Ah U = F
	// modified equation approach :
	// (D_h + dt/2 Sh) U^{n+1} = 2 D_h U^n - (D_h - dt/2 Sh) U^{n-1}  + dt^2 ( F^n + Kh U^n)
	//                           + dt^4/12  Dh d^4 U/dt^4  + dt^4/6  Sh d^3 U/dt^3
	// with
	// dU/dt = (3 U^n - 4 U^{n-1} + U^{n-2} )/ (2 dt)  (uncentered stencil)
	// Dh d^2 U/ dt^2 =  ( F^n + Kh U^n) - Sh dU/dt
	// Dh d^3 U/ dt^3 =  ( dF^n / dt + Kh dU/dt) - Sh d^2 U/dt^2
	// Dh d^4 U/ dt^4 =  ( d^2 F^n / dt^2 + Kh d^2 U/dt^2) - Sh d^3 U/dt^3
	
	// now we compute Prod_Uh = delta t^2 (Fh_n + Ah u_n)
	sys.ApplyOperatorKh(Real_wp(1), t+dt, U2, Real_wp(0), Prod_Uh);
	// we add the source term Fh(n)
	sys.AddScalarTimeSource(Real_wp(1), t+dt, 1, Prod_Uh);
	
	Mlt(dt2, Prod_Uh);
	
	// computing deltat^2 dU/dt with the uncentered stencil 
        // dU/dt = (3 U^n - 4 U^{n-1} + U^{n-2} )/ (2 dt) 
	Copy(U0, dU0_dt); Add(T(-4), U1, dU0_dt); Add(T(3), U2, dU0_dt);
	Mlt(T(0.5*dt), dU0_dt);
	
	// computing deltat^2  d^2 U/dt^2
	Copy(Prod_Uh, U0);
	sys.ApplyOperatorSh(-Real_wp(1), t+dt, dU0_dt, Real_wp(1), U0);
	Copy(U0, dU_dt2); sys.SolveOperatorDh(dU_dt2);
	sys.SetDirichletCondition(t, 2, dU_dt2, dt2);
	
	// computing deltat^2 d^3 U/dt^3 (for this term, we could use Kh 
        // restricted to absorbing media)
	sys.ApplyOperatorKh(Real_wp(1), t+dt, dU0_dt, Real_wp(0), U0);
	sys.ApplyOperatorSh(-Real_wp(1), t+dt, dU_dt2, Real_wp(1), U0);
	sys.AddScalarTimeSource(dt2, t+dt, 2, U0);
	
	Copy(U0, dU_dt3); sys.SolveOperatorDh(dU_dt3);
	sys.SetDirichletCondition(t+dt, 3, dU_dt3, dt2);
	
	// computing deltat^2 Dh d^4 U/dt^4
	sys.ApplyOperatorKh(Real_wp(1), t+dt, dU_dt2, Real_wp(0), dU_dt4);
	sys.ApplyOperatorSh(-Real_wp(1), t+dt, dU_dt3, Real_wp(1), dU_dt4);
	sys.AddScalarTimeSource(dt2, t+dt, 3, dU_dt4);
	
	// now forming the right hand side of the equation
	// (D_h + dt/2 Sh) U^{n+1} = Prod_Uh
	Add(dt2/12, dU_dt4, Prod_Uh);
	sys.ApplyOperatorSh(dt2/6, t+dt, dU_dt3, Real_wp(1), Prod_Uh);
	
	sys.ApplyOperatorDh(Real_wp(2), t+dt, U2, Real_wp(1), Prod_Uh);
	sys.ApplyOperatorDhMinusdtSh(-Real_wp(1), t+dt, U1, Real_wp(1), Prod_Uh);
	
	// solving the final equation
	Copy(U1, U0); Copy(U2, U1);
	Copy(Prod_Uh, U2); sys.SolveOperatorDhPlusdtSh(U2);

	sys.SetDirichletCondition(t+2*dt, 0, U2);
      } // end if order == 4
  }
  

  //! returns the size used by the object in bytes
  template<class T>
  size_t ModifiedEquationIterator<T>::GetMemorySize() const
  {
    size_t taille = U0.GetMemorySize() + U1.GetMemorySize() + Ah_u1.GetMemorySize() + U2.GetMemorySize();
    taille += Prod_Uh.GetMemorySize() + dU0_dt.GetMemorySize() + dU_dt2.GetMemorySize()
      + dU_dt3.GetMemorySize() + dU_dt4.GetMemorySize();

    return taille;
  }

  
  /********************************
   * Function RunFirstOrderScheme *
   ********************************/


  //! Solves an evolution system du/dt = g(t, u, v)  dv/dt = h(t, u, v)
  /*!
    \param[in] t_begin initial time t0
    \param[in] t_end final time T
    \param[in] deltat time step
    \param[inout] U0 solution u for t = t0
    \param[inout] V0 solution v for t = t0
    \param[in] sys description of the function g(t, u)
    \param[in] scheme time scheme used    
    g(t, u, v)  and h(t, u, v) are given in object sys
   */  
  template<class T>
  void RunFirstOrderScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
			   Vector<T>& U0, Vector<T>& V0, VirtualOdeSystem<T>& sys,
			   VirtualTimeScheme<T>& scheme)
  {
    Real_wp dt, t, tn_demi; 
    int nb_max_iter;
    if (deltat != Real_wp(0))
      nb_max_iter = toInteger(ceil(abs(t_end - t_begin)/deltat ) );
    else
      return;
    
    // if t_end <= t_begin, we exit 
    if (nb_max_iter < 1)
      return;
    
    dt = deltat; t = t_begin;
    sys.GiveNumberIterations(dt, nb_max_iter);

    // first we initialize the scheme
    scheme.SetInitialConditionS(t, deltat, U0, V0, sys);
    
    // iterations in time
    // main loop
    for (int n_time = 0; n_time < nb_max_iter; n_time++)
      {
	t = t_begin + n_time*dt;
        tn_demi = t + 0.5*dt;
        
	// the iterate is provided 
        sys.GiveIterate(n_time, t, scheme.GetIterate());
        sys.GiveVectorialIterate(n_time, tn_demi, scheme.GetVectorialIterate());
        
        scheme.Advance(t, n_time, sys);
      }
    
    // end of the scheme
    t = t_begin + nb_max_iter*dt;
    tn_demi = t + 0.5*dt;
    sys.GiveFinalIterate(nb_max_iter, t, scheme.GetIterate());
    sys.GiveVectorialIterate(nb_max_iter, tn_demi, scheme.GetVectorialIterate());
    
    scheme.ClearSecond(U0, V0);
  }


  /**********************************
   * ModifiedEquationSystemIterator *
   **********************************/  
  
  
  //! default constructor
  template<class T>
  ModifiedEquationSystemIterator<T>::ModifiedEquationSystemIterator()
  {
    dt = 0; dt2 = 0;
    order = 2;
  }
  
  
  //! sets the order of time-scheme
  template<class T>
  void ModifiedEquationSystemIterator<T>::SetOrder(int r)
  {
    order = r;
    if ((order != 2) && (order != 4))
      {
        cout << "The order " << r << " is not imlemented " << endl;
        abort();
      }
  }

  
  //! initialisation of the scheme with u0 and v0
  template<class T>
  void ModifiedEquationSystemIterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                         Vector<T>& v0, VirtualOdeSystem<T>& sys)
  {
    T zero; SetComplexZero(zero);
    dt = deltat; dt2 = dt*dt;
    int N = u0.GetM(), Nvec = v0.GetM();
    if (order == 2 )
      {
	Un.Reallocate(N); Vn.Reallocate(Nvec);
	Copy(u0, Un); Copy(v0, Vn); u0.Clear();
	Prod_Uh.Reallocate(N); Prod_Vh.Reallocate(Nvec);	
	Prod_Uh.Fill(zero); Prod_Vh.Fill(zero);
	
	// initialization of V^{1/2} from V^0 
	// basic second order initialization
	// V^{1/2} = V^0 + dt/2 dV/dt + dt^2/8 d^2 V/dt^2
	
	// Prod_Vh = dV/dt
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, Un, Real_wp(0), Prod_Vh);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, Vn, Real_wp(1), Prod_Vh);
	sys.SolveOperatorBh(Prod_Vh);
        
	// Prod_Uh = dU/dt
	sys.ApplyOperatorRhScalar(Real_wp(1), t, Vn, Real_wp(0), Prod_Uh);
	sys.ApplyOperatorSh(-Real_wp(1), t, Un, Real_wp(1), Prod_Uh);  
	sys.SolveOperatorDh(Prod_Uh);
        
	// v0 = d^2 V/dt^2
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, Prod_Uh, Real_wp(0), v0);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, Prod_Vh, Real_wp(1), v0);
	sys.SolveOperatorBh(v0);
        
	Add(dt/2, Prod_Vh, Vn);
	Add(dt2/8, v0, Vn);
	
        v0.Clear();
        
        // non-homogeneous Dirichlet condition
        sys.SetDirichletCondition(t, 0, Un);

      }
    else if (order == 4)
      {
	Un.Reallocate(N); Vn.Reallocate(v0.GetM());
	Copy(u0, Un); Copy(v0, Vn); u0.Clear(); v0.Clear();
	
	Prod_Uh.Reallocate(N); Prod_Vh.Reallocate(Nvec); //dU0_dt.Reallocate(N);
	Prod_Uh.Fill(zero); Prod_Vh.Fill(zero);
        dU_dt.Reallocate(N); dU_dt2.Reallocate(N); dU_dt3.Reallocate(N);
        dV_dt.Reallocate(Nvec); dV_dt2.Reallocate(Nvec); dV_dt3.Reallocate(Nvec);
        Un_m1.Reallocate(N); Vn_m05.Reallocate(Nvec);
        Un_demi.Reallocate(N); Vn_p1.Reallocate(Nvec);
        dU_dt.Fill(zero); dU_dt2.Fill(zero); dU_dt3.Fill(zero);
        dV_dt.Fill(zero); dV_dt2.Fill(zero); dV_dt3.Fill(zero);
        Un_m1.Fill(zero); Vn_m05.Fill(zero); Un_demi.Fill(zero); Vn_p1.Fill(zero);

	// initialization of V^{1/2} from V^0 
	// basic fourth order initialization
	// V^{1/2} = V^0 + dt/2 dV/dt + dt^2/8 d^2 V/dt^2 
	//         + dt^3/48 d^3 V/dt^3 + dt^4/384 d^4 V/dt^4
	// similar stuff for V^{-1/2} and U^{-1}
	
	// dV/dt
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, Un, Real_wp(0), dV_dt);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, Vn, Real_wp(1), dV_dt);
	sys.SolveOperatorBh(dV_dt);
        
	// dU/dt
	sys.ApplyOperatorRhScalar(Real_wp(1), t, Vn, Real_wp(0), dU_dt);
	sys.ApplyOperatorSh(-Real_wp(1), t, Un, Real_wp(1), dU_dt);  
	sys.SolveOperatorDh(dU_dt);
        
	// d^2 V/dt^2
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, dU_dt, Real_wp(0), dV_dt2);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, dV_dt, Real_wp(1), dV_dt2);
	sys.SolveOperatorBh(dV_dt2);
	
	// d^2 U/dt^2
	sys.ApplyOperatorRhScalar(Real_wp(1), t, dV_dt, Real_wp(0), dU_dt2);
	sys.ApplyOperatorSh(-Real_wp(1), t, dU_dt, Real_wp(1), dU_dt2);  
	sys.SolveOperatorDh(dU_dt2);
	
	// d^3 V/dt^3
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, dU_dt2, Real_wp(0), dV_dt3);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, dV_dt2, Real_wp(1), dV_dt3);
	sys.SolveOperatorBh(dV_dt3);
	
	// d^3 U/dt^3
	sys.ApplyOperatorRhScalar(Real_wp(1), t, dV_dt2, Real_wp(0), dU_dt3);
	sys.ApplyOperatorSh(-Real_wp(1), t, dU_dt2, Real_wp(1), dU_dt3);  
	sys.SolveOperatorDh(dU_dt3);
	
	// Prod_Vh = d^4 V/dt^4
	sys.ApplyOperatorRhVectorial(Real_wp(1), t, dU_dt3, Real_wp(0), Prod_Vh);
	sys.ApplyOperatorShVectorial(-Real_wp(1), t, dV_dt3, Real_wp(1), Prod_Vh);
	sys.SolveOperatorBh(Prod_Vh);
	
	// Prod_Uh = d^4 U/dt^4
	sys.ApplyOperatorRhScalar(Real_wp(1), t, dV_dt3, Real_wp(0), Prod_Uh);
	sys.ApplyOperatorSh(-Real_wp(1), t, dU_dt3, Real_wp(1), Prod_Uh);  
	sys.SolveOperatorDh(Prod_Uh);
        
	// forming V(t^1/2), V(t^{-1/2})
	Real_wp coef1 = dt/2, coef2 = dt2/8, coef3 = coef2*dt/6, coef4 = coef3*dt/8;
	for (int i = 0; i < Vn.GetM(); i++)
	  {
	    Vn_m05(i) = Vn(i) - coef1*dV_dt(i) + coef2*dV_dt2(i) 
	      - coef3*dV_dt3(i) + coef4*Prod_Vh(i);
	    Vn(i) += coef1*dV_dt(i) + coef2*dV_dt2(i) + coef3*dV_dt3(i) + coef4*Prod_Vh(i);
	  }
	
	// and U^{-1}
	coef1 = dt; coef2 = dt2/2; coef3 = coef2*dt/3; coef4 = coef3*dt/4;
	for (int i = 0; i < Un.GetM(); i++)
	  Un_m1(i) = Un(i) - coef1*dU_dt(i) + coef2*dU_dt2(i) 
	    - coef3*dU_dt3(i) + coef4*Prod_Uh(i);
	
        sys.SetDirichletCondition(t-dt, 0, Un_m1);
        sys.SetDirichletCondition(t, 0, Un);
      }
  }
  
  
  //! clears vectors stored in the scheme
  template<class T>
  void ModifiedEquationSystemIterator<T>::Clear()
  {
    Un.Clear(); Vn.Clear(); Prod_Vh.Clear();
    Prod_Uh.Clear();
    dU_dt.Clear(); dU_dt2.Clear(); dU_dt3.Clear();
    dV_dt.Clear(); dV_dt2.Clear(); dV_dt3.Clear();
    Un_m1.Clear(); Vn_m05.Clear(); Vn_p1.Clear(); Un_demi.Clear();
  }


  //! clears vectors stored in the scheme
  template<class T>
  void ModifiedEquationSystemIterator<T>::ClearSecond(Vector<T>& U, Vector<T>& V)
  {
    Prod_Vh.Clear(); Prod_Uh.Clear();
    dU_dt.Clear(); dU_dt2.Clear(); dU_dt3.Clear();
    dV_dt.Clear(); dV_dt2.Clear(); dV_dt3.Clear();
    Un_m1.Clear(); Vn_m05.Clear(); Vn_p1.Clear(); Un_demi.Clear();
    
    U = Un; V = Vn;
    Un.Clear(); Vn.Clear(); 
  }


  //! returns stability function of the scheme
  template<class T>
  UnivariatePolynomial<Real_wp> ModifiedEquationSystemIterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the CFL of the scheme
  template<class T>
  Real_wp ModifiedEquationSystemIterator<T>::GetCFL(bool real_p) const
  {
    if (order == 2)
      return 2.0;
    else if (order == 4)
      return 4.0*sqrt(Real_wp(2));

    return Real_wp(0);
  }


  //! returns u(t^n)
  template<class T>
  Vector<T>& ModifiedEquationSystemIterator<T>::GetIterate()
  {
    return Un;
  }


  //! returns u(t^n)
  template<class T>
  const Vector<T>& ModifiedEquationSystemIterator<T>::GetIterate() const
  {
    return Un;
  }


  //! returns v(t^n)
  template<class T>
  Vector<T>& ModifiedEquationSystemIterator<T>::GetVectorialIterate()
  {
    return Vn;
  }
  

  //! changes time step during the simulation
  template<class T>
  void ModifiedEquationSystemIterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int ModifiedEquationSystemIterator<T>::GetNumberOfIterates() const
  {
    if (order == 4)
      return 2;
    
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& ModifiedEquationSystemIterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Un;
    else
      return Un_m1;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int ModifiedEquationSystemIterator<T>::GetNumberOfVectorialIterates() const
  {
    if (order == 4)
      return 2;
    
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& ModifiedEquationSystemIterator<T>::GetVectorialIterate(int k)
  {
    if (k == 0)
      return Vn;
    else
      return Vn_m05;
  }

  
  //! computation of u(t^n+1) and v(t^(n+3/2))
  template<class T>
  void ModifiedEquationSystemIterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    Real_wp one(1), zero(0);
    if (order == 2)
      {
	// we are solving the following system
	// D_h d U /dt  + Sh U =  F_U + Rh^S V
	// B_h d V /dt  + ShV V = F_V +  Rh^V U 
        
        // with the following scheme
	// (Dh + dt/2 Sh) U^{n+1} = (Dh - dt/2 Sh) U^n + dt (F_U^{n+1/2} + Rh^S V^{n+1/2} )
	// (Bh + dt/2 Sh^V) V^{n+3/2} = (Bh - dt/2 Sh^V) V^{n+1/2} + dt (F_V^{n+1} + Rh^V U^{n+1} )
	Real_wp tn_demi = t + 0.5*dt;
	
	// scalar equation in U
        
        //string suffix = to_str(rank_processor) + "N" + to_str(nb_processors) + ".dat";
        //sys.var_harmonic.mesh_num.GlobDofNumber_Subdomain.Write("numGlob"+suffix);
        //Un.Write("Un"+suffix); Vn.Write("Vn"+suffix);
	
	sys.ApplyOperatorRhScalar(one, tn_demi, Vn, zero, Prod_Uh);
        sys.AddScalarTimeSource(one, tn_demi, 0, Prod_Uh);
        sys.ApplyOperatorDhMinusdtSh(one, tn_demi, Un, dt, Prod_Uh);
        
        // overwriting U^n by U^{n+1}
	Un.Fill(0);
        sys.SolveOperatorDhPlusdtSh(1.0, Prod_Uh, Un);
        sys.SetDirichletCondition(t+dt, 0, Un);
        //Un.Write("UnNext"+suffix);
        
	// vectorial equation in V
	sys.ApplyOperatorRhVectorial(dt, t+dt, Un, zero, Prod_Vh);
        sys.AddVectorialTimeSource(dt, t+dt, 0, Prod_Vh);
        
	// we compute V^{n+3/2}
	sys.ApplyOperatorBhMinusdtSh(one, t+dt, Vn, one, Prod_Vh);
        
	// overwriting V^{n+1/2} by V^{n+3/2}
	sys.SolveOperatorBhPlusdtSh(Prod_Vh);
        Copy(Prod_Vh, Vn);
        //Vn.Write("VnNext"+suffix);
        //int test_input; cout << "we wait" << endl; cin >> test_input;
      }
    else if (order == 4)
      {
	// (Dh + dt/2 Sh) U^{n+1} = (Dh - dt/2 Sh) U^n + dt (F_U^{n+1/2} + Rh^S V^{n+1/2} )  
        //                           + dt^3/24 Dh d^3 U /dt^3  + dt^3/8  Sh d^2 U/dt^2   
	// (Bh + dt/2 Sh^V) V^{n+3/2} = (Bh - dt/2 Sh^V) V^{n+1/2} + dt (F_V^{n+1} + Rh^V U^{n+1} )
        //                              + dt^3/24 Bh d^3 V /dt^3  + dt^3/8  Sh^V d^2 V/dt^2   
	Real_wp tn_demi = t + 0.5*dt;

	// scalar equation in U
	sys.ApplyOperatorRhScalar(one, tn_demi, Vn, zero, Prod_Uh);
	sys.AddScalarTimeSource(one, tn_demi, 0, Prod_Uh);
        
        // evaluation of U^n+1/2 and V^n+1
        Real_wp half = Real_wp(1)/2, three_half = Real_wp(3)/2;
        for (int i = 0; i < Un.GetM(); i++)
          Un_demi(i) = three_half*Un(i) - half*Un_m1(i);
        
        for (int i = 0; i < Vn.GetM(); i++)
          Vn_p1(i) = three_half*Vn(i) - half*Vn_m05(i);
        
	// computing dU/dt
	Copy(Prod_Uh, dU_dt);
        sys.ApplyOperatorSh(-one, tn_demi, Un_demi, one, dU_dt);
        sys.SolveOperatorDh(dU_dt);
	sys.SetDirichletCondition(t, 1, dU_dt);
	
        // then dV/dt
        sys.ApplyOperatorRhVectorial(one, tn_demi, Un_demi, zero, dV_dt);
        sys.AddVectorialTimeSource(one, tn_demi, 0, dV_dt);
        sys.ApplyOperatorShVectorial(-one, tn_demi, Vn, one, dV_dt);
        sys.SolveOperatorBh(dV_dt);
        	
        // then d^2 U/dt^2
        sys.ApplyOperatorRhScalar(one, tn_demi, dV_dt, zero, dU_dt2);
        sys.ApplyOperatorSh(-one, tn_demi, dU_dt, one, dU_dt2);
        sys.AddScalarTimeSource(one, tn_demi, 1, dU_dt2);
        sys.SolveOperatorDh(dU_dt2);
	sys.SetDirichletCondition(t, 2, dU_dt2);
        	
        // d^2 V/dt^2
        sys.ApplyOperatorRhVectorial(one, tn_demi, dU_dt, zero, dV_dt2);
        sys.AddVectorialTimeSource(one, tn_demi, 1, dV_dt2);
        sys.ApplyOperatorShVectorial(-one, tn_demi, dV_dt, one, dV_dt2);
        sys.SolveOperatorBh(dV_dt2);
		
        // and finally Dh d^3 U/dt^3
        sys.ApplyOperatorRhScalar(one, tn_demi, dV_dt2, zero, dU_dt3);
        sys.ApplyOperatorSh(-one, tn_demi, dU_dt2, one, dU_dt3);
        sys.AddScalarTimeSource(one, tn_demi, 2, dU_dt3);

        // adding terms in Prod_Uh 
        Add(dt2/24, dU_dt3, Prod_Uh);
        sys.ApplyOperatorSh(dt2/8, tn_demi, dU_dt2, one, Prod_Uh);

	// overwriting U^n by U^{n+1}
	sys.ApplyOperatorDhMinusdtSh(one, tn_demi, Un, dt, Prod_Uh);
	Copy(Un, Un_m1);
        Un.Fill(0); sys.SolveOperatorDhPlusdtSh(1.0, Prod_Uh, Un);
	sys.SetDirichletCondition(t+dt, 0, Un);
		  
        // vectorial equation in V
	sys.ApplyOperatorRhVectorial(one, t+dt, Un, zero, Prod_Vh);
 	sys.AddVectorialTimeSource(one, t+dt, 0, Prod_Vh);
        
        sys.ApplyOperatorRhScalar(one, t+dt,Vn_p1, zero, dU_dt);
        sys.AddScalarTimeSource(one, t+dt, 0, dU_dt);
        sys.ApplyOperatorSh(-one, t+dt, Un, one, dU_dt);
        sys.SolveOperatorDh(dU_dt);
	sys.SetDirichletCondition(t, 1, dU_dt);
	
        // evaluating dV_dt
        Copy(Prod_Vh, dV_dt);
        sys.ApplyOperatorShVectorial(-one, t+dt, Vn_p1, one, dV_dt);
        sys.SolveOperatorBh(dV_dt);
		
        // evaluating d^2 U/dt^2
        sys.ApplyOperatorRhScalar(one, t+dt, dV_dt, zero, dU_dt2);
        sys.ApplyOperatorSh(-one, t+dt, dU_dt, one, dU_dt2);
        sys.AddScalarTimeSource(one, t+dt, 1, dU_dt2);
        sys.SolveOperatorDh(dU_dt2);
	sys.SetDirichletCondition(t, 2, dU_dt2);
	
        // evaluating d^2 V/dt^2
        sys.ApplyOperatorRhVectorial(one, t+dt, dU_dt, zero, dV_dt2);
        sys.AddVectorialTimeSource(one, t+dt, 1, dV_dt2);
        sys.ApplyOperatorShVectorial(-one, t+dt, dV_dt, one, dV_dt2);
        sys.SolveOperatorBh(dV_dt2);
	
        // and finally Bh d^3 V/dt^3
        sys.ApplyOperatorRhVectorial(one, t+dt, dU_dt2, zero, dV_dt3);
        sys.AddVectorialTimeSource(one, t+dt, 2, dV_dt3);
        sys.ApplyOperatorShVectorial(-one, t+dt, dV_dt2, one, dV_dt3);        
        
        // adding terms to Prod_Vh
        Add(dt2/24, dV_dt3, Prod_Vh);
        sys.ApplyOperatorShVectorial(dt2/8, t+dt, dV_dt2, one, Prod_Vh);
        
	// overwriting V^{n+1/2} by V^{n+3/2}        
	sys.ApplyOperatorBhMinusdtSh(one, t+dt, Vn, dt, Prod_Vh);
	sys.SolveOperatorBhPlusdtSh(Prod_Vh);
        Copy(Vn, Vn_m05);
        Copy(Prod_Vh, Vn);
      }
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t ModifiedEquationSystemIterator<T>::GetMemorySize() const
  {
    size_t taille = Un.GetMemorySize() + Vn.GetMemorySize() + Prod_Uh.GetMemorySize() + Prod_Vh.GetMemorySize();
    taille += dU_dt.GetMemorySize() + dU_dt2.GetMemorySize() + dU_dt3.GetMemorySize()
      + dV_dt.GetMemorySize() + dV_dt2.GetMemorySize() + dV_dt3.GetMemorySize()
      + Un_m1.GetMemorySize() + Vn_m05.GetMemorySize() + Vn_p1.GetMemorySize() + Un_demi.GetMemorySize();
    
    return taille;
  }


  /************************************
   * OptimalModifiedEquation_Iterator *
   ************************************/
  
  //! default constructor
  template<class T>
  OptimalModifiedEquation_Iterator<T>::OptimalModifiedEquation_Iterator()
  {
    order = 0; nb_additional_terms = 0; nb_rk_step_first = 10;
    fourth_order = true; order_rk_first = 4;
    dt = 0.0;
    coef_xsi = 0; xsi = 1; alpha1 = 0; alpha2 = 0;
    unstable_algo = false;
  }
  
  
  //! computation of coefficients for the scheme for a given order
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::SetOrder(int r, int m, bool horner, int nb_steps, int order_rk)
  {
    order = r;
    order_rk_first = order_rk;
    nb_additional_terms = m;
    nb_rk_step_first = nb_steps;
    unstable_algo = horner;
    fourth_order = false;
    if (r == 4)
      fourth_order = true;
    
    r = (r/2-1) + m;
    xsi = 1.0;
    if (fourth_order)
      xsi = 0.5*Real_wp(2.0*r+1)/sqrt(Real_wp(r*(r+1)));

    coef_xsi = xsi/Real_wp(2.0*r+1);
    alpha1 = 1.0/(24.0*xsi)*square(2.0*r+1);
    alpha2 = 1.0/3.0*(2.0*r+1);

    if (!unstable_algo)
      {
	alpha.Reallocate(r+1);
        if (order == 2)
          {
            alpha.Fill(-Real_wp(2)*coef_xsi);
            alpha(0) = -coef_xsi;
          }
	else if (order == 4)
	  {
	    // analytical Chebychev recurrence
            // implemented directly without array alpha
	  }
	else if (order == 6)
	  {
	    switch (m)
	      {
	      case 1 :
		alpha(0) = to_num<Real_wp>("4.7712478328541499384411147216535e-1");
		alpha(1) = to_num<Real_wp>("4.5382816380944056326045516101987e-1");
		alpha(2) = to_num<Real_wp>("2.4751846721128200684731662082238e-1");
		alpha(3) = to_num<Real_wp>("-1.784714143061375639518832540076e-1");
		break;
	      case 2 :
		alpha(0) = to_num<Real_wp>("4.6738115615071999290320108046655e-1");
		alpha(1) = to_num<Real_wp>("2.9597462301379930942349843359471e-1");
		alpha(2) = to_num<Real_wp>("3.3503175738600421868974778169782e-1");
		alpha(3) = to_num<Real_wp>("-2.0555024643695872407917120146859e-1");
		alpha(4) = to_num<Real_wp>("1.0716270988643520306272390570951e-1");
		break;
	      case 3 :
		alpha(0) = to_num<Real_wp>("2.9224243417512268335833091361302e-1");
		alpha(1) = to_num<Real_wp>("-8.9202434237559788648773165170728e-2");
		alpha(2) = to_num<Real_wp>("2.5251625981886731370112791735254e-1");
		alpha(3) = to_num<Real_wp>("3.4343688888586223777327935286367e-1");
		alpha(4) = to_num<Real_wp>("3.0032088850707324428825380348675e-1");
		alpha(5) = to_num<Real_wp>("-9.9314037149365690472218822145252e-2");
		break;
	      case 4 :
		alpha(0) = to_num<Real_wp>("2.2395900163536791802764161656827e-1");
		alpha(1) = to_num<Real_wp>("-6.7077881335753207763257732786072e-2");
		alpha(2) = to_num<Real_wp>("4.2332742112632627124180774942438e-1");
		alpha(3) = to_num<Real_wp>("3.7400940590564672348484238416183e-1");
		alpha(4) = to_num<Real_wp>("-7.7731741233199199524913475117231e-2");
		alpha(5) = to_num<Real_wp>("2.2182970785947539794704999768374e-1");
		alpha(6) = to_num<Real_wp>("-9.8315913957863903413170539934915e-2");
		break;
	      case 6 :
		alpha(0) = to_num<Real_wp>("2.4559951476696645100384548595996e-2");
		alpha(1) = to_num<Real_wp>("3.4144475348043145844925090247724e-1");
		alpha(2) = to_num<Real_wp>("2.3549926038050137716179350461615e-1");
		alpha(3) = to_num<Real_wp>("1.7748891054646703008883902106167e-1");
		alpha(4) = to_num<Real_wp>("-4.3380710685044432890891686513673e-2");
		alpha(5) = to_num<Real_wp>("2.406360789615018759357036617304e-1");
		alpha(6) = to_num<Real_wp>("-6.4802046070944728826202442142175e-2");
		alpha(7) = to_num<Real_wp>("1.7934699963699066878640931443166e-1");
		alpha(8) = to_num<Real_wp>("-9.0793197726599893805286824257263e-2");
		break;
	      case 8 :
		alpha(0) = to_num<Real_wp>("1.2040751270970788175644703898331e-1");
		alpha(1) = to_num<Real_wp>("4.131853203377171099835627363804e-2");
		alpha(2) = to_num<Real_wp>("-6.3003256866610149349743908391762e-2");
		alpha(3) = to_num<Real_wp>("4.0423344693037111594045242438108e-1");
		alpha(4) = to_num<Real_wp>("5.9060000199996495358978394180229e-2");
		alpha(5) = to_num<Real_wp>("-8.3628333355102973186323966566207e-2");
		alpha(6) = to_num<Real_wp>("1.4838299129345249126267657344487e-1");
		alpha(7) = to_num<Real_wp>("8.6719073936178474062333667076248e-2");
		alpha(8) = to_num<Real_wp>("3.563242275308527739553481801568e-1");
		alpha(9) = to_num<Real_wp>("-1.1496197567363796489381937698259e-1");
		alpha(10) = to_num<Real_wp>("4.5147781261020144095294700079972e-2");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}
	      }
	  }
	else if (order == 8)
	  {
	    switch (m)
	      {
	      case 1 :
		alpha(0) = to_num<Real_wp>("4.6546286152609731656995458076557e-1");
		alpha(1) = to_num<Real_wp>("2.9978566238717635832838526911618e-1");
		alpha(2) = to_num<Real_wp>("3.3212469181088812254494359669628e-1");
		alpha(3) = to_num<Real_wp>("-2.0559715225706462868930882931558e-1");
		alpha(4) = to_num<Real_wp>("1.0822393653290283124602538273755e-1");
		break;
	      case 2 :
		alpha(0) = to_num<Real_wp>("4.6287506420559339479337877646677e-1");
		alpha(1) = to_num<Real_wp>("-1.179716051534695734029427926526e-1");
		alpha(2) = to_num<Real_wp>("1.0852272488644990318956747531917e-1");
		alpha(3) = to_num<Real_wp>("2.5353914373431698586695354141969e-1");
		alpha(4) = to_num<Real_wp>("4.2668101101691484558964707019528e-1");
		alpha(5) = to_num<Real_wp>("-1.3364633868980555603660407074831e-1");
		break;
	      case 3 :
		alpha(0) = to_num<Real_wp>("1.1962653454866653477755542300625");
		alpha(1) = to_num<Real_wp>("-8.4938276551507434897333257558147e-2");
		alpha(2) = to_num<Real_wp>("1.5310904131990127015454948776179e-1");
		alpha(3) = to_num<Real_wp>("-3.8752549214423999512448248727239e-1");
		alpha(4) = to_num<Real_wp>("2.0077771354624367917930840468844e-1");
		alpha(5) = to_num<Real_wp>("-1.3511895028115730956350175806574e-1");
		alpha(6) = to_num<Real_wp>("5.7430618624094442475905380383523e-2");
		break;
	      case 4 :
		alpha(0) = to_num<Real_wp>("1.2110872818888856902686463537983");
		alpha(1) = to_num<Real_wp>("1.0707418399638751215208993193975e-1");
		alpha(2) = to_num<Real_wp>("-3.8591608428736696126798796068688e-1");
		alpha(3) = to_num<Real_wp>("6.7114370953512942113827492140198e-2");
		alpha(4) = to_num<Real_wp>("-1.6322953560033474932904955068175e-1");
		alpha(5) = to_num<Real_wp>("9.3871549143045677481338396570823e-2");
		alpha(6) = to_num<Real_wp>("1.4373212439051477670888584986589e-1");
		alpha(7) = to_num<Real_wp>("-7.3733890484644888127750512946343e-2");
		break;
	      case 6 :
		alpha(0) = to_num<Real_wp>("5.0421637341732516002754059189441e-1");
		alpha(1) = to_num<Real_wp>("9.6610575123860088832387942711815e-2");
		alpha(2) = to_num<Real_wp>("-6.294144807805286443180114053434e-2");
		alpha(3) = to_num<Real_wp>("-1.4706726906080055932703459823596e-1");
		alpha(4) = to_num<Real_wp>("1.1134775922019986189964901464941e-1");
		alpha(5) = to_num<Real_wp>("-5.9306566432592223965031987957005e-2");
		alpha(6) = to_num<Real_wp>("6.7300479163882318850003640803486e-1");
		alpha(7) = to_num<Real_wp>("-7.9933708410202110452798661958135e-2");
		alpha(8) = to_num<Real_wp>("-8.9853587963009016904396142038864e-2");
		alpha(9) = to_num<Real_wp>("5.392308054444847582144857343381e-2");
		break;
	      case 8 :
		alpha(0) = to_num<Real_wp>("3.324741606787252324833730495948e-1");
		alpha(1) = to_num<Real_wp>("-3.2329681338824666045105288384506e-2");
		alpha(2) = to_num<Real_wp>("1.046660382347564631433982428804e-1");
		alpha(3) = to_num<Real_wp>("2.0623524271911872315906908841745e-1");
		alpha(4) = to_num<Real_wp>("-5.9162994978463740639142663026909e-2");
		alpha(5) = to_num<Real_wp>("3.3366440974162990686321110225295e-2");
		alpha(6) = to_num<Real_wp>("1.8217077992193741916340704776251e-1");
		alpha(7) = to_num<Real_wp>("2.4912814707402698752086659456889e-1");
		alpha(8) = to_num<Real_wp>("-5.0683462243717944358717551415408e-2");
		alpha(9) = to_num<Real_wp>("1.1538203417979971659882798266512e-1");
		alpha(10) = to_num<Real_wp>("-1.3521426749514085792380439122683e-1");
		alpha(11) = to_num<Real_wp>("5.3967562273619676211506777939185e-2");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}		
	      }
	  }
	else if (order == 10)
	  {
	    switch (m)
	      {
	      case 1 :
		alpha(0) = to_num<Real_wp>("5.5462054282148278034367069010781e-1");
		alpha(1) = to_num<Real_wp>("4.5363382341602965850743180707892e-1");
		alpha(2) = to_num<Real_wp>("1.4178440193067762573900598300572e-1");
		alpha(3) = to_num<Real_wp>("-2.0979055655164342167501261941812e-1");
		alpha(4) = to_num<Real_wp>("1.5840540082698138121279345227151e-1");
		alpha(5) = to_num<Real_wp>("-9.8653612443528024127889313045831e-2");
		break;
	      case 2 :
		alpha(0) = to_num<Real_wp>("4.3760135822253414895771360927587e-1");
		alpha(1) = to_num<Real_wp>("-8.3631254886620617208720550189317e-2");
		alpha(2) = to_num<Real_wp>("2.9877478091430204142717369460918e-1");
		alpha(3) = to_num<Real_wp>("4.0487172367966396967827560433562e-1");
		alpha(4) = to_num<Real_wp>("-9.8145353288666396520247973806684e-2");
		alpha(5) = to_num<Real_wp>("1.5675429068734770812302608945921e-1");
		alpha(6) = to_num<Real_wp>("-1.1622554532856085445722047368389e-1");
		break;
	      case 3 :
		alpha(0) = to_num<Real_wp>("5.9557460646229699802429897608598e-2");
		alpha(1) = to_num<Real_wp>("-4.9986760394208784616160888188979e-2");
		alpha(2) = to_num<Real_wp>("5.729341901557947077785286038705e-1");
		alpha(3) = to_num<Real_wp>("-5.6777384181713962416624236877106e-2");
		alpha(4) = to_num<Real_wp>("1.582972225538703635436514677692e-1");
		alpha(5) = to_num<Real_wp>("-1.4678186817581440509558824062881e-1");
		alpha(6) = to_num<Real_wp>("5.9915306017002097575561189318235e-1");
		alpha(7) = to_num<Real_wp>("-1.3639592077417859475184849673575e-1");
		break;
	      case 4 :
		alpha(0) = to_num<Real_wp>("6.6833858919114633214837360668518e-2");
		alpha(1) = to_num<Real_wp>("-5.1969249648370063292483204907096e-2");
		alpha(2) = to_num<Real_wp>("3.0601916678616665991902243165809e-1");
		alpha(3) = to_num<Real_wp>("2.4713675471060738872583744147318e-1");
		alpha(4) = to_num<Real_wp>("-6.3333549389545541232300213496116e-2");
		alpha(5) = to_num<Real_wp>("2.7465944538321828515677345332194e-1");
		alpha(6) = to_num<Real_wp>("2.6800103014397551773078051573676e-1");
		alpha(7) = to_num<Real_wp>("7.1967340546913629768925666219499e-2");
		alpha(8) = to_num<Real_wp>("-1.1931479745208050999139345067477e-1");
		break;
	      case 6 :
		alpha(0) = to_num<Real_wp>("4.2598135471478686970806284089629e-2");
		alpha(1) = to_num<Real_wp>("-3.7823767249716001418989228033099e-2");
		alpha(2) = to_num<Real_wp>("4.2053173977281746605937521886161e-1");
		alpha(3) = to_num<Real_wp>("-4.4043175596230162458537358850622e-2");
		alpha(4) = to_num<Real_wp>("1.5024336202668584129606480298488e-1");
		alpha(5) = to_num<Real_wp>("2.3123279143084998432008913172275e-1");
		alpha(6) = to_num<Real_wp>("3.2396528584918277121379384850532e-1");
		alpha(7) = to_num<Real_wp>("-1.4694794337104915004404981693123e-1");
		alpha(8) = to_num<Real_wp>("1.0665713758152413987001839383556e-1");
		alpha(9) = to_num<Real_wp>("-8.4033837253175126290524236890858e-2");
		alpha(10) = to_num<Real_wp>("3.7620271337631550481952960706065e-2");
		break;
	      case 8 :
		alpha(0) = to_num<Real_wp>("3.1950635567016969381315851945034e-2");
		alpha(1) = to_num<Real_wp>("-2.9845572536789141857221170812021e-2");
		alpha(2) = to_num<Real_wp>("4.976486688903574382840742977027e-1");
		alpha(3) = to_num<Real_wp>("7.3429268558378318772773274232592e-2");
		alpha(4) = to_num<Real_wp>("1.5992240140708422186927481210951e-1");
		alpha(5) = to_num<Real_wp>("-4.773782106274187797644124991281e-2");
		alpha(6) = to_num<Real_wp>("3.6893988180036987766849460487039e-1");
		alpha(7) = to_num<Real_wp>("-4.3486807806047725694972757807044e-2");
		alpha(8) = to_num<Real_wp>("1.145269578790083225186833290496e-1");
		alpha(9) = to_num<Real_wp>("-2.0859728364568566301845372526954e-1");
		alpha(10) = to_num<Real_wp>("6.6082874370738628504100490392772e-2");
		alpha(11) = to_num<Real_wp>("6.1948439986443482492760627742939e-2");
		alpha(12) = to_num<Real_wp>("-4.478164340813285094438838424412e-2");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}				
	      }
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
	
	return;
      }
    
    if (m == 0)
      {
	// original modified equation approach
	alpha.Reallocate(r);
	for (int i = 0; i < r; i++)
	  alpha(i) = Real_wp(1) / ((2*i+2)*(2*i+3));
      }
    else
      {
	alpha.Reallocate(r);
	alpha.Fill(Real_wp(0));
	if (order == 2)
	  {
            UnivariatePolynomial<Real_wp> P0, P1, Pkp1, one, X;
            one.SetOrder(0); one(0) = 1.0;
            X.SetOrder(1); X(0) = 0; X(1) = 1.0;
            P0 = one;
            P1 = coef_xsi*X;
            for (int k = 1; k <= 2*r; k++)
              {
                Pkp1 = 2.0*coef_xsi*X*P1 + P0;
                P0 = P1;
                P1 = Pkp1;
              }

            P1 *= 1.0/xsi;
            
            for (int i = 1; i <= m; i++)
              alpha(i-1) = abs(P1(2*i+1)/P1(2*i-1));
	  }
	else if (order == 4)
	  {
	    alpha(0) = Real_wp(1) / 6;
	    switch (m)
	      {
	      case 1:
		alpha(1) = to_num<Real_wp>("3.321004576763916860463822226342e-2");
		break;
	      case 2:
		alpha(1) = to_num<Real_wp>("4.158244457517423840322624805706e-2");
		alpha(2) = to_num<Real_wp>("1.1872990171559114891715044843538e-2");
		break;
	      case 3 :
		alpha(1) = to_num<Real_wp>("4.4943502384555571422315308924359e-2");
		alpha(2) = to_num<Real_wp>("1.6639074834934280452710693679794e-2");
		alpha(3) = to_num<Real_wp>("5.5453793449274486946536454924283e-3");
		break;
	      case 4 :
		alpha(1) = to_num<Real_wp>("4.6626932320931937472972230017335e-2");
		alpha(2) = to_num<Real_wp>("1.9026247910506906481516170021695e-2");
		alpha(3) = to_num<Real_wp>("8.3229924993947598097285644437148e-3");
		alpha(4) = to_num<Real_wp>("3.0263534600476272841965437518111e-3");
		break;
	      case 6 :
		alpha(1) = to_num<Real_wp>("4.8191943618304255596853603305586e-2");
		alpha(2) = to_num<Real_wp>("2.1245530381153312569725871327494e-2");
		alpha(3) = to_num<Real_wp>("1.0905335913402067629471543043222e-2");
		alpha(4) = to_num<Real_wp>("5.8400155191334256083058008516651e-3");
		alpha(5) = to_num<Real_wp>("2.974022872138068593879179764621e-3");
		alpha(6) = to_num<Real_wp>("1.1895940982090829425971281755595e-3");
		break;
	      case 8 :
		alpha(1) = to_num<Real_wp>("4.8874695743076886955469346949026e-2");
		alpha(2) = to_num<Real_wp>("2.221372847702251933739013483606e-2");
		alpha(3) = to_num<Real_wp>("1.2031951418545319963535534138818e-2");
		alpha(4) = to_num<Real_wp>("7.0675704675870093194342425402095e-3");
		alpha(5) = to_num<Real_wp>("4.2715555267765642141718876889805e-3");
		alpha(6) = to_num<Real_wp>("2.5385043858632303054492068210675e-3");
		alpha(7) = to_num<Real_wp>("1.3882371885239739148360580180699e-3");
		alpha(8) = to_num<Real_wp>("5.8451872173811675900956759097036e-4");
		break;
	      case 12 : 
		alpha(1) = to_num<Real_wp>("4.944340972775587336110460631068e-2");
		alpha(2) = to_num<Real_wp>("2.3020218749892244697073093711375e-2");
		alpha(3) = to_num<Real_wp>("1.2970413057543734799072267784427e-2");
		alpha(4) = to_num<Real_wp>("8.0901236915297244078326413103336e-3");
		alpha(5) = to_num<Real_wp>("5.3524072466504807599696384464078e-3");
		alpha(6) = to_num<Real_wp>("3.6621590478873963807991409676975e-3");
		alpha(7) = to_num<Real_wp>("2.5446552662828974436802557260194e-3");
		alpha(8) = to_num<Real_wp>("1.7668231602569338295595914022556e-3");
		alpha(9) = to_num<Real_wp>("1.2032737674017194925872356701586e-3");
		alpha(10) = to_num<Real_wp>("7.8164120345531003131240476104291e-4");
		alpha(11) = to_num<Real_wp>("4.5776635977950915353726130736503e-4");
		alpha(12) = to_num<Real_wp>("2.0345158811963579484300198539577e-4");
		break;
	      case 16 :
		alpha(1) = to_num<Real_wp>("4.9668927377134292489247830027491e-2");
		alpha(2) = to_num<Real_wp>("2.3340026194186144849632831799587e-2");
		alpha(3) = to_num<Real_wp>("1.3342556002207686972976566354119e-2");
		alpha(4) = to_num<Real_wp>("8.4956153097689212388873312107415e-3");
		alpha(5) = to_num<Real_wp>("5.7810188271835453493094402578353e-3");
		alpha(6) = to_num<Real_wp>("4.1077454227849979763730129673038e-3");
		alpha(7) = to_num<Real_wp>("3.0032349561734667412506489498295e-3");
		alpha(8) = to_num<Real_wp>("2.2356688465421158977620715038804e-3");
		alpha(9) = to_num<Real_wp>("1.6804354385888896192329038478951e-3");
		alpha(10) = to_num<Real_wp>("1.2656763343448935230473847519948e-3");
		alpha(11) = to_num<Real_wp>("9.475778445535902188763559077273e-4");
		alpha(12) = to_num<Real_wp>("6.9818559656253441027793642875987e-4");
		alpha(13) = to_num<Real_wp>("4.9897899555258415290675794355397e-4");
		alpha(14) = to_num<Real_wp>("3.3729102674178403168397897314643e-4");
		alpha(15) = to_num<Real_wp>("2.0421913538767605223508207146765e-4");
		alpha(16) = to_num<Real_wp>("9.3357303408879366532419621879073e-5");
		break;
	      case 20 :
		alpha(1) = to_num<Real_wp>("4.9780708675943487142644716153112e-2");
		alpha(2) = to_num<Real_wp>("2.349854413988161841100977340361e-2");
		alpha(3) = to_num<Real_wp>("1.3527015620857234776279350869376e-2");
		alpha(4) = to_num<Real_wp>("8.6966053574155796348273907530784e-3");
		alpha(5) = to_num<Real_wp>("5.9934691290661245750819971087893e-3");
		alpha(6) = to_num<Real_wp>("4.3286099107941627807890185193932e-3");
		alpha(7) = to_num<Real_wp>("3.2305400615827473202644049452685e-3");
		alpha(8) = to_num<Real_wp>("2.4680626694154883695729349108834e-3");
		alpha(9) = to_num<Real_wp>("1.9169513803404858280525112940855e-3");
		alpha(10) = to_num<Real_wp>("1.5055993516132142129793427471339e-3");
		alpha(11) = to_num<Real_wp>("1.1903641146897232026462337618762e-3");
		alpha(12) = to_num<Real_wp>("9.434118861078125272064245770311e-4");
		alpha(13) = to_num<Real_wp>("7.4630945798462117758193593538233e-4");
		alpha(14) = to_num<Real_wp>("5.8645469252317198953432017871164e-4");
		alpha(15) = to_num<Real_wp>("4.5499421195307188078369668838993e-4");
		alpha(16) = to_num<Real_wp>("3.4555995368591384944569167920759e-4");
		alpha(17) = to_num<Real_wp>("2.5347614906693613550536407811197e-4");
		alpha(18) = to_num<Real_wp>("1.7524665405605959898003616390338e-4");
		alpha(19) = to_num<Real_wp>("1.082148014961314507287816236285e-4");
		alpha(20) = to_num<Real_wp>("5.0332462844052142276297138139476e-5");
		break;
	      case 30 :
		alpha(1) = to_num<Real_wp>("4.9897865747404548689329358076311e-2");
		alpha(2) = to_num<Real_wp>("2.3664685844033886312712522326636e-2");
		alpha(3) = to_num<Real_wp>("1.3720347138285461333148280558674e-2");
		alpha(4) = to_num<Real_wp>("8.9072627449107070571865521396862e-3");
		alpha(5) = to_num<Real_wp>("6.2161382653402640090052983972563e-3");
		alpha(6) = to_num<Real_wp>("4.5600981499835073063941533853816e-3");
		alpha(7) = to_num<Real_wp>("3.4687788600044238273841863559801e-3");
		alpha(8) = to_num<Real_wp>("2.7116350670823221338496627916537e-3");
		alpha(9) = to_num<Real_wp>("2.1648442607248166020974987771879e-3");
		alpha(10) = to_num<Real_wp>("1.7570632596313571181564960876404e-3");
		alpha(11) = to_num<Real_wp>("1.4448290562106052396017510871034e-3");
		alpha(12) = to_num<Real_wp>("1.2004342597597285127712006117581e-3");
		alpha(13) = to_num<Real_wp>("1.005537254983693799893465872319e-3");
		alpha(14) = to_num<Real_wp>("8.4760390337596068122587533165793e-4");
		alpha(15) = to_num<Real_wp>("7.1783237083814077836868629569554e-4");
		alpha(16) = to_num<Real_wp>("6.0989437691059053912441210113984e-4");
		alpha(17) = to_num<Real_wp>("5.1914534739942991752732320905173e-4");
		alpha(18) = to_num<Real_wp>("4.4211394373668342992202814124363e-4");
		alpha(19) = to_num<Real_wp>("3.7616347129722225315963878935953e-4");
		alpha(20) = to_num<Real_wp>("3.1926206351242938327024177614656e-4");
		alpha(21) = to_num<Real_wp>("2.698234571674193107454091813788e-4");
		alpha(22) = to_num<Real_wp>("2.265946197008383028434715991204e-4");
		alpha(23) = to_num<Real_wp>("1.885751080055427264434558006883e-4");
		alpha(24) = to_num<Real_wp>("1.54958315102334234608630659972e-4");
		alpha(25) = to_num<Real_wp>("1.2508806818661371688383872583175e-4");
		alpha(26) = to_num<Real_wp>("9.8426158586739406846611754086613e-5");
		alpha(27) = to_num<Real_wp>("7.4527765832083746276320110993243e-5");
		alpha(28) = to_num<Real_wp>("5.3022655856918603041030361394682e-5");
		alpha(29) = to_num<Real_wp>("3.3600653079360952789041157645651e-5");
		alpha(30) = to_num<Real_wp>("1.6000310853176845081284981921182e-5");
		break;
	      case 40 :
                alpha(1) = to_num<Real_wp>("4.994116209311143011533160530653929833321707668960385798712111619e-2");
                alpha(2) = to_num<Real_wp>("2.372608494529513297777392027107929546553733263935631885049713326e-2");
                alpha(3) = to_num<Real_wp>("1.37917946004011613761330637521946067519635774636411754141214761e-2");
                alpha(4) = to_num<Real_wp>("8.985113242820817459285317204622726983067800075835621950763293096e-3");
                alpha(5) = to_num<Real_wp>("6.298427893183013014844512407581986538609203470074163970059748138e-3");
                alpha(6) = to_num<Real_wp>("4.645647017064986934369832146598624035096602987143443814255068702e-3");
                alpha(7) = to_num<Real_wp>("3.556822503034715953366536120813828906954008329639597345824386661e-3");
                alpha(8) = to_num<Real_wp>("2.801649825410798757444174056146560926011244690365610060846055764e-3");
                alpha(9) = to_num<Real_wp>("2.256455720508941344919329158930903564818359478714590028382176259e-3");
                alpha(10) = to_num<Real_wp>("1.8499944474505592669940201856561169076728420318635860967910932e-3");
                alpha(11) = to_num<Real_wp>("1.538869321430919564347226146433307490702663188174173516447002691e-3");
                alpha(12) = to_num<Real_wp>("1.295419662309283539633405960718499574617296039282527980778986683e-3");
                alpha(13) = to_num<Real_wp>("1.101337704366042743063678495531848710754180777770511845828358479e-3");
                alpha(14) = to_num<Real_wp>("9.44114439367844398675582575426856627569286138761463835425414767e-4");
                alpha(15) = to_num<Real_wp>("8.149670820361973728885208188571971336807688533248823044941016207e-4");
                alpha(16) = to_num<Real_wp>("7.075820540532531277386899939769789839360938296098863125111366868e-4");
                alpha(17) = to_num<Real_wp>("6.173263096487276646023772383126662052930842261553558121964955156e-4");
                alpha(18) = to_num<Real_wp>("5.407376775196103337394834989833745999452337699243985119082705641e-4");
                alpha(19) = to_num<Real_wp>("4.751868442115938289221732614929432590311034035391839054873435163e-4");
                alpha(20) = to_num<Real_wp>("4.186479530284257547168338461731556728750792161647514043098269644e-4");
                alpha(21) = to_num<Real_wp>("3.695396836557823527995546437248650854868071350231650438756614337e-4");
                alpha(22) = to_num<Real_wp>("3.266131061846749571870216998389819455289531711578493385317745996e-4");
                alpha(23) = to_num<Real_wp>("2.888712111789442542317079014240891092926011819344698768223366318e-4");
                alpha(24) = to_num<Real_wp>("2.555102873650601153292054084446476893057211318177229421690970086e-4");
                alpha(25) = to_num<Real_wp>("2.258766210645531057240784262463458072427588303753229734838744547e-4");
                alpha(26) = to_num<Real_wp>("1.99434105441976358204926357712551450233099335241743959633979065e-4");
                alpha(27) = to_num<Real_wp>("1.757397271878298686795592728158440991241990258208105673914214768e-4");
                alpha(28) = to_num<Real_wp>("1.544248146152602594047655338277954825624326063984804241295531726e-4");
                alpha(29) = to_num<Real_wp>("1.351805498266901560608784608095379522941377497405074748513993091e-4");
                alpha(30) = to_num<Real_wp>("1.177466716249935058190152259880929857506643920837244530545407163e-4");
                alpha(31) = to_num<Real_wp>("1.019025905182544356554898511128843533263418666133974891202958035e-4");
                alpha(32) = to_num<Real_wp>("8.746034461302998486967411400736058116307191985329227998177386405e-5");
                alpha(33) = to_num<Real_wp>("7.42589729945426753375210700671350093143166919866374770821878375e-5");
                alpha(34) = to_num<Real_wp>("6.215998968842688018928875323197415892761251830658845506587695327e-5");
                alpha(35) = to_num<Real_wp>("5.104371884097577332337915419036716332725158613642494024923030461e-5");
                alpha(36) = to_num<Real_wp>("4.08063087731506625031396682464318824085756611529171114299789111e-5");
                alpha(37) = to_num<Real_wp>("3.135728487939624298896237464649825055743954922158656655050232696e-5");
                alpha(38) = to_num<Real_wp>("2.261753302099693185605180262008486716365289383809723283420113694e-5");
                alpha(39) = to_num<Real_wp>("1.451762897448778924426505272993143913183330853625839213723230973e-5");
                alpha(40) = to_num<Real_wp>("6.996447683615974576136753170729964886548898165889048997508577005e-6");
		break;
	      case 50 :
                alpha(1) = to_num<Real_wp>("4.996179499065407733551511382082532627418155878310683779100155028e-2");
                alpha(2) = to_num<Real_wp>("2.375534473776872726178066349043351421597373024851511281443212889e-2");
                alpha(3) = to_num<Real_wp>("1.382584297241458856640326777342531942727656978311591049148233633e-2");
                alpha(4) = to_num<Real_wp>("9.022213007936955075917291174307468493576299244858969855429547284e-3");
                alpha(5) = to_num<Real_wp>("6.337643145400172323417899416652895580330632699447554184229580114e-3");
                alpha(6) = to_num<Real_wp>("4.686415474529860069463063343208781074673654365734508941682229771e-3");
                alpha(7) = to_num<Real_wp>("3.598779857259949432627650138087931050075606841716323938175649805e-3");
                alpha(8) = to_num<Real_wp>("2.844546523416512883368675358398188771659587248801213212265640518e-3");
                alpha(9) = to_num<Real_wp>("2.300113333570260952815815747200334980404593110428483672177419373e-3");
                alpha(10) = to_num<Real_wp>("1.894280982547854928458369194571026657964837087520086977006796134e-3");
                alpha(11) = to_num<Real_wp>("1.583684391997580807698242866119624810478275914960702733594317951e-3");
                alpha(12) = to_num<Real_wp>("1.340685141930870631310626181764768881836625070394539707792606489e-3");
                alpha(13) = to_num<Real_wp>("1.146991597812668765528195111466206169490445117269222819678222779e-3");
                alpha(14) = to_num<Real_wp>("9.901067273663130303494091801915692041157446094828618144671700651e-4");
                alpha(15) = to_num<Real_wp>("8.612568231222811405541130531409312203230053245077663634244700073e-4");
                alpha(16) = to_num<Real_wp>("7.541353131399718269904470221393649242772340934932974182716823334e-4");
                alpha(17) = to_num<Real_wp>("6.641146455747462253924922818634444935661383881709681332011172889e-4");
                alpha(18) = to_num<Real_wp>("5.877370178271301299435587779230252968437854167514319027238911903e-4");
                alpha(19) = to_num<Real_wp>("5.223766339742214284897325334496685131127720352116229812289392628e-4");
                alpha(20) = to_num<Real_wp>("4.660105013520804032401792452173418529297815474633268886292696671e-4");
                alpha(21) = to_num<Real_wp>("4.170596552080425397654074131367191447593551418644123786496859597e-4");
                alpha(22) = to_num<Real_wp>("3.742771207704223649604103765157071228691332640379271748190034475e-4");
                alpha(23) = to_num<Real_wp>("3.366675249340799025783836344612137954352070356061240810713715477e-4");
                alpha(24) = to_num<Real_wp>("3.034285363865862763046716427439308183359425602990594518851129127e-4");
                alpha(25) = to_num<Real_wp>("2.739076133565533145404100682779465502972004760309897492817808236e-4");
                alpha(26) = to_num<Real_wp>("2.475696506474475750388979744869136082847869464934661446822375324e-4");
                alpha(27) = to_num<Real_wp>("2.239724961599355811942166757807622240102634203545516253977482625e-4");
                alpha(28) = to_num<Real_wp>("2.027482227652454507249608233714380125996125953561798230834051765e-4");
                alpha(29) = to_num<Real_wp>("1.835886595757735491076075688064061284158968409752257856146626279e-4");
                alpha(30) = to_num<Real_wp>("1.662341103249266907547084105447537558074896071369108694556860884e-4");
                alpha(31) = to_num<Real_wp>("1.504644809874938665081745910549332766910454580438383456198309478e-4");
                alpha(32) = to_num<Real_wp>("1.360922460310849032032819454081685038200455667038533666002761198e-4");
                alpha(33) = to_num<Real_wp>("1.229568303546862298618766709109281460418889761976153552817670173e-4");
                alpha(34) = to_num<Real_wp>("1.109200903635436316147392450598812752715664416845817169318611655e-4");
                alpha(35) = to_num<Real_wp>("9.986265509425214527587819829275390139401134003012934297442271244e-5");
                alpha(36) = to_num<Real_wp>("8.968094526321730773247829482584406278863391521679691911114341733e-5");
                alpha(37) = to_num<Real_wp>("8.028473038219130996825159302357646892360763229329425346241959275e-5");
                alpha(38) = to_num<Real_wp>("7.159511572833729056834834991164641316106840957000804865853285546e-5");
                alpha(39) = to_num<Real_wp>("6.354287483987043151446160512496627097407623848145413840012539033e-5");
                alpha(40) = to_num<Real_wp>("5.606706137501478896368887078474678633694109181256559806322685649e-5");
                alpha(41) = to_num<Real_wp>("4.91138480918796125787261943871505746125851253025798760463802648e-5");
                alpha(42) = to_num<Real_wp>("4.263555144701064488031933028393207384862664655852084381725623042e-5");
                alpha(43) = to_num<Real_wp>("3.658980865194494311717671723315108452960062138456597256748892432e-5");
                alpha(44) = to_num<Real_wp>("3.093888054646052754898832052393629961940583014522939362612151078e-5");
                alpha(45) = to_num<Real_wp>("2.564905877287368234259646339456793129506526452294316139802737626e-5");
                alpha(46) = to_num<Real_wp>("2.06901597885785522158835601104031510911142710746387678270848633e-5");
                alpha(47) = to_num<Real_wp>("1.60350914759844889335538270278842650315358644525252951578979154e-5");
                alpha(48) = to_num<Real_wp>("1.16594806836993629999923859392457657372192070482532486582865888e-5");
                alpha(49) = to_num<Real_wp>("7.541352100451665490606301028824423033348913679028548551012853686e-6");
                alpha(50) = to_num<Real_wp>("3.660850531562567611745159650747594060239038579502867261306501637e-6");
		break;
	      default : 
		{
		  cout << "not implemented" << endl;
		  abort();
		}
	      }
	  }
	else if (order == 6)
	  {
	    alpha(0) = Real_wp(1)/6;
	    alpha(1) = Real_wp(1)/20;
	    switch (m)
	      {
	      case 1:
		alpha(2) = to_num<Real_wp>("2.3011697377029043407137268362218e-2");
		break;
	      case 2 : 
		alpha(2) = to_num<Real_wp>("2.3671834611574532146983775017284e-2");
		alpha(3) = to_num<Real_wp>("1.1303808929497150379726017497471e-2");
		break;
	      case 3 :
		alpha(2) = to_num<Real_wp>("2.3770702828833944092974290244293e-2");
		alpha(3) = to_num<Real_wp>("1.313630315137605811295969934376e-2");
		alpha(4) = to_num<Real_wp>("5.9789879979290942693573093044385e-3");
		break;
	      case 4 :
		alpha(2) = to_num<Real_wp>("2.3795210096341713441124048517593e-2");
		alpha(3) = to_num<Real_wp>("1.3606835297959850475749430092625e-2");
		alpha(4) = to_num<Real_wp>("7.8177651910433930556330179079954e-3");
		alpha(5) = to_num<Real_wp>("3.4417621011747959026436267649733e-3");
		break;
	      case 6 :
		alpha(2) = to_num<Real_wp>("2.3806423899657473179717002532468e-2");
		alpha(3) = to_num<Real_wp>("1.3826850651801212149815920844056e-2");
		alpha(4) = to_num<Real_wp>("8.7857525130798407591525144806219e-3");
		alpha(5) = to_num<Real_wp>("5.5896081026259280280347599429867e-3");
		alpha(6) = to_num<Real_wp>("3.2257203621165462043285414376186e-3");
		alpha(7) = to_num<Real_wp>("1.3968629806319902542872602116545e-3");
		break;
	      case 8 :
		alpha(2) = to_num<Real_wp>("2.3808550345974640565329651301955e-2");
		alpha(3) = to_num<Real_wp>("1.3869271149763264265678942352333e-2");
		alpha(4) = to_num<Real_wp>("8.990575501281959986095813575754e-3");
		alpha(5) = to_num<Real_wp>("6.1203427158164662728947766762799e-3");
		alpha(6) = to_num<Real_wp>("4.166486953103846109769755755135e-3");
		alpha(7) = to_num<Real_wp>("2.7050601994944288248611996066655e-3");
		alpha(8) = to_num<Real_wp>("1.5709872453023789509447478171578e-3");
		alpha(9) = to_num<Real_wp>("6.8708250082360026290934737504511e-4");
		break;
	      case 12 :
		alpha(2) = to_num<Real_wp>("2.380934811643195557947136952088e-2");
		alpha(3) = to_num<Real_wp>("1.388532687500807191734851605841e-2");
		alpha(4) = to_num<Real_wp>("9.0720676114924947354212825493769e-3");
		alpha(5) = to_num<Real_wp>("6.3521648942908270969571271283195e-3");
		alpha(6) = to_num<Real_wp>("4.6309524600302645239964961466522e-3");
		alpha(7) = to_num<Real_wp>("3.4378414397986296713323712886961e-3");
		alpha(8) = to_num<Real_wp>("2.5517462631038664639183464506957e-3");
		alpha(9) = to_num<Real_wp>("1.8639530455392648336703483515325e-3");
		alpha(10) = to_num<Real_wp>("1.3175042333166788643137695001446e-3");
		alpha(11) = to_num<Real_wp>("8.7857566048384580997357014483491e-4");
		alpha(12) = to_num<Real_wp>("5.2386081512707588795775112412014e-4");
		alpha(13) = to_num<Real_wp>("2.3564908501275929557513095842806e-4");
		break;
	      case 16 :
		alpha(2) = to_num<Real_wp>("2.3809473872691212010168064368753e-2");
		alpha(3) = to_num<Real_wp>("1.3887873913014871223059162341031e-2");
		alpha(4) = to_num<Real_wp>("9.0854645850400557656995310587251e-3");
		alpha(5) = to_num<Real_wp>("6.3930033220156899142492151241719e-3");
		alpha(6) = to_num<Real_wp>("4.7213570125571350470987841030048e-3");
		alpha(7) = to_num<Real_wp>("3.5984284946865465788135773062227e-3");
		alpha(8) = to_num<Real_wp>("2.7942189883394158145779546897036e-3");
		alpha(9) = to_num<Real_wp>("2.1881699512965537227953491905061e-3");
		alpha(10) = to_num<Real_wp>("1.714010925603092669038770523698e-3");
		alpha(11) = to_num<Real_wp>("1.3336086692140833723631292224068e-3");
		alpha(12) = to_num<Real_wp>("1.0235708831494438381818511985283e-3");
		alpha(13) = to_num<Real_wp>("7.6832362129155115524010448509009e-4");
		alpha(14) = to_num<Real_wp>("5.5665202304951710508091463165356e-4");
		alpha(15) = to_num<Real_wp>("3.8002010759017268691530016271549e-4");
		alpha(16) = to_num<Real_wp>("2.3173347868547866539081531871396e-4");
		alpha(17) = to_num<Real_wp>("1.0647442820707665885704181008165e-4");
		break;
	      case 20 :
		alpha(2) = to_num<Real_wp>("2.3809505334863033744021129008969e-2");
		alpha(3) = to_num<Real_wp>("1.3888512894191842595305916498786e-2");
		alpha(4) = to_num<Real_wp>("9.0888774792547729897583305799134e-3");
		alpha(5) = to_num<Real_wp>("6.4037262035104822599277701278365e-3");
		alpha(6) = to_num<Real_wp>("4.7462171124625421527731303453253e-3");
		alpha(7) = to_num<Real_wp>("3.6453600060525291778882383033238e-3");
		alpha(8) = to_num<Real_wp>("2.8702955639037260064246130525525e-3");
		alpha(9) = to_num<Real_wp>("2.2977465760068952178055852819084e-3");
		alpha(10) = to_num<Real_wp>("1.8578800573403187977108671027101e-3");
		alpha(11) = to_num<Real_wp>("1.5094157116917289680582647735856e-3");
		alpha(12) = to_num<Real_wp>("1.2270059862417344693969358807742e-3");
		alpha(13) = to_num<Real_wp>("9.943924409887095773051833600099e-4");
		alpha(14) = to_num<Real_wp>("8.005873234616707444654962037276e-4");
		alpha(15) = to_num<Real_wp>("6.377461857667847351627277913974e-4");
		alpha(16) = to_num<Real_wp>("4.9999403774550792747228498018872e-4");
		alpha(17) = to_num<Real_wp>("3.8277396457762694629321056881951e-4");
		alpha(18) = to_num<Real_wp>("2.824717049565680539636208353317e-4");
		alpha(19) = to_num<Real_wp>("1.9618321624477360812677857609978e-4");
		alpha(20) = to_num<Real_wp>("1.2155772934274885978071022662889e-4");
		alpha(21) = to_num<Real_wp>("5.6683448818455260955878114780877e-5");
		break;
	      case 30 :
		alpha(2) = to_num<Real_wp>("2.3809520867821512898252648916764e-2");
		alpha(3) = to_num<Real_wp>("1.3888828931019110502740902664794e-2");
		alpha(4) = to_num<Real_wp>("9.0905824733090999410027502592486e-3");
		alpha(5) = to_num<Real_wp>("6.4091897976824895193325767463739e-3");
		alpha(6) = to_num<Real_wp>("4.7592788967150568799562917841547e-3");
		alpha(7) = to_num<Real_wp>("3.6710834665011855592872136944028e-3");
		alpha(8) = to_num<Real_wp>("2.9142667430636658251181987992293e-3");
		alpha(9) = to_num<Real_wp>("2.3650880588988427563914821462591e-3");
		alpha(10) = to_num<Real_wp>("1.9523108597809635666495375878683e-3");
		alpha(11) = to_num<Real_wp>("1.6326857798041110470025863861354e-3");
		alpha(12) = to_num<Real_wp>("1.3788512312735149394143554307739e-3");
		alpha(13) = to_num<Real_wp>("1.1729219754490690663052062507646e-3");
		alpha(14) = to_num<Real_wp>("1.0028792939716610777256012604795e-3");
		alpha(15) = to_num<Real_wp>("8.6043513838370241357627256778613e-4");
		alpha(16) = to_num<Real_wp>("7.3972159322311341516775875396156e-4");
		alpha(17) = to_num<Real_wp>("6.3646828970794518297008205943675e-4");
		alpha(18) = to_num<Real_wp>("5.4747978545411630986846629914146e-4");
		alpha(19) = to_num<Real_wp>("4.7030109380198551809332422698024e-4");
		alpha(20) = to_num<Real_wp>("4.0300178067938011750846349870611e-4");
		alpha(21) = to_num<Real_wp>("3.4403454159222865589613649758404e-4");
		alpha(22) = to_num<Real_wp>("2.9214041749819617499860602055681e-4");
		alpha(23) = to_num<Real_wp>("2.4628330898713936411666497920746e-4");
		alpha(24) = to_num<Real_wp>("2.0560314521353369595379406684145e-4");
		alpha(25) = to_num<Real_wp>("1.6938123011728817592847274059075e-4");
		alpha(26) = to_num<Real_wp>("1.3701381549876624073365499614919e-4");
		alpha(27) = to_num<Real_wp>("1.0799145447230932869708116027642e-4");
		alpha(28) = to_num<Real_wp>("8.1882576087967479274258779755178e-5");
		alpha(29) = to_num<Real_wp>("5.8320247578894908191604740455116e-5");
		alpha(30) = to_num<Real_wp>("3.6991408080138315231466907118753e-5");
		alpha(31) = to_num<Real_wp>("1.7628055804869344769131988784907e-5");
		break;
	      case 40 :
                alpha(2) = to_num<Real_wp>("2.380952302434591865965470071428502318006228436408212531912021243e-2");
                alpha(3) = to_num<Real_wp>("1.388887287568579920951877013521200819690778294007301465838315258e-2");
                alpha(4) = to_num<Real_wp>("9.090821570309453435384865450735585370712419377911493771786370687e-3");
                alpha(5) = to_num<Real_wp>("6.409968765673327459083641698598586039261140194147109449866201719e-3");
                alpha(6) = to_num<Real_wp>("4.761189624416853965441660783145783174132871982241147947356255796e-3");
                alpha(7) = to_num<Real_wp>("3.674983300499606997493557759788448278992049919789696905140376281e-3");
                alpha(8) = to_num<Real_wp>("2.921247778113989407496359358889332232370367452891682668578338349e-3");
                alpha(9) = to_num<Real_wp>("2.376394050205437588175122895517816261490897847030593606125663385e-3");
                alpha(10) = to_num<Real_wp>("1.969211177544332643288643513275892356991509786851960786658834983e-3");
                alpha(11) = to_num<Real_wp>("1.656329761141787852001300045978929002748841264930193311616359419e-3");
                alpha(12) = to_num<Real_wp>("1.410134854530365839183615765129359003667636557598491365824505213e-3");
                alpha(13) = to_num<Real_wp>("1.21239662786449710620497370952538098739966382630825783652915151e-3");
                alpha(14) = to_num<Real_wp>("1.05071950141043272800528986486926225788337720850819208382171815e-3");
                alpha(15) = to_num<Real_wp>("9.164643611349105472958404682054312440404475389426727401064266421e-4");
                alpha(16) = to_num<Real_wp>("8.034805089287953340159795354486507470122875703895883854450915704e-4");
                alpha(17) = to_num<Real_wp>("7.07303454917275601160493676046842990204732823999347950017005412e-4");
                alpha(18) = to_num<Real_wp>("6.246323042090537456092991836510013958558862976140710570600228924e-4");
                alpha(19) = to_num<Real_wp>("5.529815672309463823229431101104506151629370961081931812304561191e-4");
                alpha(20) = to_num<Real_wp>("4.904454035478075795112898945928185643851113241852164251373601954e-4");
                alpha(21) = to_num<Real_wp>("4.355361439762628228266129404167553302862079775480622796152987378e-4");
                alpha(22) = to_num<Real_wp>("3.870726692730516369762078276200261018990137349508468928617062016e-4");
                alpha(23) = to_num<Real_wp>("3.441025443991654358001711244762691918131086842869336017663124533e-4");
                alpha(24) = to_num<Real_wp>("3.058471046802693562430557845053115675919578194953887705853541629e-4");
                alpha(25) = to_num<Real_wp>("2.716621941010035619411461956018037621582197108028239041315197222e-4");
                alpha(26) = to_num<Real_wp>("2.410096236950923101767551704784082957094164193117269143702367151e-4");
                alpha(27) = to_num<Real_wp>("2.134360276090505273065297918728795644776578426160380070181354324e-4");
                alpha(28) = to_num<Real_wp>("1.88556884237510501114043142383754491748211658743976640984971042e-4");
                alpha(29) = to_num<Real_wp>("1.660442014050243701591129918026270802496082592918816353109615553e-4");
                alpha(30) = to_num<Real_wp>("1.456168513317136525635886019448259109003243353757041744443503013e-4");
                alpha(31) = to_num<Real_wp>("1.27032862965615656976488318275796233726908771974010101712556041e-4");
                alpha(32) = to_num<Real_wp>("1.100831916630908520210416914613171929401921624103516248664739075e-4");
                alpha(33) = to_num<Real_wp>("9.45866268620223811674315805276794671459557574883429940115650808e-5");
                alpha(34) = to_num<Real_wp>("8.03855924236442408683680925253534059913589571228334530168672966e-5");
                alpha(35) = to_num<Real_wp>("6.734265812696208183057807222583678464279764810975268640204975796e-5");
                alpha(36) = to_num<Real_wp>("5.533762496502419708502933466502722663420306963453278654953959132e-5");
                alpha(37) = to_num<Real_wp>("4.426507817779421288465961965070102369193091238307355613394911785e-5");
                alpha(38) = to_num<Real_wp>("3.40323246715198137796375814634176208197791124670259289848368758e-5");
                alpha(39) = to_num<Real_wp>("2.455764837246703883634467318358141418153388273963046906274424321e-5");
                alpha(40) = to_num<Real_wp>("1.576882991770437184955308782067272662158688034885550782636296463e-5");
                alpha(41) = to_num<Real_wp>("7.601887056904971622968360285591191837907817351000712664257795507e-6");
		break;
	      case 50 :
                alpha(2) = to_num<Real_wp>("2.380952352958075732711079261492891020694170893227839152480341849e-2");
                alpha(3) = to_num<Real_wp>("1.388888317783629727774974990552466086189811954126368976508926655e-2");
                alpha(4) = to_num<Real_wp>("9.090877823421086882722785983627480285501292098879759642352199146e-3");
                alpha(5) = to_num<Real_wp>("6.410153307613924897773839168351478516555298294357957463142677198e-3");
                alpha(6) = to_num<Real_wp>("4.761647138000470447290450029199100876134086967995505052049714646e-3");
                alpha(7) = to_num<Real_wp>("3.675931072360931208198544870707083192152607797932835311857562377e-3");
                alpha(8) = to_num<Real_wp>("2.922977732830275731408377714008825928371809563292946358875050192e-3");
                alpha(9) = to_num<Real_wp>("2.379264836335488645727624430657160454504098849238740565975612011e-3");
                alpha(10) = to_num<Real_wp>("1.973629764745397644281054344397263254572453329078440742406008233e-3");
                alpha(11) = to_num<Real_wp>("1.662723796671743248885318537891813523971433474379831478283403333e-3");
                alpha(12) = to_num<Real_wp>("1.418919075567510377483655723351703894355959840092688300593609403e-3");
                alpha(13) = to_num<Real_wp>("1.22393829000231399389907104684899860217520891945691858302282146e-3");
                alpha(14) = to_num<Real_wp>("1.065308384093167395365496401689365446427679786638401097912058112e-3");
                alpha(15) = to_num<Real_wp>("9.3429208647493415533077317271552528906304157567123213091321027e-4");
                alpha(16) = to_num<Real_wp>("8.246318535968312835423133176460803521739396334249461325788231692e-4");
                alpha(17) = to_num<Real_wp>("7.317597326276777158749635137702821687644435010047484750112077319e-4");
                alpha(18) = to_num<Real_wp>("6.522845987605955123643709952119334582862632131950583521758944256e-4");
                alpha(19) = to_num<Real_wp>("5.836502254869514597136643227522517780157029154903694206785217334e-4");
                alpha(20) = to_num<Real_wp>("5.239019327471547855765390330361244761085984813082139630087627049e-4");
                alpha(21) = to_num<Real_wp>("4.715242731171399979241634123485537221019160275356420674423697754e-4");
                alpha(22) = to_num<Real_wp>("4.253264057284633766695977670610687610025355444647910011304243285e-4");
                alpha(23) = to_num<Real_wp>("3.843601770910438595476584123241931548767567179107014822907962843e-4");
                alpha(24) = to_num<Real_wp>("3.478610092837036676526954668005271622218062125114420148101301339e-4");
                alpha(25) = to_num<Real_wp>("3.152048824068325969763233808893616053370768580097581974275223462e-4");
                alpha(26) = to_num<Real_wp>("2.85876766371436300304964851784441740690022756033990596563475507e-4");
                alpha(27) = to_num<Real_wp>("2.59447243870641189584449790013261256618343682156046288911152665e-4");
                alpha(28) = to_num<Real_wp>("2.355550210423023863404559561753784315083991051934395131451939276e-4");
                alpha(29) = to_num<Real_wp>("2.138936912414524591893941022614694387211882334527996420818830848e-4");
                alpha(30) = to_num<Real_wp>("1.942015903658437476264992811321739504692386499596458628331785333e-4");
                alpha(31) = to_num<Real_wp>("1.762539175889601506395714312908664217721236715837709043825946609e-4");
                alpha(32) = to_num<Real_wp>("1.598565328959587261991154307648086101635525950522508208459298928e-4");
                alpha(33) = to_num<Real_wp>("1.448410105724541591045600428543541098455708608985664306172144183e-4");
                alpha(34) = to_num<Real_wp>("1.310606459585436342897090925202350634032291342196365221845782877e-4");
                alpha(35) = to_num<Real_wp>("1.18387195919646085096297423642487946530369156622389759016896827e-4");
                alpha(36) = to_num<Real_wp>("1.067081920534606969518089519071594272369892621404125315062565916e-4");
                alpha(37) = to_num<Real_wp>("9.592470707517853419210300819997613494264817551591726895368904992e-5");
                alpha(38) = to_num<Real_wp>("8.594948432257199293770028557191452944921718632290548849044800598e-5");
                alpha(39) = to_num<Real_wp>("7.670536153507678105583742628167606498633857511992123700588412258e-5");
                alpha(40) = to_num<Real_wp>("6.812393549993727995594744658710420651614410592870814136966093297e-5");
                alpha(41) = to_num<Real_wp>("6.014442555283771960064592456764301984770971071841626738353692095e-5");
                alpha(42) = to_num<Real_wp>("5.271270245721206802145967575657126541729915563099868758908785486e-5");
                alpha(43) = to_num<Real_wp>("4.578045568156732544019440482473867137385416939465982684079924271e-5");
                alpha(44) = to_num<Real_wp>("3.930447711098175015636582209805519895896545205169084134053515128e-5");
                alpha(45) = to_num<Real_wp>("3.324604315922343022364012876694854514566471554670547024180790753e-5");
                alpha(46) = to_num<Real_wp>("2.757038036693016100035916805888417616101878769275365167894895211e-5");
                alpha(47) = to_num<Real_wp>("2.224620207484304031454707161691884618382688005563343971766646118e-5");
                alpha(48) = to_num<Real_wp>("1.724530579057341649910546649508180872772980014239630338437987125e-5");
                alpha(49) = to_num<Real_wp>("1.254222252648798441571680811735122635330980948747507555685278227e-5");
                alpha(50) = to_num<Real_wp>("8.113910752303153114577456290631009790333310214122414692263905435e-6");
                alpha(51) = to_num<Real_wp>("3.93948873742696610806644046619641931899794743536002006349509225e-6");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}
	      }
	  }
	else if (order == 8)
	  {
	    alpha(0) = Real_wp(1) / 6;
	    alpha(1) = Real_wp(1) / 20;
	    alpha(2) = Real_wp(1) / 42;
	    switch (m)
	      {
	      case 1 :
		alpha(3) = to_num<Real_wp>("1.1513847175056756080917762052169e-2");
		break;
	      case 2 :
		alpha(3) = to_num<Real_wp>("1.3222447684161570426849902052954e-2");
		alpha(4) = to_num<Real_wp>("6.0449069748149949351825212939304e-3");
		break;
	      case 3 :
		alpha(3) = to_num<Real_wp>("1.3643769902713123547107711225411e-2");
		alpha(4) = to_num<Real_wp>("7.8607878623178028223836177901558e-3");
		alpha(5) = to_num<Real_wp>("3.4659096284103431332868420571128e-3");
		break;
	      case 4 :
		alpha(3) = to_num<Real_wp>("1.3781690789753742769593017430471e-2");
		alpha(4) = to_num<Real_wp>("8.5224887613035640207559602554457e-3");
		alpha(5) = to_num<Real_wp>("4.926458222869362742343693072044e-3");
		alpha(6) = to_num<Real_wp>("2.1394199384971023268172870289296e-3");
		break;
	      case 6 :
		alpha(3) = to_num<Real_wp>("1.3860247628094751256373177379183e-2");
		alpha(4) = to_num<Real_wp>("8.9298425176728946712298743489426e-3");
		alpha(5) = to_num<Real_wp>("5.9447410912673720688382415398576e-3");
		alpha(6) = to_num<Real_wp>("3.8326242338863140864706900141723e-3");
		alpha(7) = to_num<Real_wp>("2.2164935082237868781683980381764e-3");
		alpha(8) = to_num<Real_wp>("9.6333400865042246178206050297396e-4");
		break;
	      case 8 :
		alpha(3) = to_num<Real_wp>("1.3878703036034084274506546055696e-2");
		alpha(4) = to_num<Real_wp>("9.0319311818930614338794147985114e-3");
		alpha(5) = to_num<Real_wp>("6.2303787602604084528196494532573e-3");
		alpha(6) = to_num<Real_wp>("4.3756706660057363935374794880686e-3");
		alpha(7) = to_num<Real_wp>("3.0193196084923164474157440476702e-3");
		alpha(8) = to_num<Real_wp>("1.9735396854736650941016478916547e-3");
		alpha(9) = to_num<Real_wp>("1.1531463100761187014974550925766e-3");
		alpha(10) = to_num<Real_wp>("5.0797618312497558724965586957245e-4");
		break;
	      case 12 :
		alpha(3) = to_num<Real_wp>("1.388677215388488587372089269688e-2");
		alpha(4) = to_num<Real_wp>("9.0783144456827633698495311090641e-3");
		alpha(5) = to_num<Real_wp>("6.3698129297808748348137317536736e-3");
		alpha(6) = to_num<Real_wp>("4.6682743487766084799222200215844e-3");
		alpha(7) = to_num<Real_wp>("3.5017195115571212511631844298331e-3");
		alpha(8) = to_num<Real_wp>("2.6449655466004855974019709008603e-3");
		alpha(9) = to_num<Real_wp>("1.9847633953582002038340844167761e-3");
		alpha(10) = to_num<Real_wp>("1.4612733889461765832322827737143e-3");
		alpha(11) = to_num<Real_wp>("1.0398987306542755440437675781188e-3");
		alpha(12) = to_num<Real_wp>("6.9793503145264972730105749815623e-4");
		alpha(13) = to_num<Real_wp>("4.1876649698307044270526741525613e-4");
		alpha(14) = to_num<Real_wp>("1.8951076252703209920665889505778e-4");
		break;
	      case 16 :
		alpha(3) = to_num<Real_wp>("1.3888238583520640668318201423914e-2");
		alpha(4) = to_num<Real_wp>("9.0869925208788247685350495375532e-3");
		alpha(5) = to_num<Real_wp>("6.3973854899971635089041432821863e-3");
		alpha(6) = to_num<Real_wp>("4.7310411179653036989287959963169e-3");
		alpha(7) = to_num<Real_wp>("3.6161126443443360119034974050256e-3");
		alpha(8) = to_num<Real_wp>("2.8220980101880568243063662453007e-3");
		alpha(9) = to_num<Real_wp>("2.2273182138378769749685958325591e-3");
		alpha(10) = to_num<Real_wp>("1.764221096346891734468552083213e-3");
		alpha(11) = to_num<Real_wp>("1.3936889436161913857789758238446e-3");
		alpha(12) = to_num<Real_wp>("1.0918430905956956831994903939189e-3");
		alpha(13) = to_num<Real_wp>("8.4305801766163119636026959618868e-4");
		alpha(14) = to_num<Real_wp>("6.3632382003616800017950656771628e-4");
		alpha(15) = to_num<Real_wp>("4.6340324906773887636183942172673e-4");
		alpha(16) = to_num<Real_wp>("3.17895304498812710771458045327e-4");
		alpha(17) = to_num<Real_wp>("1.9472780957554215058215991682633e-4");
		alpha(18) = to_num<Real_wp>("8.9846760431952766725559492128052e-5");
		break;
	      case 20 :
		alpha(3) = to_num<Real_wp>("1.3888636146366966952242378144951e-2");
		alpha(4) = to_num<Real_wp>("9.0893769011309352245836442996147e-3");
		alpha(5) = to_num<Real_wp>("6.4051581115105011522256546412663e-3");
		alpha(6) = to_num<Real_wp>("4.7494432915615761284862392168781e-3");
		alpha(7) = to_num<Real_wp>("3.6514546031026460806805838070469e-3");
		alpha(8) = to_num<Real_wp>("2.8803463023218605738186983845469e-3");
		alpha(9) = to_num<Real_wp>("2.31261510809207179875495285175e-3");
		alpha(10) = to_num<Real_wp>("1.878019530956835190565029323458e-3");
		alpha(11) = to_num<Real_wp>("1.5348175146878549768683045039843e-3");
		alpha(12) = to_num<Real_wp>("1.2572768929119320241371032063264e-3");
		alpha(13) = to_num<Real_wp>("1.0289037135497032574546136228267e-3");
		alpha(14) = to_num<Real_wp>("8.3862603896449192325476396884189e-4");
		alpha(15) = to_num<Real_wp>("6.7862589242346500061188043957262e-4");
		alpha(16) = to_num<Real_wp>("5.4311516007311661273668813628366e-4");
		alpha(17) = to_num<Real_wp>("4.2764440747399966553724638971411e-4");
		alpha(18) = to_num<Real_wp>("3.2870340366737930144150307536314e-4");
		alpha(19) = to_num<Real_wp>("2.4347760195105562676651656377841e-4");
		alpha(20) = to_num<Real_wp>("1.6968823259449201023265722031592e-4");
		alpha(21) = to_num<Real_wp>("1.0547914026011770924615580885405e-4");
		alpha(22) = to_num<Real_wp>("4.9331930601417170622762178557163e-5");
		break;
	      case 30 :
		alpha(3) = to_num<Real_wp>("1.3888845770785215860332454498218e-2");
		alpha(4) = to_num<Real_wp>("9.0906457011430625305289592769815e-3");
		alpha(5) = to_num<Real_wp>("6.4093673033131115048368534549294e-3");
		alpha(6) = to_num<Real_wp>("4.7596830966307095393562313405198e-3");
		alpha(7) = to_num<Real_wp>("3.6718723560534619652440179442121e-3");
		alpha(8) = to_num<Real_wp>("2.9156347372203142371204739154394e-3");
		alpha(9) = to_num<Real_wp>("2.3672466523151249086542215611365e-3");
		alpha(10) = to_num<Real_wp>("1.955462252763584214500763685729e-3");
		alpha(11) = to_num<Real_wp>("1.636995513464608157768565730315e-3");
		alpha(12) = to_num<Real_wp>("1.384426875704909002052524899757e-3");
		alpha(13) = to_num<Real_wp>("1.1798030928797243678008284321331e-3");
		alpha(14) = to_num<Real_wp>("1.0110399962461136303457577503557e-3");
		alpha(15) = to_num<Real_wp>("8.6979657892595763816053027680953e-4");
		alpha(16) = to_num<Real_wp>("7.5016950498601951892870589179709e-4");
		alpha(17) = to_num<Real_wp>("6.4787072014951915118576579343808e-4");
		alpha(18) = to_num<Real_wp>("5.5970174300089676513971991343332e-4");
		alpha(19) = to_num<Real_wp>("4.8321475055349733265841890759245e-4");
		alpha(20) = to_num<Real_wp>("4.1649235831385828213757210156193e-4");
		alpha(21) = to_num<Real_wp>("3.5800278517576758410782894791272e-4");
		alpha(22) = to_num<Real_wp>("3.0650275963804745011057127020358e-4");
		alpha(23) = to_num<Real_wp>("2.6097069235266199543111311567563e-4");
		alpha(24) = to_num<Real_wp>("2.2055921540819962537597885375699e-4");
		alpha(25) = to_num<Real_wp>("1.8456035773998871613777521213961e-4");
		alpha(26) = to_num<Real_wp>("1.5237920763556399115073657297206e-4");
		alpha(27) = to_num<Real_wp>("1.2351348023284425704406770067221e-4");
		alpha(28) = to_num<Real_wp>("9.7537347674718873086453849796118e-5");
		alpha(29) = to_num<Real_wp>("7.4088452628695960266185964409405e-5");
		alpha(30) = to_num<Real_wp>("5.2857367214080719225543484995485e-5");
		alpha(31) = to_num<Real_wp>("3.3578971568260933529877095066048e-5");
		alpha(32) = to_num<Real_wp>("1.6025363134599779765912731942608e-5");
		break;
	      case 40 :
                alpha(3) = to_num<Real_wp>("1.388887695884160160524670699336172069321799209151847558390076358e-2");
                alpha(4) = to_num<Real_wp>("9.090835981668384723915862911909540834635851982266911687642868678e-3");
                alpha(5) = to_num<Real_wp>("6.410008135304827026467687674166398225158358672289231174682415717e-3");
                alpha(6) = to_num<Real_wp>("4.761278685101799314521296837711038406172626877505553799704777512e-3");
                alpha(7) = to_num<Real_wp>("3.6751582879405290342156014358764875340416683356544847694779885e-3");
                alpha(8) = to_num<Real_wp>("2.921556309532232154122693849090507004348861869477471820349269333e-3");
                alpha(9) = to_num<Real_wp>("2.376893209571357043413068527173363879112196223844653231365463336e-3");
                alpha(10) = to_num<Real_wp>("1.96996380293396318641978638026217589781022097364409719126759656e-3");
                alpha(11) = to_num<Real_wp>("1.657399304222413974488100183499242223451432567629323885004425798e-3");
                alpha(12) = to_num<Real_wp>("1.411579568427026401362531374190892737602876735067725926451970088e-3");
                alpha(13) = to_num<Real_wp>("1.214264090163279650153612881678433136657520144661848828508244073e-3");
                alpha(14) = to_num<Real_wp>("1.053042466455823425934735709358784591682628061297216205913375962e-3");
                alpha(15) = to_num<Real_wp>("9.192586319710809710673752015217407919964947423797237653940230832e-4");
                alpha(16) = to_num<Real_wp>("8.06745025280204298846395660242661153889976940635081136574029497e-4");
                alpha(17) = to_num<Real_wp>("7.110222918085132871939110173874224469956745741146589395391894528e-4");
                alpha(18) = to_num<Real_wp>("6.287779342494350996276661520472733183118703626889515630039808342e-4");
                alpha(19) = to_num<Real_wp>("5.575186347515981292654452841605495089409105280014487330964183106e-4");
                alpha(20) = to_num<Real_wp>("4.953343615521453109770034994933325892610530368610434592454778817e-4");
                alpha(21) = to_num<Real_wp>("4.407363226938112618887145523329912556082593762989516220381814183e-4");
                alpha(22) = to_num<Real_wp>("3.925445799995539252973827395984856562221282426908410161337697487e-4");
                alpha(23) = to_num<Real_wp>("3.498094090719338448849895982911970267585997381484917127301610311e-4");
                alpha(24) = to_num<Real_wp>("3.117557163869342143596481727720486752245991058967876816040048538e-4");
                alpha(25) = to_num<Real_wp>("2.777432626189039626569213822253950991465362414458996907720689836e-4");
                alpha(26) = to_num<Real_wp>("2.472377629940307799340509559184727667113158437052062969183878866e-4");
                alpha(27) = to_num<Real_wp>("2.197895199598030213869726536767774469298363469870651069810457096e-4");
                alpha(28) = to_num<Real_wp>("1.950173241144763294322862706323912583833901232301488857441268878e-4");
                alpha(29) = to_num<Real_wp>("1.725960914506440803520225837200185724924923659280779813568892738e-4");
                alpha(30) = to_num<Real_wp>("1.52247196845792821153601441288561226518038021366410354543438069e-4");
                alpha(31) = to_num<Real_wp>("1.337307919695357162556810336557712097869791488267083369366039072e-4");
                alpha(32) = to_num<Real_wp>("1.168396140964601601266405638559802785412297978708282184008482903e-4");
                alpha(33) = to_num<Real_wp>("1.013939377345364383207389466117564584084683774001887371420220646e-4");
                alpha(34) = to_num<Real_wp>("8.723741850260908586496759897924997813076690294771118336982671025e-5");
                alpha(35) = to_num<Real_wp>("7.42336448932630930014687973433730852166040254247020856783267388e-5");
                alpha(36) = to_num<Real_wp>("6.226325926958020680637849495092928082215994388842027880829750879e-5");
                alpha(37) = to_num<Real_wp>("5.122154167013999765559282192981667125393146933815394720738295483e-5");
                alpha(38) = to_num<Real_wp>("4.101637324243933709391955539894840233306699291018460682198659515e-5");
                alpha(39) = to_num<Real_wp>("3.156651329039809767294102760758433094787369847225276162694107246e-5");
                alpha(40) = to_num<Real_wp>("2.280013688315954302273946885426192513762259017620297217801433601e-5");
                alpha(41) = to_num<Real_wp>("1.465358995774604045762650975645233742324647569226894229694667702e-5");
                alpha(42) = to_num<Real_wp>("7.070326673884508889708332448264486576566226030156565276463681173e-6");
		break;
	      case 50 :
                alpha(3) = to_num<Real_wp>("1.388888454059991309155979904279535407649062306051872845398810604e-2");
                alpha(4) = to_num<Real_wp>("9.090882398876619444338237645633195429524589465235831709428205408e-3");
                alpha(5) = to_num<Real_wp>("6.410165480013272939945912962171401991126950789368141939850486927e-3");
                alpha(6) = to_num<Real_wp>("4.761674368190432506110877348168900727868316325339199621593491647e-3");
                alpha(7) = to_num<Real_wp>("3.675984493604248800345031706159483175475608724218684546247907563e-3");
                alpha(8) = to_num<Real_wp>("2.923072420958827403522934193103033525937509672851971306521444638e-3");
                alpha(9) = to_num<Real_wp>("2.379419671518877833559112382106157619581896500721834672802970407e-3");
                alpha(10) = to_num<Real_wp>("1.973866833970715561953546511864711379360534675773130175379344528e-3");
                alpha(11) = to_num<Real_wp>("1.663067339632296253201419362576688085290471983545549594408231903e-3");
                alpha(12) = to_num<Real_wp>("1.419394048375554100654818830444236295520276987089102414773297383e-3");
                alpha(13) = to_num<Real_wp>("1.224568696982017699672934121141008693612498794649105782722031408e-3");
                alpha(14) = to_num<Real_wp>("1.066115586643734558427773634834155473331177306635762041303868261e-3");
                alpha(15) = to_num<Real_wp>("9.352933186062267739580434146894497351185052979381982528321425458e-4");
                alpha(16) = to_num<Real_wp>("8.258391448162422070333308454115580018356248467489910548052394648e-4");
                alpha(17) = to_num<Real_wp>("7.331793650550234418517984776361503491259572396953289643604155851e-4");
                alpha(18) = to_num<Real_wp>("6.539171268052484640802278976271572540143719455568398502999666678e-4");
                alpha(19) = to_num<Real_wp>("5.854909908153148522787911786323639879434236471330099054609002243e-4");
                alpha(20) = to_num<Real_wp>("5.259419387204524193706937768300097801160218891622601482423219862e-4");
                alpha(21) = to_num<Real_wp>("4.737512535912488822220335294946659654033947748738041017671180264e-4");
                alpha(22) = to_num<Real_wp>("4.277259313873010515971567565614151745406071912164593901519746807e-4");
                alpha(23) = to_num<Real_wp>("3.869166767164695576055894287384079594770059187567908059278924299e-4");
                alpha(24) = to_num<Real_wp>("3.505586285712711379960835412197931038170658416805806204604625349e-4");
                alpha(25) = to_num<Real_wp>("3.180281463794318021454886602892122798737665461695209770851223579e-4");
                alpha(26) = to_num<Real_wp>("2.888110454072443945079215564474479207565578776740932903110289901e-4");
                alpha(27) = to_num<Real_wp>("2.624790452300000750656511057458995428019713706421755730338721076e-4");
                alpha(28) = to_num<Real_wp>("2.386721385438359876529273236162069150600358842417306343392590402e-4");
                alpha(29) = to_num<Real_wp>("2.170852481738490932500549612756927679810145478257422870971812815e-4");
                alpha(30) = to_num<Real_wp>("1.97458007940433885750264315663200713496771227424370917449865184e-4");
                alpha(31) = to_num<Real_wp>("1.795668358695294360813813933466339299398510456064933119058874162e-4");
                alpha(32) = to_num<Real_wp>("1.632187050388818249108115895558481785282032305172067702187952218e-4");
                alpha(33) = to_num<Real_wp>("1.48246185477596954060274451433515463243099989878129892532195032e-4");
                alpha(34) = to_num<Real_wp>("1.345034496001275132600497339745763079761615286960461192548996443e-4");
                alpha(35) = to_num<Real_wp>("1.218630178496218901383222700084369901526700496596004134836492297e-4");
                alpha(36) = to_num<Real_wp>("1.102130807884640440010398535816580072234117922986551646999551234e-4");
                alpha(37) = to_num<Real_wp>("9.945527613552347206178110046853859770769786374440265711908635724e-5");
                alpha(38) = to_num<Real_wp>("8.950282940579326426806582179159808538096029848188615183036085603e-5");
                alpha(39) = to_num<Real_wp>("8.027898850684107743732207100410052586876646041771824293210782078e-5");
                alpha(40) = to_num<Real_wp>("7.171569842913061143354230836991390815153193781545550128518393736e-5");
                alpha(41) = to_num<Real_wp>("6.375247379417008284980146140411712968114774902507646414272360832e-5");
                alpha(42) = to_num<Real_wp>("5.633543571150063062433369258000435814565178252078264788692856405e-5");
                alpha(43) = to_num<Real_wp>("4.941648598255556603738608726868981955320109746662522021439012171e-5");
                alpha(44) = to_num<Real_wp>("4.295259675815923661095749179045566197402464701134412898187812196e-5");
                alpha(45) = to_num<Real_wp>("3.690519771226749693419415480927173850570152278674461878701609282e-5");
                alpha(46) = to_num<Real_wp>("3.123964592240020959703982198277117900442506740588685377883831665e-5");
                alpha(47) = to_num<Real_wp>("2.592476614957364821276375387087086257633107454299217172773644393e-5");
                alpha(48) = to_num<Real_wp>("2.093245123322462575660579738234524961271137416447527004747270118e-5");
                alpha(49) = to_num<Real_wp>("1.623731396617938919293830391221151575428363296824521276675106538e-5");
                alpha(50) = to_num<Real_wp>("1.181638317022868257069586069724437412663160538714881091306148309e-5");
                alpha(51) = to_num<Real_wp>("7.648837813925198552501916830119745899033832070506403336803854819e-6");
                alpha(52) = to_num<Real_wp>("3.715773946429793399247838789107899995319988934227815379918492425e-6");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}	       
	      }
	  }
	else if (order == 10)
	  {
	    alpha(0) = Real_wp(1) / 6;
	    alpha(1) = Real_wp(1) / 20;
	    alpha(2) = Real_wp(1) / 42;
	    alpha(3) = Real_wp(1) / 72;
	    
	    switch (m)
	      {
	      case 1 :
		alpha(4) = to_num<Real_wp>("8.9487953793969117109360032567148e-3");
		break;
	      case 2 :
		alpha(4) = to_num<Real_wp>("9.0770056937720791822140631304913e-3");
		alpha(5) = to_num<Real_wp>("5.7245403968713726005592787448706e-3");
		break;
	      case 3 :
		alpha(4) = to_num<Real_wp>("9.0882867216563711904053685904659e-3");
		alpha(5) = to_num<Real_wp>("6.2764396383854641620591837659885e-3");
		alpha(6) = to_num<Real_wp>("3.6119671385414798004013998240676e-3");
		break;
	      case 4 :
		alpha(4) = to_num<Real_wp>("9.0902370171045884709168537732197e-3");
		alpha(5) = to_num<Real_wp>("6.3753270093796891352712756346268e-3");
		alpha(6) = to_num<Real_wp>("4.4250110429503202339569763211363e-3");
		alpha(7) = to_num<Real_wp>("2.3411085538962653495788454648497e-3");
		break;
	      case 6 :
		alpha(4) = to_num<Real_wp>("9.0908307152701187049282739153277e-3");
		alpha(5) = to_num<Real_wp>("6.4061084431019129689930573128023e-3");
		alpha(6) = to_num<Real_wp>("4.7172407606122558237025453119287e-3");
		alpha(7) = to_num<Real_wp>("3.4550331294953030289396673656237e-3");
		alpha(8) = to_num<Real_wp>("2.2766195870925816880408519441128e-3");
		alpha(9) = to_num<Real_wp>("1.0918297102652750106184104573804e-3");
		break;
	      case 8 :
		alpha(4) = to_num<Real_wp>("9.0908943942613077306456466911964e-3");
		alpha(5) = to_num<Real_wp>("6.4094719260678431057861072244046e-3");
		alpha(6) = to_num<Real_wp>("4.753028962334446153069122987192e-3");
		alpha(7) = to_num<Real_wp>("3.6277053528658318826485503686442e-3");
		alpha(8) = to_num<Real_wp>("2.7573499570294287287598263347376e-3");
		alpha(9) = to_num<Real_wp>("1.9822141700574696348011494658325e-3");
		alpha(10) = to_num<Real_wp>("1.2493064770177848466860078551778e-3");
		alpha(11) = to_num<Real_wp>("5.7847571555118805632762696678576e-4");
		break;
	      case 12 :
		alpha(4) = to_num<Real_wp>("9.0909079243372493302471220429473e-3");
		alpha(5) = to_num<Real_wp>("6.4101936602010178725809194487835e-3");
		alpha(6) = to_num<Real_wp>("4.761163654483022801872920073166e-3");
		alpha(7) = to_num<Real_wp>("3.6720272828680719763570187487212e-3");
		alpha(8) = to_num<Real_wp>("2.9065278632869116837161868628551e-3");
		alpha(9) = to_num<Real_wp>("2.3305386786210983224733838886137e-3");
		alpha(10) = to_num<Real_wp>("1.8619447136129781207459745357933e-3");
		alpha(11) = to_num<Real_wp>("1.4537191988571819857123279328166e-3");
		alpha(12) = to_num<Real_wp>("1.0860005487046372271933124810253e-3");
		alpha(13) = to_num<Real_wp>("7.555966671727067838400752780563e-4");
		alpha(14) = to_num<Real_wp>("4.6457622528297486448008371803363e-4");
		alpha(15) = to_num<Real_wp>("2.135102829085184689153089843468e-4");
		break;
	      case 16 :
		alpha(4) = to_num<Real_wp>("9.0909089159431088528101481096065e-3");
		alpha(5) = to_num<Real_wp>("6.4102469675440318254393357532294e-3");
		alpha(6) = to_num<Real_wp>("4.761791208578347460115407389561e-3");
		alpha(7) = to_num<Real_wp>("3.6757650337904358638131353343586e-3");
		alpha(8) = to_num<Real_wp>("2.9210442015831499367396579007086e-3");
		alpha(9) = to_num<Real_wp>("2.3717673489225248001403610690835e-3");
		alpha(10) = to_num<Real_wp>("1.953145413031450473833585527559e-3");
		alpha(11) = to_num<Real_wp>("1.6178676314021742810904600723697e-3");
		alpha(12) = to_num<Real_wp>("1.3358597926787645167194583259827e-3");
		alpha(13) = to_num<Real_wp>("1.0895059623843965517583010137701e-3");
		alpha(14) = to_num<Real_wp>("8.7011191475581348060094308393588e-4");
		alpha(15) = to_num<Real_wp>("6.7425444344124083106702234252641e-4");
		alpha(16) = to_num<Real_wp>("5.0062341868213795617212953501734e-4");
		alpha(17) = to_num<Real_wp>("3.4812394456359215693091055685304e-4");
		alpha(18) = to_num<Real_wp>("2.1520489950347744638715786879238e-4");
		alpha(19) = to_num<Real_wp>("9.9885290707325320078209138314919e-5");
		break;
	      case 20:
		alpha(4) = to_num<Real_wp>("9.0909090524997271738327826165832e-3");
		alpha(5) = to_num<Real_wp>("6.4102543337055427541980708409505e-3");
		alpha(6) = to_num<Real_wp>("4.7618795538441633074365252680546e-3");
		alpha(7) = to_num<Real_wp>("3.6763110605957059793874313906548e-3");
		alpha(8) = to_num<Real_wp>("2.9232942560263499495347540767675e-3");
		alpha(9) = to_num<Real_wp>("2.378725953900710448623127224021e-3");
		alpha(10) = to_num<Real_wp>("1.9703617330323309100533802537043e-3");
		alpha(11) = to_num<Real_wp>("1.6532915333217447582949874249047e-3");
		alpha(12) = to_num<Real_wp>("1.398169323005611423966764552438e-3");
		alpha(13) = to_num<Real_wp>("1.1854389790710817489821283805356e-3");
		alpha(14) = to_num<Real_wp>("1.0023812388170133488179321004504e-3");
		alpha(15) = to_num<Real_wp>("8.4126217265438977577009441082373e-4");
		alpha(16) = to_num<Real_wp>("6.9776945082323528400030954166713e-4");
		alpha(17) = to_num<Real_wp>("5.6960220962767122358207158951823e-4");
		alpha(18) = to_num<Real_wp>("4.5538473870882362559005116064199e-4");
		alpha(19) = to_num<Real_wp>("3.5401460747772348755664040795971e-4");
		alpha(20) = to_num<Real_wp>("2.6438398420991984916978747837108e-4");
		alpha(21) = to_num<Real_wp>("1.8532651519057011517418771799004e-4");
		alpha(22) = to_num<Real_wp>("1.1565829562828747121224946127055e-4");
		alpha(23) = to_num<Real_wp>("5.4235874128731300367789556190715e-5");
		break;
	      case 30 :
		alpha(4) = to_num<Real_wp>("9.0909090886449049500217768788615e-3");
		alpha(5) = to_num<Real_wp>("6.4102562876019102325809106057472e-3");
		alpha(6) = to_num<Real_wp>("4.761903257081513236588000587873e-3");
		alpha(7) = to_num<Real_wp>("3.676460869345854578705750132732e-3");
		alpha(8) = to_num<Real_wp>("2.9239337127188038892307497459784e-3");
		alpha(9) = to_num<Real_wp>("2.3808061209811899061925956193738e-3");
		alpha(10) = to_num<Real_wp>("1.9758720562696231223751624142823e-3");
		alpha(11) = to_num<Real_wp>("1.665662868702428153760897367514e-3");
		alpha(12) = to_num<Real_wp>("1.4223354414549132453474286151729e-3");
		alpha(13) = to_num<Real_wp>("1.227303319926643425936047386104e-3");
		alpha(14) = to_num<Real_wp>("1.0677229898117236753827720180519e-3");
		alpha(15) = to_num<Real_wp>("9.3449051226218499987461862331364e-4");
		alpha(16) = to_num<Real_wp>("8.2106898145214918341509563396151e-4");
		alpha(17) = to_num<Real_wp>("7.227724528827584330253016812976e-4");
		alpha(18) = to_num<Real_wp>("6.3628870663065578844251225619811e-4");
		alpha(19) = to_num<Real_wp>("5.5932047182966527639406339840173e-4");
		alpha(20) = to_num<Real_wp>("4.902922064064930563517361200545e-4");
		alpha(21) = to_num<Real_wp>("4.2811043141266661737937204101815e-4");
		alpha(22) = to_num<Real_wp>("3.7198061913810806676218221641058e-4");
		alpha(23) = to_num<Real_wp>("3.2128102946646357108011943708105e-4");
		alpha(24) = to_num<Real_wp>("2.7548576512112104742600134281671e-4");
		alpha(25) = to_num<Real_wp>("2.3412400567905214870888005601324e-4");
		alpha(26) = to_num<Real_wp>("1.9676198468118532326351066457029e-4");
		alpha(27) = to_num<Real_wp>("1.6299709316674852915934714878463e-4");
		alpha(28) = to_num<Real_wp>("1.3245718172154248153577855009624e-4");
		alpha(29) = to_num<Real_wp>("1.0480124766269375639768560560659e-4");
		alpha(30) = to_num<Real_wp>("7.9719786397340367280606835188404e-5");
		alpha(31) = to_num<Real_wp>("5.693427007553190534194292235442e-5");
		alpha(32) = to_num<Real_wp>("3.6195779860683920915057937929702e-5");
		alpha(33) = to_num<Real_wp>("1.7283023664252610691428078176689e-5");
		break;
	      case 40 :
                alpha(4) = to_num<Real_wp>("9.090909090618301343189078912151172689759826675775209923407339437e-3");
                alpha(5) = to_num<Real_wp>("6.410256394490991338088032081599415982621587809571574899718918609e-3");
                alpha(6) = to_num<Real_wp>("4.761904567644927215238247498271365500802888223508903254744414227e-3");
                alpha(7) = to_num<Real_wp>("3.676469323296412025717048053552554036463892276267041365984132822e-3");
                alpha(8) = to_num<Real_wp>("2.923970955378915800386055519712903812338966922196128427246582363e-3");
                alpha(9) = to_num<Real_wp>("2.380932773073584926181284057034173992367415474997238591558309296e-3");
                alpha(10) = to_num<Real_wp>("1.976228023441466647462493224225189223554482250382207003178611857e-3");
                alpha(11) = to_num<Real_wp>("1.666525060927011699079965340650229815481730395465876158039695114e-3");
                alpha(12) = to_num<Real_wp>("1.4241849275210626726148762833574718865264689533224826770865259e-3");
                alpha(13) = to_num<Real_wp>("1.230883104735152849043905031223743780407838215226372461125160938e-3");
                alpha(14) = to_num<Real_wp>("1.07405905090666954848547938243817918174075034152233654796680104e-3");
                alpha(15) = to_num<Real_wp>("9.448496271861447350385957986018836969571581794318864997044055698e-4");
                alpha(16) = to_num<Real_wp>("8.368425556008883879167208495552847409162253944339633371488338905e-4");
                alpha(17) = to_num<Real_wp>("7.453032935098330178570089938159552917021936405889654213353957342e-4");
                alpha(18) = to_num<Real_wp>("6.666843727960812746110982962330243425578414575900780216898238009e-4");
                alpha(19) = to_num<Real_wp>("5.983061367108831271652557659922612072227089548073060637986209805e-4");
                alpha(20) = to_num<Real_wp>("5.381411892845878160920125491297713507136603885157166572994986521e-4");
                alpha(21) = to_num<Real_wp>("4.846606148823068601428424666585505255975007513514961861957873842e-4");
                alpha(22) = to_num<Real_wp>("4.367169622050421888827872446481491238295674525470870275344011412e-4");
                alpha(23) = to_num<Real_wp>("3.934506430199105969561506623990229309544733426579886949677963552e-4");
                alpha(24) = to_num<Real_wp>("3.542138411066083880581373773553460119197014060313451199866137496e-4");
                alpha(25) = to_num<Real_wp>("3.185098106629847154482080194847114961506899022647167861674974135e-4");
                alpha(26) = to_num<Real_wp>("2.859465050532332014489419067256233272842147336058777271876660316e-4");
                alpha(27) = to_num<Real_wp>("2.562030918838027778820316917913460071744178107293281473118231172e-4");
                alpha(28) = to_num<Real_wp>("2.29007206581856217209572586989821049614464838463196798494712743e-4");
                alpha(29) = to_num<Real_wp>("2.041204110477757420783325651320964564442584759730266874408323616e-4");
                alpha(30) = to_num<Real_wp>("1.81329403058693489263388258084544375582966524141622565498015075e-4");
                alpha(31) = to_num<Real_wp>("1.604409299087960452896119122621994643492039429113430647624604822e-4");
                alpha(32) = to_num<Real_wp>("1.412788916608672291723533301609968849940049435239532161952332841e-4");
                alpha(33) = to_num<Real_wp>("1.236826212852983488655662824032341663268142373446882251788098416e-4");
                alpha(34) = to_num<Real_wp>("1.075057258370143810122689886092512113121725460942367690452566866e-4");
                alpha(35) = to_num<Real_wp>("9.261515020232316512748207294993403889561335967885405477051532892e-5");
                alpha(36) = to_num<Real_wp>("7.889030024525778993265615334366385358112902345259090527064118968e-5");
                alpha(37) = to_num<Real_wp>("6.622216285380423220264099439501460642877649581424615381841050379e-5");
                alpha(38) = to_num<Real_wp>("5.451241216180427739284777879312279356624923117951244986153691577e-5");
                alpha(39) = to_num<Real_wp>("4.367251377705980706650733709908382445747804983456416066280919057e-5");
                alpha(40) = to_num<Real_wp>("3.362284548837555909446526510296038456764383297474127582306937721e-5");
                alpha(41) = to_num<Real_wp>("2.429185186666753534769324279826114331063184448416564694178228126e-5");
                alpha(42) = to_num<Real_wp>("1.561524608519041243160116815814572537228164403754998397310368537e-5");
                alpha(43) = to_num<Real_wp>("7.535267583986648545040087511758183501006964642275423600092988109e-6");
		break;
	      case 50 :
                alpha(4) = to_num<Real_wp>("9.090909090850944748947832053126606158063983107413948124439140668e-3");
                alpha(5) = to_num<Real_wp>("6.410256407102627723148187850664372308898848935464722998040863433e-3");
                alpha(6) = to_num<Real_wp>("4.761904722957470554512687263678459306874403713640203302715954635e-3");
                alpha(7) = to_num<Real_wp>("3.67647033355927625279139852988624633300671739540563647093702113e-3");
                alpha(8) = to_num<Real_wp>("2.923975462855588549212265183345451326018765855020005243531622241e-3");
                alpha(9) = to_num<Real_wp>("2.38094837363757788316184880920505833968699140221421490866216038e-3");
                alpha(10) = to_num<Real_wp>("1.976272895075786466380743033665225125440927866199199896645352723e-3");
                alpha(11) = to_num<Real_wp>("1.666636986086980115334626170074278873251767355350491555560193614e-3");
                alpha(12) = to_num<Real_wp>("1.424433934732583424309605513788416548326748919188705129645478258e-3");
                alpha(13) = to_num<Real_wp>("1.231386893861986581388062904347982388214134550802331506689424955e-3");
                alpha(14) = to_num<Real_wp>("1.074998895072159294445113498842552180850810605564160639617730945e-3");
                alpha(15) = to_num<Real_wp>("9.464829034751085609055789689509249357365275587147880036208383994e-4");
                alpha(16) = to_num<Real_wp>("8.395070681986809758776206037539421412531717830752841658081271621e-4");
                alpha(17) = to_num<Real_wp>("7.494088769052271938473873608515315166785335162712346574355858772e-4");
                alpha(18) = to_num<Real_wp>("6.72689426867267278703738980436184075937264916551386732199159748e-4");
                alpha(19) = to_num<Real_wp>("6.066805613395911921030452055763444352678097752898478389456186961e-4");
                alpha(20) = to_num<Real_wp>("5.493214941623722974739493494802946586747366886298584303238456691e-4");
                alpha(21) = to_num<Real_wp>("4.990059960293087492313202552610369808870058296512960099805649e-4");
                alpha(22) = to_num<Real_wp>("4.544753844490093409491154766193904966877761906360955683472341846e-4");
                alpha(23) = to_num<Real_wp>("4.147412138142699983297653351216762957402811875030829001224427565e-4");
                alpha(24) = to_num<Real_wp>("3.790272081298609643647165478180597531796253836251214274784420046e-4");
                alpha(25) = to_num<Real_wp>("3.467238719254267709470781301283678276279820449332656792054971569e-4");
                alpha(26) = to_num<Real_wp>("3.173519217680673110214969349944650296568330836379930580524214576e-4");
                alpha(27) = to_num<Real_wp>("2.905324505047312679727224982128417420806056977419195151941955968e-4");
                alpha(28) = to_num<Real_wp>("2.659627377704435649619138909220305224530251995245928889385455961e-4");
                alpha(29) = to_num<Real_wp>("2.433970501321952675582825367809089225291187836991157036019914685e-4");
                alpha(30) = to_num<Real_wp>("2.226318577661380782543789304310601453071255644191961703641405909e-4");
                alpha(31) = to_num<Real_wp>("2.034948355628991956372484634484373692689459868575777079570273874e-4");
                alpha(32) = to_num<Real_wp>("1.858369537073643695248991744870501178906934731899063880099804187e-4");
                alpha(33) = to_num<Real_wp>("1.695269594542721977311787190277533295497915463978658185785824712e-4");
                alpha(34) = to_num<Real_wp>("1.544476143282139156039880524048194763155376669590338723211602709e-4");
                alpha(35) = to_num<Real_wp>("1.404931567348702796901494281682362735682239522313253481440644097e-4");
                alpha(36) = to_num<Real_wp>("1.275675805408895108944306436136165545377477192142816416251515202e-4");
                alpha(37) = to_num<Real_wp>("1.155834336784291971053208379923794897313166843973430416053548621e-4");
                alpha(38) = to_num<Real_wp>("1.044609353178895225461882139520836486149514817950609050848780328e-4");
                alpha(39) = to_num<Real_wp>("9.412728203793097419595159541155153544782629408186778045667011437e-5");
                alpha(40) = to_num<Real_wp>("8.451606430217137369676502513446939452729414304168190310338744602e-5");
                alpha(41) = to_num<Real_wp>("7.556674837511039976451949027821142205401030921377149633995803572e-5");
                alpha(42) = to_num<Real_wp>("6.722419999412566224664416470303660568937032174564669915511554299e-5");
                alpha(43) = to_num<Real_wp>("5.943823857234959468538495186458214260233350531043763998646961244e-5");
                alpha(44) = to_num<Real_wp>("5.216321749433766827508012558714633036147204955071447399826699406e-5");
                alpha(45) = to_num<Real_wp>("4.535762937613114414811680320925434013629837376272091778572712748e-5");
                alpha(46) = to_num<Real_wp>("3.898373646962071356116309122640245679632535544675503318438169258e-5");
                alpha(47) = to_num<Real_wp>("3.300722663736985323906666628228047071093011769850405249997757152e-5");
                alpha(48) = to_num<Real_wp>("2.73968950744124620072033630355892934364127042140501727931574536e-5");
                alpha(49) = to_num<Real_wp>("2.212435152904083924341338416901502013799193296267395166048753292e-5");
                alpha(50) = to_num<Real_wp>("1.716375234316818501971996387158894255041764527077565963995038303e-5");
                alpha(51) = to_num<Real_wp>("1.249155627141042067739187306834432694333119370528384009100859114e-5");
                alpha(52) = to_num<Real_wp>("8.086302772452241674216226043777879032943690547005116025686863985e-6");
                alpha(53) = to_num<Real_wp>("3.928411295903893530637443855120873172111325198884330702931079338e-6");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}
	      }
	  }
	else if (order == 12)
	  {
	    alpha(0) = Real_wp(1) / 6;
	    alpha(1) = Real_wp(1) / 20;
	    alpha(2) = Real_wp(1) / 42;
	    alpha(3) = Real_wp(1) / 72;
	    alpha(4) = Real_wp(1) / 110;

	    switch (m)
	      {
	      case 1 :
		alpha(5) = to_num<Real_wp>("5.7847497421331866358935581782571e-3");
		break;
	      case 2 :
		alpha(5) = to_num<Real_wp>("6.2925920631742580957112168596888e-3");
		alpha(6) = to_num<Real_wp>("3.6372356636689892039118463761479e-3");
		break;
	      case 3 :
		alpha(5) = to_num<Real_wp>("6.380151698636695182401884833588e-3");
		alpha(6) = to_num<Real_wp>("4.4371593819682448714239715809101e-3");
		alpha(7) = to_num<Real_wp>("2.3524421852418974185425948927227e-3");
		break;
	      case 4 :
		alpha(5) = to_num<Real_wp>("6.4007351918368898203185100176413e-3");
		alpha(6) = to_num<Real_wp>("4.6517382613361953954912148217469e-3");
		alpha(7) = to_num<Real_wp>("3.1608463394914218139651201435745e-3");
		alpha(8) = to_num<Real_wp>("1.5762708577383673317199104950308e-3");
		break;
	      case 6 :
		alpha(5) = to_num<Real_wp>("6.408807755637427330436533080735e-3");
		alpha(6) = to_num<Real_wp>("4.7438455924366508263452702849322e-3");
		alpha(7) = to_num<Real_wp>("3.5781979531171383563829568631203e-3");
		alpha(8) = to_num<Real_wp>("2.605429745153102500242687094762e-3");
		alpha(9) = to_num<Real_wp>("1.6740968264337878241778053288548e-3");
		alpha(10) = to_num<Real_wp>("7.8507464176719977997781364030904e-4");
		break;
	      case 8 :
		alpha(5) = to_num<Real_wp>("6.4099358549437697384355629381837e-3");
		alpha(6) = to_num<Real_wp>("4.7577573954256927755643947033966e-3");
		alpha(7) = to_num<Real_wp>("3.6521558039329316419714352546247e-3");
		alpha(8) = to_num<Real_wp>("2.835456818297022692268795021029e-3");
		alpha(9) = to_num<Real_wp>("2.153417436099474390569933290856e-3");
		alpha(10) = to_num<Real_wp>("1.5297102860088721554490431974251e-3");
		alpha(11) = to_num<Real_wp>("9.524093020660170789558662069906e-4");
		alpha(12) = to_num<Real_wp>("4.3826660868530519408969668531571e-4");
		break;
	      case 12 :
		alpha(5) = to_num<Real_wp>("6.4102254745545576841912719248559e-3");
		alpha(6) = to_num<Real_wp>("4.7614902592759905095702639280205e-3");
		alpha(7) = to_num<Real_wp>("3.6738681817541435852582996758924e-3");
		alpha(8) = to_num<Real_wp>("2.9134133673890815819369295702955e-3");
		alpha(9) = to_num<Real_wp>("2.3493842759901509100212377186318e-3");
		alpha(10) = to_num<Real_wp>("1.9019178141507144033124816458224e-3");
		alpha(11) = to_num<Real_wp>("1.5223449702529469394590798629016e-3");
		alpha(12) = to_num<Real_wp>("1.1854590606315241325751264629338e-3");
		alpha(13) = to_num<Real_wp>("8.8253958508962137929772500270624e-4");
		alpha(14) = to_num<Real_wp>("6.1284458321125502466395715469435e-4");
		alpha(15) = to_num<Real_wp>("3.7688988233967289375513123998507e-4");
		alpha(16) = to_num<Real_wp>("1.735617510536722285729689796945e-4");
		break;
	      case 16 :
		alpha(5) = to_num<Real_wp>("6.410251221847467568524169379705e-3");
		alpha(6) = to_num<Real_wp>("4.7618341387284295934429613853698e-3");
		alpha(7) = to_num<Real_wp>("3.6760136556524036533822268290093e-3");
		alpha(8) = to_num<Real_wp>("2.9220321757422055850305203084402e-3");
		alpha(9) = to_num<Real_wp>("2.3747374061002720661202101554484e-3");
		alpha(10) = to_num<Real_wp>("1.9602912804912865635470121525489e-3");
		alpha(11) = to_num<Real_wp>("1.6321281689271859966990274325229e-3");
		alpha(12) = to_num<Real_wp>("1.3601105924012392528769917090859e-3");
		alpha(13) = to_num<Real_wp>("1.1255334477536699190099673663312e-3");
		alpha(14) = to_num<Real_wp>("9.1806450394018846658130001736829e-4");
		alpha(15) = to_num<Real_wp>("7.3287345960545362993249446537273e-4");
		alpha(16) = to_num<Real_wp>("5.6792801280284953832390710535732e-4");
		alpha(17) = to_num<Real_wp>("4.2205590198438395215609394862461e-4");
		alpha(18) = to_num<Real_wp>("2.9399529158207406812687741111595e-4");
		alpha(19) = to_num<Real_wp>("1.821668976069834901888409833588e-4");
		alpha(20) = to_num<Real_wp>("8.4778577305326917484215336618087e-5");
		break;
	      case 20 :
		alpha(5) = to_num<Real_wp>("6.4102551858199305362806619905462e-3");
		alpha(6) = to_num<Real_wp>("4.7618879518297351055515351505832e-3");
		alpha(7) = to_num<Real_wp>("3.6763600628684070405842104421884e-3");
		alpha(8) = to_num<Real_wp>("2.9234944059346705288163950003017e-3");
		alpha(9) = to_num<Real_wp>("2.3793557656834260539958401318495e-3");
		alpha(10) = to_num<Real_wp>("1.9719779578083330674210795919304e-3");
		alpha(11) = to_num<Real_wp>("1.6567963994886236284380917116823e-3");
		alpha(12) = to_num<Real_wp>("1.404747983987460693826623529034e-3");
		alpha(13) = to_num<Real_wp>("1.1963259178495698834280753227131e-3");
		alpha(14) = to_num<Real_wp>("1.0185315684589281058902934634766e-3");
		alpha(15) = to_num<Real_wp>("8.6309696135070663225345837479102e-4");
		alpha(16) = to_num<Real_wp>("7.2512555641672980421233759631122e-4");
		alpha(17) = to_num<Real_wp>("6.0187646067859029692733811482846e-4");
		alpha(18) = to_num<Real_wp>("4.9175815887374934436641034926507e-4");
		alpha(19) = to_num<Real_wp>("3.936459965443840673980862154827e-4");
		alpha(20) = to_num<Real_wp>("3.0652577818590496661232183478118e-4");
		alpha(21) = to_num<Real_wp>("2.2936949001080401201099268853576e-4");
		alpha(22) = to_num<Real_wp>("1.611317515350306607255572989572e-4");
		alpha(23) = to_num<Real_wp>("1.0078756121175010630860866964639e-4");
		alpha(24) = to_num<Real_wp>("4.7370862948428357364018606531735e-5");
		break;
	      case 30 :
                alpha(5) = to_num<Real_wp>("6.410256330220243208090910546833302408219747959677117640950335529e-3");
                alpha(6) = to_num<Real_wp>("4.761903652075985975270475782767941435612128326468477942966781811e-3");
                alpha(7) = to_num<Real_wp>("3.676463154878692274680397211706177744651511073169156760353494475e-3");
                alpha(8) = to_num<Real_wp>("2.923943248788025469663147601346424193404850614630087908099537122e-3");
                alpha(9) = to_num<Real_wp>("2.380837481708255300326310013018412898575002089996307894882959436e-3");
                alpha(10) = to_num<Real_wp>("1.975958023763252792840838899427735822412414159405648542784942302e-3");
                alpha(11) = to_num<Real_wp>("1.665866580396155700085374818370560611702600387886298102565585359e-3");
                alpha(12) = to_num<Real_wp>("1.422763060008473635196056317851272433004536967343148284131082578e-3");
                alpha(13) = to_num<Real_wp>("1.228112129836488893166787546337798588024529248050376890338690166e-3");
                alpha(14) = to_num<Real_wp>("1.069118591968516984221429591433752716659353136009629518970626117e-3");
                alpha(15) = to_num<Real_wp>("9.367086485492972909911430288758534503049136594782574529789134666e-4");
                alpha(16) = to_num<Real_wp>("8.243431363187731848017022702334527939867577913822892316759101899e-4");
                alpha(17) = to_num<Real_wp>("7.27295563460982765123509488610529977502826428150354122175432441e-4");
                alpha(18) = to_num<Real_wp>("6.421815140721896856181641009317530860128256673167413264965469949e-4");
                alpha(19) = to_num<Real_wp>("5.666168222728210208514092131665094908730319998125636259383103372e-4");
                alpha(20) = to_num<Real_wp>("4.989440901834857886732680583667308070448034967231507854763910596e-4");
                alpha(21) = to_num<Real_wp>("4.380077031278544184376806192487600553704274581922218412462433766e-4");
                alpha(22) = to_num<Real_wp>("3.829765855510274801516945151806341426998984540050648411168050732e-4");
                alpha(23) = to_num<Real_wp>("3.332161164779480034149647222349647563458250814841937778524276448e-4");
                alpha(24) = to_num<Real_wp>("2.882049968211488255372512387787020476012276037361612113325886649e-4");
                alpha(25) = to_num<Real_wp>("2.47487451840417750854367892654082040059261808858728746209099353e-4");
                alpha(26) = to_num<Real_wp>("2.106493810915520864843473818770344358133112668229512458787775995e-4");
                alpha(27) = to_num<Real_wp>("1.773085114814384142970649807386196850597766239941041845724015299e-4");
                alpha(28) = to_num<Real_wp>("1.471114644137449128995836320908670466130861525250313167960799322e-4");
                alpha(29) = to_num<Real_wp>("1.197334511312244835381123700549767909101397857802602199808751403e-4");
                alpha(30) = to_num<Real_wp>("9.487839904536567606747107693212912299268120948727188746888541546e-5");
                alpha(31) = to_num<Real_wp>("7.227861111697812821644383683746362843263697344375911544445870652e-5");
                alpha(32) = to_num<Real_wp>("5.169375176576182889703555998553830218026683045285800546520133108e-5");
                alpha(33) = to_num<Real_wp>("3.290925882220152374612411128434043686141035888312204017159692954e-5");
                alpha(34) = to_num<Real_wp>("1.573437638753801726154625454517815358498758908167173629683776004e-5");
		break;
	      case 40 :
                alpha(5) = to_num<Real_wp>("6.410256399398110903361640994245513583912021817877942491025372488e-3");
                alpha(6) = to_num<Real_wp>("4.761904610703437877570303768920881978483094122304867349808631046e-3");
                alpha(7) = to_num<Real_wp>("3.676469567702216515998660511357345591001993888741797375026615758e-3");
                alpha(8) = to_num<Real_wp>("2.923971974897010799063722246524639381377966820629475439986368687e-3");
                alpha(9) = to_num<Real_wp>("2.38093616658814959133961201847753322162929660801133713928188656e-3");
                alpha(10) = to_num<Real_wp>("1.97623753300645954506780389071940568903586451129749850101236194e-3");
                alpha(11) = to_num<Real_wp>("1.666548315795813638218966297886409562176378624408447565788700473e-3");
                alpha(12) = to_num<Real_wp>("1.424235796847421825946714960742475038357436929275485236714163868e-3");
                alpha(13) = to_num<Real_wp>("1.230984401598711828723739182926703265059303993821003719775625728e-3");
                alpha(14) = to_num<Real_wp>("1.074245024880069666295620834588038120591426119776969686372793792e-3");
                alpha(15) = to_num<Real_wp>("9.451674085777489654424497109202524766644247271835706263680036943e-4");
                alpha(16) = to_num<Real_wp>("8.373516257592810466395147937737324331841592215477597342182636935e-4");
                alpha(17) = to_num<Real_wp>("7.460723091238380165918428334220494155733110818263907920763555571e-4");
                alpha(18) = to_num<Real_wp>("6.677853206462038703877471892327211092050567168206910676919567229e-4");
                alpha(19) = to_num<Real_wp>("5.998066812067465564575281130021124656198437995356769594497607506e-4");
                alpha(20) = to_num<Real_wp>("5.400968316329349446950984801638777958720486537831751231949534799e-4");
                alpha(21) = to_num<Real_wp>("4.871085099675614873384760454820859850376606834993260126611768313e-4");
                alpha(22) = to_num<Real_wp>("4.39672770280943583160771665855350031250199839844512959472047159e-4");
                alpha(23) = to_num<Real_wp>("3.969088356015046640237165992213333426502731085081532967689693683e-4");
                alpha(24) = to_num<Real_wp>("3.581508901738816851128600820378178330374852476416667094172186219e-4");
                alpha(25) = to_num<Real_wp>("3.228890768101051257480503132453050722626175081800086815406858146e-4");
                alpha(26) = to_num<Real_wp>("2.90723545945096506639559273825694227958642105152579008385626514e-4");
                alpha(27) = to_num<Real_wp>("2.613304112929951981492492785263146852794900690962982381900939481e-4");
                alpha(28) = to_num<Real_wp>("2.344379155811224246318145737021000779219094787003112218329073383e-4");
                alpha(29) = to_num<Real_wp>("2.098106762268162514038195823891830378635343466645787552546991411e-4");
                alpha(30) = to_num<Real_wp>("1.872398179421057487286301672139768559515912531956499431718091005e-4");
                alpha(31) = to_num<Real_wp>("1.665370621849734543006734290519952869938854990138811916921139984e-4");
                alpha(32) = to_num<Real_wp>("1.475312724262895561256637844698684863548137809608440979857206558e-4");
                alpha(33) = to_num<Real_wp>("1.300664018358833011510617578411153531202166122667547902320212417e-4");
                alpha(34) = to_num<Real_wp>("1.140001689868197754541755109437899537205944241855664273387920619e-4");
                alpha(35) = to_num<Real_wp>("9.920306738130189843003335752457292732029862574840455759745967598e-5");
                alpha(36) = to_num<Real_wp>("8.555750130964242304205134376279890359347392765547268339973767696e-5");
                alpha(37) = to_num<Real_wp>("7.295695407173868525730804102989470778348547268805773588860127223e-5");
                alpha(38) = to_num<Real_wp>("6.130515737509419054409141952563510559305554755879576499790091018e-5");
                alpha(39) = to_num<Real_wp>("5.051526150063329951609524864891641155210771624038408002703142183e-5");
                alpha(40) = to_num<Real_wp>("4.050901807665717464688940104515959298392323174503107717274634019e-5");
                alpha(41) = to_num<Real_wp>("3.121598980005820482690505145697752437234211255960956325583267116e-5");
                alpha(42) = to_num<Real_wp>("2.257279945960875068927502714093533389872253036741540221878390742e-5");
                alpha(43) = to_num<Real_wp>("1.452242706695565937923320957363410267132312936364930529168945449e-5");
                alpha(44) = to_num<Real_wp>("7.01356025540966938544468665725022968235683832751915839372189003e-6");
		break;
	      case 50 :
                alpha(5) = to_num<Real_wp>("6.41025640801022058760914934402816498805547983467145776703155092e-3");
                alpha(6) = to_num<Real_wp>("4.761904730557585830693696097183947211343402596394970101559518626e-3");
                alpha(7) = to_num<Real_wp>("3.676470375805363970306216415790678084541079598629495986719920417e-3");
                alpha(8) = to_num<Real_wp>("2.923975637935062016015813546167668000387980843373533673460615899e-3");
                alpha(9) = to_num<Real_wp>("2.380948957665775369763160050941903030843120290656680622784275507e-3");
                alpha(10) = to_num<Real_wp>("1.976274545552122783940888657263542817857152727250134097883678284e-3");
                alpha(11) = to_num<Real_wp>("1.666641078121245984290406240792210896314058510162303150974627796e-3");
                alpha(12) = to_num<Real_wp>("1.424443056273142531279417593152893367720934306236193445604668428e-3");
                alpha(13) = to_num<Real_wp>("1.231405499850504096099930512621915855259074017463501616460294959e-3");
                alpha(14) = to_num<Real_wp>("1.075034077289196199035448741388141148114348855438419986793799912e-3");
                alpha(15) = to_num<Real_wp>("9.465451765519377268341862537762702426347925552974276664718300557e-4");
                alpha(16) = to_num<Real_wp>("8.396110139696943161964630972587151644930131120072159669510049255e-4");
                alpha(17) = to_num<Real_wp>("7.495734481818099661028032383702453495340755837697501358904842974e-4");
                alpha(18) = to_num<Real_wp>("6.729377119863743723307230115722732869316896155152220809587551299e-4");
                alpha(19) = to_num<Real_wp>("6.070388728619203187105231257493214780706596237113186458689865725e-4");
                alpha(20) = to_num<Real_wp>("5.498177724912436223144779777153927681212257956648708664250089151e-4");
                alpha(21) = to_num<Real_wp>("4.996676873020633042464509096509096861061592712210139506494516559e-4");
                alpha(22) = to_num<Real_wp>("4.553271102115028841015976429384740429377927512796829061661027262e-4");
                alpha(23) = to_num<Real_wp>("4.158026365933482437939718683819530910239247100085106563177135663e-4");
                alpha(24) = to_num<Real_wp>("3.803114712319489979515280087767268635800115035699284325087494974e-4");
                alpha(25) = to_num<Real_wp>("3.482368572661431533445468263373915310114389963686099152000826072e-4");
                alpha(26) = to_num<Real_wp>("3.190923753129138091956685494319425678598514835829793835509130643e-4");
                alpha(27) = to_num<Real_wp>("2.924928406201776447054471088273516444636134859108293882229907983e-4");
                alpha(28) = to_num<Real_wp>("2.68130593567275325878894087667132279318028702683520889011175447e-4");
                alpha(29) = to_num<Real_wp>("2.457564953873926952053317401388982508139009498880169949289957997e-4");
                alpha(30) = to_num<Real_wp>("2.25165096722433017627647300608813636902161110372219359574680565e-4");
                alpha(31) = to_num<Real_wp>("2.061834264828913304787564182286261757082123348220713813401267276e-4");
                alpha(32) = to_num<Real_wp>("1.886627911499792688821737506978548577737580955665807394220422845e-4");
                alpha(33) = to_num<Real_wp>("1.724729549248528562431557205875981930722292098438132224449277718e-4");
                alpha(34) = to_num<Real_wp>("1.574981088943388786138311756633931524894207394256020744847519779e-4");
                alpha(35) = to_num<Real_wp>("1.436341199606207938364968817934127494854727724124964205165832333e-4");
                alpha(36) = to_num<Real_wp>("1.307866539509715124264735046808562931283828941508923642555724832e-4");
                alpha(37) = to_num<Real_wp>("1.188698709111378696189671744365012154019564029507497649699330932e-4");
                alpha(38) = to_num<Real_wp>("1.078054807712311291584231235826770163964544609722159917538555568e-4");
                alpha(39) = to_num<Real_wp>("9.752201880290448921569098626084143609585907485496517017221564527e-5");
                alpha(40) = to_num<Real_wp>("8.795425244068155495436656740242517987846371788274274097713949529e-5");
                alpha(41) = to_num<Real_wp>("7.904266687279295241856221404933341555016717918271146676475679107e-5");
                alpha(42) = to_num<Real_wp>("7.073300004226324154989338816597775859508008137072783840579544197e-5");
                alpha(43) = to_num<Real_wp>("6.297581192175339752479723937615784544182452392121131759699815706e-5");
                alpha(44) = to_num<Real_wp>("5.572608108813494198354074687433273490355051808250134878131614629e-5");
                alpha(45) = to_num<Real_wp>("4.894282593274479728580598551675000614574359574084033824859383024e-5");
                alpha(46) = to_num<Real_wp>("4.258874983994911580184915669397581073898918527000822141093991014e-5");
                alpha(47) = to_num<Real_wp>("3.662991036342811969370577705432257404163210642836031248939503689e-5");
                alpha(48) = to_num<Real_wp>("3.103541245973454971626627846524143976313152549826910916602230893e-5");
                alpha(49) = to_num<Real_wp>("2.577712557708198141319517148259594580161346610627861122758260885e-5");
                alpha(50) = to_num<Real_wp>("2.082942405207699308445871915675676103656788585604407137858731464e-5");
                alpha(51) = to_num<Real_wp>("1.616894994359160255106782012686734875899799683841668360006926357e-5");
                alpha(52) = to_num<Real_wp>("1.177439717638453040305130258771114708251825578786419971271145093e-5");
                alpha(53) = to_num<Real_wp>("7.626315689471183498290177395011958962758069189372861682420678352e-6");
                alpha(54) = to_num<Real_wp>("3.706934180934565169791260943108023924981783805033919713089060276e-6");
		break;
	      default :
		{
		  cout << "not implemented" << endl;
		  abort();
		}		
	      }
	  }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }
    
  }
  
  
  //! initialisation of the scheme with the initial condition U^0
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::
  SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
		      VirtualOdeSystem<T>& sys)
  {
    if (order_rk_first <= 4)
      {
        // using low-storage RK scheme as a first step
        LowStorageRK_Iterator<T> rk_scheme;
        rk_scheme.SetOrder(4);
        Unm1 = u0;
        rk_scheme.SetInitialCondition(t, deltat/nb_rk_step_first, u0, sys);
        for (int i = 0; i < nb_rk_step_first; i++)
          rk_scheme.Advance(t+Real_wp(i)*deltat/nb_rk_step_first, i, sys);
        
        Y = rk_scheme.GetIterate();
        rk_scheme.Clear();
      }
    else
      {
        // for higher orders, a Runge-Kutta scheme is used
        RungeKutta_Iterator<T> rk_scheme;
        rk_scheme.SetOrder(order_rk_first);
        Unm1 = u0;
        rk_scheme.SetInitialCondition(t, deltat/nb_rk_step_first, u0, sys);
        for (int i = 0; i < nb_rk_step_first; i++)
          rk_scheme.Advance(t+Real_wp(i)*deltat/nb_rk_step_first, i, sys);
        
        Y = rk_scheme.GetIterate();
        rk_scheme.Clear();
      }
    
    // allocating other vectors
    Q0 = Y;
    Q1 = Y;
    Qkp1 = Y;
    dt = deltat;
    sys.SetDirichletCondition(t, 0, Unm1);
    sys.SetDirichletCondition(t+dt, 0, Y);
  }
  

  //! memory is released
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory is released
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y.Clear();
    Q0.Clear(); Q1.Clear(); Qkp1.Clear();
    Yfinal = Unm1; Unm1.Clear();
  }
   

  //! returns stability function of the scheme
  template<class T>
  UnivariatePolynomial<Real_wp> OptimalModifiedEquation_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! U^n+1 is computed from U^n
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::
  Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {    
    T one; SetComplexOne(one);
    if (unstable_algo)
      {
        // Horner algorithm
        Copy(Y, Q0);
	
	Real_wp coef;
	for (int k = alpha.GetM()-1; k >= 0; k--)
	  {
	    sys.EvaluateFunction(t+dt, Q0, Q1);
	    sys.EvaluateFunction(t+dt, Q1, Q0);
	
	    coef = alpha(k)*dt*dt;
	    for (int i = 0; i < Q0.GetM(); i++)
	      Q0(i) = coef*Q0(i) + Y(i);
	  }
	
	sys.EvaluateFunction(t+dt, Q0, Q1);
	T un_next;
	coef = Real_wp(2)*dt;
	for (int i = 0; i < Q0.GetM(); i++)
	  {
	    un_next = Unm1(i) + coef*Q1(i);
	    Unm1(i) = Y(i);
	    Y(i) = un_next;
	  }
      }
    else if (fourth_order)
      {
        // Q0 = Un
        // Q1 = c dt Ah Unp
        Copy(Y, Q0);
        sys.EvaluateFunction(t+dt, Y, Q1);
        sys.SetDirichletCondition(t+dt, 0, Q1);
        Mlt(coef_xsi*dt, Q1);
        
        // Chebychev reccurence
	int r = (order/2-1) + nb_additional_terms;
        for (int k = 1; k <= 2*r; k++)
          {
            // Q_k+1 = 2 c dt Ah Q_k  + Q_k-1
            if (k%2 == 0)
              {
                if ((fourth_order) && (k == 2*r))
                  {
                    sys.EvaluateFunction(t+dt, Q1, Qkp1, false, true);
                    sys.AddPrimitiveTimeSource(alpha2, t+1.5*dt, 0, Qkp1);
                    sys.AddPrimitiveTimeSource(alpha2, t+0.5*dt, 0, Qkp1);
                    sys.AddPrimitiveTimeSource(-2.0*alpha2, t+dt, 0, Qkp1);
                    sys.SolveMassMatrix(Qkp1);
                  }
                else
                  sys.EvaluateFunction(t+dt, Q1, Qkp1, true, true);

                sys.SetDirichletCondition(t+dt, 0, Qkp1);
              }
            else
              {
                if ((fourth_order) && (k == 2*r-1))
                  {
                    sys.EvaluateFunction(t, Q1, Qkp1, false, false);
                    sys.AddPrimitiveTimeSource(alpha1, t+1.5*dt, 0, Qkp1);
                    sys.AddPrimitiveTimeSource(-alpha1, t+0.5*dt, 0, Qkp1);
                    sys.SolveMassMatrix(Qkp1);
                  }
                else
                  sys.EvaluateFunction(t+dt, Q1, Qkp1, true, false);

                sys.SetDirichletCondition(t+dt, 0, Qkp1);
              }
            
            Mlt(2.0*coef_xsi*dt, Qkp1);
            Add(one, Q0, Qkp1);
            Copy(Q1, Q0);
            Copy(Qkp1, Q1);
          }
        
        Mlt(2.0/xsi, Q1);
        Add(one, Unm1, Q1);
        Copy(Y, Unm1);
        Copy(Q1, Y);
        sys.SetDirichletCondition(t+2.0*dt, 0, Y);
      }
    else
      {
	// using three-term recurrence with coefficients alpha
	for (int i = 0; i < Y.GetM(); i++)
	  Q0(i) = -Y(i);
	
	sys.EvaluateFunction(t+dt, Y, Q1);
        sys.SetDirichletCondition(t+dt, 0, Q1);
	Mlt(alpha(0)*dt, Q1);
	
	//DISP(Q0); DISP(Q1);
	Real_wp coef;
	for (int k = 1; k < alpha.GetM(); k++)
	  {
	    sys.EvaluateFunction(t+dt, Q1, Qkp1);
	    sys.SetDirichletCondition(t+dt, 0, Qkp1);
	    Real_wp coef = -alpha(k)*dt;
	    for (int i = 0; i < Y.GetM(); i++)
	      Qkp1(i) = coef*Qkp1(i) + Q0(i);
	    
	    Copy(Q1, Q0); Copy(Qkp1, Q1);
	    
	    sys.EvaluateFunction(t+dt, Q1, Qkp1);
	    sys.SetDirichletCondition(t+dt, 0, Qkp1);
	    for (int i = 0; i < Y.GetM(); i++)
	      Qkp1(i) = coef*Qkp1(i) + Q0(i);
	    
	    Copy(Q1, Q0); Copy(Qkp1, Q1);	    
	  }
	
	T un_next; 
	coef = Real_wp(2);
	for (int i = 0; i < Q0.GetM(); i++)
	  {
	    un_next = Unm1(i) - coef*Q1(i);
	    Unm1(i) = Y(i);
	    Y(i) = un_next;
	  }
	
      }
  }


  //! changes time step during the simulation
  template<class T>
  void OptimalModifiedEquation_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented " << endl;
    abort();
  }

    
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& OptimalModifiedEquation_Iterator<T>::GetIterate()
  {
    return Unm1;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& OptimalModifiedEquation_Iterator<T>::GetIterate() const
  {
    return Unm1;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int OptimalModifiedEquation_Iterator<T>::GetNumberOfIterates() const
  {
    return 2;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& OptimalModifiedEquation_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else
      return Unm1;
  }

  
  //! returns the size used by the object in bytes
  template<class T>
  size_t OptimalModifiedEquation_Iterator<T>::GetMemorySize() const
  {
    return 0;
  }
  
}    

#define MONTJOIE_FILE_MODIFIED_EQUATION_CXX
#endif
