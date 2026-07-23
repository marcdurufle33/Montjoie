#ifndef MONTJOIE_FILE_TIME_SCHEMES_CXX

namespace Montjoie
{

  /********************
   * VirtualOdeSystem *
   ********************/
  
  
  //! provides the number of iterations that will be performed
  template<class T>
  void VirtualOdeSystem<T>
  ::GiveNumberIterations(const Real_wp& dt, int nb_max_iter)
  {
  }
    
  
  //! provides the iterate U^n
  template<class T>
  void VirtualOdeSystem<T>::GiveIterate(int n, const Real_wp& t, Vector<T>& Un)
  {
  }
  
  
  //! provides the final iterate U^N
  template<class T>
  void VirtualOdeSystem<T>::GiveFinalIterate(int n, const Real_wp&, Vector<T>&)
  {
  }

  
  //! provides the vectorial iterate V^N (for split formulation)
  template<class T>
  void VirtualOdeSystem<T>::GiveVectorialIterate(int n, const Real_wp&, Vector<T>&)
  {
  }

  
  //! evaluation of main function f in system y' = f(t, y)
  template<class T>
  void VirtualOdeSystem<T>::EvaluateFunction(const Real_wp& t, const Vector<T>& y, Vector<T>& f,
					     bool invert_mass, bool source)
  {
    cout << "EvaluateFunction not present in this ode system" << endl;
    abort();
  }
  

  //! evaluation of main function f in system y'' = f(t, y, y')
  template<class T>
  void VirtualOdeSystem<T>
  ::EvaluateFunctionS(const Real_wp& t, const Vector<T>& y, const Vector<T>& yp, Vector<T>& f,
		      bool invert_mass, bool source)
  {
    cout << "EvaluateFunctionS not present in this ode system" << endl;
    abort();
  }  
  
  
  //! Dirichlet condition
  template<class T>
  void VirtualOdeSystem<T>::SetDirichletCondition(const Real_wp& t, int n, Vector<T>&, Real_wp alpha)
  {
  }


  //! Dirichlet condition
  template<class T>
  void VirtualOdeSystem<T>::SetDirichletConditionSource(const Real_wp& t, int n, Vector<T>&, Real_wp alpha)
  {
  }


  //! Evaluates n-th derivative of f
  template<class T>
  void VirtualOdeSystem<T>::EvaluateDerivativeFunction(const Real_wp& t, int n,
						       const Vector<T>& y, Vector<T>& df,
						       bool invert_mass, bool source)
  {
    cout << "EvaluateDerivativeFunction not present in this ode system" << endl;
    abort();
  }

  
  //! overwrites x by the solution of M x = x where M is the mass matrix
  template<class T>
  void VirtualOdeSystem<T>::SolveMassMatrix(Vector<T>& X)
  {
    cout << "SolveMassMatrix not implemented for this system of ode" << endl;
    abort();
  }


  //! computes y = beta y + alpha M x where M is the mass matrix
  template<class T>
  void VirtualOdeSystem<T>::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t, 
					    const Vector<T>& x, const Real_wp& beta, Vector<T>& y)
  {
    cout << "ApplyMassMatrix not implemented for this system of ode" << endl;
    abort();
  }

  
  //! returns the final time of the source
  template<class T>
  Real_wp VirtualOdeSystem<T>::GetFinalTimeSource() const
  {
    cout << "GetFinalTimeSource not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! fills F with the space source only
  template<class T>
  void VirtualOdeSystem<T>::FillSource(Vector<T>& F) const
  {
    cout << "FillSource not implemented for this system of ode" << endl;
    abort();
  }
  

  //! sets pulse = g^(n)(t0) where g is the time source only
  /*!
    \param[in] t0 time where g^(n) is evaluated
    \param[in] n derivative number (if n=0, g is required)
    \param[out] pulse value g^(n)(t0)
   */
  template<class T>
  void VirtualOdeSystem<T>::SourceOnlyTime(const Real_wp& t0, int n, T& pulse)
  {
    cout << "SourceOnlyTime not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! Yh = beta Yh + alpha Kh Uh
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
		    const Real_wp& beta, Vector<T>& Yh)
  {
    cout << "ApplyOperatorKh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! computes Y = beta Y + alpha D_h U where D_h is the scalar mass matrix
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& U,
		    const Real_wp& beta, Vector<T>& Y)
  {
    cout << "ApplyOperatorDh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! overwrites Y with the solution X of D_h X = Y where D_h is the scalar mass matrix
  template<class T>
  void VirtualOdeSystem<T>::SolveOperatorDh(Vector<T>& Y)
  {
    cout << "SolveOperatorDh not implemented for this system of ode" << endl;
    abort();
  }

  
  //! Computes Prod_Uh = beta Prod_Uh + alpha (Dh - dt/2 Sh) Uh
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
			     const Real_wp& beta, Vector<T>& Prod_Uh)
  {
    cout << "ApplyOperatorDhMinusdtSh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! overwrites Y with the solution X of (Dh + dt/2 Sh) X = Y
  template<class T>
  void VirtualOdeSystem<T>::SolveOperatorDhPlusdtSh(Vector<T>& Y)
  {
    cout << "SolveOperatorDhPlusdtSh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! computes Y = Y + alpha (Dh + dt/2 Sh)^{-1} X
  template<class T>
  void VirtualOdeSystem<T>::SolveOperatorDhPlusdtSh(const Real_wp&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "SolveOperatorDhPlusdtSh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! Computes Prod_Uh = beta Prod_Uh + alpha Sh Uh
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
		    const Real_wp& beta, Vector<T>& Prod_Uh)
  {
    cout << "ApplyOperatorSh not implemented for this system of ode" << endl;
    abort();
  }
  
    
  //! Computes Vh = beta Vh + alpha RhV U
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
			     const Real_wp& beta, Vector<T>& Vh, bool extrapol)
  {
    cout << "ApplyOperatorRhVectorial not implemented for this system of ode" << endl;
    abort();
  }

  
  //! Computes Vh = beta Vh + alpha RhS U
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
			  const Real_wp& beta, Vector<T>& Vh, bool extrapol)
  {
    cout << "ApplyOperatorRhScalar not implemented for this system of ode" << endl;
    abort();
  }

  
  //! Computes Prod_Uh = beta Prod_Uh + alpha ShV Uh
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const Vector<T>& Uh,
			     const Real_wp& beta, Vector<T>& Prod_Uh)
  {
    cout << "ApplyOperatorShVectorial not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! overwrites Y with the solution X of Bh X = Y
  template<class T>
  void VirtualOdeSystem<T>::SolveOperatorBh(Vector<T>& Y)
  {
    cout << "SolveOperatorBh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! computes Vh = beta Vh + alpha (Bh - dt/2 ShV) Uh
  template<class T>
  void VirtualOdeSystem<T>
  ::ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const Vector<T>& Uh,
			     const Real_wp& beta, Vector<T>& Vh)
  {
    cout << "ApplyOperatorBhMinusdtSh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! overwrites Y with the solution X of (Bh + dt/2 ShV) X = Y
  template<class T>
  void VirtualOdeSystem<T>::SolveOperatorBhPlusdtSh(Vector<T>& Y)
  {
    cout << "SolveOperatorBhPlusdtSh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! computes b_src = b_src + alpha F^(n)(t)
  template<class T>
  void VirtualOdeSystem<T>::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
						   int nb_deriv, Vector<T>& b_src)
  {
    cout << "AddPrimitiveTimeSource not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! computes Y = Y + alpha F^(n)(t) where F is the scalar source
  template<class T>
  void VirtualOdeSystem<T>
  ::AddScalarTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, Vector<T>& Y)
  {
    cout << "AddScalarTimeSource not implemented for this system of ode" << endl;
    abort();
  }
  

  //! computes Y = Y + alpha F^(n)(t) where F is the vectorial source  
  template<class T>
  void VirtualOdeSystem<T>
  ::AddVectorialTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, Vector<T>& Y)
  {
    cout << "AddVectorialTimeSource not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! solves the system Op(num) Un_next = Un
  template<class T>
  void VirtualOdeSystem<T>::SolveSplitOperator(const Real_wp& t, const Real_wp& dt,
					       const Vector<T>& Un, Vector<T>& Un_next, int num)
  {
    cout << "SolveSplitOperator not implemented for this system of ode" << endl;
    abort();
  }

  
  //! Factorizes alpha M_h + beta S_h + gamma K_h for a single set of coefficients
  template<class T>
  void VirtualOdeSystem<T>::FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta,
							   const Real_wp& gamma)
  {
    cout << "FactorizeOperatorDhPlusGammaKh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! Solves (alpha M_h + beta S_h + gamma K_h) Y = b for a single set of coefficients
  template<class T>
  void VirtualOdeSystem<T>
  ::SolveOperatorDhPlusGammaKh(const Real_wp& t, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "SolveOperatorDhPlusGammaKh not implemented for this system of ode" << endl;
    abort();
  }


  //! Factorizes alpha M_h + beta S_h + gamma K_h for a set of real coefficients
  template<class T>
  void VirtualOdeSystem<T>
  ::FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
			  const VectReal_wp& gamma)    
  {
    cout << "FactorizeOperatorReal not implemented for this system of ode" << endl;
    abort();
  }


  //! Solves (alpha M_h + beta S_h + gamma K_h) Y = b for a set of real coefficients
  template<class T>
  void VirtualOdeSystem<T>
  ::SolveOperatorReal(const Real_wp& t, const Vector<T>& X, Vector<T>& Y,
		      int num_system)
  {
    cout << "SolveOperatorReal not implemented for this system of ode" << endl;
    abort();
  }
  

  //! Factorizes the system k_i = f(y^n + dt \sum a_ij k_j)
  template<class T>
  void VirtualOdeSystem<T>
  ::FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& AgaussRK, const Real_wp& dt,
				 const Vector<T>& Y, const Vector<T>& ki)
  {
    cout << "FactorizeOperatorDhPlusAijKh not implemented for this system of ode" << endl;
    abort();
  }
  
  
  //! Solves the system k_i = f(y^n + dt \sum a_ij k_j)
  template<class T>
  void VirtualOdeSystem<T>
  ::SolveOperatorDhPlusAijKh(const VectReal_wp& t, const Matrix<Real_wp>& A,
			     const Real_wp& dt, const Vector<T>& Y, Vector<T>& ki)
  {
    cout << "SolveOperatorDhPlusAijKh not implemented for this system of ode" << endl;
    abort();
  }


  //! Factorizes alpha M_h + beta S_h + gamma K_h for a set of complex coefficients
  template<class T>
  void VirtualOdeSystem<T>
  ::FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
			     const VectComplex_wp& gamma)    
  {
    cout << "FactorizeOperatorComplex not implemented for this system of ode" << endl;
    abort();
  }


  //! Solves (alpha M_h + beta S_h + gamma K_h) Y = b for a set of complex coefficients
  template<class T>
  void VirtualOdeSystem<T>
  ::SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			 int num_system)
  {
    cout << "SolveOperatorComplex not implemented for this system of ode" << endl;
    abort();
  }
    
  
  /*********************
   * VirtualTimeScheme *
   *********************/
  

  //! destructor
  template<class T>
  VirtualTimeScheme<T>::~VirtualTimeScheme()
  {
  }
  

  //! initializes the time scheme with an initial condition U0 and time step dt
  template<class T>
  void VirtualTimeScheme<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& dt,
			Vector<T>& U0, VirtualOdeSystem<T>& sys)
  {
    cout << "Function SetInitialCondition not implemented for this time scheme" << endl;
    abort();
  }
  

  //! initializes the time scheme for the ode y'' = f(y, y') with two initial conditions y(0) and y'(0)
  template<class T>
  void VirtualTimeScheme<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& dt,
                         Vector<T>& U0, Vector<T>& dU0, VirtualOdeSystem<T>& sys)
  {
    cout << "Function SetInitialCondition not implemented for this time scheme" << endl;
    abort();
  }


  //! not used
  template<class T>
  void VirtualTimeScheme<T>::SetParametersOutputEnergy(const Vector<string>& parameters)
  {
    cout << "No energy implemented for this scheme" << endl;
    abort();
  }


  //! not used 
  template<class T>
  void VirtualTimeScheme<T>::SetEnergyFile(const Vector<string>& parameters) 
  {
    cout << "No energy implemented for this scheme" << endl;
    abort();
  }


  //! returns the vectorial unknown (for staggered schemes)
  template<class T>
  Vector<T>& VirtualTimeScheme<T>::GetVectorialIterate()
  {
    cout << "Function GetVectorialIterate not present for this time scheme" << endl;
    abort();
  }


  //! returns the number of vectorial unknowns needed to advance the scheme
  template<class T>
  int VirtualTimeScheme<T>::GetNumberOfVectorialIterates() const
  {
    cout << "Function GetNumberOfVectorialIterates not present for this time scheme" << endl;
    abort();
  }


  //! returns the vectorial unknown k
  template<class T>
  Vector<T>& VirtualTimeScheme<T>::GetVectorialIterate(int k)
  {
    cout << "Function GetVectorialIterate not present for this time scheme" << endl;
    abort();
  }


  //! clears the scheme and returns the final iterate in Un
  template<class T>
  void VirtualTimeScheme<T>::ClearFirst(Vector<T>& Un)
  {
    cout << "Function Clear not implemented in this time scheme" << endl;
    abort();
  }
  

  //! clears the scheme (for ode y'' = f(y, y)) and returns the final iterates y and y'
  template<class T>
  void VirtualTimeScheme<T>::ClearSecond(Vector<T>& Un, Vector<T>& dUn)
  {
    cout << "Function Clear not implemented in this time scheme" << endl;
    abort();
  }


  //! returns the CFL number associated with the scheme
  template<class T>
  Real_wp VirtualTimeScheme<T>::GetCFL(bool real_p) const
  {
    cout << "We don't know the cfl of this scheme" << endl;
    abort();
    
    return Real_wp(0);
  }

  template<class T>
  UnivariatePolynomial<Real_wp> VirtualTimeScheme<T>::GetStabilityFunction() const
  {
    UnivariatePolynomial<Real_wp> R;
    R.SetOrder(0); R(0) = 1.0;
    return R;
  }


  //! returns the memory used by the object
  template<class T>
  size_t VirtualTimeScheme<T>::GetMemorySize() const
  {
    return 0;
  }

  
  //! fills memory usage
  template<class T>
  void VirtualTimeScheme<T>::GetMemoryUsed(map<string, size_t>& var) const
  {
  }

  
  /******************
   * All_TimeScheme *
   ******************/
  

  //! default constructor
  template<class T>
  All_TimeScheme<T>::All_TimeScheme()
  {
    // default scheme
    type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE;
    order_time_scheme = 4;
    nb_iterates = 1;
    
    LowStorageRK_Iterator<T>* rk_scheme;
    rk_scheme = new LowStorageRK_Iterator<T>();
    rk_scheme->SetOrder(order_time_scheme);
    time_scheme = rk_scheme;
  }  
  
  
  //! releases memory used by the object
  template<class T>
  All_TimeScheme<T>::~All_TimeScheme()
  {
    Clear();
  }
  

  //! returns the size used by the object in bytes
  template<class T>
  size_t All_TimeScheme<T>::GetMemorySize() const
  {
    return time_scheme->GetMemorySize();
  }
  

  template<class T>
  void All_TimeScheme<T>::GetMemoryUsed(map<string, size_t>& var) const
  {
    time_scheme->GetMemoryUsed(var);
  }
  
  
  //! returns the actual time scheme contained in the class
  template<class T>
  VirtualTimeScheme<T>& All_TimeScheme<T>::GetTimeScheme() const
  {
    CheckPresenceScheme();
    
    return *time_scheme;
  }


  //! returns the type of time scheme selected
  template<class T>
  int All_TimeScheme<T>::GetTimeSchemeType() const
  {
    return type_time_scheme;
  }
  

  //! returns the order of the selected time scheme
  template<class T>
  int All_TimeScheme<T>::GetOrder() const
  {
    return order_time_scheme;
  }
  
  
  //! sets time scheme to use with a line of the datafile
  template<class T>
  void All_TimeScheme<T>::SetInputData(const string& keyword, const Vector<string>& parameters)
  {
    if (keyword == "OrderTimeScheme")
      {

        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of All_TimeScheme" << endl;
            cout << "OrderTimeScheme needs more parameters, for instance :" << endl;
            cout << "OrderTimeScheme = r type_scheme" << endl;
            cout << "Current parameters are " << endl << parameters << endl;
            abort();
          }
        
        order_time_scheme = to_num<int>(parameters(0));
        Vector<string> param(parameters.GetM() - 1);
        for (int i = 1; i < parameters.GetM(); i++)
          param(i-1) = parameters(i);
        
        SetTimeScheme(order_time_scheme, param);
      }

    else if (keyword == "ParametersOutputEnergy")
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of All_TimeScheme" << endl;
            cout << "ParametersOuputEnergy needs more parameters, for instance :" << endl;
            cout << "ParametersOuputEnergy = t_begin t_end dt buffer_size" << endl;
            cout << "Current parameters are " << endl << parameters << endl;
            abort();
          }
        
	time_scheme->SetParametersOutputEnergy(parameters);
      }
    else if (keyword == "FileOutputEnergy")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of All_TimeScheme" << endl;
            cout << "FileOuputEnergy needs more parameters, for instance :" << endl;
            cout << "FileOuputEnergy = file_name" << endl;
            cout << "Current parameters are " << endl << parameters << endl;
            abort();
          }
       
	time_scheme->SetEnergyFile(parameters);
      }
  }

  
  //! is the time scheme implicit ?
  template<class T>
  bool All_TimeScheme<T>::IsImplicitScheme() const
  {
    switch (type_time_scheme)
      {
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
      case TimeSchemeEnum::SDIRK :
      case TimeSchemeEnum::LINEAR_SDIRK :
      case TimeSchemeEnum::DIRK :
      case TimeSchemeEnum::CRANK_NICOLSON :
      case TimeSchemeEnum::ESDIRK :
      case TimeSchemeEnum::THETA_SCHEME :
      case TimeSchemeEnum::PREDICTOR :
      case TimeSchemeEnum::ADAMS_IMPLICIT :
      case TimeSchemeEnum::MILNE_SIMPSON :
      case TimeSchemeEnum::BDF :
      case TimeSchemeEnum::PADE_SCHEME :
      case TimeSchemeEnum::IMPLICIT_SYMMETRIC_MULTISTEP :
        return true;
      }
    
    return false;
  }
  

  //! is the time scheme adapted to a first-order formulation ?
  template<class T>
  bool All_TimeScheme<T>::FirstOrderScheme() const
  {
    switch (type_time_scheme)
      {
      case TimeSchemeEnum::THETA_SCHEME :
      case TimeSchemeEnum::LEAP_FROG_ORDER2 :
      case TimeSchemeEnum::RUNGE_KUTTA_NYSTROM :
      case TimeSchemeEnum::OPTIMAL_LEAP_FROG :
      case TimeSchemeEnum::SYMMETRIC_MULTISTEP :
      case TimeSchemeEnum::IMPLICIT_SYMMETRIC_MULTISTEP :
        return false;
      }
    
    return true;
  }


  //! is the time scheme with staggered times (t^n and t^{n+1/2}) ?
  template<class T>
  bool All_TimeScheme<T>::StaggeredScheme() const
  {
    switch (type_time_scheme)
      {
      case TimeSchemeEnum::LEAP_FROG_SYSTEM :
        return true;
      }
    
    return false;
  }


  template<class T>
  bool All_TimeScheme<T>::LocalTimeStepping() const
  {
    switch (type_time_scheme)
      {
      case TimeSchemeEnum::LOCAL_IMPLICIT_RK :
      case TimeSchemeEnum::LOCAL_IMEX_RK :
      case TimeSchemeEnum::LOCAL_IMPERIALE_SCHEME :
      case TimeSchemeEnum::LOCAL_PIPERNO_SCHEME :
        return true;
      }
    
    return false;
  }
  
  
  //! setting the time scheme with parameters of the data file
  /*!
    \param[in] r order of the scheme
    \param[in] param description of the scheme
    param is a list of strings, the first string names the time scheme
    the other strings are additional parameters (for example theta).    
   */
  template<class T>
  void All_TimeScheme<T>::SetTimeScheme(int r, const Vector<string>& param)
  {
    // previous scheme is removed if present
    this->Clear();
    
    if (param.GetM() <= 0)
      {
        cout << "Provide a non-empty vector" << endl;
        abort();
      }

    if (param(0) == "LEAP_FROG_TRUE")
      {
	Nystrom_Iterator<T>* nystrom = new Nystrom_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::LEAP_FROG;
        // leap frog on first order formulation
        // without split unknowns (u, v)
        // it is equivalent to a Nystrom scheme with order 1
        nystrom->SetOrder(1);
        nb_iterates = nystrom->GetNumberOfIterates();
	time_scheme = nystrom;
      }
    else if (param(0) == "LEAP_FROG_EXACT")
      {
        type_time_scheme = TimeSchemeEnum::LEAP_FROG_EXACT;
        // not implemented in this general interface
      }
    else if (param(0) == "LEAP_FROG")
      {
	ModifiedEquationIterator<T>* second_order_leap_frog
	  = new ModifiedEquationIterator<T>();
	
        type_time_scheme = TimeSchemeEnum::LEAP_FROG_ORDER2;
        // centered scheme for second-order formulation
        second_order_leap_frog->SetOrder(r);
        nb_iterates = second_order_leap_frog->GetNumberOfIterates();
	time_scheme = second_order_leap_frog;
      }
    else if (param(0) == "LEAP_FROG_PML")
      {
	ModifiedEquationSystemIterator<T>* first_order_leap_frog
	  = new ModifiedEquationSystemIterator<T>();
	
        type_time_scheme = TimeSchemeEnum::LEAP_FROG_SYSTEM;
        // leap frog on first order split formulation (with two unknowns u, v)
        // and a staggered discretization (u evaluate at t^n, v at t^{n+1/2})
        first_order_leap_frog->SetOrder(r);
        nb_iterates = first_order_leap_frog->GetNumberOfIterates();
	time_scheme = first_order_leap_frog;
      }
    else if (param(0) == "OPTIMAL_LEAP_FROG")
      {
	OptimalModifiedEquationSecond_Iterator<T>* opt_leap_frog_centered
	  = new OptimalModifiedEquationSecond_Iterator<T>();
	
        type_time_scheme = TimeSchemeEnum::OPTIMAL_LEAP_FROG;
        // modified equation on centered scheme for second-order formulation
        // with additional terms in order to obtain an optimal cfl
        opt_leap_frog_centered->SetOrder(r);
        nb_iterates = opt_leap_frog_centered->GetNumberOfIterates();
	time_scheme = opt_leap_frog_centered;
      }
    else if (param(0) == "OPTIMAL_LEAP_FROG_PML")
      {
        if (param.GetM() < 4)
          {
            cout << "Not enough parameters in SetTimeScheme" << endl;
            cout << "For example OrderTimeScheme = 2 OPTIMAL_LEAP_FROG_PML 10 STABLE 10" << endl;
            cout << "Current parameters are " << endl << param << endl;
            abort();
          }
        
	OptimalModifiedEquation_Iterator<T>* opt_leap_frog_pml
	  = new OptimalModifiedEquation_Iterator<T>();
	
        type_time_scheme = TimeSchemeEnum::OPTIMAL_LEAP_FROG_PML;
        // modified equation on leap frog scheme (as in LEAP_FROG_TRUE)
        // with additional terms in order to obtain an optimal cfl
        int m = to_num<int>(param(1));
        bool unstable_algo = true;
        if (param(2) == "STABLE")
          unstable_algo = false;
        
	int nb_rk_steps = to_num<int>(param(3));
        int order_rk = 4;
        if (param.GetM() > 4)
          order_rk = to_num<int>(param(4));
        
        opt_leap_frog_pml->SetOrder(r, m, unstable_algo, nb_rk_steps, order_rk);
        nb_iterates = opt_leap_frog_pml->GetNumberOfIterates();
	time_scheme = opt_leap_frog_pml;
      }
    else if (param(0) == "SPLIT_EXACT")
      {
	SplitScheme_Iterator<T>* split_scheme = new SplitScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::SPLIT_EXACT;
        // schemes based on a splitting, each step is assumed to be solved exactly
        split_scheme->SetOrder(r);
        nb_iterates = split_scheme->GetNumberOfIterates();
	time_scheme = split_scheme;
      }
    else if (param(0) == "SYMPLECTIC")
      {
        type_time_scheme = TimeSchemeEnum::SYMPLECTIC;
        // not implemented at all
	abort();
      }
    else if (param(0) == "NYSTROM")
      {
	Nystrom_Iterator<T>* nystrom = new Nystrom_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::NYSTROM;
        // Nystrom scheme (see Hairer's book, vol I)
        if (param.GetM() >= 2)
          nystrom->SetOrder(r, nystrom->MILNE_PREDICTOR);
        else
          nystrom->SetOrder(r);

        nb_iterates = nystrom->GetNumberOfIterates();
	time_scheme = nystrom;
      }
    else if (param(0) == "RUNGE_KUTTA_NYSTROM")
      {
	RungeKuttaNystrom_Iterator<T>* rk_nystrom = new RungeKuttaNystrom_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_NYSTROM;
        // explicit Runge-Kutta Nystrom
        if (param.GetM() >= 2)
          {
            // case y'' = f(t, y, y') => f depends on y'
            rk_nystrom->SetOrder(r, true);
          }
        else
          {
            // case y'' = f(t, y)
            rk_nystrom->SetOrder(r);
          }
        
        nb_iterates = rk_nystrom->GetNumberOfIterates();
	time_scheme = rk_nystrom;
      }
    else if (param(0) == "SYMMETRIC_MULTISTEP")
      {
	SymmetricMultistep_Iterator<T>* sym_multistep_scheme
	  = new SymmetricMultistep_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::SYMMETRIC_MULTISTEP;
        if (r >= 2)
          {
            VectReal_wp theta((r-2)/2);
            
            if (param.GetM() == 1)
              {
                // using acceptable choice (cfl is correct)
                if (r == 4)
                  {
                    theta(0) = pi_wp/2;
                  }
                else if (r == 6)
                  {
                    theta(0) = pi_wp/2;
                    theta(1) = 2*pi_wp/3;
                  }
                else if (r == 8)
                  {
                    theta(0) = pi_wp/6;
                    theta(1) = pi_wp/2;
                    theta(2) = 2*pi_wp/3;
                  }
                else
                  {
                    cout << "Default choice not available" << endl;
                    abort();
                  }
              }
            else
              {
                if (param.GetM() <= theta.GetM())
                  {
                    cout << "Not enough parameters for SYMMETRIC_MULTISTEP"<< endl;
                    cout << "You need to specify theta parameters" << endl;
                    cout << "For example OrderTimeScheme = 4 SYMMETRIC_MULTISTEP theta_0" << endl;
                    cout << "Current parameters are " << param << endl;
                    abort();
                  }
                
                for (int k = 0; k < theta.GetM(); k++)
                  theta(k) = to_num<Real_wp>(param(k+1));
              }

            // Symmetric Multistep
            sym_multistep_scheme->SetOrder(r, theta);
            nb_iterates = sym_multistep_scheme->GetNumberOfIterates();
          }
	
	time_scheme = sym_multistep_scheme;
      }
    else if (param(0) == "IMPLICIT_SYMMETRIC_MULTISTEP")
      {
	ImplicitSymmetricMultistep_Iterator<T>* sym_multistep_scheme
	  = new ImplicitSymmetricMultistep_Iterator<T>();

        type_time_scheme = TimeSchemeEnum::IMPLICIT_SYMMETRIC_MULTISTEP;
        if (r >= 2)
          {
            VectReal_wp theta((r-2)/2);
            Real_wp phi(to_num<Real_wp>(param(theta.GetM()+1)));
            for (int k = 0; k < theta.GetM(); k++)
              theta(k) = to_num<Real_wp>(param(k+1));
            sym_multistep_scheme->SetOrder(r, theta, phi);
            nb_iterates = sym_multistep_scheme->GetNumberOfIterates();
          }

	time_scheme = sym_multistep_scheme;
      }
    else if (param(0) == "ADAMS_BASHFORTH")
      {
	AdamsBashforth_Moulton_Iterator<T>* AB_scheme
	  = new AdamsBashforth_Moulton_Iterator<T>();
	
        type_time_scheme = TimeSchemeEnum::ADAMS_BASHFORTH;
        // Adams-Bashforth schemes
        AB_scheme->SetOrder(r);
        nb_iterates = AB_scheme->GetNumberOfIterates();
	time_scheme = AB_scheme;
      }
    else if (param(0) == "ADAMS_BASHFORTH_MOULTON")
      {
	AdamsBashforth_Moulton_Iterator<T>* AB_scheme
	  = new AdamsBashforth_Moulton_Iterator<T>();

        type_time_scheme = TimeSchemeEnum::ADAMS_BASHFORTH;
        // Adams-Bashforth schemes with Moulton correction
        AB_scheme->SetOrder(r, AB_scheme->BASHFORTH_MOULTON);
        nb_iterates = AB_scheme->GetNumberOfIterates();
	time_scheme = AB_scheme;
      }
    else if (param(0) == "MULTISTEP_BUTCHER")
      {
	MultiStepButcher_Iterator<T>* butcher_scheme
	  = new MultiStepButcher_Iterator<T>();
	
        type_time_scheme = TimeSchemeEnum::MULTISTEP_BUTCHER;
        // multistep schemes as described in Butcher's book
        // schemes that involve high order derivatives of f
        if (param.GetM() >= 2)
          butcher_scheme->SetOrder(r, butcher_scheme->PECE);
        else
          butcher_scheme->SetOrder(r);

        nb_iterates = butcher_scheme->GetNumberOfIterates();
	time_scheme = butcher_scheme;
      }
    else if (param(0) == "TALEZER")
      {
        if (param.GetM() <= 2)
          {
            cout 
              << "Talezer schemes needs two parameters : R and the order of integration" << endl;
            cout << "Current provided parameters are : " << endl << param << endl;
            abort();
          }
	
	Talezer_Iterator<T>* talezer_scheme = new Talezer_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::TALEZER;
        talezer_scheme->SetSpectralRadius(to_num<Real_wp>(param(1)));
        talezer_scheme->SetOrderIntegration(to_num<int>(param(2)));        
        // SetOrder of this scheme is called in SetInitialCondition
        // since it needs the time step 
        nb_iterates = talezer_scheme->GetNumberOfIterates();
	time_scheme = talezer_scheme;
      }
    else if (param(0) == "TAYLOR_SERIES")
      {
	TaylorSeries_Iterator<T>* taylor_scheme = new TaylorSeries_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::TAYLOR_SERIES;
        // scheme based on Taylor expansion (using high order derivatives of f)
        taylor_scheme->SetOrder(r);
        nb_iterates = taylor_scheme->GetNumberOfIterates();
	time_scheme = taylor_scheme;
      }
    else if (param(0) == "RUNGE_KUTTA_LOW")
      {
	LowStorageRK_Iterator<T>* LowRK_scheme = new LowStorageRK_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE;
        // low-storage Runge-Kutta schemes
        if (param.GetM() >= 2)
          {
            // alternative sets of coefficients
            int nb_steps = to_num<int>(param(1));
            LowRK_scheme->SetOrder(r, nb_steps);
          }
        else
          LowRK_scheme->SetOrder(r);

        nb_iterates = LowRK_scheme->GetNumberOfIterates();
	time_scheme = LowRK_scheme;
      }
    else if (param(0) == "RUNGE_KUTTA")
      {
	RungeKutta_Iterator<T>* RK_scheme = new RungeKutta_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA;
        // explicit Runge-Kutta schemes
        // default choice : coefficients found in Butcher's book
        int type = RK_scheme->BUTCHER;
        if (param.GetM() >= 2)
          {
            // other sets of coefficients
            if (param(1) == "FEHLBERG")
              type = RK_scheme->FEHLBERG;
            else if (param(1) == "DOPRI")
              type = RK_scheme->DOPRI;
          }
        
        RK_scheme->SetOrder(r, type);
        nb_iterates = RK_scheme->GetNumberOfIterates();
	time_scheme = RK_scheme;
      }
    else if (param(0) == "LINEAR_RUNGE_KUTTA")
      {
	LinearRungeKutta_Iterator<T>* RK_scheme = new LinearRungeKutta_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::LINEAR_RUNGE_KUTTA;
	// Number of additional stages
	int extraS_l = 0;
        if (param.GetM() >= 2)
          extraS_l = to_num<int>(param(1));

        bool stable_source = false;
        if ((param.GetM() >= 3) && (param(2) == "SourceStable"))
          stable_source = true;
        
        RK_scheme->SetOrder(r, extraS_l, true, stable_source);
        nb_iterates = RK_scheme->GetNumberOfIterates();
	time_scheme = RK_scheme;
      }
    else if (param(0) == "RUNGE_KUTTA_EXACT")
      {
        type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_EXACT;
        // not implemented in this general interface
      }
    else if (param(0) == "GAUSS_RUNGE_KUTTA")
      {
	GaussRungeKutta_Iterator<T>* gauss_rk_scheme = new GaussRungeKutta_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::GAUSS_RUNGE_KUTTA;
        // implicit Runge-Kutta schemes based on Gauss-Legendre formulae
        gauss_rk_scheme->SetOrder(r);
        nb_iterates = gauss_rk_scheme->GetNumberOfIterates();
	time_scheme = gauss_rk_scheme;
      }
    else if (param(0) == "SDIRK")
      {
	SdirkScheme_Iterator<T>* sdirk_scheme = new SdirkScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::SDIRK;
        // Singly diagonally implicit Runge-Kutta schemes
        sdirk_scheme->SetOrder(r);
        nb_iterates = sdirk_scheme->GetNumberOfIterates();
	time_scheme = sdirk_scheme;
      }
    else if (param(0) == "LINEAR_SDIRK")
      {
	LinearSdirkScheme_Iterator<T>* linear_sdirk_scheme = new LinearSdirkScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::LINEAR_SDIRK;
        // Number of additional stages
	int extraS_l;
	if(param.GetM() >= 2)
	  extraS_l = to_num<int>(param(1));
	else 
	  extraS_l = 0;
        
        int type_algo = linear_sdirk_scheme->POLYNOMIAL;
        if (param.GetM() >= 3)
          {
            if (param(2) == "Stable")
              type_algo = linear_sdirk_scheme->STABLE_WEIGHTS;
            else if (param(2) == "SourceStablePoly")
              type_algo = linear_sdirk_scheme->POLY_ALTERNATIVE;
            else if (param(2) == "SourceStable")
              type_algo = linear_sdirk_scheme->STABLE_ALTERNATIVE;
          }

	// Linear Singly diagonally implicit Runge-Kutta schemes
        linear_sdirk_scheme->SetOrder(r, extraS_l, type_algo);
        nb_iterates = linear_sdirk_scheme->GetNumberOfIterates();
	time_scheme = linear_sdirk_scheme;
      }
    else if (param(0) == "DIRK")
      {
	DirkScheme_Iterator<T>* dirk_scheme = new DirkScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::DIRK;
        //Diagonally implicit Runge-Kutta schemes
        dirk_scheme->SetOrder(r);
        nb_iterates = dirk_scheme->GetNumberOfIterates();
	time_scheme = dirk_scheme;
      }
    else if (param(0) == "CRANK_NICOLSON")
      {
	GaussRungeKutta_Iterator<T>* gauss_rk_scheme = new GaussRungeKutta_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::CRANK_NICOLSON;
        // Crank-Nicolson scheme is a particular case of Gauss-Runge-Kutta schemes
        // for order = 2
        gauss_rk_scheme->SetOrder(2);
        nb_iterates = gauss_rk_scheme->GetNumberOfIterates();
	time_scheme = gauss_rk_scheme;
      }
    else if (param(0) == "ESDIRK")
      {
	SdirkScheme_Iterator<T>* sdirk_scheme = new SdirkScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::ESDIRK;
        // same as sdirk scheme
        sdirk_scheme->SetOrder(r);
        nb_iterates = sdirk_scheme->GetNumberOfIterates();
	time_scheme = sdirk_scheme;
      }
    else if (param(0) == "THETA_SCHEME")
      {
	TetaScheme_Iterator<T>* teta_scheme = new TetaScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::THETA_SCHEME;
        if (param.GetM() >= 2)
          teta_scheme->SetTheta(to_num<Real_wp>(param(1)));
        else
          teta_scheme->SetTheta(0.25);        
        
        nb_iterates = teta_scheme->GetNumberOfIterates();
	time_scheme = teta_scheme;
      }
    else if (param(0) == "PREDICTOR")
      {
        type_time_scheme = TimeSchemeEnum::PREDICTOR;
        // not implemented in this general interface
      }
    else if (param(0) == "ADAMS_IMPLICIT")
      {
	AdamsImplicit_Iterator<T>* implicit_adams = new AdamsImplicit_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::ADAMS_IMPLICIT;
        // implicit adams method
        implicit_adams->SetOrder(r);
        nb_iterates = implicit_adams->GetNumberOfIterates();
	time_scheme = implicit_adams;
      }
    else if (param(0) == "MILNE_SIMPSON")
      {
	MilneSimpson_Iterator<T>* milne_simpson = new MilneSimpson_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::MILNE_SIMPSON;
        // implicit Milne-Simpson method
        milne_simpson->SetOrder(r);
        nb_iterates = milne_simpson->GetNumberOfIterates();
	time_scheme = milne_simpson;
      }
    else if (param(0) == "BDF")
      {
	BackwardDifferentiation_Iterator<T>* bdf_scheme = new BackwardDifferentiation_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::BDF;
        // BDF scheme
        bdf_scheme->SetOrder(r);
        nb_iterates = bdf_scheme->GetNumberOfIterates();
	time_scheme = bdf_scheme;
      }
    else if (param(0) == "PADE_SCHEME")
      {
	PadeScheme_Iterator<T>* pade_scheme = new PadeScheme_Iterator<T>();
        type_time_scheme = TimeSchemeEnum::PADE_SCHEME;
        bool stable = false;
        if (param.GetM() >= 2)
          if (param(1) == "Stable")
            stable = true;
        
        // Pade schemes
        pade_scheme->SetOrder(r, stable);
        nb_iterates = pade_scheme->GetNumberOfIterates();
	time_scheme = pade_scheme;
      }
    else if (param(0) == "LOCAL_IMPLICIT_RK")
      {
	LocalTimeSteppingRK<T>* scheme = new LocalTimeSteppingRK<T>();
	type_time_scheme = TimeSchemeEnum::LOCAL_IMPLICIT_RK;

	scheme->SetInputData(r, param);
	nb_iterates = scheme->GetNumberOfIterates();
	time_scheme = scheme;
      }
    else if (param(0) == "LOCAL_IMEX_RK")
      {
	LocalImexRK<T>* scheme = new LocalImexRK<T>();
	type_time_scheme = TimeSchemeEnum::LOCAL_IMEX_RK;

	scheme->SetInputData(r, param);
	nb_iterates = scheme->GetNumberOfIterates();
	time_scheme = scheme;
      }
    else if (param(0) == "LOCAL_IMPERIALE")
      {
	// Local time stepping scheme based on Imperiale's idea
	//LocalImperialeScheme_Iterator<Real_wp> local_time_scheme; 

	abort();
        //time_scheme = TimeSchemeEnum::LOCAL_IMPERIALE_SCHEME;
        //if (param.GetM() > 1)
	//local_time_scheme.SetInputData(param);
      }
    else if (param(0) == "LOCAL_PIPERNO")
      {
	// Local time stepping (symplectic approach of Piperno)
	// LocalPipernoScheme_Iterator<Real_wp> local_piperno_scheme;
	
	abort();
        //time_scheme = TimeSchemeEnum::LOCAL_PIPERNO_SCHEME;
        //if (param.GetM() > 1)
	//local_piperno_scheme.SetInputData(param);
      }
    else
      {
        cout << "Unknown time scheme : " << param(0) << endl;
        abort();
      }
  }

  
  //! Checks if a time scheme is present
  template<class T>
  void All_TimeScheme<T>::CheckPresenceScheme() const
  {
    if (time_scheme == NULL)
      {
	cout << "Time scheme not present" << endl;
	cout << "Did you write OrderTimeScheme in your data file?" << endl;
	abort();
      }
  }


  //! changes time step during the simulation
  template<class T>
  void All_TimeScheme<T>::ChangeTimeStep(const Real_wp& dt)
  {
    CheckPresenceScheme();    
    time_scheme->ChangeTimeStep(dt);
  }
  
  
  //! returns CFL of the time scheme used
  template<class T>
  Real_wp All_TimeScheme<T>::GetCflScheme(bool real_p) const
  {
    // the CFL is computed on the imaginary axis for first-order formulation
    // and real axis on second-order formulation

    // First-order formulation : dU/dt = A U
    // if A is assumed to have only purely imaginary eigenvalues
    // explicit time schemes have a stability condition (also known as CFL) :
    //  dt <= alpha / |Lambda_{max}(A)| 
    // where Lambda_{max}(A) is the maximum eigenvalue (in modulus) of A
    // This function returns the coefficient alpha
    
    // Second-order formulation : d^2 U/dt^2 + A U = 0
    // if A is assumed to have only real positive eigenvalues
    // explicit time schemes have a stability condition (aka CFL)
    // dt <= alpha / sqrt(Lambda_{max}(A))
    // where Lambda_{max}(A) is the maximum eigenvalue (in modulus) of A
    // This function returns the coefficient alpha

    CheckPresenceScheme();    
    return time_scheme->GetCFL(real_p);
  }
  
  
  //! returns the stability function associated with the time scheme
  template<class T>
  UnivariatePolynomial<Real_wp> All_TimeScheme<T>::GetStabilityFunction() const
  {
    CheckPresenceScheme();
    return time_scheme->GetStabilityFunction();
  }
  

  //! gives initial condition to the selected scheme
  /*!
    \param[in] t0 initial time
    \param[in] dt time step
    \param[inout] U0 on input, initial condition (cleared on output)
    \param[inout] sys system of equations
   */
  template<class T>
  void All_TimeScheme<T>
  ::SetInitialCondition(const Real_wp& t0, const Real_wp& dt, Vector<T>& U0,
                        VirtualOdeSystem<T>& sys)
  {
    CheckPresenceScheme();
    time_scheme->SetInitialCondition(t0, dt, U0, sys);
  }

  
  //! gives initial condition to the selected scheme
  /*!
    \param[in] t0 initial time
    \param[in] dt time step
    \param[inout] U0 on input, initial condition (cleared on output)
    \param[inout] V0 on input additional initial condition
    \param[inout] sys system of equations
    V0 can the initial condition for du/dt (for second-order formulation)
    or the initial condition for v (for split first-order formulation)
   */
  template<class T>
  void All_TimeScheme<T>
  ::SetInitialCondition(const Real_wp& t0, const Real_wp& dt, Vector<T>& U0,
                        Vector<T>& V0, VirtualOdeSystem<T>& sys)
  {
    CheckPresenceScheme();
    time_scheme->SetInitialConditionS(t0, dt, U0, V0, sys);
  }
  
  
  
  //! computes u^(n+1) from u^n with the selected scheme
  template<class T>
  void All_TimeScheme<T>::Advance(const Real_wp& z, int n, VirtualOdeSystem<T>& sys)
  {    
    CheckPresenceScheme();
    time_scheme->Advance(z, n, sys);
  }
  
  
  //! returns the current iterate u^n
  template<class T>
  Vector<T>& All_TimeScheme<T>::GetIterate()
  {
    return time_scheme->GetIterate();
  }

  
  //! returns the current iterate u^n
  template<class T>
  const Vector<T>& All_TimeScheme<T>::GetIterate() const
  {
    return time_scheme->GetIterate();
  }

  
  //! returns the unknown v
  template<class T>
  Vector<T>& All_TimeScheme<T>::GetVectorialIterate()
  {
    return time_scheme->GetVectorialIterate();
  }


  //! returns the number of iterates stored by the current time scheme
  template<class T>
  int All_TimeScheme<T>::GetNumberOfIterates() const
  {
    return nb_iterates;
  }
  
  
  //! returns the k- th iterate
  template<class T>
  Vector<T>& All_TimeScheme<T>::GetIterate(int k)
  {
    return time_scheme->GetIterate(k);
  }
  
  
  //! releases memory used by current time scheme
  template<class T>
  void All_TimeScheme<T>::Clear()
  {
    if (time_scheme != NULL)
      {
	time_scheme->Clear();
	delete time_scheme;
	time_scheme = NULL;
      }
  }
  
  
  //! releases memory used by current time scheme
  template<class T>
  void All_TimeScheme<T>::Clear(Vector<T>& U)
  {
    if (time_scheme != NULL)
      time_scheme->ClearFirst(U);
  }
  

  //! releases memory used by current time scheme
  template<class T>
  void All_TimeScheme<T>::Clear(Vector<T>& U, Vector<T>& V)
  {
    if (time_scheme != NULL)
      time_scheme->ClearSecond(U, V);
  }
  
}

#define MONTJOIE_FILE_TIME_SCHEMES_CXX
#endif  
