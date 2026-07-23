#ifndef MONTJOIE_FILE_TIME_SPECTRAL_SCHEME1D_CXX

namespace Montjoie
{
  
  //! default constructor
  template<class TypeEquation>
  TimeSpectralScheme1D<TypeEquation>::TimeSpectralScheme1D()
  {
    type_scheme = TimeSchemeEnum::RUNGE_KUTTA;
    order_time_scheme = 4;
    deltat = 0.0;
    advance_in_real_space = true;
    //advance_in_real_space = false;
  }
  
  
  //! modification of parameters of the scheme with a line of the data file
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "UnknownToAdvance")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of TimeSpectralScheme1D" << endl;
	    cout << "UnknownToAdvance needs more parameters, for instance :" << endl;
	    cout << "UnknownToAdvance = Spectral" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "Spectral")
          advance_in_real_space = false;
        else
          advance_in_real_space = true;
      }
  }
  
  
  //! sets the number of points and interval of computation
  /*!
    This method sets the number of points on which the solution is computed
    and the computational interval [a, b]
   */
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::SetInterval(const Real_wp& a, const Real_wp& b, int nelt, int r)
  {
    int N = nelt*(r+1);
    xmin = a;
    xmax = b;
    nb_points_fft = N;
    
    Real_wp dx = (b-a) / N;
    Xdof.Reallocate(N);
    for (int i = 0; i < N; i++)
      Xdof(i) = xmin + i*dx;
    
    fft.Init(N);
    
    nu_fft.Reallocate(N/2+1);
    //Real_wp bminus = b - dx;
    Real_wp bminus = b;
    for (int i = 0; i < N/2+1; i++)
      nu_fft(i) = 2.0*pi_wp*Real_wp(i)/(bminus-a);     
  }
  
  
  //! initialisation of coefficients of the equation
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>::SetOrder(int r)
  {
    var_eq.ComputePhysicalProperty(fft, nu_fft, xmin, xmax, nb_points_fft);
  }

  
  //! Computation of the interpolation grid
  /*!
    \param[in] a first extremity of the interval
    \param[in] b second extremity of the interval
    \param[in] N number of points on the interpolation grid
    Each point of the interpolation grid is located on the computational interval
   */
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::InitGrid(const Real_wp& a, const Real_wp& b, int N)
  {
    Mesh<Dimension1> mesh;
    mesh.CreateRegularMesh(xmin, xmax, nb_points_fft+1, 1);
    grid.Init(a, b, N);
    grid.LocalizePoints(mesh);
  }
  
  
  //! computation of coefficients for different time schemes
  /*!
    \param[in] type_time_scheme time scheme used for the evolution
    \param[in] order order of the time scheme
    \param[in] dt time step
    \param[in] Tf final time
    Initialization of time schemes with the initial condition
   */
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::InitTimeScheme(int type_scheme_, int order_, const Real_wp& dt, const Real_wp& Tf)
  {
    type_scheme = type_scheme_;
    order_time_scheme = order_;
    int Nvol = nb_points_fft;
    int Nvol_c = nu_fft.GetM();
    int nodl = Nvol*TypeEquation::nb_unknowns;
    int nodl_c = Nvol_c*TypeEquation::nb_unknowns;
    VectReal_wp Y0(nodl), x(Nvol);
    VectComplex_wp Y0_hat(nodl_c), y(Nvol_c);
    Y0.Fill(0); Y0_hat.Fill(0);
    this->var_eq.ComputeInitialCondition(this->xmin, this->xmax, this->Xdof, Y0);
    
    deltat = dt;
    if ( (type_scheme == TimeSchemeEnum::LEAP_FROG_EXACT)
         || (type_scheme == TimeSchemeEnum::RUNGE_KUTTA_EXACT) )
      this->var_eq.InitSinus(deltat, nu_fft);
    
    if (advance_in_real_space)
      {
	switch (type_scheme)
	  {
	  case TimeSchemeEnum::LEAP_FROG :
	  case TimeSchemeEnum::LEAP_FROG_EXACT :
	    {
	      // basic RK-scheme to get U^1 from U^0
	      Real_wp t = 0.0;
	      Ur_time.Reallocate(nodl);
	      Ur_next.Reallocate(nodl);
	      Ur_prev.Reallocate(nodl);
	      EvaluateFunction(t, Y0, Ur_time);
	      Copy(Ur_time, Ur_next);
	      Copy(Y0, Ur_prev);
	      Add(0.5*dt, Ur_time, Ur_prev);
	      
	      EvaluateFunction(t+0.5*dt, Ur_prev, Ur_time);
	      Add(Real_wp(2.0), Ur_time, Ur_next);
	      Copy(Y0, Ur_prev);
	      Add(0.5*dt, Ur_time, Ur_prev);
	      
	      EvaluateFunction(t+0.5*dt, Ur_prev, Ur_time);
	      Add(Real_wp(2.0), Ur_time, Ur_next);
	      Copy(Y0, Ur_prev);
	      Add(dt, Ur_time, Ur_prev);
	      
	      EvaluateFunction(t+dt, Ur_prev, Ur_time);
	      Add(Real_wp(1.0), Ur_time, Ur_next);
	      
	      for (int i = 0; i < nodl; i++)
		{
		  Ur_prev(i) = Y0(i);
		  Ur_time(i) = Y0(i) + dt/6.0*Ur_next(i);
		}
	    }
	    break;	    	  
	  case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
            {
              Ur_half.Reallocate(nodl); Ur_half.Fill(0);
              rhs.Reallocate(nodl); rhs.Fill(0);
              Ur_time.Reallocate(nodl); Ur_time.Fill(0);
              AhUn.Reallocate(order_time_scheme/2*nodl); AhUn.Fill(0);
              
              fd_scheme.var_eq.Copy(this->var_eq);
              fd_scheme.SetInterval(xmin, xmax, nb_points_fft);
              
              fd_scheme.SetOrder(2);
              fd_scheme.InitTimeScheme(type_scheme, order_time_scheme, dt, Tf);
              
              GaussRK_scheme.SetOrder(order_time_scheme);
	      GaussRK_scheme.SetInitialCondition(0.0, dt, Y0, *this); 
            }
            break;
	  case TimeSchemeEnum::RUNGE_KUTTA :
	  case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	    {
	      RK_scheme_real.SetOrder(order_time_scheme);
	      RK_scheme_real.SetInitialCondition(0.0, dt, Y0, *this);
	    }
	    break;
	  case TimeSchemeEnum::ADAMS_BASHFORTH :
	    {
	      AB_scheme_real.SetOrder(order_time_scheme);
	      AB_scheme_real.SetInitialCondition(0.0, dt, Y0, *this);
	    }
	    break;
	  case TimeSchemeEnum::TALEZER :
	    {
	      talezer_scheme_real.SetOrder(order_time_scheme, dt);
	      talezer_scheme_real.SetInitialCondition(0.0, dt, Y0, *this);
	    }
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	    {
	      LowRK_scheme_real.SetOrder(order_time_scheme);
	      LowRK_scheme_real.SetInitialCondition(0.0, dt, Y0, *this);
	    }
	    break;
          default :
            {
              cout << "Unknown time scheme" << endl;
              abort();
            }
	  }
      }
    else
      {
	// computation of initial condition in spectral domain
	for (int k = 0; k < TypeEquation::nb_unknowns; k++)
	  {
	    for (int i = 0; i < Nvol; i++)
	      x(i) = Y0(i+k*Nvol);
	    
	    fft.ApplyForward(x, y);
	    
	    for (int i = 0; i < Nvol_c; i++)
	      Y0_hat(i + Nvol_c*k) = y(i);
	  }
        
	// initialization of time schemes
	switch (type_scheme)
	  {
	  case TimeSchemeEnum::LEAP_FROG :
	  case TimeSchemeEnum::LEAP_FROG_EXACT :
	    {
	      // basic RK-scheme to get U^1 from U^0
	      Real_wp t = 0.0;
	      Un_time.Reallocate(nodl_c);
	      Un_next.Reallocate(nodl_c);
	      Un_prev.Reallocate(nodl_c);
	      EvaluateFunction(t, Y0_hat, Un_time);
	      Copy(Un_time, Un_next);
	      Copy(Y0_hat, Un_prev);
	      Add(0.5*dt, Un_time, Un_prev);
	      
	      EvaluateFunction(t+0.5*dt, Un_prev, Un_time);
	      Add(Complex_wp(2.0), Un_time, Un_next);
	      Copy(Y0_hat, Un_prev);
	      Add(0.5*dt, Un_time, Un_prev);
	      
	      EvaluateFunction(t+0.5*dt, Un_prev, Un_time);
	      Add(Complex_wp(2.0), Un_time, Un_next);
	      Copy(Y0_hat, Un_prev);
	      Add(dt, Un_time, Un_prev);
	      
	      EvaluateFunction(t+dt, Un_prev, Un_time);
	      Add(Complex_wp(1.0), Un_time, Un_next);
	      
	      for (int i = 0; i < nodl_c; i++)
		{
		  Un_prev(i) = Y0_hat(i);
		  Un_time(i) = Y0_hat(i) + dt/6.0*Un_next(i);
		}
	    }
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA :
	  case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	    {
	      RK_scheme.SetOrder(order_time_scheme);
	      RK_scheme.SetInitialCondition(0.0, dt, Y0_hat, *this);
	    }
	    break;
	  case TimeSchemeEnum::ADAMS_BASHFORTH :
	    {
	      AB_scheme.SetOrder(order_time_scheme);
	      AB_scheme.SetInitialCondition(0.0, dt, Y0_hat, *this);
	    }
	    break;
	  case TimeSchemeEnum::TALEZER :
	    {
	      talezer_scheme.SetOrder(order_time_scheme, dt);
	      talezer_scheme.SetInitialCondition(0.0, dt, Y0_hat, *this);
	    }
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	    {
	      LowRK_scheme.SetOrder(order_time_scheme);
	      LowRK_scheme.SetInitialCondition(0.0, dt, Y0_hat, *this);
	    }
	    break;
          default :
            {
              cout << "Unknown time scheme" << endl;
              abort();
            }
	  }
      }
  }
  

  //! factorisation of M - dt a_ij df/dx_j(y^n + dt \sum a_ij k_j)
  template<class TypeEquation> void TimeSpectralScheme1D<TypeEquation>::
  FactorizeOperatorDhPlusAijKh(const Matrix<Real_wp>& A, const Real_wp& dt,
                               const VectReal_wp& Y, const VectReal_wp& ki)
  {
    fd_scheme.FactorizeOperatorDhPlusAijKh(A, dt, Y, ki);
    return;
    
    /*fd_scheme.mat_sparse.WriteText("matDiff.dat");
    
    VectReal_wp kx(ki), fplus(ki), fminus(ki);
    Real_wp h = 1e-6;
    diff_matrix.Reallocate(ki.GetM(), ki.GetM());
    for (int jcol = 0; jcol < ki.GetM(); jcol++)
      {
        kx(jcol) += h;
        
        for (int i = 0; i < A.GetM(); i++)
          {
            for (int k = 0; k < Y.GetM(); k++)
              {
                Ur_half(k) = Y(k);
                Ur_time(k) = kx(k*A.GetM() + i);
              }
            
            for (int j = 0; j < A.GetM(); j++)
              for (int k = 0; k < Y.GetM(); k++)
                Ur_half(k) += dt*A(i, j)*kx(k*A.GetM() + j);
            
            this->var_eq.ApplyMass(Ur_time, fft, nu_fft, rhs);
            this->var_eq.AddFunction(-1.0, Ur_half, fft, nu_fft, rhs, false);
            
            for (int k = 0; k < Y.GetM(); k++)
              fplus(k*A.GetM() + i) = rhs(k);
          }
        
        kx(jcol) -= 2.0*h;
        
        for (int i = 0; i < A.GetM(); i++)
          {
            for (int k = 0; k < Y.GetM(); k++)
              {
                Ur_half(k) = Y(k);
                Ur_time(k) = kx(k*A.GetM() + i);
              }
            
            for (int j = 0; j < A.GetM(); j++)
              for (int k = 0; k < Y.GetM(); k++)
                Ur_half(k) += dt*A(i, j)*kx(k*A.GetM() + j);
            
            this->var_eq.ApplyMass(Ur_time, fft, nu_fft, rhs);
            this->var_eq.AddFunction(-1.0, Ur_half, fft, nu_fft, rhs, false);
            
            for (int k = 0; k < Y.GetM(); k++)
              fminus(k*A.GetM() + i) = rhs(k);
          }
        
        for (int irow = 0; irow < ki.GetM(); irow++)
          {
            Real_wp df = (fplus(irow) - fminus(irow)) / (2.0*h);
            diff_matrix(irow, jcol) = df;
          }
        
        kx(jcol) += h;
      }   
    
    diff_matrix.WriteText("matSp.dat");
    GetLU(diff_matrix, pivot); */
  }
  
  
  //! resolution of M ki - f(y^n + dt \sum a_ij k_j) = 0
  template<class TypeEquation>  
  void TimeSpectralScheme1D<TypeEquation>::
  SolveOperatorDhPlusAijKh(const VectReal_wp& tsub,
                           const Matrix<Real_wp>& A, const Real_wp& dt,
                           const VectReal_wp& Y, VectReal_wp& ki)
  {
    // using Newton's method to solve non-linear system
    //Real_wp test = 1e30, test_prec = 2e30;
    Real_wp test = 1e30;
    int nb_iter = 0;
    Real_wp threshold = 100.0*epsilon_machine, norme_init(1);
    ki.Fill(0);
    // on continue tant que || g(Un_next) || > threshold 
    // et tant que la suite des residus est bien decroissante
    //while ((test > threshold)&&(test < 0.9*test_prec))
    while (test > threshold)
      {
        // calcul fonctionelle
        for (int i = 0; i < A.GetM(); i++)
          {
            for (int k = 0; k < Y.GetM(); k++)
              {
                Ur_half(k) = Y(k);
                Ur_time(k) = ki(k*A.GetM() + i);
              }
            
            for (int j = 0; j < A.GetM(); j++)
              for (int k = 0; k < Y.GetM(); k++)
                Ur_half(k) += dt*A(i, j)*ki(k*A.GetM() + j);
            
            this->var_eq.ApplyMass(Ur_time, fft, nu_fft, rhs);
            this->var_eq.AddFunction(-1.0, Ur_half, fft, nu_fft, rhs, false);
            
            for (int k = 0; k < Y.GetM(); k++)
              AhUn(k*A.GetM() + i) = rhs(k);
          }
        
        //test_prec = test;
        if (nb_iter == 0)
          {
            norme_init = Norm2(AhUn);
            test = 1.0;
          }
        else
          test = Norm2(AhUn)/norme_init;
        
        cout << "residu at iteration " << nb_iter << " = " << test << endl;
        //cout << "ratio = " << test_prec/test << endl;
        if (fd_scheme.DifferentialToBeComputed())
          {
            FactorizeOperatorDhPlusAijKh(A, dt, Y, ki);
            fd_scheme.SetDifferentialToCompute(false);
            fd_scheme.IncrementDifferential();
          }
        
        //SolveLU(diff_matrix, pivot, AhUn);
        fd_scheme.SolveGaussJacobian(AhUn, A.GetM());
        
        // Newton iterate : X^n+1 = X^n - DG^-1 G(X^n)
        for (int i = 0; i < ki.GetM(); i++)
          ki(i) -= AhUn(i);
        
        nb_iter++;
      }
    
    exit(0);
    // if the number of iterations is large, we consider that the 
    // jacobian should be recomputed
    if (nb_iter > 10)
      fd_scheme.SetDifferentialToCompute(true);
  }
  
  
  //! Evaluation of g(t, u) where u is solution of du/dt = g(t, u)
  template<class TypeEquation>  
  void TimeSpectralScheme1D<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                               const VectComplex_wp& Xn,
                               VectComplex_wp& ProdXn, bool invert_mass, bool source)
  {
    var_eq.EvaluateFunction(tn, Xn, fft, nu_fft, ProdXn);
  }
  
  
  //! Evaluation of g(t, u) where u is solution of du/dt = g(t, u)
  template<class TypeEquation>  
  void TimeSpectralScheme1D<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn, const VectComplex_wp& Xn,
                     VectComplex_wp& ProdXn, bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, 0, Xn, ProdXn);
  }


  //! Evaluation of g(t, u) where u is solution of du/dt = g(t, u)
  template<class TypeEquation>  
  void TimeSpectralScheme1D<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                               const VectReal_wp& Xn,
                               VectReal_wp& ProdXn, bool invert_mass, bool source)
  {
    var_eq.EvaluateFunction(tn, Xn, fft, nu_fft, ProdXn);
  }
  
  
  //! Evaluation of g(t, u) where u is solution of du/dt = g(t, u)
  template<class TypeEquation>  
  void TimeSpectralScheme1D<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn, const VectReal_wp& Xn,
                     VectReal_wp& ProdXn, bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, 0, Xn, ProdXn);
  }

  
  //! advancing scheme in time
  /*!
    \param[in] t time step
    \param[in] n iteration number
    This methods computes u^n+1 from u^n
   */
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::Advance(const Real_wp& t, int n)
  {
    if (advance_in_real_space)
      {
	switch (type_scheme)
	  {
	  case TimeSchemeEnum::LEAP_FROG :
	  case TimeSchemeEnum::LEAP_FROG_EXACT :
	    {
	      EvaluateFunction(t, Ur_time, Ur_next);
	      Real_wp coef = 2.0*deltat;
	      for (int i = 0; i < Ur_next.GetM(); i++)
		Ur_next(i) = Ur_prev(i) + coef*Ur_next(i);
	      
	      Copy(Ur_time, Ur_prev);
	      Copy(Ur_next, Ur_time);
	    }
	    break;
	  case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
            GaussRK_scheme.Advance(t, n, *this);
            break;
	  case TimeSchemeEnum::RUNGE_KUTTA :
	  case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	    RK_scheme_real.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::ADAMS_BASHFORTH :
	    AB_scheme_real.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::TALEZER :
	    talezer_scheme_real.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	    LowRK_scheme_real.Advance(t, n, *this);
	    break;
	  }
      }
    else
      {
	switch (type_scheme)
	  {
	  case TimeSchemeEnum::LEAP_FROG :
	  case TimeSchemeEnum::LEAP_FROG_EXACT :
	    {
	      EvaluateFunction(t, Un_time, Un_next);
	      Real_wp coef = 2.0*deltat;
	      for (int i = 0; i < Un_next.GetM(); i++)
		Un_next(i) = Un_prev(i) + coef*Un_next(i);
	      
	      Copy(Un_time, Un_prev);
	      Copy(Un_next, Un_time);
	    }
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA :
	  case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	    RK_scheme.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::ADAMS_BASHFORTH :
	    AB_scheme.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::TALEZER :
	    talezer_scheme.Advance(t, n, *this);
	    break;
	  case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	    LowRK_scheme.Advance(t, n, *this);
	    break;
	  }
      }
  }
  
  
  //! returns the iterate in spectral domain
  template<class TypeEquation>  
  const VectComplex_wp& TimeSpectralScheme1D<TypeEquation>::GetIterate() const
  {
    switch (type_scheme)
      {
      case TimeSchemeEnum::LEAP_FROG :
      case TimeSchemeEnum::LEAP_FROG_EXACT :
	return Un_prev;
      case TimeSchemeEnum::RUNGE_KUTTA :
      case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	return RK_scheme.GetIterate();
      case TimeSchemeEnum::ADAMS_BASHFORTH :
        return AB_scheme.GetIterate();
      case TimeSchemeEnum::TALEZER :
        return talezer_scheme.GetIterate();
      case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	return LowRK_scheme.GetIterate();
      }
    
    return Un_prev;
  }


  //! returns the iterate in real domain
  template<class TypeEquation>  
  const VectReal_wp& TimeSpectralScheme1D<TypeEquation>::GetIterateReal() const
  {
    switch (type_scheme)
      {
      case TimeSchemeEnum::LEAP_FROG :
      case TimeSchemeEnum::LEAP_FROG_EXACT :
	return Ur_prev;
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	return GaussRK_scheme.GetIterate();
      case TimeSchemeEnum::RUNGE_KUTTA :
      case TimeSchemeEnum::RUNGE_KUTTA_EXACT :
	return RK_scheme_real.GetIterate();
      case TimeSchemeEnum::ADAMS_BASHFORTH :
        return AB_scheme_real.GetIterate();
      case TimeSchemeEnum::TALEZER :
        return talezer_scheme_real.GetIterate();
      case TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE :
	return LowRK_scheme_real.GetIterate();
      }
    
    return Ur_prev;
  }
  
  
  //! computes the interpolation of u on points of the interpolation grid
  template<class TypeEquation>
  void TimeSpectralScheme1D<TypeEquation>
  ::GetInterpolateUn(int m, VectReal_wp& ZnInterp)
  {    
    ZnInterp.Reallocate(grid.GetNbPointsGrid());
    if (grid.GetNbPointsGrid() == nb_points_fft+1)
      {
	if (advance_in_real_space)
	  for (int i = 0; i < nb_points_fft; i++)
	    ZnInterp(i) = this->GetIterateReal()(i);
	else
	  fft.ApplyInverse(this->GetIterate(), ZnInterp);
	
	// periodic condition
	ZnInterp(nb_points_fft) = ZnInterp(0);
      }
    else
      {
	VectReal_wp Zn(nb_points_fft);
	if (advance_in_real_space)
	  Copy(this->GetIterateReal(), Zn);
	else
	  fft.ApplyInverse(this->GetIterate(), Zn);
	
	// use of a fourth-order local interpolation
	Real_wp x;
	TinyVector<Real_wp, 5> phi;
	int N = nb_points_fft;
	for (int i = 0; i < grid.GetNbPointsGrid(); i++)
	  {
	    int j = grid.GetElementNumber(i);
	    x = 0.5 + 0.25*grid.GetLocalCoordinate(i);
	    phi(0) = 32.0/3*(x-0.25)*(x-0.5)*(x-0.75)*(x-1.0);
	    phi(1) = -128.0/3*x*(x-0.5)*(x-0.75)*(x-1.0);
	    phi(2) = 64.0*x*(x-0.25)*(x-0.75)*(x-1.0);
	    phi(3) = -128.0/3*x*(x-0.25)*(x-0.5)*(x-1.0);
	    phi(4) = 32.0/3*x*(x-0.25)*(x-0.5)*(x-0.75);
	    
	    int jm2 = (j+N-2)%N;
	    int jm1 = (j+N-1)%N;
	    int jp1 = (j+1)%N;
	    int jp2 = (j+2)%N;
	    //DISP(i); 
	    //DISP(j); DISP(x); DISP(jm1); DISP(jm2); DISP(jp1); DISP(jp2);
	    
	    ZnInterp(i) = phi(0)*Zn(jm2) + phi(1)*Zn(jm1) + phi(2)*Zn(j) + 
	      phi(3)*Zn(jp1) + phi(4)*Zn(jp2); 
	  }
      }
  }

}

#define MONTJOIE_FILE_TIME_SPECTRAL_SCHEME1D_HXX
#endif

