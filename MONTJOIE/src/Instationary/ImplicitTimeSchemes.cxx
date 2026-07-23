#ifndef MONTJOIE_FILE_IMPLICIT_TIME_SCHEMES_CXX

namespace Montjoie
{
  
  /***********************
   * TetaScheme_Iterator *
   ***********************/

  
  //! default constructor
  template<class T>
  TetaScheme_Iterator<T>::TetaScheme_Iterator()
  {
    teta = 0.25; dt = 0; dt2 = 0;
    fourth_order_init = false;
  }
    
  
  template<class T>
  void TetaScheme_Iterator<T>::SetTheta(const Real_wp& theta)
  {
    teta = theta; 
  }
  
  
  template<class T>
  Real_wp TetaScheme_Iterator<T>::GetTheta() const
  {
    return teta;
  }


  template<class T>
  void TetaScheme_Iterator<T>::SetFourthOrderInitialScheme()
  {
    fourth_order_init = true;
  }

  
  //! initialisation of theta-scheme with initial conditions u0 and du0/dt
  /*!
    \param[in] t initial time
    \param[in] deltat time step
    \param[inout] u0 on input initial vector u(t=0), on output this vector is cleared
    \param[inout] du0_dt on input initial vector du/dt(t=0), on output this vector is cleared
    \param[inout] sys linear system considered
    we consider the following equation
    Dh d^2 u/dt^2 + Sh du/dt - Kh u = F(t)
    discretized by the following theta-scheme
    Dh (U^{n+1} - 2 U^n + U^{n-1})/dt^2 + Sh (U^{n+1} - U^{n-1}) / (2 dt)
    - Kh (\theta U^{n-1} + (1-2 \theta) U^n + \theta U^{n+1}) \, = \, F^n
    This method computes U^1 from U^0 and dU^0/dt with an explicit scheme
  */
  template<class T>
  void TetaScheme_Iterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                         Vector<T>& du0_dt, VirtualOdeSystem<T>& sys)
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    dt = deltat; dt2 = dt*dt;
    int N = u0.GetM();
    U0 = u0; u0.Clear();
    U1.Reallocate(N); U1.Fill(zero);
    Ah_u1.Reallocate(N); Ah_u1.Fill(zero);
    Ah_u0.Reallocate(N); Ah_u0.Fill(zero);
    Prod_Uh.Reallocate(N); Prod_Uh.Fill(zero);

    sys.SetDirichletCondition(t, 0, U0);
    
    // first step, we use basic second order stencil to compute U^1
    // (we assume that potential sources are null near t = initial_time)
    // basic second order initialization :
    // U^1 = U^0 + dt dU/dt(0) + dt^2/2 d^2 U/dt^2(0)
    // computing D_h d^2 U /dt^2 = Kh U - Sh dU/dt
    sys.ApplyOperatorKh(Real_wp(1), t, U0, Real_wp(0), U1);
    Copy(U1, Ah_u1); 
    sys.ApplyOperatorSh(-Real_wp(1), t, du0_dt, Real_wp(1), U1);
    sys.SolveOperatorDh(U1);
    
    // not working for non-linear system
    if (fourth_order_init)
      {	    
	// evaluating d^3 U/dt^3 in Ah_u0
	sys.ApplyOperatorKh(Real_wp(1), t, du0_dt, Real_wp(0), Ah_u0);
	sys.ApplyOperatorSh(-Real_wp(1), t, U1, Real_wp(1), Ah_u0);
	sys.SolveOperatorDh(Ah_u0);
	
	// evaluating d^4 U/dt^4 in Prod_Uh
	sys.ApplyOperatorKh(Real_wp(1), t, U1, Real_wp(0), Prod_Uh);
	sys.ApplyOperatorSh(-Real_wp(1), t, U1, Real_wp(1), Prod_Uh);
	sys.SolveOperatorDh(Prod_Uh);
	
	// then forming U1
	Mlt(dt/4, Prod_Uh); Add(one, Ah_u0, Prod_Uh);
	Mlt(dt/3, Prod_Uh); Add(one, U1, Prod_Uh);
	Mlt(dt/2, Prod_Uh); Add(one, du0_dt, Prod_Uh);
	Mlt(dt, Prod_Uh); Add(one, U0, Prod_Uh);
	Copy(Prod_Uh, U1);
      }
    else
      {
	// then forming U1
	Mlt(T(0.5*dt), U1); Add(one, du0_dt, U1);
	Mlt(dt, U1); Add(one, U0, U1);	
      }
    
    sys.SetDirichletCondition(t+dt, 0, U1);
     
    // we can clear dU0_dt
    du0_dt.Clear();
    
    // factorization of linear system
    sys.FactorizeOperatorDhPlusGammaKh(1.0, 0.5*dt, dt*dt*teta);
  }
  

  //! releasing memory used by arrays for the time scheme
  template<class T>
  void TetaScheme_Iterator<T>::Clear()
  {
    Vector<T> Ufinal, dUfinal; ClearSecond(Ufinal, dUfinal);
  }
  
  
  //! releasing memory used by arrays for the time scheme
  template<class T>
  void TetaScheme_Iterator<T>::ClearSecond(Vector<T>& Ufinal, Vector<T>& dUfinal)
  {
    Ah_u1.Clear(); Ah_u0.Clear(); Prod_Uh.Clear(); U1.Clear();
    Ufinal = U0;     U0.Clear(); dUfinal.Clear();
  }
  
  
  //! computation of U^n+1 from U^n
  template<class T>
  void TetaScheme_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    T zero, one, two;
    SetComplexZero(zero);
    SetComplexOne(one);
    two = Real_wp(2)*one;
    
    // scheme
    // D_h (U^{n+1} - 2 U^n + U^{n-1}) / dt2 
    // - K_h (teta U^{n+1} + (1 - 2 teta) U^n + teta U^{n-1}) 
    // + \Sigma_h (U^{n+1} - U^{n-1}) / (2 dt) = F^n    
    // Ah_u0 = Ah U^{n-1}
    Copy(Ah_u1, Ah_u0);
       
    // computing Prod_Uh = 2 Dh U^n - (Dh - dt/2 Sh) U^{n-1} + teta dt^2 Kh U^n-1 
    //  + (1-2teta) Kh U^n + dt^2 F^n
    sys.ApplyOperatorKh(Real_wp(1), t+dt, U1, Real_wp(0), Ah_u1);
    sys.ApplyOperatorDh(Real_wp(2), t+dt, U1, Real_wp(0), Prod_Uh);
    sys.ApplyOperatorDhMinusdtSh(-Real_wp(1), t, U0, Real_wp(1), Prod_Uh);
    sys.AddScalarTimeSource(dt2, t+dt, 1, Prod_Uh);

    Add(teta*dt2, Ah_u0, Prod_Uh); Add((one-two*teta)*dt2, Ah_u1, Prod_Uh);
    
    // solving (Dh + dt/2 Sh - teta dt^2 Kh) U^{n+1} = Prod_Uh
    Copy(U1, U0);
    sys.SetDirichletConditionSource(t+2*dt, 0, Prod_Uh);
    sys.SolveOperatorDhPlusGammaKh(t+2*dt, Prod_Uh, U1); 

    // int test_input; cout << "waiting" << endl; cin >> test_input;
  }
  

  //! changes time step during the simulation
  template<class T>
  void TetaScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }
  
  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& TetaScheme_Iterator<T>::GetIterate()
  {
    return U0;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& TetaScheme_Iterator<T>::GetIterate() const
  {
    return U0;
  }
  

  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int TetaScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 3;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& TetaScheme_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return U1;
    else if (k == 1)
      return U0;
    else
      return Ah_u1;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t TetaScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = U1.GetMemorySize() + Ah_u1.GetMemorySize() + Ah_u0.GetMemorySize()
      + Prod_Uh.GetMemorySize() + U0.GetMemorySize();
    
    return taille;
  }

  
  /************************
   * SdirkScheme_Iterator *
   ************************/


  //! default constructor
  template<class T>
  SdirkScheme_Iterator<T>::SdirkScheme_Iterator()
  {
    SetOrder(2, false, 0);
  }
  
  
  //! computation of coefficients for SDIRK scheme of order r
  template<class T>
  void SdirkScheme_Iterator<T>::SetOrder(int r, bool esdirk_scheme, int type)
  {
    Real_wp one(1), zero(0);
    gamma = zero;
    alternative_type = type;
    explicit_first_stage = esdirk_scheme;    
    explicit_first_stage = false;
    switch (r)
      {
      case 2 :
	{
	  // Crank-Nicolson scheme
	  gamma = one/2;
	  A.Reallocate(1, 1);
	  A(0, 0) = one/2;
	  c.Reallocate(1); c(0) = one/2;
	  b.Reallocate(1); b(0) = one;
	}
	break;
      case 3 :
	{
          // basic SDIRK method of order 3 (in volume 1 of Hairer)
          gamma = (Real_wp(3) + sqrt(Real_wp(3)))/6;
          A.Reallocate(2, 2);
          A(0, 0) = gamma;
          A(1, 1) = gamma;
          A(1, 0) = one - 2*gamma;
          b.Reallocate(2); b(0) = one/2; b(1) = one/2;
          c.Reallocate(2); c(0) = gamma; c(1) = one-gamma;
          break;          
	}
	break;
      case 4 :
	{
          // Crouzeix-Raviart, 1980 (cited in Hairer's book)
          gamma = one/sqrt(Real_wp(3))*cos(pi_wp/18) + one/2;
          Real_wp delta = one/(6*square(2*gamma-one));
          A.Reallocate(3, 3);
          A(0, 0) = gamma;
          A(1, 0) = one/2-gamma;
          A(1, 1) = gamma;
          A(2, 0) = 2*gamma;
          A(2, 1) = one-4*gamma;
          A(2, 2) = gamma;
          c.Reallocate(3); b.Reallocate(3);
          b(0) = delta; b(1) = one-2*delta; b(2) = delta;
          c(0) = gamma; c(1) = one/2; c(2) = one-gamma;
	  
	  /*/ My opt sdirk in range [0, 2pi]
	    gamma = Real_wp(0.1289);
	    A.Reallocate(3, 3);
	    A(0, 0) = gamma;
	    A(1, 0) = Real_wp(0.3711);
	    A(1, 1) = gamma;
	    A(2, 0) = Real_wp(0.2578);
	    A(2, 1) = Real_wp(0.4845);
	    A(2, 2) = gamma;
	    c.Reallocate(3); b.Reallocate(3);
	    b(0) = Real_wp(0.3025); b(1) = Real_wp(0.3949); b(2) = Real_wp(0.3025);
	    c(0) = gamma; c(1) = A(1,0) + gamma; c(2) = A(2,0) + A(2,1) + gamma;*/

	  
	  /*/L. M. Skvortsov SDIRK54
	    //Real_wp six = Real_wp(6);
	    //Real_wp sqrt_six = sqrt(six);
	    gamma = one/4;
	    A.Reallocate(5, 5);
	    A(0, 0) = gamma;
	    A(1, 0) = -gamma;
	    A(1, 1) = gamma;
	    A(2, 0) = one/8;
	    A(2, 1) = one/8;
	    A(2, 2) = gamma;
	    A(3, 0) = Real_wp(-3)/2;
	    A(3, 1) = Real_wp(3)/4;
	    A(3, 2) = Real_wp(3)/2;
	    A(3, 3) = gamma;
	    A(4, 0) = zero;
	    A(4, 1) = one/6;
	    A(4, 2) = Real_wp(2)/3;
	    A(4, 3) = -one/12;
	    A(4, 4) = gamma;

	    c.Reallocate(5); b.Reallocate(5);
	    c(0) = gamma; c(1) = zero; c(2) = one/2;
	    c(3) = one; c(4) = one;

	    b(0) = zero; b(1) = one/6; b(2) = Real_wp(2) / 3;
	    b(3) = -one/12; b(4) = gamma;*/
	  
	}
	break;
      case 5 :
	{
	  //Cooper and Sayfy SDIRK55
	  Real_wp six = Real_wp(6);
	  Real_wp sqrt_six = sqrt(six);
	  gamma = ( six - sqrt_six) / Real_wp(10);
	  A.Reallocate(5, 5);
	  A(0, 0) = gamma;
	  A(1, 0) = (-six + 5*sqrt_six) / 14;
	  A(1, 1) = gamma;
	  A(2, 0) = (Real_wp(888) + 607*sqrt_six) / 2850;
	  A(2, 1) = (Real_wp(126) - 161*sqrt_six) / 1425;
	  A(2, 2) = gamma;
	  A(3, 0) = (Real_wp(3153) - 3082*sqrt_six) / 14250;
	  A(3, 1) = (Real_wp(3213) + 1148*sqrt_six) / 28500;
	  A(3, 2) = (Real_wp(-267) + 88*sqrt_six) / 500;
	  A(3, 3) = gamma;
	  A(4, 0) = (Real_wp(-32583) + 14638*sqrt_six) / 71250;
	  A(4, 1) = (Real_wp(-17199) + 364*sqrt_six) / 142500;
	  A(4, 2) = (Real_wp(1329) - 544*sqrt_six) / 2500;
	  A(4, 3) = (Real_wp(-96) + 131*sqrt_six) / 625;
	  A(4, 4) = gamma;

	  c.Reallocate(5); b.Reallocate(5);
	  c(0) = gamma; c(1) = (six + 9*sqrt_six)/35; c(2) = one;
	  c(3) = (Real_wp(4) - sqrt_six) / 10;
	  c(4) = (Real_wp(4) + sqrt_six) / 10;

	  b(0) = 0; b(1) = 0; b(2) = one / 9;
	  b(3) = (Real_wp(16) - sqrt_six) / 36;
	  b(4) = (Real_wp(16) + sqrt_six) / 36;
	}
	break;
      case 6 :
	{
	}
	break;
      default :
	{
	  cout << "Order " << r << " Not implemented " <<endl;
	  abort();
	}
      }
  }
  
  
  //! initialisation of the scheme with u^0
  template<class T>
  void SdirkScheme_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                        Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();   
    sys.SetDirichletCondition(t, 0, Y);
    Fn.Reallocate(N); Fn.Fill(0);
    Un_tmp.Reallocate(N); Un_tmp.Fill(0);
    int nb_stages = A.GetM();
    Zn.Reallocate(nb_stages);
    for (int i = 0; i < nb_stages; i++)
      {
	Zn(i).Reallocate(N);
	Zn(i).Fill(0);
      }

    // we ask the user to factorize the matrix
    // M_h + dt gamma S_h + dt gamma K_h
    sys.FactorizeOperatorDhPlusGammaKh(1.0, dt*gamma, dt*gamma);
  }
  

  //! clearing arrays used by time scheme
  template<class T>
  void SdirkScheme_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! clearing arrays used by time scheme
  template<class T>
  void SdirkScheme_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Zn.Clear();
    Fn.Clear(); Un_tmp.Clear();
    Yfinal = Y; Y.Clear();
  }
  
  
  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void SdirkScheme_Iterator<T>::Advance(const Real_wp& t, int n_time,
					VirtualOdeSystem<T>& sys)
  {
    Real_wp tcurrent;
    int ilow = 0;
    if (explicit_first_stage)
      {
	ilow = 1;
	// computation of Fn = f(U^n, t^n)
	sys.EvaluateFunction(t, Y, Zn(0));
      }

    // computation of all intermediary unknowns
    for (int i = ilow; i < A.GetM(); i++)
      {
	// current time
	tcurrent = t + c(i)*dt;
	
	// computation of Un_tmp = U^n + delta t \sum_{j=1}^{i-1} a_ij z_j 
	Copy(Y, Un_tmp);
	for (int j = 0; j < i; j++)
	  Add(dt*A(i, j), Zn(j), Un_tmp);
        
	// computation of f(Un_tmp, tcurrent)
	sys.EvaluateFunction(tcurrent, Un_tmp, Fn, false);

	// then, we solve system z_i = f(Un_tmp + gamma dt z_i, t^n + c_i dt)
	sys.SolveOperatorDhPlusGammaKh(tcurrent, Fn, Zn(i));
      }
    
    // Once we have all the intermediary unknowns, we deduce
    // U^n+1 = U^n + dt \sum b_i z_i
    for (int i = 0; i < b.GetM(); i++)
      Add(dt*b(i), Zn(i), Y);

    //sys.SetDirichletCondition(t+dt, 0, Y);
    
    //int test_input; cout << "waiting... " << endl; cin >> test_input;
  }


  //! changes time step during the simulation
  template<class T>
  void SdirkScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& SdirkScheme_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& SdirkScheme_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int SdirkScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& SdirkScheme_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }
  

  //! returns the size used by the object in bytes
  template<class T>
  size_t SdirkScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = A.GetMemorySize() + b.GetMemorySize() + c.GetMemorySize();
    taille += Y.GetMemorySize() + Fn.GetMemorySize() + Un_tmp.GetMemorySize();
    taille += Seldon::GetMemorySize(Zn);
    return taille;
  }

  
  /**************************
   * SdirkNonLinearEquation *
   **************************/
  
  
  void SdirkNonLinearEquation
  ::Init(const UnivariatePolynomial<Real_wp>& P, const UnivariatePolynomial<Real_wp>& Q,
         const Real_wp& gam_, int order, Matrix<Real_wp>& A, VectReal_wp& b, VectReal_wp& c)
  {
    gamma = gam_;    
    Numer = P;
    r = order;
    s = Numer.GetOrder();

    //DISP(r); DISP(s);

    int nb_cond_ordre = Numer.GetM();

    if (r >= 3)
      nb_cond_ordre++;

    if (r >= 4)
      nb_cond_ordre += 3;

    if (r >= 5)
      nb_cond_ordre += 3;

    if (r >= 6)
      nb_cond_ordre += 4;

    if (r >= 7)
      nb_cond_ordre += 5;

    if (r >= 8)
      nb_cond_ordre += 6;

    int nb_unknowns = s*(s+1)/2;
    //DISP(nb_unknowns); DISP(nb_cond_ordre);
    
    this->m_ = nb_cond_ordre;
    this->n_ = nb_unknowns;

    VectReal_wp x(nb_unknowns); x.Zero();
    int nb = 0;
    for (int i = 0; i < s; i++)
      for (int j = 0; j < i; j++)
        x(nb++) = A(i, j);
    
    for (int i = 0; i < s; i++)
      x(nb++) = b(i);
    
    SolveLeastSquaresLvm(*this, x, 10*epsilon_machine, -1, 20000, 2, 1e-2);

    FillAB(x, A, b, c);
    //DISP(A); DISP(b); DISP(c);
    
    /* ofstream file_out("coefs.txt");
    file_out.precision(15);

    for (int i = 0; i < x.GetM(); i++)
      file_out << "x(" << i << ") = " << x(i) << ";\n";

      file_out.close(); */
  }
  
  
  void SdirkNonLinearEquation
  ::FillAB(const VectReal_wp& x, Matrix<Real_wp>& A, VectReal_wp& b, VectReal_wp& c)
  {
    c.Reallocate(s);
    A.Reallocate(s, s);
    A.Zero();
    c.Zero();

    for (int i = 0; i < s; i++)
      {
	c(i) = gamma;
	A(i, i) = gamma;
      }

    int nb = 0;
    for (int i = 0; i < s; i++)
      for (int j = 0; j < i; j++)
	{
	  A(i, j) = x(nb++);
	  c(i) += A(i, j);
	}
    
    b.Reallocate(s);
    for (int i = 0; i < s; i++)
      b(i) = x(nb++);
  }

  void SdirkNonLinearEquation
  ::GetStability(const Matrix<Real_wp>& A, const VectReal_wp& b,
                 UnivariatePolynomial<Real_wp>& N)
  {
    // monome = (1- gamma z) 
    UnivariatePolynomial<Real_wp> monome, z;
    monome.SetOrder(1);
    monome(0) = Real_wp(1); monome(1) = -gamma;
    
    z.SetOrder(1);
    z(0) = Real_wp(0); z(1) = Real_wp(1);

    // vec_mono(i) = (1 - gamma z)^i
    Vector<UnivariatePolynomial<Real_wp> > vec_mono(s+1);
    vec_mono(0).SetOrder(0); vec_mono(0)(0) = Real_wp(1);
    for (int i = 0; i < s; i++)
      vec_mono(i+1) = vec_mono(i)*monome;

    // on calcule le numerateur de la fonction de stabilite
    N = vec_mono(s);
    Vector<UnivariatePolynomial<Real_wp> > ki(s);
    for (int i = 0; i < s; i++)
      {
	ki(i) = vec_mono(i);
	for (int j = 0; j < i; j++)
	  ki(i) += A(i, j)*ki(j)*vec_mono(i-1-j);
	
	ki(i) *= z;
	N += b(i)*ki(i)*vec_mono(s-1-i);
      }
  }

  void SdirkNonLinearEquation::FindInitGuess(VectReal_wp& x)
  {
  }

  void SdirkNonLinearEquation::EvaluateF(const VectReal_wp& x, VectReal_wp& f)
  {
    Matrix<Real_wp> A; VectReal_wp b, c;
    
    FillAB(x, A, b, c);

    UnivariatePolynomial<Real_wp> N;
    GetStability(A, b, N);

    Matrix<Real_wp> powC(c.GetM(), r);
    powC.Fill(Real_wp(1));
    for (int i = 0; i < r-1; i++)
      for (int j = 0; j < c.GetM(); j++)
	powC(j, i+1) = powC(j, i) * c(j);
    
    f.Reallocate(this->m_);
    f.Zero();
    
    for (int i = 0; i <= s; i++)
      f(i) = Numer(i) - N(i);

    int nb = s+1;
    if (r >= 3)
      {
	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 2);

	f(nb) -= Real_wp(1)/3;
	nb++;
      }

    if (r >= 4)
      {
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*A(i, j)*powC(j, 2);

	f(nb) -= Real_wp(1)/12;
	nb++;

	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 3);

	f(nb) -= Real_wp(1)/4;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*c(i)*A(i, j)*c(j);
	
	f(nb) -= Real_wp(1)/8;
	nb++;
      }
    
    if (r >= 5)
      {
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      f(nb) += b(i)*A(i, j)*A(j, k)*powC(k, 2);

	f(nb) -= Real_wp(1)/60;
	nb++;
	
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*A(i, j)*powC(j, 3);

	f(nb) -= Real_wp(1)/20;
	nb++;

	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 4);

	f(nb) -= Real_wp(1)/5;
	nb++;
      }

    if (r >= 6)
      {
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
		f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*powC(l, 2);
	
	f(nb) -= Real_wp(1)/360;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      f(nb) += b(i)*A(i, j)*A(j, k)*powC(k, 3);

	f(nb) -= Real_wp(1)/120;
	nb++;
	
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*A(i, j)*powC(j, 4);

	f(nb) -= Real_wp(1)/30;
	nb++;

	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 5);

	f(nb) -= Real_wp(1)/6;
	nb++;
      }

    if (r >= 7)
      {
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
                for (int m = 0; m <= l; m++)
                  f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*A(l, m)*powC(m, 2);
	
	f(nb) -= Real_wp(1)/2520;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
		f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*powC(l, 3);
	
	f(nb) -= Real_wp(1)/840;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      f(nb) += b(i)*A(i, j)*A(j, k)*powC(k, 4);

	f(nb) -= Real_wp(1)/210;
	nb++;
	
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*A(i, j)*powC(j, 5);

	f(nb) -= Real_wp(1)/42;
	nb++;

	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 6);

	f(nb) -= Real_wp(1)/7;
	nb++;
      }

    if (r >= 8)
      {
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
                for (int m = 0; m <= l; m++)
                  for (int n = 0; n <= m; n++)
                    f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*A(l, m)*A(m, n)*powC(n, 2);
	
	f(nb) -= Real_wp(1)/20160;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
                for (int m = 0; m <= l; m++)
                  f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*A(l, m)*powC(m, 3);
	
	f(nb) -= Real_wp(1)/6720;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      for (int l = 0; l <= k; l++)
		f(nb) += b(i)*A(i, j)*A(j, k)*A(k, l)*powC(l, 4);
	
	f(nb) -= Real_wp(1)/1680;
	nb++;

	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    for (int k = 0; k <= j; k++)
	      f(nb) += b(i)*A(i, j)*A(j, k)*powC(k, 5);

	f(nb) -= Real_wp(1)/336;
	nb++;
	
	for (int i = 0; i < s; i++)
	  for (int j = 0; j <= i; j++)
	    f(nb) += b(i)*A(i, j)*powC(j, 6);

	f(nb) -= Real_wp(1)/56;
	nb++;

	for (int i = 0; i < s; i++)
	  f(nb) += b(i)*powC(i, 7);

	f(nb) -= Real_wp(1)/8;
	nb++;
      }
  }

  void SdirkNonLinearEquation
  ::EvaluateJacobian(const VectReal_wp& x, VectReal_wp& f, Matrix<Real_wp, General, ColMajor>& df)
  {
    Real_wp h = pow(epsilon_machine, Real_wp(1)/3);
    EvaluateF(x, f);
    
    df.Reallocate(this->m_, this->n_);
    for (int j = 0; j < this->n_; j++)
      {
	VectReal_wp xm(x), fp, fm;
	xm(j) += h;

	this->EvaluateF(xm, fp);

	xm(j) -= Real_wp(2)*h;
	this->EvaluateF(xm, fm);

	for (int i = 0; i < this->m_; i++)
	  df(i, j) = (fp(i) - fm(i)) / (Real_wp(2)*h);
      }
  }
  
  
  /******************************
   * LinearSdirkScheme_Iterator *
   ******************************/


  //! default constructor
  template<class T>
  LinearSdirkScheme_Iterator<T>::LinearSdirkScheme_Iterator()
  {
    SetOrder(2, 0);
    stable_algo = false;
    sdirk_algo = false;
    lobatto_points = true;
  }


  //! constructor with the order
  template<class T>
  LinearSdirkScheme_Iterator<T>::LinearSdirkScheme_Iterator(int r, int extraS, int stable_)
  {
    alternative_source = true;
    stable_algo = false;
    sdirk_algo = false;
    SetOrder(r, extraS, stable_);
  }

  
  // case where no additional stages are added
  // We assume we know optimal gamma with DenomR = (1-gamma)^{s-1}
  // order is the scheme order
  template<class T>
  void LinearSdirkScheme_Iterator<T>::FindSdirkMinimal(int order, Real_wp gamma)
  {
    UnivariatePolynomial<Real_wp> Numer, Denom;     
    
    // on genere les coefs C_n^k directement
    //IVect row;
    //GetCnk(order, row);
    
    Real_wp invFacto(1);
    Numer.SetOrder(order);
    Numer(0) = Real_wp(1);
    for (int i = 0; i <= order-1; i++)
      {      
	// NumerR : developpement limite de l'exponentielle
	Numer(i+1) = abs(invFacto);
	invFacto *= Real_wp(1)/(i+2);
      }
    
    // vrai denominateur avec gamma
    //Denom = 1 - gamma z
    Denom.SetOrder(1);
    Denom(0) = Real_wp(1);
    Denom(1) = -gamma;
    // DenomR = (1-gamma)^s
    DenomR = Pow(Denom,order-1);
    
    // vrai numerateur en enlevant les zeros
    // dus aux constantes alpha, beta, etc
    this->NumerR = Numer*DenomR;
    this->NumerR.ResizeOrder(order-1);
  }

  
  // on calcule la fonction de stabilite R = Numer/Denom
  // en fonction de s (la methode est d'ordre s+1) et des parametres
  // param(0) = gamma
  // param(1:) = les constantes alpha, beta, etc
  template<class T>
  void LinearSdirkScheme_Iterator<T>::GetPolynomial_Sdirk(const VectReal_wp& param, int s,
                                                          UnivariatePolynomial<Real_wp>& Numer,
                                                          UnivariatePolynomial<Real_wp>& Denom,
                                                          UnivariatePolynomial<Real_wp>& Pol)
  {
    // premiere parametre : gamma
    this->gamma = param(0);
    int p = param.GetM() + s;
  
    // Pol = 1 - gamma z
    Pol.SetOrder(1); 
    Pol(0) = Real_wp(1); 
    Pol(1) = -this->gamma;
    
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
    //Numer.WriteText("NumerS.dat");
  }


  template<class T>
  void LinearSdirkScheme_Iterator<T>
  ::AddParamAlternativeSource(int i, UnivariatePolynomial<Real_wp>& Q)
  {
    int s = this->order-1;
    int p = NumerR.GetOrder();
    VectReal_wp param(p);
    param.Zero();
    switch (s)
      {
      case 7:
        if (p == 10)
          {
            switch(i)
              {
              case 1:
                param(0) = 5.44372738272868e-06;
                param(1) = -6.82082252911936e-06;
                param(2) = 3.91642147148413e-06;
                break;
              case 2:
                param(0) = 9.68812003968254e-08;
                param(1) = -1.18780947914244e-07;
                param(2) = 7.06910034807584e-08;
                break;
              case 3:
                param(0) = 7.1264251608397e-09;
                param(1) = -5.76968672352194e-10;
                param(2) = 5.64717829915781e-10;
                break;                  
              }
          }
        break;
      }
    
    for (int i = s; i < p; i++)
      Q(i+1) += param(i-s);
  }
  
  
  //! Compute ci and mat_wi use to evaluate source function for Linear Sdirk Scheme
  template<class T>
  void LinearSdirkScheme_Iterator<T>::ComputeCoefForLinearSdirkSourceFunction()
  {
    // numerator is factorized for a "stable" algorithm
    // ie we list the real and complex roots of the numerators
    VectReal_wp R, Ri;
    SolvePolynomialEquation(NumerR, R, Ri);
    
    real_roots.Clear(); complex_roots.Clear();
    int n = 0;
    while (n < R.GetM())
      {
        if (Ri(n) == Real_wp(0))
          {
            // real root
            real_roots.PushBack(R(n));
            n++;
          }
        else
          {
            // complex conjugate roots
            complex_roots.PushBack(Complex_wp(R(n), Ri(n)));
            n += 2;
          }
      }
    
    Sort(real_roots);
    Sort(complex_roots);
    
    // number of terms for the right hand side phi (each term involves 
    // a different power of (dt A)^{r-1} )
    int nb_terms = this->order-1;
    if (NumerR.GetOrder() != nb_terms)
      nb_terms = NumerR.GetOrder();
    
    Matrix<Real_wp> A;
    // if you want to see high-order error terms
    int nb_terms_add = nb_terms+1;
    A.Reallocate(nb_terms_add, nb_terms_add);
    A.Zero();
    
    // calcul de 1/k! et 1 / 2^k k!
    VectReal_wp invFacto(2*nb_terms_add+1), CoefC(2*nb_terms_add+1);
    invFacto(0) = Real_wp(1); CoefC(0) = Real_wp(1);
    for (int k = 1; k <= 2*nb_terms_add; k++)
      {
        invFacto(k) = invFacto(k-1) / Real_wp(k);
        CoefC(k) = CoefC(k-1) / Real_wp(2*k);
      }
    
    // terms with X^{0} (to see the error constants for the homogeneous case
    if (!alternative_source)
      {
        VectReal_wp err(nb_terms_add);
        err.Zero();
        for (int r = 0; r < nb_terms_add; r++)
          {
            Real_wp sum(0);
            for (int i = 0; i <= min(r, NumerR.GetOrder()); i++)
              {
                Real_wp signK(1);
                if ((r-i)%2 == 1)
                  signK = -signK;
                
                sum += (this->DenomR(i) - signK*this->NumerR(i))*CoefC(r-i);
              }
            
            err(r) = sum;
          }
        
        //cout << "Errors for the homogeneous case" << endl << err << endl;
      }
    
    // computation of coefficients alpha_j^r for the right hand side phi
    // phi = dt \sum_{r=1}^\infty (dt A)^{r-1} \sum_{j=1}^\infty \alpha_j^r dt^{j-1} F^{j-1}
    // where \alpha_j^r = \sum_{i=0}^{min(r-1, m)} (D_i - (-1)^k N_i) / 2^k k!
    // where k = r+j-i-1
    //cout << endl << "Coefficients alpha for inhomogeneous case" << endl;
    if (!alternative_source)
      for (int r = 1; r <= nb_terms_add; r++)
        {
          VectReal_wp err(nb_terms_add);
          err.Zero();
          for (int j = 1; j <= nb_terms_add; j++)
            {
              for (int i = 0; i <= min(r-1, this->NumerR.GetOrder()); i++)
                {
                  int k = r + j-i - 1;
                  Real_wp signK(1);
                  if (k%2 == 1)
                    signK = -signK;
                  
                  err(j-1) += (this->DenomR(i) - signK*this->NumerR(i))*CoefC(k);
                }
              
              A(r-1, j-1) = err(j-1);
            }
          
          //cout << "Terms with (dt A)^" << r-1 << endl << err << endl;
        }
    
    Matrix<Real_wp> VDM; VectReal_wp rhs(nb_terms+1);
    Vector<int> pivot;
    VectReal_wp weights;

    if (alternative_source)
      {
        nb_terms = NumerR.GetOrder() + 1;
        if (!lobatto_points)
          nb_terms--;
        
        int rquad = (this->order+1)/2;
        //DISP(this->order); DISP(rquad);
        
        // dans ce cas, on prefere Gauss-Lobatto
        // car les deux extremites sont immediates a traiter
        if (lobatto_points)
          ComputeGaussLobatto(this->ci, weights, rquad);
        else
          ComputeGaussLegendre(this->ci, weights, rquad-1);
        
        //DISP(this->ci); DISP(weights);

        A.Reallocate(nb_terms, this->ci.GetM());
        A.Zero();

        int p = NumerR.GetOrder(); //DISP(this->order); DISP(nb_terms); DISP(p);
        invFacto.Fill(1.0);
        UnivariatePolynomial<Real_wp> Pol, Q, P;
        for (int i = 0; i < this->ci.GetM(); i++)
          {
            // Pol = 1 - gamma z
            Pol.SetOrder(1); 
            Pol(0) = Real_wp(1); 
            Pol(1) = -this->gamma;
            
            if (abs(this->ci(i)) <= 1e-12)
              {
                Q = this->NumerR;
              }
            else if (abs(this->ci(i)-1.0) <= 1e-12)
              {
                Q = Pow(Pol, p);
              }
            else
              {
                // produit du DL de l'exponentielle par (1+gamma z)^(m+l)
                P.SetOrder(this->order-1);
                P(0) = 1.0;
                Real_wp xsi = (1.0-this->ci(i));
                for (int j = 1; j <= this->order-1; j++)
                  P(j) = P(j-1) / Real_wp(j)*xsi;

                Q = P*Pow(Pol, p);
                this->AddParamAlternativeSource(i, Q);

                // on tronque 
                Q.ResizeOrder(nb_terms-1);
              }
            
            //DISP(this->ci(i)); DISP(Q);
            
            // on forme la matrice A
            for (int r = 1; r <= nb_terms; r++)
              A(r-1, i) = weights(i)*Q(r-1);
          }

        //DISP(this->ci); DISP(weights); DISP(this->gamma);
        //A.WriteText("A.dat");
      }
    else
      {
        // calcul des ci : on prend les points de Gauss-Legendre
        ComputeGaussLegendre(this->ci, weights, nb_terms);
        //ComputeGaussLobatto(this->ci, weights, nb_terms);

        // Vandermonde Matrix for ci
        // les factorielles sont mises dans le second membre            
        VDM.Reallocate(this->ci.GetM(), this->ci.GetM());  
        for (int i = 0; i < this->ci.GetM(); i++)
          for (int j = 0; j < this->ci.GetM(); j++)
            VDM(i, j) = pow(this->ci(j) - Real_wp(0.5), i);

        GetLU(VDM, pivot);        
      }
    
    int orderN = nb_terms;
    if (stable_algo)
      {        
        orderN = NumerR.GetOrder()+1;
        if (!lobatto_points)
          orderN--;
        
        // Pol = 1 - gamma z
        UnivariatePolynomial<Real_wp> Pol;
        Pol.SetOrder(1); 
        Pol(0) = Real_wp(1); 
        Pol(1) = -this->gamma;
        
        // series (1 - z / lambda_n) (1 - z / lambda_{n-1}) .. (1 - z/lambda_2)
        UnivariatePolynomial<Real_wp> Q, P, z;
        z.SetOrder(1); z(0) = Real_wp(0); z(1) = Real_wp(1);
        Q.SetOrder(0); Q(0) = Real_wp(1);
        Vector<UnivariatePolynomial<Real_wp> > PolNumer(orderN+1);
        int nb = 0;
        PolNumer(nb) = Q; nb++;
        for (int i = real_roots.GetM()-1; i >= 0; i--)
          {
            P.SetOrder(1); P(0) = Real_wp(1); P(1) = -Real_wp(1)/real_roots(i);
            Q *= P; 
            PolNumer(nb) = Q; nb++;
          }
        
        for (int i = complex_roots.GetM()-1; i >= 0; i--)
          {
            PolNumer(nb) = z*Q; nb++;
            P.SetOrder(2);
            P(0) = Real_wp(1);
            P(1) = -Real_wp(2)*realpart(Real_wp(1) / complex_roots(i));
            P(2) = Real_wp(1)/absSquare(complex_roots(i));
            Q *= P;
            PolNumer(nb) = Q; nb++;
          }
        
        // polynomials needed to expand the source are enumerated
        Vector<UnivariatePolynomial<Real_wp> > pol_source(orderN);
        UnivariatePolynomial<Real_wp> monoGamma;
        monoGamma.SetOrder(0); monoGamma(0) = Real_wp(1);
        int num = 0;
        for (int i = 0; i < complex_roots.GetM(); i++)
          {
            pol_source(num) = PolNumer(nb-2-num)*monoGamma;
            num++;
            pol_source(num) = PolNumer(nb-2-num)*monoGamma;
            num++;
            
            monoGamma *= Pol*Pol;
          }
       
        for (int i = 0; i < real_roots.GetM(); i++)
          {
            pol_source(num) = PolNumer(nb-2-num)*monoGamma;
            num++;

            monoGamma *= Pol;
          }

        if (lobatto_points)
          pol_source(num) = monoGamma;
        
        // inverse of coefficients, to express x^i as a combination of pol_source
        Matrix<Real_wp> coefBase;
        coefBase.Reallocate(orderN, orderN);
        coefBase.Zero();
        for (int j = 0; j < orderN; j++)
          for (int i = 0; i < pol_source(j).GetM(); i++)
            coefBase(i, j) = pol_source(j)(i);
        
        GetInverse(coefBase);
        polA.Reallocate(orderN, A.GetN());
        polA.Zero();
        for (int powA = 0; powA < nb_terms; powA++)
          {
            rhs.Zero();	
            
            for (int j = 1; j <= A.GetN(); j++)
              rhs(j-1) = A(powA, j-1);
            
            for (int j = 0; j < A.GetN(); j++)
              for (int k = 0; k < orderN; k++)
                polA(k, j) += coefBase(k, powA)*rhs(j);
          }
      }
    
    this->mat_wi.Reallocate(orderN, this->ci.GetM());

    // loop over rows of mat_wi
    for (int powA = 0; powA < orderN; powA++)
      {
	rhs.Zero();	
        
        if (stable_algo)
          for (int j = 1; j <= polA.GetN(); j++)
            rhs(j-1) = polA(powA, j-1) / invFacto(j-1);
        else
          for (int j = 1; j <= A.GetN(); j++)
            rhs(j-1) = A(powA, j-1) / invFacto(j-1);
        
        if (alternative_source)
          SetRow(rhs, powA, this->mat_wi);
        else
          {
            SolveLU(VDM, pivot, rhs);
            SetRow(rhs, powA, this->mat_wi);
          }
        
        // on evalue le reste
        /* for (int  i = 0; i < nb_terms_add; i++)
          {
            Real_wp coef(0);
            for (int j = 0; j < this->ci.GetM(); j++)
              coef += pow(this->ci(j) - Real_wp(0.5), i)*invFacto(i)*rhs(j);
            
            DISP(i); DISP(coef-A(powA, i)); 
            } */
      }

    //DISP(this->real_roots); DISP(this->complex_roots);
    //DISP(this->ci); DISP(this->mat_wi);
    //this->mat_wi.WriteText("w.dat");
  }


  //! computation of coefficients for Linear SDIRK scheme of order s+1
  template<class T>
  void LinearSdirkScheme_Iterator<T>::SetOrder(int order, int extraS, int type_algo)
  {
    stable_algo = false;
    sdirk_algo = false;
    alternative_source = false;
    if (type_algo == STABLE_WEIGHTS)
      {
        stable_algo = true;
      }
    else if (type_algo == POLY_ALTERNATIVE)
      {
        alternative_source = true;
        
      }
    else if (type_algo == STABLE_ALTERNATIVE)
      {
        alternative_source = true;
        stable_algo = true;
      }

    if (type_algo == SDIRK)
      sdirk_algo = true;
    
    Real_wp one(1), zero(0);
    gamma = zero;
    this->order = order;
    UnivariatePolynomial<Real_wp> Pol;

    VectReal_wp param(extraS);
    param.Fill(one);
	  
    switch (order)
      {
      case 2 :
	{
	  switch(extraS)
	    {
	    case 0: { this->gamma = 0.5;}
	      break;
	    default :
	      {
		cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
		abort();
	      }
	    }
	}
	break;
      case 3 : 
	{
	  switch(extraS)
            {
	    case 0: { this->gamma = to_num<Real_wp>("0.78867513459481288225457439025098");}
	      break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }          
	}
	break;
      case 4 :
	{
	  switch(extraS)
            {
	    case 0: { this->gamma = to_num<Real_wp>("1.06857902130162880641883397596");}
	      break;
            case 1: { param(0) = 0.394337567297407;}
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 6 :
	{
	  switch(extraS)
            {
	    case 0 : { this->gamma = to_num<Real_wp>("0.47326839125829532445558852540261");}
	      break;
            case 1: { param(0) = 0.284064638011799;}
              break;
            case 2: {param(0) = 0.204071; param(1) = 1.9839430662e-4;} 
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 8:
	{
	  switch(extraS)
            {
            case 1: { param(0) = 0.217049743094304;}
              break;
	    case 2: {param(0) = 0.16689; param(1) = 2.9259251764e-6;} 
              break;
	      
            case 3: 
              {
                // valeurs quasi-optimales
                //param(0) = 0.136339 ; param(1) = 2.766997072e-06; param(2) = -3.464371455249e-06;
                
                // valeurs plus sures (dans la zone de stabilite)
                param(0) = 0.136339 ; param(1) = 2.767416226e-6; param(2) = -3.464398093e-6;
              } 
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 10:
	{
	  switch(extraS)
            {
            case 2: {param(0) = 0.141940; param(1) = 2.2982637210e-8;} 
              break;
	    case 3 :
              {
                // valeurs quasi optimales
                //param(0) = 0.151706; param(1) = 2.459142766e-08; param(2) = -4.3141010528e-08;

                // valeurs plus sures (dans la zone de stabilite)
                param(0) = 0.151706; param(1) = 2.459114959e-8; param(2) = -4.3140917546e-8;
              }
	      break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 12:
	{
	  switch(extraS)
            {
	    case 3 : 
              {
                // parametres quasi-optimaux
                //param(0) = 0.132571; param(1) = 1.64451085e-10; param(2) = -2.89888832092e-10;
                
                // valeurs plus sures (dans la zone de stabilite)
                param(0) = 0.132572; param(1) = 1.644515143e-10; param(2) = -2.89891484131e-10;
              }
	      break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_SDIRK number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      default :
	{
	  cout << "Order " << order << " Not implemented " <<endl;
	  abort();
	}
      }

    // Find stability function R =NumerR/DenomR  of linear sdirk of order s+1 from param
    if(extraS==0)
      {
	this->FindSdirkMinimal(order,this->gamma);
      }
    else
      {
	this->GetPolynomial_Sdirk(param, this->order-1, NumerR, DenomR, Pol);
      }
    
    if (sdirk_algo)
      {
        A_rk.Reallocate(order + extraS-1, order + extraS-1); 
        A_rk.Zero();
        for (int i = 0; i < A_rk.GetM(); i++)
          A_rk(i, i) = gamma;
        
        b_rk.Reallocate(order + extraS-1); b_rk.Zero();
        c_rk.Reallocate(order + extraS-1); c_rk.Zero();
        switch (order)
          {
          case 2 :
            switch(extraS)
              {
              case 0:
                b_rk(0) = Real_wp(1);
                c_rk(0) = Real_wp(0.5);
                break;
              }
            break;
          case 3 : 
            switch(extraS)
              {
	    case 0:
              A_rk(1, 0) = Real_wp(1) - 2*gamma;
              b_rk(0) = one/2; b_rk(1) = one/2;
              c_rk(0) = gamma; c_rk(1) = one-gamma;
	      break;
              }          
            break;
          case 4 :
            switch(extraS)
              {
              case 0:
                {
                  Real_wp delta = one/(6*square(2*gamma-one));
                  A_rk(1, 0) = one/2-gamma; A_rk(2, 0) = 2*gamma; A_rk(2, 1) = one-4*gamma;
                  b_rk(0) = delta; b_rk(1) = one-2*delta; b_rk(2) = delta;
                  c_rk(0) = gamma; c_rk(1) = one/2; c_rk(2) = one-gamma;
                }
                break;
              case 1:
                A_rk(1, 0) = to_num<Real_wp>("-0.2886751345948147762484678898");
                A_rk(2, 0) = to_num<Real_wp>("-0.1201707048685739554617350983");
                A_rk(2, 1) = to_num<Real_wp>("0.4927599105238210712717179457");
                A_rk(3, 0) = to_num<Real_wp>("0.5862179118900888194222121688");
                A_rk(3, 1) = to_num<Real_wp>("0.210931547279092506040574263");
                A_rk(3, 2) = to_num<Real_wp>("-0.3487361007262912198383726178");
                b_rk(0) = to_num<Real_wp>("0.3759622071506074656024816985");
                b_rk(1) = to_num<Real_wp>("0.2298047916130732484017089428");
                b_rk(2) = to_num<Real_wp>("0.06301398783510507853815187951");
                b_rk(3) = to_num<Real_wp>("0.3312190134012142074576574791");
                c_rk(0) = to_num<Real_wp>("0.3943375672974069878762293229");
                c_rk(1) = to_num<Real_wp>("0.1056624327025922116277614331");
                c_rk(2) = to_num<Real_wp>("0.7669267729526541036862121703");
                c_rk(3) = to_num<Real_wp>("0.8427509257402970935006431368");
                break;
              }
            break;
          case 6 :
            switch(extraS)
              {
              case 0 :
                break;
              case 1:
                break;
              case 2:
                A_rk(1, 0) = to_num<Real_wp>("-0.166045901047610158434879897");
                A_rk(2, 0) = to_num<Real_wp>("0.4096575353468192508516437702");
                A_rk(2, 1) = to_num<Real_wp>("-0.280017680327768296686800165");
                A_rk(3, 0) = to_num<Real_wp>("0.2957916453358002917992326443");
                A_rk(3, 1) = to_num<Real_wp>("0.1713748967266998304285877486");
                A_rk(3, 2) = to_num<Real_wp>("-0.3378514714085876401203255972");
                A_rk(4, 0) = to_num<Real_wp>("0.4364137446395151231562725723");
                A_rk(4, 1) = to_num<Real_wp>("0.01281537228266837107077785072");
                A_rk(4, 2) = to_num<Real_wp>("0.1849418647187038009172024896");
                A_rk(4, 3) = to_num<Real_wp>("-0.3234234927228179618725047744");
                A_rk(5, 0) = to_num<Real_wp>("-0.03698714673820067516178143592");
                A_rk(5, 1) = to_num<Real_wp>("0.2045321965187032072186450404");
                A_rk(5, 2) = to_num<Real_wp>("0.1190866732359654531201458836");
                A_rk(5, 3) = to_num<Real_wp>("-0.03512223693158210584554581939");
                A_rk(5, 4) = to_num<Real_wp>("0.5687727338603634533809266257");
                A_rk(6, 0) = to_num<Real_wp>("-0.07537345700323794294469495066");
                A_rk(6, 1) = to_num<Real_wp>("0.1864431937998533267508718057");
                A_rk(6, 2) = to_num<Real_wp>("0.1814994252778522182386200051");
                A_rk(6, 3) = to_num<Real_wp>("0.2486163764378303204173848473");
                A_rk(6, 4) = to_num<Real_wp>("0.05648477197606868917634557114");
                A_rk(6, 5) = to_num<Real_wp>("-0.03062566373857441623608633552");
                b_rk(0) = to_num<Real_wp>("0.0833948370374948088659322941");
                b_rk(1) = to_num<Real_wp>("0.1193469581907871075164042306");
                b_rk(2) = to_num<Real_wp>("0.1059284090283850903841566058");
                b_rk(3) = to_num<Real_wp>("0.1928545304954600586890754047");
                b_rk(4) = to_num<Real_wp>("0.07637102566629057372233113448");
                b_rk(5) = to_num<Real_wp>("0.05576375704868218367317241305");
                b_rk(6) = to_num<Real_wp>("0.3663404825329001771489279173");
                c_rk(0) = to_num<Real_wp>("0.2039400000000000101607611214");
                c_rk(1) = to_num<Real_wp>("0.03789409895238985172588122432");
                c_rk(2) = to_num<Real_wp>("0.3335798550190509643256047266");
                c_rk(3) = to_num<Real_wp>("0.3332550706539124922682559171");
                c_rk(4) = to_num<Real_wp>("0.5146874889180693434325092597");
                c_rk(5) = to_num<Real_wp>("1.024222219945249342873151416");
                c_rk(6) = to_num<Real_wp>("0.7709846467497922055632020644");
                break;
              }
            break;
          case 8:
            switch(extraS)
              {
              case 1:
                A_rk(1, 0) = to_num<Real_wp>("-0.1796128810692800908174243662");
                A_rk(2, 0) = to_num<Real_wp>("-0.8328495931145723339945206342");
                A_rk(2, 1) = to_num<Real_wp>("1.031105136398725003032049804");
                A_rk(3, 0) = to_num<Real_wp>("-0.1325672897658342356835593468");
                A_rk(3, 1) = to_num<Real_wp>("-0.1483967394890516001442034053");
                A_rk(3, 2) = to_num<Real_wp>("0.9955028933330461834845539904");
                A_rk(4, 0) = to_num<Real_wp>("0.2620862983501705389158956236");
                A_rk(4, 1) = to_num<Real_wp>("-0.2775708792908560351564389476");
                A_rk(4, 2) = to_num<Real_wp>("0.6731341054500475074875676381");
                A_rk(4, 3) = to_num<Real_wp>("-0.01962368612446165917572132688");
                A_rk(5, 0) = to_num<Real_wp>("-0.1835083484276407591089726347");
                A_rk(5, 1) = to_num<Real_wp>("0.4569793640889019458315131778");
                A_rk(5, 2) = to_num<Real_wp>("0.04093539207370280483632473435");
                A_rk(5, 3) = to_num<Real_wp>("0.09174967540211776267529535633");
                A_rk(5, 4) = to_num<Real_wp>("-0.1324347156284396595801915333");
                A_rk(6, 0) = to_num<Real_wp>("-0.4839461037980150520728672588");
                A_rk(6, 1) = to_num<Real_wp>("0.2527173012727905006575453537");
                A_rk(6, 2) = to_num<Real_wp>("0.563161587333493452147665423");
                A_rk(6, 3) = to_num<Real_wp>("0.08056999208694265543726937033");
                A_rk(6, 4) = to_num<Real_wp>("-0.1484228223961586720450496233");
                A_rk(6, 5) = to_num<Real_wp>("0.4710568793660774232700754679");
                A_rk(7, 0) = to_num<Real_wp>("0.8427709054434147859959936593");
                A_rk(7, 1) = to_num<Real_wp>("-0.3102674519954351307067487574");
                A_rk(7, 2) = to_num<Real_wp>("-0.03417573995885483394078008582");
                A_rk(7, 3) = to_num<Real_wp>("-0.078687181753846607195136314");
                A_rk(7, 4) = to_num<Real_wp>("0.2286406237181812775167520585");
                A_rk(7, 5) = to_num<Real_wp>("0.06679423159911813965588679436");
                A_rk(7, 6) = to_num<Real_wp>("-0.1261849640312153152002533346");
                b_rk(0) = to_num<Real_wp>("0.2583465021432374273718210536");
                b_rk(1) = to_num<Real_wp>("0.1004205266062366276302700841");
                b_rk(2) = to_num<Real_wp>("-0.07984841397746895478995993257");
                b_rk(3) = to_num<Real_wp>("-0.09672915297643563735547697433");
                b_rk(4) = to_num<Real_wp>("-0.06380695096217242503149697703");
                b_rk(5) = to_num<Real_wp>("0.3811630470331293756303135862");
                b_rk(6) = to_num<Real_wp>("0.188718419871379788301765723");
                b_rk(7) = to_num<Real_wp>("0.311736022262093798242763437");
                c_rk(0) = to_num<Real_wp>("0.2170497430943039984807541032");
                c_rk(1) = to_num<Real_wp>("0.03743686202502390766332973696");
                c_rk(2) = to_num<Real_wp>("0.4153052863784566675182832729");
                c_rk(3) = to_num<Real_wp>("0.9315886071724643461375453415");
                c_rk(4) = to_num<Real_wp>("0.8550755814792043505520570904");
                c_rk(5) = to_num<Real_wp>("0.4907711106029460931347232037");
                c_rk(6) = to_num<Real_wp>("0.952186576959434305875392836");
                c_rk(7) = to_num<Real_wp>("0.8059401661156663146064681235");
                break;
              case 2:
                break;
              case 3:
                A_rk(1, 0) = to_num<Real_wp>("-0.1390654076366021986567843099");
                A_rk(2, 0) = to_num<Real_wp>("0.216007406154890626442183881");
                A_rk(2, 1) = to_num<Real_wp>("-0.07147684702610457924332899644");
                A_rk(3, 0) = to_num<Real_wp>("0.2253277199711876849965597452");
                A_rk(3, 1) = to_num<Real_wp>("0.0221600616190271465097082392");
                A_rk(3, 2) = to_num<Real_wp>("-0.3691182179503144062048907546");
                A_rk(4, 0) = to_num<Real_wp>("-0.02993396704315138394939769341");
                A_rk(4, 1) = to_num<Real_wp>("0.1728707660748781924891760372");
                A_rk(4, 2) = to_num<Real_wp>("0.02300271188763571431093385542");
                A_rk(4, 3) = to_num<Real_wp>("0.4037330591900536616571121884");
                A_rk(5, 0) = to_num<Real_wp>("0.2602360368085311804993313111");
                A_rk(5, 1) = to_num<Real_wp>("0.0980269586088953542239373001");
                A_rk(5, 2) = to_num<Real_wp>("0.06344575235946783089633034242");
                A_rk(5, 3) = to_num<Real_wp>("0.09199080775625423132563212789");
                A_rk(5, 4) = to_num<Real_wp>("0.146171012679266390072105416");
                A_rk(6, 0) = to_num<Real_wp>("0.1306870995235180256698315927");
                A_rk(6, 1) = to_num<Real_wp>("-0.01972135957966156454320418012");
                A_rk(6, 2) = to_num<Real_wp>("0.2844409312478850603947033786");
                A_rk(6, 3) = to_num<Real_wp>("-0.03317110012091200929457248685");
                A_rk(6, 4) = to_num<Real_wp>("0.1175544726498403993802063099");
                A_rk(6, 5) = to_num<Real_wp>("-0.1454930235325671236520496958");
                A_rk(7, 0) = to_num<Real_wp>("0.09465203502138632957306334712");
                A_rk(7, 1) = to_num<Real_wp>("0.1107181809653364435253052189");
                A_rk(7, 2) = to_num<Real_wp>("0.04926116827596873208383366745");
                A_rk(7, 3) = to_num<Real_wp>("-0.01201552137652603106622827866");
                A_rk(7, 4) = to_num<Real_wp>("0.09857632081838423693636408196");
                A_rk(7, 5) = to_num<Real_wp>("0.2100177535462915089277625953");
                A_rk(7, 6) = to_num<Real_wp>("0.2125341246753185229288819955");
                A_rk(8, 0) = to_num<Real_wp>("0.07822953770733603919809567538");
                A_rk(8, 1) = to_num<Real_wp>("0.005056150390280394630374148505");
                A_rk(8, 2) = to_num<Real_wp>("0.05899985920278310719856665836");
                A_rk(8, 3) = to_num<Real_wp>("0.1925818100152683195227059192");
                A_rk(8, 4) = to_num<Real_wp>("0.2009134406705905707755680409");
                A_rk(8, 5) = to_num<Real_wp>("0.1141455783455047494054309695");
                A_rk(8, 6) = to_num<Real_wp>("0.2718102944989661446723999484");
                A_rk(8, 7) = to_num<Real_wp>("-0.1123800346025666231231068593");
                A_rk(9, 0) = to_num<Real_wp>("0.07918172298026633806349301806");
                A_rk(9, 1) = to_num<Real_wp>("0.1415399502784503626755313433");
                A_rk(9, 2) = to_num<Real_wp>("0.110492226083356750811420412");
                A_rk(9, 3) = to_num<Real_wp>("0.1282655026535273824160208811");
                A_rk(9, 4) = to_num<Real_wp>("0.222254970395505940513691109");
                A_rk(9, 5) = to_num<Real_wp>("0.181496999258891444633042827");
                A_rk(9, 6) = to_num<Real_wp>("-0.02452805931085999092150967159");
                A_rk(9, 7) = to_num<Real_wp>("-0.1737016053283928162832276522");
                A_rk(9, 8) = to_num<Real_wp>("0.09069419619459767472553845391");
                b_rk(0) = to_num<Real_wp>("0.1902755904701471365410244351");
                b_rk(1) = to_num<Real_wp>("0.02338724337456116512132980698");
                b_rk(2) = to_num<Real_wp>("0.1011035970316270618692108953");
                b_rk(3) = to_num<Real_wp>("0.02127517911851249796636462503");
                b_rk(4) = to_num<Real_wp>("0.08784214730779824345175950664");
                b_rk(5) = to_num<Real_wp>("0.2315789618059805477716587305");
                b_rk(6) = to_num<Real_wp>("0.278936452252409535352519933");
                b_rk(7) = to_num<Real_wp>("0.2142249432364535352148213689");
                b_rk(8) = to_num<Real_wp>("0.14811733754936201011422021");
                b_rk(9) = to_num<Real_wp>("-0.2967414521468517334029095114");
                c_rk(0) = to_num<Real_wp>("0.1358800000000000007815970093");
                c_rk(1) = to_num<Real_wp>("-0.003185407636602197875187300543");
                c_rk(2) = to_num<Real_wp>("0.2804105591287860479804518939");
                c_rk(3) = to_num<Real_wp>("0.01424956363990042608297423915");
                c_rk(4) = to_num<Real_wp>("0.7055525701094161852894213969");
                c_rk(5) = to_num<Real_wp>("0.7957505682124149877989335069");
                c_rk(6) = to_num<Real_wp>("0.4701770201881027887365119277");
                c_rk(7) = to_num<Real_wp>("0.8996240619261597436905796369");
                c_rk(8) = to_num<Real_wp>("0.9452366362281627030616315103");
                c_rk(9) = to_num<Real_wp>("0.89157590320534308741559773");
                break;
              }
            break;
          case 10:
            switch(extraS)
              {
              case 2:
              break;
	    case 3 :
	      break;
              }
            break;
          case 12:
            switch(extraS)
              {
              case 3 :
	      break;
              }
            break;
          }
            
        if (epsilon_machine < 1e-30 )
          {
            SdirkNonLinearEquation eq;
            
            eq.Init(this->GetNumeratorStabilityFunction(),
                    this->GetDenominatorStabilityFunction(),
                    this->GetGammaCoefficient(), this->order,
                    A_rk, b_rk, c_rk);

            this->WriteCoefficients("coef"+to_str(NumerR.GetM())+".txt");
          }
        
        return;
      }

    // Compute data needed (ci and wi) to evaluate the source function
    if (stable_algo)
      {
        this->ComputeCoefForLinearSdirkSourceFunction();
        //this->WriteCoefficients("coef.txt");
      }
    else
      this->ComputeCoefForLinearSdirkSourceFunction();

    //DISP(real_roots); DISP(complex_roots); DISP(ci); DISP(polA); DISP(mat_wi); 
  }  

  
  //! Initialisation of the scheme with u^0
  template<class T>
  void LinearSdirkScheme_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                        Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();   
    sys.SetDirichletCondition(t, 0, Y);
    Fn.Reallocate(N); Fn.Zero();
    EvalF.Reallocate(ci.GetM());
    if (sdirk_algo)
      EvalF.Reallocate(c_rk.GetM());
    
    for(int i = 0; i < EvalF.GetM(); i++)
      {
        EvalF(i).Reallocate(N);
        EvalF(i).Zero();
      }

    KhUn.Reallocate(N); KhUn.Fill(0);
    
    // we ask the user to factorize the matrix
    // M_h + dt gamma S_h + dt gamma K_h
    sys.FactorizeOperatorDhPlusGammaKh(Real_wp(1), dt*gamma, dt*gamma);
    
  }
  

  //! clearing arrays used by time scheme
  template<class T>
  void LinearSdirkScheme_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! clearing arrays used by time scheme
  template<class T>
  void LinearSdirkScheme_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Fn.Clear(); KhUn.Clear();
    EvalF.Clear();
    Yfinal = Y; Y.Clear();
  }
  
  
  //! coefficients are writeen on a file
  template<class T>
  void LinearSdirkScheme_Iterator<T>::WriteCoefficients(string file_name)
  {    
    ofstream file_out(file_name.data());
    file_out.precision(28);

    if (sdirk_algo)
      {
        for (int i = 0; i < A_rk.GetM(); i++)
          for (int j = 0; j < i; j++)
            file_out << "A_rk(" << i << ", " << j << ") = to_num<Real_wp>(\"" 
                     << A_rk(i, j) << "\");\n";
        
        for (int i = 0; i < b_rk.GetM(); i++)
          file_out << "b_rk(" << i << ") = to_num<Real_wp>(\"" 
                     << b_rk(i) << "\");\n";

        for (int i = 0; i < c_rk.GetM(); i++)
          file_out << "c_rk(" << i << ") = to_num<Real_wp>(\"" 
                     << c_rk(i) << "\");\n";
          
        return;
      }

    file_out << "{\n string data(\"" << real_roots << "\");\n"
             << "istringstream flux(data); \n real_roots.ReadText(flux); \n } \n";

    file_out << "{\n string data(\"" << complex_roots << "\");\n"
             << "istringstream flux(data); \n complex_roots.ReadText(flux); \n } \n";

    file_out << "{\n string data(\"" << ci << "\");\n"
             << "istringstream flux(data); \n ci.ReadText(flux); \n } \n";

    file_out << "{\n string data(\"";
    for (int i = 0; i < mat_wi.GetM(); i++)
      {
        for (int j = 0; j < mat_wi.GetN(); j++)
          file_out << mat_wi(i, j) << ' ';
        
        file_out << "\\n";
      }

    file_out << "\");\n"
             << "istringstream flux(data); \n mat_wi.ReadText(flux); \n } \n";
    
    file_out << "{\n string data(\"";
    for (int i = 0; i < polA.GetM(); i++)
      {
        for (int j = 0; j < polA.GetN(); j++)
          file_out << polA(i, j) << ' ';
        
        file_out << "\\n";
      }

    file_out << "\");\n"
             << "istringstream flux(data); \n polA.ReadText(flux); \n } \n";
  }
  

  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void LinearSdirkScheme_Iterator<T>::Advance(const Real_wp& t, int n_time,
                                              VirtualOdeSystem<T>& sys)
  {
    if (sdirk_algo)
      {
        // computation of all intermediary unknowns
        for (int i = 0; i < A_rk.GetM(); i++)
          {
            // current time
            Real_wp tcurrent = t + c_rk(i)*dt;
            
            // computation of Un_tmp = U^n + delta t \sum_{j=1}^{i-1} a_ij z_j 
            Copy(Y, Fn);
            for (int j = 0; j < i; j++)
              Add(dt*A_rk(i, j), EvalF(j), Fn);
        
            // computation of f(Un_tmp, tcurrent)
            sys.EvaluateFunction(tcurrent, Fn, KhUn, false);

            // then, we solve system z_i = f(Un_tmp + gamma dt z_i, t^n + c_i dt)
            sys.SolveOperatorDhPlusGammaKh(tcurrent, KhUn, EvalF(i));
          }
        
        // Once we have all the intermediary unknowns, we deduce
        // U^n+1 = U^n + dt \sum b_i z_i
        for (int i = 0; i < b_rk.GetM(); i++)
          Add(dt*b_rk(i), EvalF(i), Y);

        return;
      }

    if (stable_algo)
      {
        // on evalue d'abord la source
        bool presence_source = false;
        if (t <= sys.GetFinalTimeSource())
          {
            for (int i = 0; i < this->mat_wi.GetN(); i++)
              {
                Real_wp tcurrent = t + ci(i)*dt;
                // compute y = y + alpha f(t)	
                EvalF(i).Zero();
                sys.AddPrimitiveTimeSource(Real_wp(1), tcurrent, 0, EvalF(i));
              }
            
            presence_source = true;
          }
        
        // cas ou on divise au fur et a mesure
        int num_source = 0;
        for (int k = 0; k < complex_roots.GetM(); k++)
          {
            Real_wp b = -Real_wp(2)*dt*realpart(Real_wp(1) / complex_roots(k));
            Real_wp a = dt*dt/absSquare(complex_roots(k));
            sys.EvaluateFunction(t, Y, Fn, false, false);
            Fn *= a;
            KhUn.Zero();
            if (presence_source)
              {
                for (int i = 0; i < mat_wi.GetN(); i++)
                  KhUn += dt*dt*mat_wi(num_source, i)*EvalF(i);
              }
            
            Fn += KhUn;
            sys.SolveMassMatrix(Fn);
            KhUn = b*Y + Fn;
            sys.EvaluateFunction(t, KhUn, Fn, false, false);
            
            sys.ApplyMassMatrix(Real_wp(1), t, Y, Real_wp(1), Fn); 
            KhUn.Zero();
            if (presence_source)
              {
                for (int i = 0; i < mat_wi.GetN(); i++)
                  KhUn += dt*mat_wi(num_source+1, i)*EvalF(i);
              }

            Fn += KhUn;
            sys.SolveOperatorDhPlusGammaKh(t, Fn, KhUn);
            sys.ApplyMassMatrix(Real_wp(1), t, KhUn, Real_wp(0), Fn);
            sys.SolveOperatorDhPlusGammaKh(t, Fn, KhUn);
            Y = KhUn;
            num_source += 2;
          }
        
        for (int k = 0; k < real_roots.GetM(); k++)
          {
            Real_wp b = -dt / real_roots(k);
            sys.EvaluateFunction(t, Y, Fn, false, false);
            sys.ApplyMassMatrix(Real_wp(1), t, Y, b, Fn);
            KhUn.Zero();
            if (presence_source)
              {
                for (int i = 0; i < mat_wi.GetN(); i++)
                  KhUn += dt*mat_wi(num_source, i)*EvalF(i);
              }
            
            Fn += KhUn;
            sys.SolveOperatorDhPlusGammaKh(t, Fn, Y);            
            num_source++;
          }

        if (presence_source && lobatto_points)
          {
            KhUn.Zero();
            for (int i = 0; i < mat_wi.GetN(); i++)
              KhUn += dt*mat_wi(num_source, i)*EvalF(i);
            
            sys.SolveMassMatrix(KhUn);
            Y += KhUn;
          }
        
        return;
      }

    if (t <= sys.GetFinalTimeSource())
      for(int i = 0; i < this->mat_wi.GetN(); i++)
	{
	  Real_wp tcurrent = t + ci(i)*dt;
	  // compute y = y + alpha f(t)	
	  EvalF(i).Zero();
	  sys.AddPrimitiveTimeSource(Real_wp(1), tcurrent, 0, EvalF(i));
          sys.SolveMassMatrix(EvalF(i));
	}

    // We evaluate Fn = NumerR*U^n + dt*Fn
    int last_coef = NumerR.GetM()-1;
    
    Fn = this->NumerR(last_coef)*Y;
    if (lobatto_points)
      if((t <= sys.GetFinalTimeSource()) && (last_coef < this->mat_wi.GetM()))
        {
          KhUn.Zero();
          for(int i = 0; i < mat_wi.GetN(); i++)
            KhUn += dt*mat_wi(last_coef, i)*EvalF(i);
          
          Fn += KhUn;
        }
    
    for(int k = last_coef -1; k >= 0; k--)
      {
	// we multiply by K
	sys.EvaluateFunction(t, Fn, KhUn, true, false);
	Fn = dt*KhUn + NumerR(k)*Y;
	if((t <= sys.GetFinalTimeSource()) && (k < this->mat_wi.GetM()))
	  {
	    KhUn.Zero();
	    for(int i = 0; i < mat_wi.GetN(); i++)
	      KhUn += dt*mat_wi(k, i)*EvalF(i);
	    
	    Fn += KhUn;
          }	
      }
    
    // now we have to solve DenomR*U^n+1 = Fn
    for(int k = 0; k < NumerR.GetM()-1; k++)
      {
	sys.ApplyMassMatrix(Real_wp(1), t, Fn, Real_wp(0), KhUn); 
	sys.SolveOperatorDhPlusGammaKh(t, KhUn, Fn);
      }    
    
    // updating U^n to the next value
    Y = Fn;    
    //int test_input; cout << "waiting" << endl; cin >> test_input;
  }


  //! changes time step during the simulation
  template<class T>
  void LinearSdirkScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "LinearSdirkScheme_Iterator<T>::ChangeTimeStep not implemented" << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& LinearSdirkScheme_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& LinearSdirkScheme_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int LinearSdirkScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& LinearSdirkScheme_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t LinearSdirkScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = NumerR.GetMemorySize() + DenomR.GetMemorySize() + ci.GetMemorySize();
    taille += mat_wi.GetMemorySize() + polA.GetMemorySize() + real_roots.GetMemorySize()
      + complex_roots.GetMemorySize() + A_rk.GetMemorySize() + b_rk.GetMemorySize()
      + c_rk.GetMemorySize();
    
    taille += Y.GetMemorySize() + Fn.GetMemorySize() + KhUn.GetMemorySize();
    taille += Seldon::GetMemorySize(EvalF);
    return taille;
  }


  /***********************
   * DirkScheme_Iterator *
   ***********************/
  

  //! default constructor
  template<class T>
  DirkScheme_Iterator<T>::DirkScheme_Iterator()
  {
    SetOrder(4);
  }

  
  //! computation of coefficients for DIRK scheme of order r
  template<class T>
  void DirkScheme_Iterator<T>::SetOrder(int r)
  {
    Real_wp one(1);
    switch (r)
      {
      case 2 :
	{
	  // implicite midpoint rule
	  A.Reallocate(1, 1);
	  A(0, 0) = one/2;
	  c.Reallocate(1); c(0) = one/2;
	  b.Reallocate(1); b(0) = one;
	}
	break;
      case 3 :
	{
	  cout << "DIRK order 3 not implemented yet" << endl;
	  abort();
	}
	break;
      case 4 :
	{
	  //L-D L-D Dirk schemes A. Najafi-Yazdi and L. Mongeau
	  A.Reallocate(3, 3);
          A(0, 0) = Real_wp(0.377847764031163);
          A(1, 0) = Real_wp(0.385232756462588);
          A(1, 1) = Real_wp(0.461548399939329);
          A(2, 0) = Real_wp(0.675724855841358);
          A(2, 1) = Real_wp(-0.061710969841169);
          A(2, 2) = Real_wp(0.241480233100410);
	  
          c.Reallocate(3); b.Reallocate(3);
          b(0) = Real_wp(0.750869573741408);
	  b(1) = Real_wp(-0.362218781852651);
	  b(2) = Real_wp(0.611349208111243);
	  
          c(0) = Real_wp(0.257820901066211);
	  c(1) = Real_wp(0.434296446908075);
	  c(2) = Real_wp(0.758519768667167);
	  
	}
	break;
      case 5 :
	{
	  cout << "DIRK order 5 not implemented yet" << endl;
	  abort();
	}
	break;
      case 6 :
	{
	  cout << "DIRK order 6 not implemented yet" << endl;
	  abort();
	}
	break;
      default :
	{
	  cout << "Order " << r << " Not implemented " <<endl;
	  abort();
	}
      }
  }

  
  //! initialisation of the scheme with u^0
  template<class T>
  void DirkScheme_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                        Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();   
    sys.SetDirichletCondition(t, 0, Y);
    Fn.Reallocate(N); Fn.Fill(0);
    Un_tmp.Reallocate(N); Un_tmp.Fill(0);
    int nb_stages = A.GetM();
    VectReal_wp alpha, vect_dt_gamma;
    
    alpha.Reallocate(nb_stages);
    vect_dt_gamma.Reallocate(nb_stages);
    Zn.Reallocate(nb_stages);

    alpha.Fill(1.0);
    for (int i = 0; i < nb_stages; i++)
      {
	Zn(i).Reallocate(N);
	Zn(i).Fill(0);
	vect_dt_gamma(i) = dt*A(i,i);
      }

    //DISP(vect_dt_gamma);
    // we ask the user to factorize the matrix
    // M_h + dt gamma S_h + dt gamma K_h
    sys.FactorizeOperatorReal(alpha, vect_dt_gamma, vect_dt_gamma);

    //alpha.clear();
    //vect_dt_gamma.clear();
  }
  

  //! clearing arrays used by time scheme
  template<class T>
  void DirkScheme_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! clearing arrays used by time scheme
  template<class T>
  void DirkScheme_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Zn.Clear();
    Fn.Clear(); Un_tmp.Clear();
    Yfinal = Y; Y.Clear();
  }
  
  
  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void DirkScheme_Iterator<T>::Advance(const Real_wp& t, int n_time,
				       VirtualOdeSystem<T>& sys)
  {
    Real_wp tcurrent;

    // computation of all intermediary unknowns
    for (int i = 0; i < A.GetM(); i++)
      {
	// current time
	tcurrent = t + c(i)*dt;
	
	// computation of Un_tmp = U^n + delta t \sum_{j=1}^{i-1} a_ij z_j 
	Copy(Y, Un_tmp);
	for (int j = 0; j < i; j++)
	  Add(dt*A(i, j), Zn(j), Un_tmp);
        
	// computation of f(Un_tmp, tcurrent)
	sys.EvaluateFunction(tcurrent, Un_tmp, Fn, false);

	// then, we solve the i_th system z_i = f(Un_tmp + A(i,i) dt z_i, t^n + c_i dt)
	sys.SolveOperatorReal(tcurrent, Fn, Zn(i), i);
      }
    
    // Once we have all the intermediary unknowns, we deduce
    // U^n+1 = U^n + dt \sum b_i z_i
    for (int i = 0; i < b.GetM(); i++)
      Add(dt*b(i), Zn(i), Y);

  }


  //! changes time step during the simulation
  template<class T>
  void DirkScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& DirkScheme_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& DirkScheme_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int DirkScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& DirkScheme_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t DirkScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = A.GetMemorySize() + b.GetMemorySize() + c.GetMemorySize();
    taille += Y.GetMemorySize() + Fn.GetMemorySize() + Un_tmp.GetMemorySize();
    taille += Seldon::GetMemorySize(Zn);
    return taille;
  }

  
  /***********************
   * PadeScheme_Iterator *
   ***********************/
  

  //! default constructor
  template<class T>
  PadeScheme_Iterator<T>::PadeScheme_Iterator()
  {
    stable_algo = false;
    dt = Real_wp(0);
    
    SetOrder(2);
  }


  //! constructor with the order
  template<class T>
  PadeScheme_Iterator<T>::PadeScheme_Iterator(int r, bool s)
  {
    dt = Real_wp(0);
    SetOrder(r, s);
  }

  
  //! Polynomials P and Q are computed, Q is factorized
  template<class T>
  void PadeScheme_Iterator<T>::FindPadeCoeff(int r)
  {
    Real_wp one(1);
    //
    this->order = r;
    if(this->order%2 != 0)
      {
	cout << "Scheme not compatible with Diagonal Pade scheme" << endl;
	abort();
      }
    
    // Polynum coeff
    //Numerator
    unsigned int deg_num;
    UnivariatePolynomial<Real_wp> P;
    //Denominator
    unsigned int deg_deno;
    UnivariatePolynomial<Real_wp> Q;
    // roots of Q
    //VectComplex_wp rootsQ;
    // Real and imaginary part of complex roots
    VectReal_wp rRealQ;
    VectReal_wp rImagQ;

    deg_num = order/2;
    deg_deno = order - deg_num;
  
    P.SetOrder(deg_num);
    Q.SetOrder(deg_deno);
    
    // P_(m,n) = \sum_{i=0}^n (m+n-i)!s!/(m+n)!i!(n-i)! z^i
    unsigned int n = deg_num;
    unsigned int m = deg_deno;
    for(unsigned i = 0; i <= n; i++)
      {
	P(i)=this->fac(m+n-i)*this->fac(n)/(this->fac(m+n)*this->fac(i)*this->fac(n-i));
      }
    
    // Q_(m,n) = \sum_{i=0}^m (m+n-i)!m!/(m+n)!i!(m-i)! (-z)^i
    for(unsigned i = 0; i <= m; i++)
      {
	Q(i)=this->fac(m+n-i)*this->fac(m)/(this->fac(m+n)*this->fac(i)*this->fac(m-i));
	if(i%2 != 0)
	  Q(i) *= -one; 
      }
    
    //DISP(P); DISP(Q);
    //Q.WriteText("Denom.dat");
    
    //! Compute roots of Q
    SolvePolynomialEquation(Q, rRealQ, rImagQ);
    //! Sort roots regarding the real part
    Sort(rRealQ, rImagQ);
    
    //Separate complex root from the unique real possible root.
    int k = 0;
    real_root = Real_wp(0);
    complex_root.Clear();
    while (k < rRealQ.GetM())
      {
	if (rImagQ(k) == Real_wp(0))
	  {
	    real_root = rRealQ(k);
	    k++;
	  }
	else
	  {
	    complex_root.PushBack(Complex_wp(rRealQ(k), rImagQ(k)));
	    k += 2;
	  }
      }
    
    //DISP(complex_root);
    //DISP(real_root);
    
    //Comput P-Q
    P -= Q;
    //int nb_coef = (deg_num+1)/2;
    int nb_coef = deg_num;
    if (deg_num%2 == 1)
      nb_coef++;
    
    coefPmQ.Reallocate(nb_coef);
    coefPmQ.Zero();
    for (k = 0; k < nb_coef; k++)
      //coefPmQ(k) = P(2*k+1);
      coefPmQ(k) = P(k);
    
    //DISP(coefPmQ);
  }


  //! Compute ci and mat_wi use to evaluate source function for pade Scheme 
  template<class T>
  void PadeScheme_Iterator<T>::ComputeCoefForPadeSourceFunction()
  {
    int m = this->order/2;
    Real_wp one(1);

    UnivariatePolynomial<Real_wp> P;
    int n = m;
    P.SetOrder(m);
    for(int i = 0; i <= m; i++)
      P(i) = this->fac(m+n-i)*this->fac(n)/(this->fac(m+n)*this->fac(i)*this->fac(n-i));
    
    VectReal_wp CoefC(2*m+1);
    for (int k = 0; k <= 2*m; k++)
      CoefC(k) = one/(pow(Real_wp(2), k-1) * this->fac(k));
    
    // Computation of ci using Gauss-Legendre points
    VectReal_wp weights;
    if (m == 1)
      {
	this->ci.Reallocate(1);
	this->ci(0) = Real_wp(0.5);
      }
    else
      ComputeGaussLegendre(this->ci, weights, m-1);

    //DISP(this->ci); DISP(weights);
    this->ci = this->ci - Real_wp(0.5);

    // Vandermonde matrix for ci = cj^i
    // factorials are put in the right hand side, so we dont have c_j^i / j!
    Matrix<Real_wp> VDM(this->ci.GetM(), this->ci.GetM());
    for(int i = 0; i < this->ci.GetM(); i++)
      for(int j = 0; j < this->ci.GetM(); j++)
	VDM(i, j) = pow(this->ci(j), i);
    
    // Factorize VDM
    Vector<int> pivot;
    GetLU(VDM, pivot);

    // compute matrix A
    this->mat_wi.Reallocate(m, this->ci.GetM());
    Matrix<Real_wp> A(m, m);
    for (int powL = 0; powL < m; powL++)
      {
	int r = powL+1;
	int jmin = 2;
	if (r%2 == 1)
	  jmin = 1;

	VectReal_wp vec_f(2*m);
	vec_f.Zero();
	for (int j = jmin; j <= 2*m-r; j += 2)
	  {
	    for (int i = 0; i <= r-1; i += 2)
	      vec_f(j-1) += P(i)*CoefC(r+j-i-1);
	    
	    for (int i = 1; i <= r-1; i += 2)
	      vec_f(j-1) -= P(i)*CoefC(r+j-i-1);	    
	  }
        
        for (int j = 0; j < this->ci.GetM(); j++)
          A(powL, j) = vec_f(j);
        
	//DISP(powL); DISP(vec_f);
	VectReal_wp rhs(this->ci.GetM());
	rhs.Zero();

        if (!stable_algo)
          {
            for(int i = 0; i < this->ci.GetM(); i++)
              rhs(i) = this->fac(i)*vec_f(i);
            
            SolveLU(VDM, pivot, rhs); //DISP(rhs);
            SetRow(rhs, powL, this->mat_wi);
          }

        //this->mat_wi.WriteText("mat.dat");
      }

    if (stable_algo)
      {
        UnivariatePolynomial<Real_wp> Q, P, z;
        z.SetOrder(1); z(0) = Real_wp(0); z(1) = Real_wp(1);
        Q.SetOrder(0); Q(0) = Real_wp(1);
        
        Vector<UnivariatePolynomial<Real_wp> > PolNumer(m+1);       
        int nb = 0;
        PolNumer(nb) = Q; nb++;
        
        for (int i = complex_root.GetM()-1; i >= 0; i--)
          {
            PolNumer(nb) = z*Q; nb++;
            P.SetOrder(2);
            P(0) = Real_wp(1);
            P(1) = Real_wp(2)*realpart(Real_wp(1) / complex_root(i));
            P(2) = Real_wp(1)/absSquare(complex_root(i));
            Q *= P;
            PolNumer(nb) = Q; nb++;            
          }
        
        if (m%2 == 1)
          {
            P.SetOrder(1);
            P(0) = Real_wp(1); P(1) = Real_wp(1) / real_root;
            Q *= P;
            PolNumer(nb) = Q; nb++;
          }

        // polynomials needed to expand the source are enumerated
        Vector<UnivariatePolynomial<Real_wp> > pol_source(m);
        UnivariatePolynomial<Real_wp> Denom;
        Denom.SetOrder(0); Denom(0) = Real_wp(1);
        int num = 0;
        if (m%2 == 1)
          {
            pol_source(num) = PolNumer(nb-2-num);
            num++;
            
            // real root present, we multiply by (1 - z / real_root)
            P.SetOrder(1);
            P(0) = Real_wp(1); P(1) = -Real_wp(1) / real_root;
            Denom = P;
          }
            
        // loop on complex roots
        for (int i = 0; i < complex_root.GetM(); i++)
          {
            pol_source(num) = PolNumer(nb-2-num)*Denom;
            num++;
            pol_source(num) = PolNumer(nb-2-num)*Denom;
            num++;
            
            P.SetOrder(2);
            P(0) = Real_wp(1);
            P(1) = -Real_wp(2)*realpart(Real_wp(1) / complex_root(i));
            P(2) = Real_wp(1)/absSquare(complex_root(i));

            Denom *= P;
          }
        
        // inverse of coefficients, to express x^i as a combination of pol_source
        Matrix<Real_wp> coefBase;
        coefBase.Reallocate(m, m);
        coefBase.Zero();
        for (int j = 0; j < m; j++)
          for (int i = 0; i < pol_source(j).GetM(); i++)
            coefBase(i, j) = pol_source(j)(i);
        
        GetInverse(coefBase);
        
        Matrix<Real_wp> polA;
        polA.Reallocate(m, m);
        polA.Zero();
        for (int powA = 0; powA < m; powA++)
          {
            VectReal_wp rhs(m);
            rhs.Zero();	
            for (int j = 1; j <= m; j++)
              rhs(j-1) = A(powA, j-1);
            
            for (int j = 0; j < m; j++)
              for (int k = 0; k < m; k++)
                polA(k, j) += coefBase(k, powA)*rhs(j);
          }
        
        
        for (int powA = 0; powA < m; powA++)
          {
            VectReal_wp rhs(m);
            rhs.Zero();
            
            for(int i = 0; i < this->ci.GetM(); i++)
              rhs(i) = this->fac(i)*polA(powA, i);
            
            SolveLU(VDM, pivot, rhs);
            SetRow(rhs, powA, this->mat_wi);
          }        
      }
    
    // we get back to the interval [0, 1]
    this->ci = this->ci + Real_wp(0.5);
    
  }
  
  
  //! computation of coefficients for Pade scheme of order r
  template<class T>
  void PadeScheme_Iterator<T>::SetOrder(int r, bool s)
  {
    stable_algo = s;
    
    //Real_wp one(1);
    this->FindPadeCoeff(r);

    // Compute data needed (ci and wi) to evaluate the source function
    this->ComputeCoefForPadeSourceFunction();
  }


  //! initialisation of the scheme with u^0
  template<class T>
  void PadeScheme_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                        Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();   
    sys.SetDirichletCondition(t, 0, Y);
    Fn.Reallocate(N); Fn.Fill(0);
    EvalF.Reallocate(ci.GetM());
    for(int i = 0; i < ci.GetM(); i++)
      EvalF(i).Reallocate(N);

    Un_tmp.Reallocate(N); Un_tmp.Fill(0);
    KhUn.Reallocate(N); KhUn.Fill(0);

    Freal.Reallocate(N); Freal.Zero();

    Fn_cplx.Reallocate(N); Fn_cplx.Zero();
    KhUn_cplx.Reallocate(N); KhUn_cplx.Zero();

    // we ask the user to factorize the matrix if a real root is present
    if (this->order%4 == 2)
      {
	sys.FactorizeOperatorDhPlusGammaKh(Real_wp(1), dt/real_root, dt/real_root);
      }
    
    if (this->order > 2)
      {
	// there are complex roots in this case
	VectComplex_wp coef_cplx;
	coef_cplx = dt / complex_root;
	
	VectComplex_wp coef_ones(coef_cplx.GetM());
	coef_ones.Fill(1);
	
	sys.FactorizeOperatorComplex(coef_ones, coef_cplx, coef_cplx);
      }
  }
  

  //! clearing arrays used by time scheme
  template<class T>
  void PadeScheme_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! clearing arrays used by time scheme
  template<class T>
  void PadeScheme_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    // we clear only arrays containing iterates
    Fn.Clear(); Un_tmp.Clear(); KhUn.Clear();
    Freal.Clear(); 
    Fn_cplx.Clear();
    KhUn_cplx.Clear();
    for(int i = 0; i < EvalF.GetM(); i++)
      EvalF(i).Clear();
    EvalF.Clear();
    Yfinal = Y; Y.Clear();
  }

  //! For real source function (Paper Real valued iterative methods for solving complex linear systems)
  template<>
  void PadeScheme_Iterator<Real_wp>::SolveOperatorP2(const Real_wp& t, const VectReal_wp& f,
						     const Complex_wp& a, int k, VectReal_wp& y,
						     VirtualOdeSystem<Real_wp>& sys, bool apply_mass)
  {
    // we compute v = (I + a A)^-1 f = (M + a K)^{-1} M f
    if (apply_mass)
      {
	sys.ApplyMassMatrix(Real_wp(1), t, f, Real_wp(0), Freal);
	for (int i = 0; i < f.GetM(); i++)
	  Fn_cplx(i) = Complex_wp(Freal(i), 0);
      }
    else
      {
	for (int i = 0; i < f.GetM(); i++)
	  Fn_cplx(i) = Complex_wp(f(i), 0);
      }
    
    sys.SolveOperatorComplex(t, Fn_cplx, KhUn_cplx, k);

    // y = 2 Re(b v)
    Complex_wp b = a / (a - conj(a));
    for (int i = 0; i < f.GetM(); i++)
      y(i) = 2.0*realpart(b*KhUn_cplx(i));
  }


  //! For complex source function
  template<>
  void PadeScheme_Iterator<Complex_wp>::SolveOperatorP2(const Real_wp& t, const VectComplex_wp& f,
							const Complex_wp& a, int k, VectComplex_wp& y,
							VirtualOdeSystem<Complex_wp>& sys, bool apply_mass)
  {
    // We compute Freal = (I+aA)^{-1} f.real = (M+aK)^{-1} M f.real
    // and Fimag = (I+aA)^{-1} f.imag = (M+aK)^{-1} M f.imag    
    if (apply_mass)
      {
	sys.ApplyMassMatrix(Real_wp(1), t, f, Real_wp(0), Fn_cplx);
      }
    else
      {
	for (int i = 0; i < f.GetM(); i++)
	  Fn_cplx(i) = f(i);
      }

    int N = f.GetM();
    VectComplex_wp realf, imagf, realY, imagY;
    realf.Reallocate(N); realf.Zero();
    imagf.Reallocate(N); imagf.Zero();
    
    realY.Reallocate(N); realY.Zero();
    imagY.Reallocate(N); imagY.Zero();
    
    for(int i = 0; i < N; i++)
      {
	realf(i) = Complex_wp(realpart(Fn_cplx(i)), 0);
	imagf(i) = Complex_wp(imagpart(Fn_cplx(i)), 0);
      }
    
    sys.SolveOperatorComplex(t, realf, realY, k);
    sys.SolveOperatorComplex(t, imagf, imagY, k);
   
   // y = 2 Re(b v)
   Complex_wp b = a / (a - conj(a));
   for (int i = 0; i < f.GetM(); i++)
     y(i) = Complex_wp(2.0*realpart(b*realY(i)), 2.0*realpart(b*imagY(i)));
   
   realf.Clear(); realY.Clear();
   imagf.Clear(); imagY.Clear();
   
  }
  
  
  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void PadeScheme_Iterator<T>::Advance(const Real_wp& t, int n_time,
				       VirtualOdeSystem<T>& sys)
  {
    bool presence_source = false;
    if (t <= sys.GetFinalTimeSource())
      {
        for(int i = 0; i < this->mat_wi.GetN(); i++)
          {
            Real_wp tcurrent = t + ci(i)*dt;
            // compute y = y + alpha f(t)	
            EvalF(i).Zero();
            sys.AddPrimitiveTimeSource(Real_wp(1), tcurrent, 0, EvalF(i));
          }
        presence_source = true;
      }

    if (stable_algo)
      {
        // real part if present
        int num_source = 0;
        if (this->order%4 == 2)
          {
            Fn.Zero();
            if (presence_source)
              for (int i = 0; i < mat_wi.GetN(); i++)
                Fn += dt*mat_wi(num_source, i)*EvalF(i);
            
            num_source++;
            sys.ApplyMassMatrix(Real_wp(1), t, Y, Real_wp(1), Fn);            
            sys.EvaluateFunction(t, Y, KhUn, false, false);
            Real_wp coef = dt / real_root;
            Fn += coef * KhUn;
            
            sys.SolveOperatorDhPlusGammaKh(t, Fn, KhUn);
            Y = KhUn;
          }
        
        // loop on complex roots
        for (int k = 0; k < complex_root.GetM(); k++)
          {
            Real_wp b = Real_wp(2)*dt*realpart(Real_wp(1) / complex_root(k));
            Real_wp a = dt*dt/absSquare(complex_root(k));
            sys.EvaluateFunction(t, Y, Fn, false, false);
            Fn *= a;

            KhUn.Zero();
            if (presence_source)
              for (int i = 0; i < mat_wi.GetN(); i++)
                KhUn += dt*dt*mat_wi(num_source, i)*EvalF(i);
            
            Fn += KhUn;
            sys.SolveMassMatrix(Fn);
            KhUn = b*Y + Fn;
            sys.EvaluateFunction(t, KhUn, Fn, false, false);
            
            sys.ApplyMassMatrix(Real_wp(1), t, Y, Real_wp(1), Fn); 
            KhUn.Zero();
            if (presence_source)
              {
                for (int i = 0; i < mat_wi.GetN(); i++)
                  KhUn += dt*mat_wi(num_source+1, i)*EvalF(i);
                
                Fn += KhUn;
              }
            
            SolveOperatorP2(t, Fn, -Real_wp(1)/complex_root(k), k, Y, sys, false);
            num_source += 2;
          }
        
        return;
      }

    Fn.Zero();
    Un_tmp = Y;
    
    int last_coef = coefPmQ.GetM()-1;
    // Fn = (P-Q)_{last_coef} u + dt F_{last_coef}
    if (t <= sys.GetFinalTimeSource())
      if (this->order%4 == 0)
	{
	  for(int i = 0; i < mat_wi.GetN(); i++)
	    Fn += dt*mat_wi(last_coef, i)*EvalF(i);
	  
	  sys.SolveMassMatrix(Fn);
	}

    Fn += coefPmQ(last_coef)*Y;
    
    for (int k = last_coef-1; k >= 0; k--)
      {
	// Fn = dt L Fn + dt F_i
	sys.EvaluateFunction(t, Fn, KhUn, false, false);
	Fn = dt*KhUn; 	
	if (t <= sys.GetFinalTimeSource())
	  for(int i = 0; i < mat_wi.GetN(); i++)
	    Fn += dt*mat_wi(k, i)*EvalF(i);	
	
	if (k > 0)
	  sys.SolveMassMatrix(Fn);
	
	// we add (P-Q)_{k} u
	if (k%2 == 1)
	  Fn += coefPmQ(k)*Y;	
      }
    
    // now we have to solve Q (U^n+1 - U^n) = Fn
   
    // real part if present
    if (this->order%4 == 2)
      {
	sys.SolveOperatorDhPlusGammaKh(t, Fn, KhUn);
	Fn = KhUn;
      }
   
    // complex part
    if (this->order > 2)
      {
	for (int k = 0; k < complex_root.GetM(); k++)
	  {
	    bool apply_mass = true;
	    if ((k == 0) && (this->order%4 == 0))
	      apply_mass = false;
	    
	    SolveOperatorP2(t, Fn, -Real_wp(1)/complex_root(k), k, KhUn, sys, apply_mass);
	    Fn = KhUn;
	  }
      }    
   
    // updating U^n to the next value
    Y += Fn;
    
    //int test_input; cout << "waiting" << endl; cin >> test_input;
  }


  //! changes time step during the simulation
  template<class T>
  void PadeScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& PadeScheme_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& PadeScheme_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int PadeScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& PadeScheme_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t PadeScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = mat_wi.GetMemorySize() + complex_root.GetMemorySize() + coefPmQ.GetMemorySize();
    taille += Y.GetMemorySize() + Fn.GetMemorySize() + Un_tmp.GetMemorySize() + KhUn.GetMemorySize();
    taille += Freal.GetMemorySize() + Fn_cplx.GetMemorySize() + KhUn_cplx.GetMemorySize();
    taille += Seldon::GetMemorySize(EvalF);
    return taille;
  }
  
  
  /****************************
   * GaussRungeKutta_Iterator *
   ****************************/

  
  //! default constructor
  template<class T>
  GaussRungeKutta_Iterator<T>::GaussRungeKutta_Iterator()
  {
    dt = 0.0; SetOrder(2);
  }
  
  
  //! initialisation of coefficients
  /*!
    \param[in] order order of the method (2r+2)
  */
  template<class T>
  void GaussRungeKutta_Iterator<T>::SetOrder(int order_)
  {
    // retrieving r
    int r = order_/2-1;
    order = r;
    
    // then computation of the matrix a_{i,j} and the vector b_i
    AgaussRK.Reallocate(r+1, r+1);
    if (r == 0)
      {
        AgaussRK(0, 0) = 0.5;
        BgaussRK.Reallocate(r+1);
        BgaussRK(0) = 1.0;
	CgaussRK.Reallocate(1);
	CgaussRK(0) = 0.5;
      }
    else
      {
        Globatto<Real_wp> gauss;
        gauss.ConstructQuadrature(r, gauss.QUADRATURE_GAUSS);
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            {
              Real_wp vloc = 0;
              for (int k = 0; k <= r; k++)
                vloc += gauss.Points(j)*gauss.Weights(k)
                  *gauss.EvaluatePhi(i, gauss.Points(j)*gauss.Points(k));
              
              AgaussRK(j, i) = vloc;
            }
        
        CgaussRK = gauss.Points();
        BgaussRK = gauss.Weights();
      }
  }

  
  //! initialization of the time scheme with u^0
  template<class T>
  void GaussRungeKutta_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
			VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    Y = u0; u0.Clear();
    ki.Reallocate((order+1)*Y.GetM());
    ki.Fill(0);
    
    // factorisation of linear (if linear) system
    // Dh k_i - g(y^n + dt*\sum a_{i, j} k_j)
    sys.FactorizeOperatorDhPlusAijKh(AgaussRK, dt, Y, ki);
  }


  //! releasing memory used by arrays for the time scheme
  template<class T>
  void GaussRungeKutta_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! releasing memory used by arrays for the time scheme
  template<class T>
  void GaussRungeKutta_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    ki.Clear();
    Yfinal = Y; Y.Clear();
  }
  
  
  //! computation of U^n+1 from U^n
  template<class T>
  void GaussRungeKutta_Iterator<T>
  ::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    VectReal_wp t_sub(order+1);
    for (int i = 0; i <= order; i++)
      t_sub(i) = t + CgaussRK(i)*dt;
    
    // solving the system Dh k_i - g(y^n + dt \sum a_ij k_j) = 0
    sys.SolveOperatorDhPlusAijKh(t_sub, AgaussRK, dt, Y, ki);
    
    // then computing U^n+1 = U^n + dt \sum b_i k_i
    int N = Y.GetM();
    for (int i = 0; i <= order; i++)
      for (int k = 0; k < N; k++)
        Y(k) = Y(k) + dt*BgaussRK(i)*ki(k*(order+1)+i);
  }
  

  //! changes time step during the simulation
  template<class T>
  void GaussRungeKutta_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& GaussRungeKutta_Iterator<T>::GetIterate()
  {
    return Y;
  }

  
  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& GaussRungeKutta_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */  
  template<class T>
  int GaussRungeKutta_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& GaussRungeKutta_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t GaussRungeKutta_Iterator<T>::GetMemorySize() const
  {
    size_t taille = AgaussRK.GetMemorySize() + BgaussRK.GetMemorySize() + CgaussRK.GetMemorySize();
    taille += Y.GetMemorySize() + ki.GetMemorySize();
    return taille;
  }

  
  /**************************
   * AdamsImplicit_Iterator *
   **************************/


  //! default constructor
  template<class T>
  AdamsImplicit_Iterator<T>::AdamsImplicit_Iterator()
    : AdamsBashforth_Moulton_Iterator<T>()
  { 
    this->type_algo = this->MOULTON;
  }
  
  
  //! sets the order of approximation
  template<class T>
  void AdamsImplicit_Iterator<T>::SetOrder(int order)
  {
    // implicit scheme is enforced
    AdamsBashforth_Moulton_Iterator<T>::SetOrder(order, this->MOULTON);
  }
  

  //! initialisation of the scheme with u^0
  template<class T>
  void AdamsImplicit_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
			VirtualOdeSystem<T>& sys)
  {
    AdamsBashforth_Moulton_Iterator<T>::SetInitialCondition(t, deltat, u0, sys);
    
    // deallocating unnecessary vectors
    int order = this->beta.GetM();
    this->fi(order-1).Clear();
    
    sys.FactorizeOperatorDhPlusGammaKh(1.0, this->dt*this->beta(0), this->dt*this->beta(0));
  }
  
  
  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void AdamsImplicit_Iterator<T>::Advance(const Real_wp& t, int n_time,
					  VirtualOdeSystem<T>& sys)
  {
    int order = this->beta.GetM();
    if (n_time < order-1)
      {
        // using RK4 scheme for the first iterates
        Copy(this->Y, this->Y_star);
        sys.EvaluateFunction(t, this->Y, this->fi(0));
        Add(this->dt/6, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt/2, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+0.5*this->dt, this->F_star, this->fi(0));
        Add(this->dt/3, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt/2, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+0.5*this->dt, this->F_star, this->fi(0));
        Add(this->dt/3, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+this->dt, this->F_star, this->fi(0));
        Add(this->dt/6, this->fi(0), this->Y_star);
        Copy(this->Y_star, this->Y);

	// once Y_star = y^n+1 is computed
	// f(y^n+1) is stored in the vector fi
	sys.EvaluateFunction(t+this->dt, this->Y, this->fi(order-2-n_time));
        
        return;
      }
    
    // forming right-hand side
    // F_star = y^n + dt \sum_{i=1}^order beta_i f(t^{n+1-i}, y^{n+1-i})
    Copy(this->Y, this->F_star);
    for (int i = 0; i < order-1; i++)
      Add(this->dt*this->beta(i+1), this->fi(i), this->F_star);
    
    // then solving the system (linear or non-linear)
    // (D_h y^{n+1} - dt beta(0) f(t^{n+1}, y^{n+1}) = F_star
    sys.SolveOperatorDhPlusGammaKh(t+this->dt, this->F_star, this->Y);
    
    // updating the scheme
    for (int i = order-2; i > 0; i--)
      Copy(this->fi(i-1), this->fi(i));
    
    if (order > 1)
      sys.EvaluateFunction(t+this->dt, this->Y, this->fi(0));
  }


  /*************************
   * MilneSimpson_Iterator *
   *************************/


  //! default constructor
  template<class T>
  MilneSimpson_Iterator<T>::MilneSimpson_Iterator()
    : Nystrom_Iterator<T>()
  {
    this->type_algo = this->MILNE;
  }
  
  
  //! sets the order of approximation
  template<class T>
  void MilneSimpson_Iterator<T>::SetOrder(int order)
  {
    if (order == 2)
      {
        cout << "Explicit method" << endl;
        abort();
      }
    
    // implicit scheme is enforced
    Nystrom_Iterator<T>::SetOrder(order, this->MILNE);
  }
  

  //! initialisation of the scheme with u^0
  template<class T>
  void MilneSimpson_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
			VirtualOdeSystem<T>& sys)
  {
    Nystrom_Iterator<T>::SetInitialCondition(t, deltat, u0, sys);
    
    // deallocating unnecessary vectors
    int order = this->beta.GetM();
    this->fi(order-1).Clear();
    
    sys.FactorizeOperatorDhPlusGammaKh(1.0, this->dt*this->beta(0), this->dt*this->beta(0));
  }
  
  
  //! we compute next iterate U^n+1 from U^n
  template<class T>
  void MilneSimpson_Iterator<T>::Advance(const Real_wp& t, int n_time,
					 VirtualOdeSystem<T>& sys)
  {
    int order = this->beta.GetM();
    if (n_time < order-1)
      {
        // using RK4 scheme for the first iterates
        Copy(this->Y, this->Y_star);
        sys.EvaluateFunction(t, this->Y, this->fi(0));
        Add(this->dt/6, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt/2, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+0.5*this->dt, this->F_star, this->fi(0));
        Add(this->dt/3, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt/2, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+0.5*this->dt, this->F_star, this->fi(0));
        Add(this->dt/3, this->fi(0), this->Y_star);
        
        Copy(this->Y, this->F_star); Add(this->dt, this->fi(0), this->F_star);
        sys.EvaluateFunction(t+this->dt, this->F_star, this->fi(0));
        Add(this->dt/6, this->fi(0), this->Y_star);
        Copy(this->Y, this->Ynm1);
        Copy(this->Y_star, this->Y);

	// once Y_star = y^n+1 is computed
	// f(y^n+1) is stored in the vector fi
	sys.EvaluateFunction(t+this->dt, this->Y, this->fi(order-2-n_time));
        
        return;
      }
    
    // forming right-hand side
    // F_star = y^{n-1} + dt \sum_{i=1}^order beta_i f(t^{n+1-i}, y^{n+1-i})
    Copy(this->Ynm1, this->F_star);
    for (int i = 0; i < order-1; i++)
      Add(this->dt*this->beta(i+1), this->fi(i), this->F_star);
    
    // then solving the system (linear or non-linear)
    // (D_h y^{n+1} - dt beta(0) f(t^{n+1}, y^{n+1}) = F_star
    sys.SolveOperatorDhPlusGammaKh(t+this->dt, this->F_star, this->Y_star);
    
    // updating scheme
    Copy(this->Y, this->Ynm1);
    Copy(this->Y_star, this->Y);
    
    // updating the scheme
    for (int i = order-2; i > 0; i--)
      Copy(this->fi(i-1), this->fi(i));
    
    if (order > 2)
      sys.EvaluateFunction(t+this->dt, this->Y_star, this->fi(0));
  }

  
  /************************************
   * BackwardDifferentiation_Iterator *
   ************************************/


  //! sets the order of time scheme
  template<class T>
  void BackwardDifferentiation_Iterator<T>::SetOrder(int r)
  {
    order = r;
    VectReal_wp delta(r+1);
    delta(0) = 0.0;
    for (int i = 1; i <= r; i++)
      delta(i) = Real_wp(1)/i;
    
    // computation of beta from delta
    beta.Reallocate(r+1);
    beta(0) = delta(r);
    VectReal_wp new_beta(order+1);
    for (int m = order-1; m >= 0; m--)
      {
        new_beta(0) = delta(m) + beta(0);
        for (int j = 1; j <= order-1-m; j++)
          new_beta(j) = -beta(j-1) + beta(j);
        
        new_beta(order-m) = -beta(order-1-m);
        beta = new_beta;
      }
  }
    
  
  //! inits time scheme with the initial condition
  template<class T>
  void BackwardDifferentiation_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat, Vector<T>& u0,
                        VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();   
    sys.SetDirichletCondition(t, 0, Y);
    Yn_history.Reallocate(beta.GetM()-2);
    for (int k = 0; k < beta.GetM()-2; k++)
      {
        Yn_history(k).Reallocate(N);
        Yn_history(k).Fill(0);
      }
    
    T zero; SetComplexZero(zero);
    Y_star.Reallocate(N); F_star.Reallocate(N);
    Y_star.Fill(zero); F_star.Fill(zero);
    Fy.Reallocate(N); Fy.Fill(zero);
    
    sys.FactorizeOperatorDhPlusGammaKh(beta(0), dt, dt);
  }
  
  
  //! releases memory used by time scheme
  template<class T>
  void BackwardDifferentiation_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! releases memory used by time scheme
  template<class T>
  void BackwardDifferentiation_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Yn_history.Clear();
    Yfinal = Y; Y.Clear();
  }
    
  
  //! computes y^n+1 from y^n, y^n is overwritten by y^n+1
  template<class T>
  void BackwardDifferentiation_Iterator<T>::Advance(const Real_wp& t, int n_time,
						    VirtualOdeSystem<T>& sys)
  {
    T zero; SetComplexZero(zero);
    int order = beta.GetM();
    if (n_time < order-2)
      {
	// y^n is stored in the vector Yn_history
        Copy(Y, Yn_history(order-3-n_time));
        
        // using RK4 scheme for the first iterates
        Copy(Y, Y_star);
        sys.EvaluateFunction(t, Y, Fy);
        Add(dt/6, Fy, Y_star);
        
        Copy(Y, F_star); Add(dt/2, Fy, F_star);
        sys.EvaluateFunction(t+0.5*dt, F_star, Fy);
        Add(dt/3, Fy, Y_star);
        
        Copy(Y, F_star); Add(dt/2, Fy, F_star);
        sys.EvaluateFunction(t+0.5*dt, F_star, Fy);
        Add(dt/3, Fy, Y_star);
        
        Copy(Y, F_star); Add(dt, Fy, F_star);
        sys.EvaluateFunction(t+dt, F_star, Fy);
        Add(dt/6, Fy, Y_star);
        Copy(Y_star, Y);

        return;
      }
    
    // forming right-hand side
    // F_star = -\sum_{k=1} beta_k y^{n-k}
    Copy(Y, Fy);
    Copy(Y, F_star); Mlt(-beta(1), F_star);
    for (int i = 2; i < order; i++)
      Add(-beta(i), Yn_history(i-2), F_star); 
    
    // then solving the system (linear or non-linear)
    // (beta D_h y^{n+1} - dt f(t^{n+1}, y^{n+1}) = F_star
    sys.SolveOperatorDhPlusGammaKh(t+dt, F_star, Y);
    
    // updating the scheme
    for (int i = order-3; i > 0; i--)
      Copy(Yn_history(i-1), Yn_history(i));

    if (order > 2)
      Copy(Fy, Yn_history(0));
  }
  

  //! changes time step during the simulation
  template<class T>
  void BackwardDifferentiation_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! returns current iterate y^n
  template<class T>
  Vector<T>& BackwardDifferentiation_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns current iterate y^n
  template<class T>
  const Vector<T>& BackwardDifferentiation_Iterator<T>::GetIterate() const
  {
    return Y;
  }

  
  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int BackwardDifferentiation_Iterator<T>::GetNumberOfIterates() const
  {
    return beta.GetM()-1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& BackwardDifferentiation_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else
      return Yn_history(k-1);
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t BackwardDifferentiation_Iterator<T>::GetMemorySize() const
  {
    size_t taille = beta.GetMemorySize();
    taille += Y.GetMemorySize() + Y_star.GetMemorySize() + F_star.GetMemorySize();
    taille += Fy.GetMemorySize() + Seldon::GetMemorySize(Yn_history);
    return taille;
  }

  
  /***************************************
   * ImplicitSymmetricMultistep_Iterator *
   ***************************************/

  
  //! default constructor
  template<class T>
  ImplicitSymmetricMultistep_Iterator<T>::ImplicitSymmetricMultistep_Iterator()
  {
    dt = 0.0;
    t_begin_ener = 0.0; t_end_ener = 0.0; deltat_ener = 0.0;
    size_buffer_ener = 10;
    file_output_ener = "MultistepEnergy.txt";
  }

  template<class T>
  void  ImplicitSymmetricMultistep_Iterator<T>::SetParametersOutputEnergy(const Vector<string>& parameters)
  {
    //DISP(parameters);
    t_begin_ener = to_num<Real_wp>(parameters(0));
    t_end_ener = to_num<Real_wp>(parameters(1));
    deltat_ener = to_num<Real_wp>(parameters(2));

    if (parameters.GetM() > 3)
      size_buffer_ener = to_num<int>(parameters(3));
    
  }
  
  template<class T>
  void  ImplicitSymmetricMultistep_Iterator<T>::SetEnergyFile(const Vector<string>& parameters)
  { 
   file_output_ener = parameters(0);
  }
  

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::GetAlphaBeta(VectReal_wp& A, VectReal_wp& B)
  {
    A = alpha;
    B = beta;
  }

  
  //! finds coefficients alpha and beta from parameters theta
  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::FindCoefAlphaBeta(const VectReal_wp& theta, const Real_wp& phi)
  {
    // expressing polynom rho(zeta) from parameters theta
    Real_wp one(1);
    UnivariatePolynomial<Real_wp> rho, zeta, pol;
    zeta.SetOrder(1);
    zeta(0) = 0; zeta(1) = one;
    pol.SetOrder(2);
    rho = (zeta-1)*(zeta-1);
    for (int i = 0; i < theta.GetM(); i++)
      {
	pol(0) = one; pol(1) = -2.0*cos(theta(i)); pol(2) = one;
	rho *= pol;
      }
    
    int r = theta.GetM()*2+2;
    beta.Reallocate(r+1);
    alpha.Reallocate(r+1);
    beta.Fill(0);
    
    // solving the linear system to find beta
    int N = beta.GetM()-2;
    Matrix<Real_wp> sys(N, N);
    Vector<Real_wp> row(N+2), rhs(N);
    row.Fill(one);
    for (int s = 2; s < N+2; s++)
      {
	rhs(s-2) = row(N+1)*(N+1)*(N+1);
	for (int i = 1; i < N+1; i++)
	  {
	    sys(s-2, i-1) = row(i);
	    rhs(s-2) += rho(i)*row(i)*i*i;
	  }
        
	rhs(s-2) /= s*(s-1);
        
	if (s == 2) 
	  rhs(0) -= 2.0 * phi;
	else
	  rhs(s-2) -= row(N+1) * phi;

	for (int i = 1; i < N+2; i++)
	  row(i) *= i; 
      }

    IVect pivot(N);
    GetLU(sys, pivot);    
    SolveLU(sys, pivot, rhs);
    
    // filling vectors alpha and beta
    //alpha.Reallocate(N+2);
    alpha(0) = one; alpha(N+1) = one;
    for (int i = 1; i < N+1; i++)
      alpha(i) = rho(i);
    
    beta(0) = phi; beta(N+1) = phi;
    for (int i = 1; i < N+1; i++)
      beta(i) = rhs(i-1);
  }

  //! inits the time scheme with an order and parameters theta
  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::SetOrder(int order, const VectReal_wp& theta, const Real_wp& phi)
  {
    if (order == 2)
      {
	alpha.Reallocate(3);
	beta.Reallocate(3);
	alpha(0) = 1.0; alpha(1) = -2.0; alpha(2) = 1.0;
	beta(0) = phi; beta(1) = 1.0 - 2.0*phi; beta(2) = phi;
      }
    else if (order%2 == 0)
      FindCoefAlphaBeta(theta, phi);
    else
      {
	cout << "The order must be even" << endl;
	abort();
      }
    
    if (alpha(0) != 1.0)
      {
	Real_wp coef = 1.0/alpha(0);
	Mlt(coef, alpha);
	Mlt(coef, beta);
      }
    
    rk_nystrom.SetOrder(order);
  }

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                         Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys)
  {
    int N = u0.GetM();  
    int order = alpha.GetM();
    
    output_ener.SetDoublePrecision(true);
    output_ener.Init(file_output_ener, size_buffer_ener, true);

    dt = deltat;
    yk.Reallocate(order-1);
    fk.Reallocate(order-1);
    for (int i = 0; i < yk.GetM(); i++)
      {
	yk(i).Reallocate(N);
	yk(i).Fill(0);

	fk(i).Reallocate(N);
	fk(i).Fill(0);
      }
    
    yk(0) = u0; 
    sys.EvaluateFunctionS(t, yk(0), yk(0), fk(0), false);
     
    Y_star.Reallocate(N);
    Y_star.Fill(0);
    F_star.Reallocate(N);
    F_star.Fill(0);

    // initializing RKN scheme
    rk_nystrom.SetInitialConditionS(t, deltat, u0, du0_dt, sys);
    u0.Clear(); du0_dt.Clear();

    sys.FactorizeOperatorDhPlusGammaKh(alpha(0), 0.0, dt*dt*beta(0));

  }

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::Clear()
  {
    yk.Clear();
    fk.Clear();
    Y_star.Clear();
    F_star.Clear();

    output_ener.CloseBuffer();
  }
  
  
  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::ClearSecond(Vector<T>& y0, Vector<T>& dy0)
  {
    y0 = yk(0);
    dy0.Clear();
    
    yk.Clear();
    fk.Clear();
    Y_star.Clear();
    F_star.Clear();

    output_ener.CloseBuffer();
  }

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    int order = beta.GetM();
    Real_wp zero, one; SetComplexZero(zero); SetComplexOne(one);
    if (n_time < order-2)
      {
	// first iterates are advanced with Runge-Kutta-Nystrom
	rk_nystrom.Advance(t, n_time, sys);
	Y_star = rk_nystrom.GetIterate();
	if (n_time == order-3)
	  rk_nystrom.Clear();
      }
    else
      {
	// on evalue F_star =  -\sum \alpha_k yk
	F_star.Zero();
	Y_star.Zero();
	
	for (int k = 0; k < order-1; k++)
	  Add(alpha(k+1), yk(k), Y_star);
	
	// on applique l operateur de masse a F_star
	sys.ApplyOperatorDh(-one, t, Y_star, zero, F_star);

	// on ajoute le terme source dt^2*phi*F
	sys.AddPrimitiveTimeSource(dt*dt*beta(0), t+dt, 1, F_star);
        
    	// on rajoute \sum \beta_k fk
	for (int k = 0; k < order-1; k++)
	  Add(dt*dt*beta(k+1), fk(k), F_star);

        // inverse par D_h + gamma K_h
	sys.SolveOperatorDhPlusGammaKh(t+dt, F_star, Y_star);

      }
    
    //Evaluation de l energie
    WriteEnergy(t,sys);

    // on decale les y_{n-k} stockes
    for (int k = order-3; k >= 0; k--)
      Copy(yk(k), yk(k+1));
    
    // on decale les f_{n-k} stockes
    for (int k = order-3; k >= 0; k--)
      Copy(fk(k), fk(k+1));
    
    // on calcule le nouveau f_{n+1}
    Copy(Y_star, yk(0));
    // on n'utilise pas EvaluateFunctionS car il appelle Assemble
    //sys.EvaluateFunctionS(t+dt, yk(0), yk(0), fk(0), false);
    
    sys.ApplyOperatorKh(one, t+dt, yk(0), zero, fk(0));
    sys.AddScalarTimeSource(one, t+dt, 1, fk(0));
    sys.SetDirichletCondition(t+dt, 2, fk(0));
  }

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::ChangeTimeStep(const Real_wp&)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  template<class T>
  Vector<T>& ImplicitSymmetricMultistep_Iterator<T>::GetIterate()
  {
    return yk(0);
  }
  
  
  template<class T>
  const Vector<T>& ImplicitSymmetricMultistep_Iterator<T>::GetIterate() const
  {
    return yk(0);
  }
  
    
  template<class T>
  int ImplicitSymmetricMultistep_Iterator<T>::GetNumberOfIterates() const
  {
    return alpha.GetM();
  }
  
  
  template<class T>
  Vector<T>& ImplicitSymmetricMultistep_Iterator<T>::GetIterate(int k)
  {
    return yk(k);
  }


  template<>
  Real_wp ImplicitSymmetricMultistep_Iterator<Complex_wp>
  ::ComputeEnergy(const Real_wp& t, VirtualOdeSystem<Complex_wp>& sys, TinyVector<Real_wp, 4>& energy)
  {
    cout << "Not implemented for complex numbers" << endl;
    abort();
  }


  template<>
  Real_wp ImplicitSymmetricMultistep_Iterator<Real_wp>
  ::ComputeEnergy(const Real_wp& t, VirtualOdeSystem<Real_wp>& sys, TinyVector<Real_wp, 4>& energy)
  {
    int order = beta.GetM()-1;
    Real_wp zero,one; SetComplexZero(zero); SetComplexOne(one);

    int dim=Y_star.GetSize();
    Vector<Real_wp> Tmp(dim), Tmp2(dim), Tmp3(dim);
    Tmp.Fill(zero);
    Tmp2.Fill(zero);
    Tmp3.Fill(zero);
    Real_wp energy_tot(0);
    if(order ==2)
      {
        Real_wp Coeff[4] = {1.0/dt,-1.0/dt,0.5,0.5};
        Real_wp mon_theta=  beta(0);
        Real_wp w = dt*dt*(mon_theta- 0.25);
        // Modified mass matrix for kinetic energy
        for (int i=0;i<2;i++)
          Add(Coeff[i], yk(i), Tmp);
	
        sys.ApplyOperatorKh(-w,zero,Tmp,zero,Tmp2);
        sys.ApplyOperatorDh(one,zero,Tmp,zero,Tmp3);
        Add(Real_wp(1.),Tmp2,Tmp3);
        
        energy(1) = DotProd(Tmp,Tmp3);
        
        // Potential energy    
        Tmp.Fill(zero);
        Tmp2.Fill(zero);
        for (int i=0;i<2;i++)
          Add(Coeff[2+i], yk(i), Tmp);
        sys.ApplyOperatorKh(-Real_wp(1.0),0.,Tmp,Real_wp(0.),Tmp2);
        
        energy(2) = DotProd(Tmp,Tmp2);

      }
    else
      {
        Real_wp n=Real_wp(2.), w=-Real_wp(1.)/6.;
        
        Real_wp a(alpha(4));
        Real_wp b(alpha(3)+alpha(4)),c(alpha(2)+alpha(3)+alpha(4)),d(-alpha(0)),e(alpha(4)/n),f((alpha(3)+3*alpha(4))/n),g((alpha(2)+3*alpha(3)+5*alpha(4))/n),h(alpha(0)/n);
        Real_wp Coeff[8] = {a,b,c,d,e,f,g,h};
        
        // Coeff contains the kinetic and potential coefficients
        for (int i=0;i<4;i++)
          Add(Coeff[i], yk(i), Tmp);
        Mlt(Real_wp(1.)/dt,Tmp);
        // Now Tmp contains (sum_{i=0}^3 coeff[i] * U^{n+4-i})/Deltat
        
        // Modified part of the kinetic energy 
        sys.ApplyOperatorKh(-w*dt*dt,zero,Tmp,zero,Tmp2);
        
        // Mass part of the kinetic energy
        sys.ApplyOperatorDh(one,zero,Tmp,zero,Tmp3);
        Add(Real_wp(1.),Tmp2,Tmp3);
        energy(1) = DotProd(Tmp,Tmp3);
        
        // Potential energy    
        Tmp.Fill(zero);
        Tmp2.Fill(zero);
        for (int i=0;i<4;i++)
          Add(Coeff[4+i], yk(i), Tmp);
        sys.ApplyOperatorKh(-Real_wp(1.0),0.,Tmp,Real_wp(0.),Tmp2);
        
        energy(2) = DotProd(Tmp,Tmp2);
      }
    energy_tot = energy(1)+ energy(2);
    energy(3) = energy_tot;
    energy(0) = t;
    
    return energy_tot;

  }

  template<class T>
  void ImplicitSymmetricMultistep_Iterator<T>::WriteEnergy(const Real_wp& t, VirtualOdeSystem<T>& sys)
  {
    bool write_ener(false);
    if (deltat_ener !=0)
      {
	int n(toInteger(round(t/deltat_ener)));
	if ((t >= t_begin_ener) && (t < t_end_ener) && (abs(t/deltat_ener - n) < dt*0.5))
	  write_ener = true;
    
	// Real_wp ener_tot(0.);
	TinyVector<Real_wp,4> ener_vect;
	
	if (write_ener)
	  {
	    //ener_tot = ComputeEnergy(t, sys, ener_vect);
	    output_ener.AddTinyVect(ener_vect);
	    

            /* if (rank_processor == 0){
	      cout << rank_processor << " At time : " << t << ",";
              cout << "Total energy :" << ener_tot << endl; 
              } */
	  }	
      }
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t ImplicitSymmetricMultistep_Iterator<T>::GetMemorySize() const
  {
    size_t taille = alpha.GetMemorySize() + beta.GetMemorySize();
    taille += Y_star.GetMemorySize() + F_star.GetMemorySize();
    taille += rk_nystrom.GetMemorySize() + Seldon::GetMemorySize(yk) + Seldon::GetMemorySize(fk);
    return taille;
  }

}

#define MONTJOIE_FILE_IMPLICIT_TIME_SCHEMES_CXX
#endif
