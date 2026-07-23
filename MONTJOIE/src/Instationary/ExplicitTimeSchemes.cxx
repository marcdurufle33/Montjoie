#ifndef MONTJOIE_FILE_EXPLICIT_TIME_SCHEMES_CXX

namespace Montjoie
{
  
  //! Solves an evolution system du/dt = g(t, u) with any scheme
  /*!
    \param[in] t_begin initial time t0
    \param[in] t_end final time T
    \param[in] deltat time step
    \param[inout] Y0 solution for t = t0 on input, final solution on output
    \param[in] sys description of the function g(t, u)
    \param[in] scheme time scheme used    
    
    g(t, u) is given by calling the method sys.EvaluateFunction(t, u, result)
    
    the intermediary solutions (at times t0+dt, t0+2 dt, ..., t_end)
    are provided to the object sys by calling the functions GiveIterate and GiveFinalIterate
    
    The number of time iterations is provided to the object sys
    by calling the function GiveNumberIterations
    
    If you want to provide your own system of equations (i.e. the function g) in the object sys,
    you need to define a class with the following minimal methods (example with real numbers) :
    
    class MyEdoToSolve : public VirtualOdeSystem<Real_wp>
    {
    public :
    
    void GiveNumberIterations(const Real_wp& dt, int nb_iter);
    void GiveIterate(int n, const Real_wp& t, VectReal_wp& yn);
    void GiveFinalIterate(int n, const Real_wp& t, VectReal_wp& yn);
    
    // for almost all the schemes
    // computes result = g(t, u)
    void EvaluateFunction(const Real_wp& t, const VectReal_wp& u, VectReal_wp& result);

    // if you have a Dirichlet condition : u_i = h_i(t)
    // you will have an equation du_i/dt = dh_i/dt
    // but you can force Dirichlet condition on u_i for the first iterates
    // this function should set y_i = d^n h_i/dt^n for Dirichlet unknowns
    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& y);

    // depending on the scheme, you could need to overload other functions
    // see the documentation on the desired scheme to see if it needs additional functions
    };

    If you want to provide your own time scheme in the object scheme, you need to define a 
    class with the following methods (example with real numbers)
    
    class MyScheme : public VirtualTimeScheme<Real_wp>
    {
    public :

    // provides the initial condition y^0 to the scheme and time step
    // this function should be the moment to allocate all the vectors needed to advance the scheme
    void SetInitialCondition(const Real_wp& t, const Real_wp& dt, VectReal_wp& Y, VirtualOdeSystem<Real_wp>& sys);
    
    // returns a reference to the current iterate y^n
    const VectReal_wp& GetIterate();

    // this method should compute y^{n+1} from y^n
    // usually the values are overwritten such that the storage
    // is reduced to a limited number of vectors
    void Advance(const Real_wp& t, int n, VirtualOdeSystem<Real_wp>& sys);
    
    // this method should release memory used by the scheme and set y to the last iterate y^N
    void Clear(VectReal_wp& yn);
    
    };
  */
  template<class T>
  void RunTimeScheme(const Real_wp& t_begin, const Real_wp& t_end, const Real_wp& deltat,
		     Vector<T>& Y0, VirtualOdeSystem<T>& sys, VirtualTimeScheme<T>& scheme,
                     int num_iter0)
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
    
    dt = deltat; t = t_begin + num_iter0*dt;
    sys.GiveNumberIterations(dt, nb_max_iter);

    // first we initialize the scheme
    scheme.SetInitialCondition(t, deltat, Y0, sys);
       
    // iterations in time
    // main loop
    for (int n_time = num_iter0; n_time < nb_max_iter; n_time++)
      {
	t = t_begin + n_time*dt;
	// the iterate is provided 
	sys.GiveIterate(n_time, t, scheme.GetIterate());
	
	// Y^n+1 is computed
	scheme.Advance(t, n_time, sys);
      }
    
    // end of the scheme
    t = t_begin + nb_max_iter*dt;
    sys.GiveFinalIterate(nb_max_iter, t, scheme.GetIterate());

    // clearing intermediary vectors and setting Y0 to the last solution
    scheme.ClearFirst(Y0);
  }
  
  
  /***********************
   * RungeKutta_Iterator *
   ***********************/
  
  
  //! default constructor
  template<class T>
  RungeKutta_Iterator<T>::RungeKutta_Iterator()
  {
    dt = Real_wp(0); type_algo = BUTCHER;
  }
  
  
  //! computation of Runge-Kutta coefficients
  /*!
    \param[in] order order of RK scheme
  */
  template<class T>
  void RungeKutta_Iterator<T>::SetOrder(int order, int type)
  {
    type_algo = type;
    Real_wp zero, one, two(2);
    zero = 0.0; one = 1.0;
    
    order_rk = order;
    if (order == 1)
      {
	// Euler explicit rule
	A.Reallocate(1,1); A(0, 0) = zero;
	B.Reallocate(1); B(0) = one;
      }
    else if (order == 2)
      {
	// mid-point rule
	A.Reallocate(2, 2); A.Fill(0); A(1, 0) = one/2;
	B.Reallocate(2); B(0) = zero; B(1) = one;
      }
    else if (order == 3)
      {
	A.Reallocate(3,3); A.Fill(0); 
	A(1, 0) = two/3; A(2, 0) = one/3; A(2, 1) = one/3;
	B.Reallocate(3); B(0) = one/4; B(1) = zero; B(2) = Real_wp(3)/4;
      }
    else if (order == 4)
      {
	// classical fourth-order Runge-Kutta
	A.Reallocate(4, 4); A.Fill(0); 
	A(1, 0) = one/2; A(2, 1) = one/2; A(3, 2) = one;
	B.Reallocate(4); B(0) = one/6; B(1) = one/3;
	B(2) = one/3; B(3) = one/6;
      }
    else if (order == 5)
      {
        if (type_algo == FEHLBERG)
          {
            A.Reallocate(6, 6); A.Fill(zero);
            A(1, 0) = one/4; A(2, 1) = Real_wp(9)/32; A(3, 2) = Real_wp(7296)/2197;
            A(4, 3) = Real_wp(-845)/4104; A(5, 4) = Real_wp(-11)/40;
            A(2, 0) = Real_wp(3)/32; A(3, 1) = Real_wp(-7200)/2197;
            A(4, 2) = Real_wp(3680)/513; A(5, 3) = Real_wp(1859)/4104;
            A(3, 0) = Real_wp(1932)/2197; A(4, 1) = Real_wp(-8); A(5, 2) = Real_wp(-3544)/2565;
            A(4, 0) = Real_wp(439)/216; A(5, 1) = Real_wp(2); A(5, 0) = Real_wp(-8)/27;
            
            B.Reallocate(6);
            B(0) = Real_wp(16)/135; B(1) = zero; B(2) = Real_wp(6656)/12825;
            B(3) = Real_wp(28561)/56430; B(4) = Real_wp(-9)/50; B(5) = Real_wp(2)/55;
            
            Bembed.Reallocate(5);
            Bembed(0) = Real_wp(25)/216; Bembed(1) = zero; Bembed(2) = Real_wp(1408)/2565;
            Bembed(3) = Real_wp(2197)/4104; Bembed(4) = -one/5;
          }
        else if (type_algo == DOPRI)
          {
            A.Reallocate(7, 7); A.Fill(zero);
            A(1, 0) = one/5; A(2, 1) = Real_wp(9)/40;
            A(3, 2) = Real_wp(32)/9; A(4, 3) = Real_wp(-212)/729;
            A(5, 4) = Real_wp(-5103)/18656; A(6, 5) = Real_wp(11)/84;
            A(2, 0) = Real_wp(3)/40; A(3, 1) = Real_wp(-56)/15; A(4, 2) = Real_wp(64448)/6561;
            A(5, 3) = Real_wp(49)/176; A(6, 4) = Real_wp(-2187)/6784;
            A(3, 0) = Real_wp(44)/45; A(4, 1) = Real_wp(-25360)/2187;
            A(5, 2) = Real_wp(46732)/5247;
            A(6, 3) = Real_wp(125)/192;
            A(4, 0) = Real_wp(19372)/6561; A(5, 1) = Real_wp(-355)/33; A(6, 2) = Real_wp(500)/1113;
            A(5, 0) = Real_wp(9017)/3168; A(6, 0) = Real_wp(35)/384;
            
            B.Reallocate(7);
            B(0) = Real_wp(5179)/57600; B(1) = zero; B(2) = Real_wp(7571)/16695;
            B(3) = Real_wp(393)/640;
            B(4) = Real_wp(-92097)/339200; B(5) = Real_wp(187)/2100; B(6) = one/40;
            //B(0) = Real_wp(-71)/57600; B(1) = zero; B(2) = Real_wp(71)/16695;
            //B(3) = Real_wp(-71)/1920;
            //B(4) = Real_wp(17253)/339200; B(5) = Real_wp(-22)/525; B(6) = one/40;
            
            Bembed.Reallocate(6);
            Bembed(0) = Real_wp(35)/384; Bembed(1) = zero; Bembed(2) = Real_wp(500)/1113; 
            Bembed(3) = Real_wp(125)/192; Bembed(4) = Real_wp(-2187)/6784;
            Bembed(5) = Real_wp(11)/84;
            
            B = Bembed;
          }
        else
          {
            // coefficients given in Butcher's book (2008, p99)
            A.Reallocate(6, 6); A.Fill(zero); 
            A(1, 0) = one/4; A(2, 1) = one/8;
            A(3, 2) = one/2; A(4, 3) = Real_wp(9)/16; A(5, 4) = Real_wp(8)/7;
            A(2, 0) = one/8; A(4, 2) = Real_wp(3)/8; A(5, 3) = Real_wp(-12)/7;
            A(4, 1) = Real_wp(-3)/8; A(5, 2) = Real_wp(6)/7; A(4, 0) = Real_wp(3)/16;
            A(5, 1) = Real_wp(8)/7; A(5, 0) = Real_wp(-3)/7;
            
            B.Reallocate(6);
            B(0) = Real_wp(7)/90; B(1) = zero; B(2) = Real_wp(16)/45;
            B(3) = two/15; B(4) = Real_wp(16)/45; B(5) = Real_wp(7)/90;
          }
      }
    else if (order == 6)
      {
        // scheme in 2008 Butcher's book, p195
	A.Reallocate(7, 7); A.Fill(zero); 
	A(1, 0) = one/3; A(2, 1) = two/3; A(3, 2) = -one/12; 
	A(4, 3) = Real_wp(15)/8; A(5, 4) = one/10; A(6, 5) = Real_wp(80)/39;
	A(2, 0) = zero; A(3, 1) = one/3; A(4, 2) = Real_wp(35)/48;
	A(5, 3) = one/2; A(6, 4) = Real_wp(32)/195; A(3, 0) = one/12;
	A(4, 1) = Real_wp(-55)/24; A(5, 2) = -one/8; A(6, 3) = Real_wp(-118)/39;
	A(4, 0) = Real_wp(25)/48; A(5, 1) = Real_wp(-11)/24; A(6, 2) = Real_wp(43)/156;
	A(5, 0) = Real_wp(3)/20; A(6, 1) = Real_wp(33)/13; A(6, 0) = Real_wp(-261)/260;
	
	B.Reallocate(7);
	B(0) = Real_wp(13)/200; B(1) = zero; 
	B(2) = Real_wp(11)/40; B(3) = Real_wp(11)/40; B(4) = Real_wp(4)/25; 
	B(5) = Real_wp(4)/25; B(6) = Real_wp(13)/200;
      }
    else if (order == 7)
      {
        // scheme in 2008 Butcher's book, p196
        A.Reallocate(9, 9); A.Fill(zero);
        A(1, 0) = one/6; A(2, 1) = one/3; A(3, 2) = Real_wp(3)/8; A(4, 3) = Real_wp(-56)/1331;
        A(5, 4) = Real_wp(10648)/1701; A(6, 5) = Real_wp(1053)/2401; A(7, 6) = Real_wp(49)/1144;
        A(8, 7) = Real_wp(21)/16;
        A(2, 0) = zero; A(3, 1) = zero; A(4, 2) = Real_wp(150)/1331; A(5, 3) = Real_wp(4024)/1701;
        A(6, 4) = Real_wp(-51909)/16807; A(7, 5) = Real_wp(-405)/2464;
        A(8, 6) = Real_wp(1029)/1408;
        A(3, 0) = one/8; A(4, 1) = zero; A(5, 2) = Real_wp(-170)/27;
        A(6, 3) = Real_wp(-19176)/16807;
        A(7, 4) = Real_wp(-1815)/20384; A(8, 5) = Real_wp(-729)/512;
        A(4, 0) = Real_wp(148)/1331; A(5, 1) = zero; A(6, 2) = Real_wp(1242)/343;
        A(7, 3) = Real_wp(96)/539; A(8, 4) = Real_wp(29403)/3584;
        A(5, 0) = Real_wp(-404)/243; A(6, 1) = zero; A(7, 2) = zero; A(8, 3) = Real_wp(32)/7;
        A(6, 0) = Real_wp(2466)/2401; A(7, 1) = zero; A(8, 2) = Real_wp(-195)/22;
        A(7, 0) = Real_wp(5)/154; A(8, 1) = zero; A(8, 0) = Real_wp(-113)/32;
        
        B.Reallocate(9);
        B(0) = zero; B(1) = zero; B(2) = zero; B(3) = Real_wp(32)/105;
        B(4) = Real_wp(1771561)/6289920; B(5) = Real_wp(243)/2560;
        B(6) = Real_wp(16807)/74880; B(7) = Real_wp(77)/1440; B(8) = Real_wp(11)/270;
      }
    else if (order == 8)
      {
        if (type_algo == FEHLBERG)
          {
            A.Reallocate(13, 13); A.Fill(zero);
            A(1, 0) = Real_wp(2)/27; A(2, 1) = one/12;
            A(3, 2) = one/8; A(4, 3) = Real_wp(25)/16; A(5, 4) = one/5;
            A(6, 5) = Real_wp(125)/54; A(7, 6) = Real_wp(13)/900; A(8, 7) = Real_wp(3);
            A(9, 8) = -one/12; A(10, 9) = Real_wp(18)/41; A(12, 11) = one;
            A(2, 0) = one/36; A(4, 2) =  Real_wp(-25)/16; A(5, 3) = one/4;
            A(6, 4) = Real_wp(-65)/27; A(7, 5) = Real_wp(-2)/9;
            A(8, 6) = Real_wp(67)/90; A(9, 7) = Real_wp(17)/6;
            A(10, 8) = Real_wp(45)/164; A(11, 9) = Real_wp(6)/41;
            A(3, 0) = one/24; A(6, 3) = Real_wp(125)/108;
            A(7, 4) = Real_wp(61)/225; A(8, 5) = Real_wp(-107)/9;
            A(9, 6) = Real_wp(-19)/60; A(10, 7) = Real_wp(45)/82;
            A(11, 8) = Real_wp(3)/41; A(12, 9) = Real_wp(12)/41;
            A(4, 0) = Real_wp(5)/12; A(8, 4) = Real_wp(704)/45; A(9, 5) = Real_wp(311)/54;
            A(10, 6) = Real_wp(2133)/4100; A(11, 7) = Real_wp(-3)/41; A(12, 8) = Real_wp(33)/164;
            A(5, 0) = one/20; A(8, 3) = Real_wp(-53)/6; A(9, 4) = Real_wp(-976)/135;
            A(10, 5) = Real_wp(-301)/82;
            A(11, 6) = Real_wp(-3)/205; A(12, 7) = Real_wp(51)/82;
            A(6, 0) = Real_wp(-25)/108; A(9, 3) = Real_wp(23)/108; A(10, 4) = Real_wp(4496)/1025;
            A(11, 5) = Real_wp(-6)/41; A(12, 6) = Real_wp(2193)/4100;
            A(7, 0) = Real_wp(31)/300; A(10, 3) = Real_wp(-341)/164; A(12, 5) = Real_wp(-289)/82;
            A(8, 0) = Real_wp(2); A(12, 4) = Real_wp(4496)/1025;
            A(9, 0) = Real_wp(-91)/108; A(12, 3) = Real_wp(-341)/164;
            A(10, 0) = Real_wp(2383)/4100; A(11, 0) = Real_wp(3)/205;
            A(12, 0) = Real_wp(-1777)/4100;
            
            B.Reallocate(13); B.Fill(zero);
            B(5) = Real_wp(34)/105; B(6) = Real_wp(9)/35; B(7) = Real_wp(9)/35;
            B(8) = Real_wp(9)/280; B(9) = B(8);
            B(11) = Real_wp(41)/840; B(12) = B(11);
            
            Bembed.Reallocate(11); Bembed.Fill(zero);
            Bembed(0) = B(12); Bembed(5) = B(5); Bembed(6) = B(6); Bembed(7) = B(7);
            Bembed(8) = B(8);
            Bembed(9) = B(9); Bembed(10) = B(12);
          }
        else if (type_algo == DOPRI)
          {
            Real_wp sqrt6 = sqrt(Real_wp(6));
            C.Reallocate(12);
            C(0) = zero; C(6) = one/4; C(7) = Real_wp(4)/13; C(9) = Real_wp(3)/5;
            C(10) = Real_wp(6)/7; C(11) = one;
            C(5) = Real_wp(4)/3*C(6); C(4) = (Real_wp(6)+sqrt6)/10*C(5);
            C(3) = (Real_wp(6)-sqrt6)/10*C(5);
            C(2) = Real_wp(2)/3*C(3); C(1) = Real_wp(2)/3*C(2);
            Real_wp s1 = C(5) + C(6) + C(7), s2 = C(5)*C(6) + C(5)*C(7) + C(6)*C(7),
              s3 = C(5)*C(6)*C(7);
            Real_wp denom = (Real_wp(6) - 21*s1 + 35*s2 - 42*s3 + 21*s1*s1 + 98*s2*s2 
                             + 735*s3*s3 - 84*s1*s2 + 168*s1*s3 - 490*s2*s3);
            Real_wp numer = 3*s1 - 28*s2 + 189*s3 + 14*s1*s2 - 168*s1*s3 + 98*s2*s3;
            C(8) = numer / (2*denom);
            
            Matrix<Real_wp> VDM(8, 8); Real_wp ci_q, ci; VectReal_wp rhs(8);
            for (int i = 0; i < 8; i++)
              {
                ci = C(i);
                if (i > 0)
                  ci = C(i+4);
                
                ci_q = one;
                for (int q = 0; q < 8; q++)
                  {
                    VDM(q, i) = ci_q;
                    ci_q *= ci;
                  }
                
                rhs(i) = one/(i+1);
              }
            
            IVect pivot_VDM(8);
            GetLU(VDM, pivot_VDM);
            SolveLU(VDM, pivot_VDM, rhs);
            
            B.Reallocate(12);
            B(0) = rhs(0); B(1) = zero; B(2) = zero; B(3) = zero; B(4) = zero;
            for (int i = 1; i < 8; i++)
              B(i+4) = rhs(i);

            VectReal_wp C2(12), C3(12), C4(12), C5(12);
            for (int i = 0; i < 12; i++)
              {
                C2(i) = C(i)*C(i); C3(i) = C2(i)*C(i);
                C4(i) = C(i)*C3(i); C5(i) = C4(i)*C(i);
              }
            
            A.Reallocate(12, 12); A.Fill(zero);
            A(1, 0) = C(1);
            A(2, 1) = C2(2)/(C(1)*2);
            A(2, 0) = C(2) - A(2, 1);
            A(3, 2) = C2(3)/(C(2)*2); 
            A(3, 0) = C(3) - A(3, 2);
            
            TinyMatrix<Real_wp, General, 2, 2> vdm2; TinyVector<Real_wp, 2> sol2, source2;
            vdm2(0, 0) = C(2); vdm2(0, 1) = C(3); vdm2(1, 0) = C2(2); vdm2(1, 1) = C2(3);
            source2(0) = C2(4)/2; source2(1) = C3(4)/3; GetInverse(vdm2);
            Mlt(vdm2, source2, sol2);
            A(4, 2) = sol2(0); A(4, 3) = sol2(1); A(4, 0) = C(4) - A(4, 2) - A(4, 3);
            
            vdm2(0, 0) = C(3); vdm2(0, 1) = C(4); vdm2(1, 0) = C2(3); vdm2(1, 1) = C2(4);
            source2(0) = C2(5)/2; source2(1) = C3(5)/3; GetInverse(vdm2);
            Mlt(vdm2, source2, sol2);
            A(5, 3) = sol2(0); A(5, 4) = sol2(1); A(5, 0) = C(5) - A(5, 3) - A(5, 4);
            
            TinyMatrix<Real_wp, General, 3, 3> vdm3; TinyVector<Real_wp, 3> sol3, source3;
            vdm3(0, 0) = C(3); vdm3(0, 1) = C(4); vdm3(0, 2) = C(5);
            vdm3(1, 0) = C2(3); vdm3(1, 1) = C2(4); vdm3(1, 2) = C2(5);
            vdm3(2, 0) = C3(3); vdm3(2, 1) = C3(4); vdm3(2, 2) = C3(5);
            source3(0) = C2(6)/2; source3(1) = C3(6)/3; source3(2) = C4(6)/4;
            GetInverse(vdm3); Mlt(vdm3, source3, sol3);
            A(6, 3) = sol3(0); A(6, 4) = sol3(1); A(6, 5) = sol3(2);
            A(6, 0) = C(6) - A(6, 3) - A(6, 4) - A(6, 5);
            
            TinyMatrix<Real_wp, General, 4, 4> vdm4; TinyVector<Real_wp, 4> sol4, source4;
            vdm4(0, 0) = C(3); vdm4(0, 1) = C(4); vdm4(0, 2) = C(5); vdm4(0, 3) = C(6);
            vdm4(1, 0) = C2(3); vdm4(1, 1) = C2(4); vdm4(1, 2) = C2(5); vdm4(1, 3) = C2(6);
            vdm4(2, 0) = C3(3); vdm4(2, 1) = C3(4); vdm4(2, 2) = C3(5); vdm4(2, 3) = C3(6);
            vdm4(3, 0) = C4(3); vdm4(3, 1) = C4(4); vdm4(3, 2) = C4(5); vdm4(3, 3) = C4(6);
            source4(0) = C2(7)/2; source4(1) = C3(7)/3; source4(2) = C4(7)/4; source4(3) = C5(7)/5;
            GetInverse(vdm4); Mlt(vdm4, source4, sol4);
            A(7, 3) = sol4(0); A(7, 4) = sol4(1); A(7, 5) = sol4(2); A(7, 6) = sol4(3);
            A(7, 0) = C(7) - A(7, 3) - A(7, 4) - A(7, 5) - A(7, 6);
            
            A(11, 10) = B(10)*(one - C(10))/B(11);
            VectReal_wp e(12);
            e.Fill(zero); e(10) = B(11)*A(11, 10) - B(10)*(one - C2(10))/2;
            
            VDM.Reallocate(6, 6); rhs.Reallocate(6);
            for (int i = 0; i < 6; i++)
              {
                ci = C(i);
                if (i > 0)
                  ci = C(i+4);
                
                ci_q = one;
                for (int q = 0; q < 6; q++)
                  {
                    VDM(q, i) = ci_q;
                    ci_q *= ci;
                  }                
              }
            
            ci_q = one; ci = C(10);
            for (int q = 0; q < 6; q++)
              {
                rhs(q) = -e(10)*ci_q;
                ci_q *= ci;
              }
            
            pivot_VDM.Reallocate(6);
            GetLU(VDM, pivot_VDM);
            SolveLU(VDM, pivot_VDM, rhs);            
            e(0) = rhs(0); e(5) = rhs(1); e(6) = rhs(2); e(7) = rhs(3);
            e(8) = rhs(4); e(9) = rhs(5);
            
            vdm2(0, 0) = B(10)*C(10); vdm2(0, 1) = B(11)*C(11);
            source2(0) = e(9) + B(9)/2*(one-C2(9));
            vdm2(1, 0) = B(10); vdm2(1, 1) = B(11); source2(1) = B(9)*(one-C(9));
            GetInverse(vdm2); Mlt(vdm2, source2, sol2);
            A(10, 9) = sol2(0); A(11, 9) = sol2(1);
            
            rhs.Reallocate(23); rhs.Fill(zero);
            rhs(0) = -3.36089; rhs(1) = -8.682193e-1; rhs(2) = 2.7592e1;
            rhs(3) = 2.0154e1; rhs(4) = -4.348988e1;
            rhs(5) = -2.48811; rhs(6) = -5.902908e-1; rhs(7) = 2.123005e1; rhs(8) = 1.52792e1;
            rhs(9) = -3.328821e1; rhs(10) = -2.0331202e-2;
            rhs(11) = 5.18637; rhs(12) = 1.091437; rhs(13) = -8.149787; rhs(14) = -1.8520065e1;
            rhs(15) = 2.273948e1;  rhs(16) = 2.4936055;
            rhs(17) = -1.053449e1; rhs(18) = -2.00087205; rhs(19) = -1.795893e1;
            rhs(20) = 2.794888e1; rhs(21) = -2.85899; rhs(22) = -8.87285;

            NewtonSolver<Real_wp> newton_solver;
            newton_solver.SetAbsoluteSchemeNorm(Real_wp(1));
            newton_solver.ForceReevaluationJacobian(false);
            newton_solver.EnableBacktracking(false);
            newton_solver.SetPrintLevel(0);
            newton_solver.SetStoppingCriterion(1e3*epsilon_machine, 1e5*epsilon_machine);
            newton_solver.SetMaxNumberOfIterations(40);
            
            newton_solver.Init(*this, rhs);
            newton_solver.Solve(*this, rhs);
            
            int n = 0;
            for (int i = 8; i < 12; i++)
              {
                int jend = min(i, 9);
                for (int j = 3; j < jend; j++)
                  A(i, j) = rhs(n++);
                
                A(i, 0) = C(i);
                for (int j = 3; j < i; j++)
                  A(i, 0) -= A(i, j);
              }
          }
        else
          {
            // Cooper-Verner coefficients
            A.Reallocate(11, 11); A.Fill(zero);
            Real_wp sqrt21 = sqrt(Real_wp(21));
            A(1, 0) = one/2; A(2, 1) = one/4; A(3, 2) = (Real_wp(21)+5*sqrt21)/49;
            A(4, 3) = (Real_wp(21)-sqrt21)/252;
            A(5, 4) = (Real_wp(63)-7*sqrt21)/80; A(6, 5) = (Real_wp(63)-13*sqrt21)/35;
            A(7, 6) = one/9;
            A(8, 7) = (Real_wp(63)+13*sqrt21)/128; A(9, 8) = (Real_wp(132)+28*sqrt21)/245;
            A(10, 9) = (Real_wp(49)-7*sqrt21)/18;
            A(2, 0) = one/4; A(3, 1) = (Real_wp(-7)-3*sqrt21)/98;
            A(4, 2) = (Real_wp(18)+4*sqrt21)/63;
            A(5, 3) = (Real_wp(-231)+14*sqrt21)/360; A(6, 4) = (Real_wp(-504)+115*sqrt21)/70;
            A(7, 5) = (Real_wp(13)-3*sqrt21)/63;
            A(8, 6) = (Real_wp(-385)-75*sqrt21)/1152; A(9, 7) = (Real_wp(-51)-11*sqrt21)/56;
            A(10, 8) = (Real_wp(28)-28*sqrt21)/45;
            A(3, 0) = one/7; A(4, 1) = zero; A(5, 2) = (Real_wp(9)+sqrt21)/36;
            A(6, 3) = (Real_wp(633)-145*sqrt21)/90;
            A(7, 4) = (Real_wp(14)-3*sqrt21)/126; A(8, 5) = Real_wp(11)/72;
            A(9, 6) = (Real_wp(515)+111*sqrt21)/504;
            A(10, 7) = (Real_wp(301)+53*sqrt21)/72;
            A(4, 0) = (Real_wp(11)+sqrt21)/84; A(5, 1) = zero;
            A(6, 2) = (Real_wp(-432)+92*sqrt21)/315;
            A(7, 3) = zero; A(8, 4) = (Real_wp(91)-21*sqrt21)/576;
            A(9, 5) = (Real_wp(-733)-147*sqrt21)/2205;
            A(10, 6) = (Real_wp(-273)-53*sqrt21)/72;
            A(5, 0) = (Real_wp(5)+sqrt21)/48; A(9, 4) = one/9;
            A(10, 5) = (Real_wp(-18)+28*sqrt21)/45;
            A(6, 0) = (Real_wp(10)-sqrt21)/42; A(10, 4) = (Real_wp(-42)+7*sqrt21)/18;
            A(7, 0) = one/14; A(8, 0) = one/32; A(9, 0) = one/14;

        
            B.Reallocate(11); B.Fill(zero);
            B(0) = one/20; B(7) = Real_wp(49)/180; B(8) = Real_wp(16)/45;
            B(9) = B(7); B(10) = one/20;
          }
      }
    else
      {
	cout << "This order is not implemented" << endl;
	abort();
      }
    
    C.Reallocate(A.GetM()); C.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
	C(i) += A(i,j);
    
    /*DISP(B); DISP(C);
      Real_wp sum0 = zero, sum1 = zero;
      for (int i = 0; i < B.GetM(); i++)
      {
      sum0 += B(i);
      sum1 += B(i)*C(i);
      }
      DISP(sum0); DISP(sum1); */
  }

  
  //! initialisation of the scheme with u0
  template<class T>
  void RungeKutta_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    Y = u0; Y_tilde.Reallocate(Y.GetM()); Y_tilde.Fill(0);
    u0.Clear(); ki.Reallocate(B.GetM());
    for (int i = 0; i < B.GetM(); i++)
      {
	ki(i).Reallocate(Y.GetM());
	ki(i).Fill(0);
      }
    
    sys.SetDirichletCondition(t, 0, Y);
  }
  

  //! memory used for intermediary vectors is released 
  template<class T>
  void RungeKutta_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory used for intermediary vectors is released 
  template<class T>
  void RungeKutta_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_tilde.Clear();
    ki.Clear();
    Yfinal = Y; Y.Clear();
  }
  
  
  //! returns the stability function of RK scheme
  template<class T>
  UnivariatePolynomial<Real_wp> RungeKutta_Iterator<T>::GetStabilityFunction() const
  {
    // PolK = polynomial associated with each ki
    Vector<UnivariatePolynomial<Real_wp> > PolK(B.GetM());
    UnivariatePolynomial<Real_wp> z, one, P;
    one.SetOrder(0); one(0) = 1.0;
    z.SetOrder(1); z(0) = 0.0; z(1) = 1.0;
    
    P = one;
    for (int i = 0; i < B.GetM(); i++)
      {
        PolK(i) = one;
        for (int j = 0; j < i; j++)
          Add(A(i, j), PolK(j), PolK(i));
        
        PolK(i) *= z;
        Add(B(i), PolK(i), P);
      }
    
    return P;
  }


  //! returns the CFL of the scheme
  template<class T>
  Real_wp RungeKutta_Iterator<T>::GetCFL(bool real_p) const
  {
    switch (order_rk)
      {
      case 3 :
	return sqrt(Real_wp(3));
      case 4 :
	// cas ou la partie imaginaire domine
	if (!real_p)
	  return 2.0*sqrt(Real_wp(2));
	else
	  // cas ou la partie reelle domaine
	  return Real_wp(2.785);
      }
    
    return Real_wp(0);
  }
  
    
  //! Runge-Kutta algorithm to compute y^n+1
  template<class T>
  void RungeKutta_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    Real_wp zero(0);
    
    for (int i = 0; i < B.GetM(); i++)
      {
        glob_chrono.Start(VirtualTimer::SCHEME);
	// we have an explicit scheme
	Copy(Y, Y_tilde);
	
	// we compute Y_tilde = y_n + h * \sum_{s=1}^R  a_{r,s} k_s
	for (int j = 0; j < i; j++)
	  if (A(i, j) != zero)
	    Add(dt*A(i,j), ki(j), Y_tilde);
	    
	Real_wp tn_tilde = t + C(i)*dt;
        glob_chrono.Stop(VirtualTimer::SCHEME);
        
	sys.EvaluateFunction(tn_tilde, Y_tilde, ki(i));
	// k_i = f (x_n + c_i * h, y_n + h * \sum_{s=1}^i  a_{i,s} k_s )
      }
    
    glob_chrono.Start(VirtualTimer::SCHEME);
    // forming y^n+1 = y^n + dt \sum b_i k_i
    for (int i = 0; i < B.GetM(); i++)
      Add(dt*B(i), ki(i), Y);

    glob_chrono.Stop(VirtualTimer::SCHEME);
    
  }
  
  
  //! changes time step during the simulation
  template<class T>
  void RungeKutta_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& RungeKutta_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& RungeKutta_Iterator<T>::GetIterate() const
  {
    return Y;
  }
  

  //! returns the size used by the object in bytes
  template<class T>
  size_t RungeKutta_Iterator<T>::GetMemorySize() const
  {
    size_t taille = B.GetMemorySize() + C.GetMemorySize() + A.GetMemorySize() + Bembed.GetMemorySize();
    taille += Y.GetMemorySize() + Y_tilde.GetMemorySize() + Seldon::GetMemorySize(ki);
    return taille;
  }
  

  //! fills the coefficients a_ij, b_i and c_i of the scheme
  template<class T>
  void RungeKutta_Iterator<T>::GetCoefficients(Matrix<Real_wp>& a, VectReal_wp& b, VectReal_wp& c) const
  {
    a = A; b = B; c = C; 
  }
  

  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int RungeKutta_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& RungeKutta_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }

  //! computes the non-linear system to solve in order to satisfy eighth-order scheme
  template<class T>
  void RungeKutta_Iterator<T>::ComputeScheme(const VectReal_wp& X, VectReal_wp& F)
  {
    int n = 0;
    for (int i = 8; i < 12; i++)
      {
        int jend = min(i, 9);
        for (int j = 3; j < jend; j++)
          A(i, j) = X(n++);
      }
    
    VectReal_wp C2(C), C3(C), C4(C), C5(C);
    for (int i = 0; i < C.GetM(); i++)
      {
        C2(i) = C(i)*C(i); C3(i) = C(i)*C2(i);
        C4(i) = C(i)*C3(i); C5(i) = C(i)*C4(i);
      }
    
    int num = 0;
    for (int i = 8; i < 12; i++)
      {
        // equation \sum a_{i, j} c_j = c_i^2/2
        F(num) = C2(i)/2;
        for (int j = 3; j < i; j++)
          F(num) -= A(i, j)*C(j);
        
        num++;

        // equation \sum a_{i, j} c_j^2 = c_i^3/3
        F(num) = C3(i)/3;
        for (int j = 3; j < i; j++)
          F(num) -= A(i, j)*C2(j);
        
        num++;
        
        // equation \sum a_{i, j} c_j^3 = c_i^4/4
        F(num) = C4(i)/4;
        for (int j = 3; j < i; j++)
          F(num) -= A(i, j)*C3(j);
        
        num++;
        
        if (i > 8)
          {
            // equation \sum a_{i, j} c_j^2 = c_i^5/5
            F(num) = C5(i)/5;
            for (int j = 3; j < i; j++)
              F(num) -= A(i, j)*C4(j);
            
            num++;                
          }
      }
    
    Real_wp one(1), zero(0);
    for (int j = 3; j < 5; j++)
      {
        // equation \sum b_i a_{i, j} = b_j (1 - c_j)
        F(num) = B(j)*(one - C(j));
        for (int i = j+1; i < 12; i++)
          F(num) -= B(i)*A(i, j);
        
        num++;
        
        // equation \sum b_i c_i a_{i, j} = 0
        F(num) = zero;
        for (int i = j+1; i < 12; i++)
          F(num) += B(i)*C(i)*A(i, j);
        
        num++;

        // equation \sum b_i c_i^2 a_{i, j} = 0
        F(num) = zero;
        for (int i = j+1; i < 12; i++)
          F(num) += B(i)*C2(i)*A(i, j);
        
        num++;

        // equation \sum b_i c_i \sum_k a_{i, k} a_{k, j} = 0
        F(num) = zero;
        for (int i = j+2; i < 12; i++)
          for (int k = j+1; k < i; k++)
            F(num) += B(i)*C(i)*A(i, k)*A(k, j);
        
        num++;
      }
  }
  
  
  //! returns the 2-norm of scheme
  template<class T>
  Real_wp RungeKutta_Iterator<T>::GetNorm2Vector(const VectReal_wp& scheme)
  {
    return Norm2(scheme);
  }
  
  
  //! computes and factorises the jacobian matrix associated the non-linear system to solve
  template<class T>
  void RungeKutta_Iterator<T>::ComputeAndFactoriseDiff(const VectReal_wp& X, const VectReal_wp& scale)
  {
    EvaluateJacobian(X, mat_jacobian_lu);  
    GetLU(mat_jacobian_lu, pivot_jacobian);
  }

  
  //! computes the non-linear system to solve in order to satisfy eighth-order conditions
  template<class T>
  void RungeKutta_Iterator<T>::EvaluateFunction(const VectReal_wp& X, VectReal_wp& F)
  {
    ComputeScheme(X, F);
  }

  
  //! computes the jacobian matrix associated the non-linear system to solve
  template<class T>
  void RungeKutta_Iterator<T>::EvaluateJacobian(const VectReal_wp& X, Matrix<Real_wp>& DF)
  {
    Real_wp one(1), two(2), three(3);
    Real_wp h = pow(epsilon_machine, one/three);
    int N = X.GetM();
    VectReal_wp f_plus(N), f_minus(N), Xtmp(X);
    DF.Reallocate(N, N);
    for (int j = 0; j < N; j++)
      {
        Xtmp(j) += h;
        ComputeScheme(Xtmp, f_plus);
        
        Xtmp(j) -= two*h;
        ComputeScheme(Xtmp, f_minus);
        
        for (int i = 0; i < N; i++)
          DF(i, j) = (f_plus(i) - f_minus(i))/(two*h);
        
        Xtmp(j) += h;
      }
  }


  //! solves the jacobian matrix associated the non-linear system to solve (X  = DF^{-1} F )
  template<class T>  
  void RungeKutta_Iterator<T>::SolveDifferential(const VectReal_wp& F, VectReal_wp& X)
  {
    Copy(F, X);
    SolveLU(mat_jacobian_lu, pivot_jacobian, X);
  }


  /**************************************
   *LinearRungeKutta
   ***************************************/
  

  //! Default constructor
  template<class T>
  LinearRungeKutta_Iterator<T>::LinearRungeKutta_Iterator()
  {
    SetOrder(2);
  }
  

  //! Compute ci and mat_wi use to evaluate source function for Linear Runge Kutta Scheme
  template<class T>
  void LinearRungeKutta_Iterator<T>::ComputeCoefForLinearRungeKuttaSourceFunction()
  {
    // int nb_terms_add = order+extraS;
    int nb_terms_add = order;
    if (stable_source)
      nb_terms_add = Pol.GetOrder() + 1;
    
    int nb_terms = nb_terms_add;
   
    // calcul de 1/k! et 1 / 2^k k!
    VectReal_wp invFacto(2*nb_terms_add+1), CoefC(2*nb_terms_add+1), CoefAlpha;
    invFacto(0) = Real_wp(1); CoefC(0) = Real_wp(1);
    for (int k = 1; k <= 2*nb_terms_add; k++)
      {
	invFacto(k) = invFacto(k-1) / Real_wp(k);
	CoefC(k) = CoefC(k-1) / Real_wp(2*k);
      }
    
    if (!stable_algo)
      {
	// basis functions associated with Gauss points
	Globatto<Real_wp> Q;
	Q.ConstructQuadrature(order-2, Globatto<Real_wp>::QUADRATURE_GAUSS);
	int p = Q.GetNbPointsQuad();	
	
	// derivatives of these basis functions
	Q.ComputeGradPhi();

	// number of terms in numerator
	int orderN = this->Pol.GetOrder();

	CoefAlpha.Reallocate(orderN);
	for (int j = 0; j < this->order; j++)
	  CoefAlpha(j) = invFacto(j+1);

	for (int j = this->order; j < orderN; j++)
	  CoefAlpha(j) = this->Pol(j+1);

	this->ci = Q.Points();
	this->mat_wi.Reallocate(orderN, p);
	this->mat_wi.Zero();
	for (int i = 0; i < p; i++)
	  this->mat_wi(0, i) = Q.Weights(i);

	VectReal_wp ValPhi0(p), qtilde(p);
	Q.ComputeValuesPhiRef(Real_wp(0), ValPhi0);
	
	const Matrix<Real_wp>& dPhi = Q.GradPhi();
	Matrix<Real_wp> Decomp(dPhi), Interm(dPhi);
	Decomp.SetIdentity();
	DerMat.Reallocate(orderN-1, p);
	for (int l = 1; l < orderN; l++)
	  {
	    Mlt(Decomp, ValPhi0, qtilde);
	    for (int i = 0; i < p; i++)
	      {
		Real_wp dphi_tilde = qtilde(i);
		DerMat(l-1, i) = dphi_tilde;
		for (int k = 1; k < orderN; k++)
		  if (l <= orderN-k)
		    this->mat_wi(k, i) += dphi_tilde*CoefAlpha(k+l-1);
	      }
	    
	    Mlt(Decomp, dPhi, Interm);
	    Decomp = Interm;
	  }

	return;
      }
    
    // terms with X^{0} (to see the error constants for the homogeneous case)
    VectReal_wp err(nb_terms_add);
    err.Zero();
    
    for (int r = 0; r < nb_terms_add; r++)
      {
	Real_wp sum = CoefC(r);
	for (int i = 0; i <= min(r, Pol.GetOrder()); i++)
	  {
	    Real_wp signK(1);
	    if ((r-i)%2 == 1)
	      signK = -signK;
	    
	    sum -= signK*Pol(i)*CoefC(r-i);
	  }
	
	err(r) = sum;
      }
    
    //cout << "Errors for the homogeneous case" << endl << err << endl;
    
    // computation of coefficients alpha_j^r for the right hand side phi
    // phi = dt \sum_{r=1}^\infty (dt A)^{r-1} \sum_{j=1}^\infty \alpha_j^r dt^{j-1} F^{j-1}
    // where \alpha_j^r = \sum_{i=0}^{min(r-1, m)} (D_i - (-1)^k N_i) / 2^k k!
    // where k = r+j-i-1
    //cout << endl << "Coefficients alpha for inhomogeneous case" << endl;
    Matrix<Real_wp> A(nb_terms_add, nb_terms_add);
    A.Zero();
    if (!stable_source)
      for (int r = 1; r <= nb_terms_add; r++)
        {
          err.Zero();
          for (int j = 1; j <= nb_terms_add+1-r; j++)
            {
              // partie denominateur
              err(j-1) = CoefC(r+j-1);
              
              // numerateur
              for (int i = 0; i <= min(r-1, Pol.GetOrder()); i++)
                {
                  int k = r + j-i - 1;
                  Real_wp signK(1);
                  if (k%2 == 1)
                    signK = -signK;
                  
                  err(j-1) -= signK*Pol(i)*CoefC(k);
                }
              
              A(r-1, j-1) = err(j-1);
            }
          
          //cout << "Terms with (dt A)^" << r-1 << endl << err << endl;
        }
    
    Matrix<Real_wp> VDM; VectReal_wp rhs(nb_terms-1);
    Vector<int> pivot;
    
    // calcul des ci : on prend les points de Gauss-Legendre
    VectReal_wp weights; 
    if (!stable_source)
      {
        ComputeGaussLegendre(this->ci, weights, nb_terms-2);
        
        // Vandermonde Matrix for ci
        // les factorielles sont mises dans le second membre    
        VDM.Reallocate(this->ci.GetM(), this->ci.GetM());  
        for (int i = 0; i < this->ci.GetM(); i++)
          for (int j = 0; j < this->ci.GetM(); j++)
            VDM(i, j) = pow(this->ci(j) - Real_wp(0.5), i);
        
        GetLU(VDM, pivot);        
      }
    else
      {
        ComputeGaussLegendre(this->ci, weights, order/2);
        for (int i = 0; i < this->ci.GetM(); i++)
          for (int r = 0; r < nb_terms_add; r++)
            A(r, i) = weights(i)*Pol(r)*pow(1.0-this->ci(i), r);
        
        //DISP(this->ci); DISP(weights);
        //A.WriteText("A.dat");
      }
    
    int orderN = Pol.GetOrder();
    
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
    
    if (stable_source)
      {
        orderN++; nb++;
      }
    
    // polynomials needed to expand the source are enumerated
    Vector<UnivariatePolynomial<Real_wp> > pol_source(orderN);
    for (int i = 0; i < orderN; i++)
      pol_source(i) = PolNumer(nb-2-i);
    
    // inverse of coefficients, to express x^i as a combination of pol_source
    Matrix<Real_wp> coefBase;
    coefBase.Reallocate(orderN, orderN);
    coefBase.Zero();
    for (int j = 0; j < orderN; j++)
      for (int i = 0; i < pol_source(j).GetM(); i++)
	coefBase(i, j) = pol_source(j)(i);
    
    GetInverse(coefBase);    
    polA.Reallocate(orderN, nb_terms-1);
    polA.Zero();
    for (int powA = 0; powA < nb_terms; powA++)
      {
	rhs.Zero();	
        
	for (int j = 1; j <= nb_terms-1; j++)
	  rhs(j-1) = A(powA, j-1);
	
	//DISP(powA); DISP(rhs);
	for (int j = 0; j <= nb_terms-2; j++)
	  for (int k = 0; k < orderN; k++)
	    polA(k, j) += coefBase(k, powA)*rhs(j);
      }

    this->mat_wi.Reallocate(orderN, this->ci.GetM());

    // loop over rows of mat_wi
    for (int powA = 0; powA < orderN; powA++)
      {
	rhs.Zero();	
        
	if (stable_source)
          for (int j = 1; j <= polA.GetN(); j++)
            rhs(j-1) = polA(powA, j-1);
        else
          {
            for (int j = 1; j <= polA.GetN(); j++)
              rhs(j-1) = polA(powA, j-1) / invFacto(j-1);
            
            SolveLU(VDM, pivot, rhs);
          }
        
        SetRow(rhs, powA, this->mat_wi);
	//DISP(powA); DISP(rhs);
      }
  }


  //! initializes a linear Runge-Kutta scheme with a given order and additional stages
  template<class T>
  void LinearRungeKutta_Iterator<T>::SetOrder(int order, int extraS, bool s, bool s2)
  {
    stable_algo = s;
    stable_source = s2;
    this->order = order;
    this->extraS = extraS;
    Pol.SetOrder(this->order + this->extraS);
    Real_wp one(1);

    
    switch (order)
      {
      case 2 :
	{
	  switch(extraS)
	    {
	    case 0: { 
	      Pol(0) =   one;
	      Pol(1) =   one;
	      Pol(2) =   4.999999999999977e-01 ;
	    }
	      break;
	    case 1: { 
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000001e-01 ; 
	      Pol(3) =   1.451277982649155e-01 ;
	    }
	      break;
	    case 2: {
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.665532314108146e-01 ; 
	      Pol(4) =   2.327815361933148e-02 ;
	    }
	      break;
	    case 3: {
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   9.999999999999999e-01 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.618342913053687e-01 ; 
	      Pol(4) =   3.289792611743811e-02 ; 
	      Pol(5) =   2.839528016518102e-03 ;
	    }
	      break;
	    case 4: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   9.999999999999999e-01 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.642981320398038e-01 ; 
	      Pol(4) =   3.657769285804588e-02 ; 
	      Pol(5) =   5.035250867609586e-03 ; 
	      Pol(6) =   3.001880509358407e-04 ;
	    }
	      break;
	    case 5: {
	      Pol(0) =   9.999999999999990e-01 ; 
	      Pol(1) =   1.000000000001165e+00 ; 
	      Pol(2) =   4.999999999988338e-01 ; 
	      Pol(3) =   1.626462249413356e-01 ; 
	      Pol(4) =   3.762678272315501e-02 ; 
	      Pol(5) =   5.996644250417070e-03 ; 
	      Pol(6) =   5.826143210213330e-04 ; 
	      Pol(7) =   2.487327304531716e-05 ;
	    }
	      break;
	    case 6: {
	      Pol(0) =   1.000000000000062e+00 ; 
	      Pol(1) =   9.999999999999879e-01 ; 
	      Pol(2) =   5.000000000000121e-01 ; 
	      Pol(3) =   1.627509585676844e-01 ; 
	      Pol(4) =   3.773348832445807e-02 ; 
	      Pol(5) =   6.387803046851333e-03 ; 
	      Pol(6) =   7.489561665296774e-04 ; 
	      Pol(7) =   5.356270766078865e-05 ; 
	      Pol(8) =   1.713109940102836e-06 ;
	    }
	      break;
	    case 7: { 
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   1.000000000000150e+00 ; 
	      Pol(2) =   4.999999999998272e-01 ; 
	      Pol(3) =   1.640094942014296e-01 ; 
	      Pol(4) =   3.840429977823329e-02 ; 
	      Pol(5) =   6.724597512047917e-03 ; 
	      Pol(6) =   8.718626803227696e-04 ; 
	      Pol(7) =   7.857554562878064e-05 ; 
	      Pol(8) =   4.327975378833797e-06 ; 
	      Pol(9) =   1.072985856243921e-07 ;
	    }
	      break;
	    case 8: {
	      Pol(0) =   1.000000000000001e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999976e-01 ; 
	      Pol(3) =   1.649990588856614e-01 ; 
	      Pol(4) =   3.927394350377206e-02 ; 
	      Pol(5) =   7.055384479248899e-03 ; 
	      Pol(6) =   9.695797812914759e-04 ; 
	      Pol(7) =   9.943224646288322e-05 ; 
	      Pol(8) =   7.129812259258231e-06 ; 
	      Pol(9) =   3.148056880771953e-07 ; 
	      Pol(10) =   6.324920988294407e-09 ;
	    }
	      break;
	    default :
	      {
		cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_RUNGE_KUTTA number_additional_stages" << endl;
		abort();
	      }
	    }
	}
	break;
      case 4 :
	{
	  switch(extraS)
            {
	    case 0: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666667e-02 ;
	    }
	      break;
            case 1: {
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000002e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666667e-02 ; 
	      Pol(5) =   4.730163010446185e-03 ;
	    }
              break;
            case 2: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999997e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666666e-02 ; 
	      Pol(5) =   6.541349497416528e-03 ; 
	      Pol(6) =   4.395282130923843e-04 ;
	    }
              break;
            case 3: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666666e-02 ; 
	      Pol(5) =   7.241999849787970e-03 ; 
	      Pol(6) =   7.614940065988191e-04 ; 
	      Pol(7) =   3.521874589831831e-05 ;
	    }
              break;
            case 4: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000000e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666667e-02 ; 
	      Pol(5) =   7.603292194142675e-03 ; 
	      Pol(6) =   9.535828377031919e-04 ; 
	      Pol(7) =   7.298469178025099e-05 ; 
	      Pol(8) =   2.500124976522895e-06 ;
	    }
              break;
            case 5: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000109e+00 ; 
	      Pol(2) =   4.999999999998921e-01 ; 
	      Pol(3) =   1.666666666666659e-01 ; 
	      Pol(4) =   4.166666666666743e-02 ; 
	      Pol(5) =   7.817918289656257e-03 ; 
	      Pol(6) =   1.075759999127459e-03 ; 
	      Pol(7) =   1.026588721744709e-04 ; 
	      Pol(8) =   6.038353896295552e-06 ; 
	      Pol(9) =   1.628169027707504e-07 ;
	    }
              break;
            case 6: {
	      Pol(0) =   1.000000000000008e+00 ; 
	      Pol(1) =   9.999999999997634e-01 ; 
	      Pol(2) =   5.000000000001779e-01 ; 
	      Pol(3) =   1.666666666670609e-01 ; 
	      Pol(4) =   4.166666666597196e-02 ; 
	      Pol(5) =   7.992535147077134e-03 ; 
	      Pol(6) =   1.180030987873825e-03 ; 
	      Pol(7) =   1.307878349087823e-04 ; 
	      Pol(8) =   1.020785594818226e-05 ; 
	      Pol(9) =   4.943966219870204e-07 ; 
	      Pol(10) =   1.097077616437946e-08 ;
	    }
              break;
            case 7: {
	      Pol(0) =   9.999999999999916e-01 ; 
	      Pol(1) =   1.000000000006203e+00 ; 
	      Pol(2) =   4.999999999938012e-01 ; 
	      Pol(3) =   1.666666666666492e-01 ; 
	      Pol(4) =   4.166666666668406e-02 ; 
	      Pol(5) =   9.619397138072583e-03 ; 
	      Pol(6) =   3.970757223041604e-03 ; 
	      Pol(7) =   1.979923031733034e-03 ; 
	      Pol(8) =   6.726632799312973e-04 ; 
	      Pol(9) =   1.385778310637994e-04 ; 
	      Pol(10) =   1.585824201586086e-05 ; 
	      Pol(11) =   7.742514686545619e-07 ;
	    }
              break;
            case 8: {
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   1.000000000000021e+00 ; 
	      Pol(2) =   4.999999999999787e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666669e-02 ; 
	      Pol(5) =   8.105487675563905e-03 ; 
	      Pol(6) =   1.249316412377197e-03 ; 
	      Pol(7) =   1.531845812394507e-04 ; 
	      Pol(8) =   1.473468121845849e-05 ; 
	      Pol(9) =   1.071860716775002e-06 ; 
	      Pol(10) =   5.510748021396615e-08 ; 
	      Pol(11) =   1.766727504578043e-09 ; 
	      Pol(12) =   2.623218531216638e-11 ;
	    }
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
	        cout << "OrderTimeScheme = order LINEAR_RUNGE_KUTTA number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 6 :
	{
	  switch(extraS)
            {
	      
            case 0: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000026e-01 ; 
	      Pol(3) =   1.666666666666798e-01 ; 
	      Pol(4) =   4.166666666669645e-02 ; 
	      Pol(5) =   8.333333333363002e-03 ; 
	      Pol(6) =   1.388888888894670e-03 ;
	    }
              break;
	    case 1 : {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   9.999999999999999e-01 ; 
	      Pol(2) =   5.000000000000001e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666669e-02 ; 
	      Pol(5) =   8.333333333333335e-03 ; 
	      Pol(6) =   1.388888888888889e-03 ; 
	      Pol(7) =   2.070461615593214e-04 ;
	    }
	      break;
            case 2: {
	      Pol(0) =   9.999999999999999e-01 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999996e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666669e-02 ; 
	      Pol(5) =   8.333333333333338e-03 ; 
	      Pol(6) =   1.388888888888888e-03 ; 
	      Pol(7) =   2.204061707466545e-04 ; 
	      Pol(8) =   1.942982735313673e-05 ;
	    } 
              break;
            case 3: {
	      Pol(0) =   9.999999999999997e-01 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999996e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666669e-02 ; 
	      Pol(5) =   8.333333333333331e-03 ; 
	      Pol(6) =   1.388888888888889e-03 ; 
	      Pol(7) =   2.073919102492977e-04 ; 
	      Pol(8) =   2.499262304459253e-05 ; 
	      Pol(9) =   1.453234258464881e-06 ;
	    } 
              break;
            case 4: {
	      Pol(0) =   9.999999999998850e-01 ; 
	      Pol(1) =   1.000000000014410e+00 ; 
	      Pol(2) =   4.999999999857789e-01 ; 
	      Pol(3) =   1.666666666583180e-01 ; 
	      Pol(4) =   4.166666667495167e-02 ; 
	      Pol(5) =   8.333333331832365e-03 ; 
	      Pol(6) =   1.388888890344787e-03 ; 
	      Pol(7) =   2.358338644436141e-04 ; 
	      Pol(8) =   4.056334413908446e-05 ; 
	      Pol(9) =   4.775871882059528e-06 ; 
	      Pol(10) =   2.442645091656458e-07 ;
	    } 
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
	        cout << "OrderTimeScheme = order LINEAR_RUNGE_KUTTA number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 8:
	{
	  switch(extraS)
            {
            case 0: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   9.999999999999998e-01 ; 
	      Pol(2) =   5.000000000000001e-01 ; 
	      Pol(3) =   1.666666666666666e-01 ; 
	      Pol(4) =   4.166666666666666e-02 ; 
	      Pol(5) =   8.333333333333328e-03 ; 
	      Pol(6) =   1.388888888888889e-03 ; 
	      Pol(7) =   1.984126984126984e-04 ; 
	      Pol(8) =   2.480158730158730e-05 ; 
	    }
              break;
	    case 1: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   5.000000000000002e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666674e-02 ; 
	      Pol(5) =   8.333333333333345e-03 ; 
	      Pol(6) =   1.388888888888891e-03 ; 
	      Pol(7) =   1.984126984126987e-04 ; 
	      Pol(8) =   2.480158730158732e-05 ; 
	      Pol(9) =   1.684112035592431e-06 ;
	    } 
              break;
	      
            case 2: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999998e-01 ; 
	      Pol(3) =   1.666666666666667e-01 ; 
	      Pol(4) =   4.166666666666669e-02 ; 
	      Pol(5) =   8.333333333333331e-03 ; 
	      Pol(6) =   1.388888888888887e-03 ; 
	      Pol(7) =   1.984126984126984e-04 ; 
	      Pol(8) =   2.480158730158725e-05 ; 
	      Pol(9) =   2.288709306973234e-06 ; 
	      Pol(10) =   9.960040692054680e-08 ;
	    } 
              break;
	      
            case 3: {
	      Pol(0) =   9.999999999997605e-01 ; 
	      Pol(1) =   1.000000000005352e+00 ; 
	      Pol(2) =   4.999999999949840e-01 ; 
	      Pol(3) =   1.666666666564997e-01 ; 
	      Pol(4) =   4.166666667634694e-02 ; 
	      Pol(5) =   8.333333344117122e-03 ; 
	      Pol(6) =   1.388888874173291e-03 ; 
	      Pol(7) =   1.984126719728717e-04 ; 
	      Pol(8) =   2.480162654363079e-05 ; 
	      Pol(9) =   2.528206540248994e-06 ; 
	      Pol(10) =   1.724423811134767e-07 ; 
	      Pol(11) =   5.449535772542617e-09 ;
	    } 
              break;
	      
            case 4: {
	      Pol(0) =   1.000000000000000e+00 ; 
	      Pol(1) =   1.000000000000000e+00 ; 
	      Pol(2) =   4.999999999999993e-01 ; 
	      Pol(3) =   1.666666666666665e-01 ; 
	      Pol(4) =   4.166666666666674e-02 ; 
	      Pol(5) =   8.333333333333345e-03 ; 
	      Pol(6) =   1.388888888888887e-03 ; 
	      Pol(7) =   1.984126984126986e-04 ; 
	      Pol(8) =   2.480158730158578e-05 ; 
	      Pol(9) =   2.638893313733145e-06 ; 
	      Pol(10) =   2.150620166601062e-07 ; 
	      Pol(11) =   1.123553506837818e-08 ; 
	      Pol(12) =   2.690758844819519e-10 ;
	    } 
              break;
	      
            case 5: {
	      Pol(0) =   9.999999999999556e-01 ; 
	      Pol(1) =   1.000000000001086e+00 ; 
	      Pol(2) =   4.999999999989788e-01 ; 
	      Pol(3) =   1.666666666645971e-01 ; 
	      Pol(4) =   4.166666666863206e-02 ; 
	      Pol(5) =   8.333333335529922e-03 ; 
	      Pol(6) =   1.388888886864611e-03 ; 
	      Pol(7) =   1.984126972643316e-04 ; 
	      Pol(8) =   2.480158854116792e-05 ; 
	      Pol(9) =   2.703333893632985e-06 ; 
	      Pol(10) =   2.435581983430564e-07 ; 
	      Pol(11) =   1.631043038503232e-08 ; 
	      Pol(12) =   6.905312067380033e-10 ; 
	      Pol(13) =   1.342332862257654e-11 ;
	    } 
              break;
	      
            case 6: {
	      Pol(0) =   9.999999999955160e-01 ; 
	      Pol(1) =   1.000000000043304e+00 ; 
	      Pol(2) =   4.999999999645489e-01 ; 
	      Pol(3) =   1.666666665722565e-01 ; 
	      Pol(4) =   4.166666674705911e-02 ; 
	      Pol(5) =   8.333333484341877e-03 ; 
	      Pol(6) =   1.388888765299902e-03 ; 
	      Pol(7) =   1.984124663131307e-04 ; 
	      Pol(8) =   2.480179708118387e-05 ; 
	      Pol(9) =   2.711246141311401e-06 ; 
	      Pol(10) =   2.500568374959440e-07 ; 
	      Pol(11) =   1.817647917892119e-08 ; 
	      Pol(12) =   9.481642471601341e-10 ; 
	      Pol(13) =   3.089127728872379e-11 ; 
	      Pol(14) =   4.655664953646905e-13 ;
	    } 
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_RUNGE_KUTTA number_additional_stages" << endl;
                abort();
              }
            }
	}
	break;
      case 10:
	{
	  switch(extraS)
            {
            case 0: {
	      Pol(0) = one; Pol(1) = one;
	      Pol(2) = Real_wp(0.500000000000001);
	      Pol(3) = Real_wp(0.166666666666670); 
	      Pol(4) = Real_wp(0.041666666666672);
	      Pol(5) = Real_wp(0.008333333333337);
	      Pol(6) = Real_wp(0.001388888888890);
	      Pol(7) = Real_wp(0.000198412698413);
	      Pol(8) = Real_wp(0.000024801587302);
	      Pol(9) = Real_wp(0.000002755731922);
	      Pol(10) = Real_wp(0.000000275573192);
	    } 
              break;
            default :
              {
                cout << "Order " << order << " with " << extraS << " additional stage(s) is not implemented" <<endl;
		cout << "Syntax:" << endl;
		cout << "OrderTimeScheme = order LINEAR_RUNGE_KUTTA number_additional_stages" << endl;
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
    
    // pour les premiers coefs, on met les valeurs exactes pour avoir le bon ordre
    // en precision quelconque
    Pol(0) = one;
    for(int k = 0; k < order; k++)
      Pol(k+1) = Pol(k)/(k+1);
    
    // numerator is factorized for a "stable" algorithm
    // ie we list the real and complex roots of the numerators
    VectReal_wp R, Ri;
    SolvePolynomialEquation(Pol, R, Ri);
    
    real_roots.Clear(); complex_roots.Clear();
    int n = 0;
    while (n < Pol.GetM()-1)
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
            n +=2;
          }
      }
    //DISP(Pol);
    //DISP(R); DISP(Ri);
    //DISP(real_roots); 
    //DISP(complex_roots);

    // source term is treated
    ComputeCoefForLinearRungeKuttaSourceFunction();
  }


  //! returns the stability function R(z)
  //! (the scheme can be written as y_{n+1} = R(\Delta t A) y^n for the homogeneous case)
  template<class T>
  UnivariatePolynomial<Real_wp> LinearRungeKutta_Iterator<T>::GetStabilityFunction() const
  {
    return Pol;
  }


  //! returns quadrature points c_i
  template<class T>
  const VectReal_wp& LinearRungeKutta_Iterator<T>::GetPoints()
  {
    return this->ci;
  }


  //! returns weights w_{i, k}
  template<class T>
  const Matrix<Real_wp>& LinearRungeKutta_Iterator<T>::GetWeights()
  {
    return this->mat_wi;
  }


  //! returns the matrix of derivatives for phi
  template<class T>
  const Matrix<Real_wp>& LinearRungeKutta_Iterator<T>::GetDerivativeMatrix()
  {
    return this->DerMat;
  }

  
  //! initialisation of the scheme with u0
  template<class T>
  void LinearRungeKutta_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    Y = u0;
    int N = Y.GetM();
    u0.Clear(); 
    Fn.Reallocate(N);
    KhUn.Reallocate(N); KhUn.Fill(0);
    
    EvalF.Reallocate(ci.GetM());
    
    for(int i = 0; i < EvalF.GetM(); i++)
      {
        EvalF(i).Reallocate(N);
        EvalF(i).Zero();
      }
    
    sys.SetDirichletCondition(t, 0, Y);
  }
  

  //! memory used for intermediary vectors is released 
  template<class T>
  void LinearRungeKutta_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory used for intermediary vectors is released 
  template<class T>
  void LinearRungeKutta_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Fn.Clear(); KhUn.Clear(); EvalF.Clear();
    Yfinal = Y; Y.Clear();
  }

  
  //! returns the CFL of the scheme
  template<class T>
  Real_wp LinearRungeKutta_Iterator<T>::GetCFL(bool real_p) const
  {
    int coef = 1;
    if (real_p)
      coef = 2;
    
    switch (order)
      {
      case 2 :
	switch(extraS)
	  {
	  case 0 :
	    return Real_wp(coef*9.454078674316406e-03);
	  case 1 : 
	    return Real_wp(coef*1.379053344726562);
	  case 2 : 
	    return Real_wp(coef*2.251435546874999);
	  case 3 : 
	    return Real_wp(coef*2.909065246582031);
	  case 4 :
	    return Real_wp(coef*3.581784667968750);
	  case 5 : 
	    return Real_wp(coef*4.265054626464844);
	  case 6 :
	    return Real_wp(coef*4.922929687499999);
	  case 7 :
	    return Real_wp(coef*5.639323425292968);
	  case 8 : 
	    return Real_wp(coef*6.311920166015623);
	  }
	break;
      case 4 :
	switch(extraS)
	  {
	  case 0 :
	    return Real_wp(coef*1.391689453125000);
	  case 1 :
	    return Real_wp(coef*2.482749938964843);
	  case 2 :
	    return Real_wp(coef*3.107010498046874);
	  case 3 :
	    return Real_wp(coef*3.961481323242187);
	  case 4 :
	    return Real_wp(coef*4.577460937499999);
	  case 5 :
	    return Real_wp(coef*5.132529602050781);
	  case 6 : 
	    return Real_wp(coef*5.744522094726562);
	  case 7 :
	    return Real_wp(coef*6.468430023193359);
	  case 8 : 
	    return Real_wp (coef*7.145903320312499);
	  }
	break;
      case 6 :
	switch(extraS)
	  {
	  case 0 : 
	    return Real_wp(coef*0.2215152740478515);
	  case 1 : 
	    return Real_wp(coef*1.946175384521484);
	  case 2 : 
	    return Real_wp(coef*2.761787109375000);
	  case 3 : 
	    return Real_wp(coef*3.555009613037109);
	  case 4 : 
	    return Real_wp(coef*4.161941528320312);
	  }
	break;
      case 8 : 
	switch(extraS)
	  {
	  case 0 : 
	    return Real_wp(coef*2.155253906250000);
	  case 1 : 
	    return Real_wp(coef*3.274285583496093);
	  case 2 : 
	    return Real_wp(coef*4.011657714843749);
	  case 3 : 
	    return Real_wp(coef*4.653706970214843);
	  case 4 : 
	    return Real_wp(coef*4.518083496093749);
	  case 5 :
	    return Real_wp(coef*6.033312377929686);
	  case 6 :
	    return Real_wp(coef*6.510294799804687);
	  }
	break;
      case 10 :
	return Real_wp(coef*0.7751884460449218);
	break;
      }

    return Real_wp(0);
  }
  
  
  //! Runge-Kutta algorithm to compute y^n+1
  template<class T>
  void LinearRungeKutta_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
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
        
        // cas ou on multiplie au fur et a mesure
        int num_source = 0;
        if (stable_source)
          {
            if (presence_source)
              {
                KhUn.Zero();
                for (int i = 0; i < mat_wi.GetN(); i++)
                  KhUn += dt*mat_wi(num_source, i)*EvalF(i);

                sys.SolveMassMatrix(KhUn);
                Y += KhUn;
              }

            num_source++;
          }
        
	for (int i = 0; i < complex_roots.GetM(); i++)
          {
            Real_wp b = -Real_wp(2)*dt*realpart(Real_wp(1) / complex_roots(i));
            Real_wp a = dt*dt/absSquare(complex_roots(i));
	    sys.EvaluateFunction(t, Y, Fn, false, false);
            Fn *= a;
	    
            KhUn.Zero();
            if (presence_source)
              {
		for (int i = 0; i < mat_wi.GetN(); i++)
		  KhUn += dt*dt*mat_wi(num_source, i)*EvalF(i);
              }
            
            Fn += KhUn;
	    // multiplication by M_h^-1
            sys.SolveMassMatrix(Fn);
            KhUn = b*Y + Fn;
            sys.EvaluateFunction(t, KhUn, Fn, false, false);
            
            KhUn.Zero();
            if (presence_source)
              {
		for (int i = 0; i < mat_wi.GetN(); i++)
		  KhUn += dt*mat_wi(num_source+1, i)*EvalF(i);
              }

            Fn += KhUn;
            sys.SolveMassMatrix(Fn);
	    KhUn = Y + Fn;
            
	    Y = KhUn;
            num_source += 2;
          }
        
        for (int i = 0; i < real_roots.GetM(); i++)
          {
            Real_wp b = -dt / real_roots(i);
	    sys.EvaluateFunction(t, Y, Fn, false, false);
            KhUn.Zero();
            if (presence_source)
              {
		for (int i = 0; i < mat_wi.GetN(); i++)
		  KhUn += dt*mat_wi(num_source, i)*EvalF(i);
              }
            
	    Fn *= b;
            Fn += KhUn;
            sys.SolveMassMatrix(Fn);
	    KhUn = Y + Fn;
	    Y = KhUn;
	    num_source++;
          }

        return;
      }

    Fn.Zero();
    Copy(Y, KhUn); 
    
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
    
    // We evaluate U_{n+1} = Pol*U^n
    int last_coef = Pol.GetM()-1;
    
    Fn = this->Pol(last_coef)*Y;
    for(int k = last_coef -1; k >= 0; k--)
      {
	sys.EvaluateFunction(t, Fn, KhUn, false, false);
	if ((k < this->mat_wi.GetM()) && presence_source)
	  for (int i = 0; i < this->mat_wi.GetN(); i++)
	    KhUn += this->mat_wi(k, i)*EvalF(i);

	sys.SolveMassMatrix(KhUn);
	Fn = dt*KhUn + Pol(k)*Y;	
      }
    
    // updating U^n to the next value
    Y = Fn;
  }
  
  
  //! changes time step during the simulation
  template<class T>
  void LinearRungeKutta_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& LinearRungeKutta_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& LinearRungeKutta_Iterator<T>::GetIterate() const
  {
    return Y;
  }  

  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int LinearRungeKutta_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& LinearRungeKutta_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t LinearRungeKutta_Iterator<T>::GetMemorySize() const
  {
    size_t taille = real_roots.GetMemorySize() + Pol.GetMemorySize() + complex_roots.GetMemorySize();
    taille += ci.GetMemorySize() + mat_wi.GetMemorySize() + polA.GetMemorySize() + DerMat.GetMemorySize();
    taille += Y.GetMemorySize() + Fn.GetMemorySize() + KhUn.GetMemorySize() + Seldon::GetMemorySize(EvalF);
    return taille;
  }
  
  
  /********************
   * Talezer_Iterator *
   ********************/
  
  
  //! default constructor
  template<class T>
  Talezer_Iterator<T>::Talezer_Iterator() : lob()
  {
    dt = 0;
    order_integration_bessel = 10;
    Rmax = 10.0;
  }

  
  //! computation of coefficients needed before using Tal-ezer schemes
  /*!
    \param[in] order number of bessel functions
    \param[in] deltat time step
  */
  template<class T>
  void Talezer_Iterator<T>::SetOrder(int order, const Real_wp& deltat)
  {
    coef.Reallocate(order+1);
    ComputeBesselFunctions(0, order, deltat*Rmax, coef);
    for (int i = 1; i < coef.GetM(); i++)
      coef(i) *= 2.0;
    
    int r = order_integration_bessel;
    lob.ConstructQuadrature(r, lob.QUADRATURE_GAUSS);
    coef_bessel.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      {
	Real_wp tau = lob.Points(i)*Rmax*deltat;
	ComputeBesselFunctions(0, order, tau, coef_bessel(i));
	Mlt(lob.Weights(i)*deltat, coef_bessel(i));
	for (int j = 1; j < coef_bessel(i).GetM(); j++)
	  coef_bessel(i)(j) *= 2.0; 
      }
    
    coef_bk.Reallocate(order + 1);
  }


  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> Talezer_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }
  

  template<class T>
  void Talezer_Iterator<T>::SetSpectralRadius(const Real_wp& R)
  {
    Rmax = R;
  }
  
  
  template<class T>
  void Talezer_Iterator<T>::SetOrderIntegration(int r)
  {
    order_integration_bessel = r; 
  }
  
  
  //! integration of time source for Tal-ezer scheme
  /*!
    \param[in] t0 integration over [t0, t0+deltat]
    \param[in] deltat time step
    \param[in] sys given problem
  */
  template<class T>
  void Talezer_Iterator<T>::
  ComputeCoefficientsSource(const Real_wp& t0, const Real_wp& deltat, VirtualOdeSystem<T>& sys)
  {
    int r = order_integration_bessel;

    coef_bk.Fill(0);
    T pulse; SetComplexZero(pulse);
    for (int j = 0; j <= r; j++)
      {
	Real_wp t = (1.0-lob.Points(j))*deltat + t0;
	sys.SourceOnlyTime(t, 0, pulse);
	// DISP(t); DISP(pulse(0));
	for (int i = 0; i < coef.GetM(); i++)
	  coef_bk(i) += coef_bessel(j)(i)*pulse;
	
      }
  }

  
  //! initialisation of the scheme with U^0
  template<class T>
  void Talezer_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM(); T zero(0);
    Y = u0; u0.Clear();
    By.Reallocate(N); Qn.Reallocate(N); Qnp1.Reallocate(N);
    Bsource.Reallocate(N); Bn.Reallocate(N); Bnp1.Reallocate(N);
    By.Fill(zero); Qn.Fill(zero); Qnp1.Fill(zero);
    Bsource.Fill(zero); Bn.Fill(zero); Bnp1.Fill(zero);
    sys.SetDirichletCondition(t, 0, Y);
  }
  

  //! intermediary vectors are cleared
  template<class T>
  void Talezer_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! intermediary vectors are cleared
  template<class T>
  void Talezer_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    By.Clear(); Qn.Clear(); Qnp1.Clear();
    Bsource.Clear(); Bn.Clear(); Bnp1.Clear();
    Yfinal = Y; Y.Clear(); 
  }
  
  
  //! Tal-ezer algorithm
  /*!
    \param[in] t_begin initial time
    \param[in] t_end final time
    \param[in] deltat time step
    \param[in,out] Y iterate
    \param[in,out] sys given problem
  */
  template<class T>
  void Talezer_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    Real_wp coef_gn = 1.0/Rmax, one(1);
    
    // first iterates
    Copy(Y, Qn);
    sys.EvaluateFunction(t, Y, Qnp1, true, false); Mlt(coef_gn, Qnp1);
    Mlt(coef(0), Y);
    Add(coef(1), Qnp1, Y);
	
    if (t < sys.GetFinalTimeSource())
      {
	ComputeCoefficientsSource(t, dt, sys);
	sys.FillSource(Bn);
	
	sys.SolveOperatorDh(Bn);
	
	// DISP(coef_bk);
	sys.EvaluateFunction(t, Bn, Bnp1, true, false); Mlt(coef_gn, Bnp1);
	Bsource = Bn; Mlt(coef_bk(0), Bsource);
	Add(coef_bk(1), Bnp1, Bsource);
      }
	
    for (int i = 2; i <= coef.GetM()-1; i++)
      {
	sys.EvaluateFunction(t, Qnp1, By, true, false);
	Mlt(2.0*coef_gn, By); Add(one, Qn, By);
	Copy(Qnp1, Qn);
	Copy(By, Qnp1);
	    
	Add(coef(i), Qnp1, Y);
	
	if (t < sys.GetFinalTimeSource())
	  {
	    sys.EvaluateFunction(t, Bnp1, By, true, false);
	    Mlt(2.0*coef_gn, By); Add(one, Bn, By);
	    Bn = Bnp1;
	    Bnp1 = By;
	    
	    Add(coef_bk(i), Bnp1, Bsource);
	  }
      }
	
    if (t < sys.GetFinalTimeSource())
      Add(one, Bsource, Y);
    
  }


  //! changes time step during the simulation
  template<class T>
  void Talezer_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented " << endl;
    abort();
  }

    
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& Talezer_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& Talezer_Iterator<T>::GetIterate() const
  {
    return Y;
  }

  
  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int Talezer_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& Talezer_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t Talezer_Iterator<T>::GetMemorySize() const
  {
    size_t taille = coef.GetMemorySize() + Seldon::GetMemorySize(coef_bessel) + coef_bk.GetMemorySize();
    taille += lob.GetMemorySize();
    taille += Y.GetMemorySize() + By.GetMemorySize() + Qn.GetMemorySize() + Qnp1.GetMemorySize()
      + Bsource.GetMemorySize() + Bn.GetMemorySize() + Bnp1.GetMemorySize();
    
    return taille;
  }


  /***********************************
   * AdamsBashforth_Moulton_Iterator *
   ***********************************/
  
  
  //! default constructor
  template<class T>
  AdamsBashforth_Moulton_Iterator<T>::AdamsBashforth_Moulton_Iterator()
  {
    dt = 0; type_algo = BASHFORTH;
  }
  

  //! computation of coefficients
  /*!
    \param[in] order order of Adams-Bashfort-Moulton scheme
  */
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>::SetOrder(int order, int type)
  {
    type_algo = type;
    
    // evaluation of coefficients gamma and gamma_star
    VectReal_wp gamma, gamma_star;
    gamma.Reallocate(order);
    gamma_star.Reallocate(order);
    gamma(0) = 1.0; gamma_star(0) = 1.0;
    for (int m = 1; m < order; m++)
      {
        gamma(m) = 1.0;
        for (int j = 0; j < m; j++)
          gamma(m) -= gamma(j)/(m+1-j);

        gamma_star(m) = 0.0;
        for (int j = 0; j < m; j++)
          gamma_star(m) -= gamma_star(j)/(m+1-j);
      }
    
    // computation of beta and beta_star from gamma, gamma_star
    beta.Reallocate(order); beta_star.Reallocate(order);
    beta(0) = gamma_star(order-1); beta_star(0) = gamma(order-1);
    VectReal_wp new_beta(order), new_beta_star(order);
    for (int m = order-2; m >= 0; m--)
      {
        new_beta(0) = gamma_star(m) + beta(0);
        for (int j = 1; j <= order-2-m; j++)
          new_beta(j) = -beta(j-1) + beta(j);
        
        new_beta(order-1-m) = -beta(order-2-m);
        beta = new_beta;

        new_beta_star(0) = gamma(m) + beta_star(0);
        for (int j = 1; j <= order-2-m; j++)
          new_beta_star(j) = -beta_star(j-1) + beta_star(j);
        
        new_beta_star(order-1-m) = -beta_star(order-2-m);
        beta_star = new_beta_star;
      }
  }
  
  
  //! initialisation of the scheme with u^0
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM(); int order = beta.GetM();
    Y = u0; Y_star.Reallocate(N); Y_star.Fill(0);
    sys.SetDirichletCondition(t, 0, Y);
    F_star.Reallocate(N); F_star.Fill(0);
    fi.Reallocate(order);
    for (int i = 0; i < order; i++)
      {
	fi(i).Reallocate(N);
	fi(i).Fill(0);
      }    
    
    if (order >= 7)
      {
        rk_scheme.SetOrder(6);
        rk_scheme.SetInitialCondition(t, deltat, u0, sys);
      }

    u0.Clear();
  }
  

  //! memory is released
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory is released
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_star.Clear();
    fi.Clear(); F_star.Clear();
    Yfinal = Y; Y.Clear(); 
    rk_scheme.Clear();
  }
  

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> AdamsBashforth_Moulton_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

    
  //! Adams-Bashforth-Moulton algorithm
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>
  ::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    int order = beta.GetM();
    
    if (n_time < order)
      {
        if (order < 7)
          {
            // use of RK4 scheme for the first iterates
            Copy(Y, Y_star);
            sys.EvaluateFunction(t, Y, fi(0));
            Add(dt/6, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt/2, fi(0), F_star);
            sys.EvaluateFunction(t+0.5*dt, F_star, fi(0));
            Add(dt/3, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt/2, fi(0), F_star);
            sys.EvaluateFunction(t+0.5*dt, F_star, fi(0));
            Add(dt/3, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt, fi(0), F_star);
            sys.EvaluateFunction(t+dt, F_star, fi(0));
            Add(dt/6, fi(0), Y_star);
            Copy(Y_star, Y);
          }
        else
          {
            rk_scheme.Advance(t, n_time, sys);
            Y = rk_scheme.GetIterate();
          }
        
	// once Y_star = y^n+1 is computed
	// f(y^n+1) is stored in the vector fi
	sys.EvaluateFunction(t+dt, Y, fi(order-1-n_time));
	
	return;
      }
    
    if (type_algo == BASHFORTH)
      {
	for (int i = 0; i < order; i++)
	  Add(dt*beta_star(i), fi(i), Y);
	
	// udpate of f(t_{n-i},y_{n-i})
	for (int i = (order-1); i > 0; i--)
	  Copy(fi(i-1), fi(i));
	
	sys.EvaluateFunction(t+dt, Y, fi(0));
      }
    else if (type_algo == MOULTON)
      {
        cout << "Scheme not implemented here" << endl;
        cout << "Use rather AdamsImplicit_Iterator" << endl;
        abort();
      }
    else
      {
	// Adams-Bashforth step
	// y*_n = y_n + h \sum_{i=1}^k  beta_star(i) fi(i)
	Copy(Y, Y_star);
	for (int i = 0; i < order; i++)
	  Add(dt*beta_star(i), fi(i), Y_star);
	
	// computation of f(t_n+1, y*_n)
	sys.EvaluateFunction(t+dt, Y_star, F_star);
	
	// Adams-Moulton step
	// y_{n+1} = y_n + h \beta_0 f(t_n+1, y*_n) + h \sum_{i=0}^k  beta(i) fi(i)
	Add(dt*beta(0), F_star, Y);
	for (int i = 0; i < order-1; i++)
	  Add(dt*beta(i+1), fi(i), Y);
	
	// udpate of f(t_{n-i},y_{n-i})
	for (int i = (order-1); i > 0; i--)
	  Copy(fi(i-1), fi(i));
	
	sys.EvaluateFunction(t+dt, Y, fi(0));
      }
  }
  

  //! changes time step during the simulation
  template<class T>
  void AdamsBashforth_Moulton_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented " << endl;
    abort();
  }

    
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& AdamsBashforth_Moulton_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& AdamsBashforth_Moulton_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int AdamsBashforth_Moulton_Iterator<T>::GetNumberOfIterates() const
  {
    return beta.GetM()+1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& AdamsBashforth_Moulton_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else
      return fi(k-1);
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t AdamsBashforth_Moulton_Iterator<T>::GetMemorySize() const
  {
    size_t taille = beta_star.GetMemorySize() + Seldon::GetMemorySize(fi) + beta.GetMemorySize();
    taille += Y.GetMemorySize() + Y_star.GetMemorySize() + F_star.GetMemorySize() + rk_scheme.GetMemorySize();
    
    return taille;
  }

  
  /*************************
   * TaylorSeries_Iterator *
   *************************/

  
  //! default constructor
  template<class T>
  TaylorSeries_Iterator<T>::TaylorSeries_Iterator()
  {
    order = 1; dt = 0;
  }


  //! constructor specifying the order
  template<class T>
  TaylorSeries_Iterator<T>::TaylorSeries_Iterator(int r)
  {
    order = r; dt = 0;    
  }
  

  template<class T>
  void TaylorSeries_Iterator<T>::SetOrder(int r)
  {
    order = r;
  }
  
  
  //! initialisation of the scheme with u^0
  template<class T>
  void TaylorSeries_Iterator<T>::
  SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
		      Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();
    Y_tilde.Reallocate(N); Fy.Reallocate(N);
    Y_tilde.Fill(0); Fy.Fill(0);
    sys.SetDirichletCondition(t, 0, Y);
    
    dtk_div_factork.Reallocate(order); dtk_div_factork(0) = dt;
    for (int i = 1; i < order; i++)
      dtk_div_factork(i) = dtk_div_factork(i-1)*dt/(i+1);
  }
  

  //! releases memory used by the scheme
  template<class T>
  void TaylorSeries_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! releases memory used by the scheme
  template<class T>
  void TaylorSeries_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_tilde.Clear(); Fy.Clear();
    Yfinal = Y; Y.Clear(); 
  }
  

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> TaylorSeries_Iterator<T>::GetStabilityFunction() const
  {
    UnivariatePolynomial<Real_wp> Pol;
    Pol.SetOrder(order);
    Pol(0) = Real_wp(1); 
    for (int i = 0; i < order; i++)
      Pol(i+1) = Pol(i) / (i+1);

    return Pol;
  }

    
  //! "Taylor-series" algorithm
  template<class T>
  void TaylorSeries_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {    
    Copy(Y, Y_tilde);
    for (int i = 0; i < order; i++)
      {
	sys.EvaluateDerivativeFunction(t, i, Y_tilde, Fy);
	Add(dtk_div_factork(i), Fy, Y); Copy(Fy, Y_tilde);
      }    
  }


  //! changes time step during the simulation
  template<class T>
  void TaylorSeries_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }
  
  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& TaylorSeries_Iterator<T>::GetIterate()
  {
    return Y;
  }
  

  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& TaylorSeries_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int TaylorSeries_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& TaylorSeries_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t TaylorSeries_Iterator<T>::GetMemorySize() const
  {
    size_t taille = dtk_div_factork.GetMemorySize();
    taille += Y.GetMemorySize() + Y_tilde.GetMemorySize() + Fy.GetMemorySize();    
    return taille;
  }

  
  /*************************
   * LowStorageRK_Iterator *
   *************************/
  
  
  //! default constructor
  template<class T>
  LowStorageRK_Iterator<T>::LowStorageRK_Iterator()
  {
    dt = 0;
  }
  
  
  //! computation of RK coefficients
  template<class T>
  void LowStorageRK_Iterator<T>::SetOrder(int order, int nb_steps)
  {
    Real_wp one(1), zero(0);
    if (order == 1)
      {
        coef_alpha.Reallocate(1);
        coef_alpha(0) = zero;
        coef_beta.Reallocate(1);
        coef_beta(0) = one;
        coef_C.Reallocate(1);
        coef_C(0) = zero;
      }
    else if (order == 2)
      {
        coef_alpha.Reallocate(2);
        coef_alpha(0) = zero; coef_alpha(1) = -one/2;
        coef_beta.Reallocate(2);
        coef_beta(0) = one/2; coef_beta(1) = one;
        coef_C.Reallocate(2);
        coef_C(0) = zero; coef_C(1) = one/2;
      }
    else if (order == 3)
      {
        if (nb_steps == 4)
          {
            coef_alpha.Reallocate(4);
            coef_alpha(0) = zero; coef_alpha(1) = -one; coef_alpha(2) = -one; coef_alpha(3) = -one;
            coef_beta.Reallocate(4);
            coef_beta(0) = one/3; coef_beta(1) = Real_wp(3)/4;
            coef_beta(2) = Real_wp(2)/3; coef_beta(3) = one/4;
            coef_C.Reallocate(4);
            coef_C(0) = zero; coef_C(1) = one/3; coef_C(2) = one/3; coef_C(3) = one;
          }
        else
          {
            coef_alpha.Reallocate(3);
            coef_alpha(0) = 0; coef_alpha(1) = -one; coef_alpha(2) = -one;
            coef_beta.Reallocate(3);
            coef_beta(0) = Real_wp(2)/3; coef_beta(1) = Real_wp(3)/4; coef_beta(2) = one/3;
            coef_C.Reallocate(3);
            coef_C(0) = zero; coef_C(1) = Real_wp(2)/3; coef_C(2) = Real_wp(2)/3;
          }
      }
    else if (order == 4)
      {
        // alpha parameter corresponds to the following linear scheme
        // y^(n+1) = (I + dt A + dt^2/2 A^2 + dt^3/6 A^3 + dt^4/24 A^4 + alpha A^5) y^n
        Real_wp alpha(0);
        if (nb_steps == 4)
          {
            // best choice for stability along imaginary axis
            // CFL along imaginary axis is equal to 3.4641, along real axis to 3.548
            alpha = Real_wp(1)/144;
            coef_alpha.Reallocate(5); coef_beta.Reallocate(5); coef_C.Reallocate(5);
            
            coef_alpha(0) = 0; coef_alpha(1) = -0.613085564532188;
            coef_alpha(2) =  -4.80974684512998;
            coef_alpha(3) =  -0.681440108994655; coef_alpha(4) = -0.323868168222252;
            
            coef_beta(0) =  0.257132922435245; coef_beta(1) = 1.19010535737294;
            coef_beta(2) = 0.343756181299216; coef_beta(3) = 0.296230012887199;
            coef_beta(4) = 0.222851128767784;
            
            coef_C(0) = 0; coef_C(1) = coef_beta(0); coef_C(2) = 0.717601864930413;
            coef_C(3) = 0.421641376245036; coef_C(4) = 0.891667466145127;
          }
        else if (nb_steps == 5)
          {
            // best choice for stability along real axis (optimum is 1.0/244.68055944...)
            // CFL along imaginary axis is equal to 3.23, along real axis to 6.036
            alpha = Real_wp(1)/244;
            coef_alpha.Reallocate(5); coef_beta.Reallocate(5); coef_C.Reallocate(5);
            
            coef_alpha(0) = 0; coef_alpha(1) = -0.417524383853296;
            coef_alpha(2) = -0.994542746678958;
            coef_alpha(3) =  -1.75243004657323; coef_alpha(4) = -2.81513792740248;
            
            coef_beta(0) = 0.109526703681785; coef_beta(1) = 0.318922544798638;
            coef_beta(2) = 0.493029456654037; coef_beta(3) = 1.75857621924023;
            coef_beta(4) = 0.135322767802179;
            
            coef_C(0) = 0; coef_C(1) = coef_beta(0); coef_C(2) = 0.295291309466447;
            coef_C(3) = 0.502710330688203; coef_C(4) = 0.964771388109944;
          }
        else
          {
            // intermediary choice taken by Carpenter/Kennedy
            // CFL along imaginary axis is equal to 3.3407, along real axis to 4.656
            alpha = Real_wp(1)/200;
            coef_alpha.Reallocate(5); coef_beta.Reallocate(5); coef_C.Reallocate(5);
            
            coef_alpha(0) = 0; coef_alpha(1) = -4.178904744998519e-01;
            coef_alpha(2) =  -1.192151694642677e+00;
            coef_alpha(3) =  -1.697784692471528e+00; coef_alpha(4) = -1.514183444257156e+00;
            
            coef_beta(0) =  1.496590219992291e-01; coef_beta(1) = 3.792103129996273e-01;
            coef_beta(2) = 8.229550293869817e-01; coef_beta(3) = 6.994504559491221e-01;
            coef_beta(4) = 1.530572479681520e-01;
            
            coef_C(0) = 0; coef_C(1) = coef_beta(0); coef_C(2) = 3.704009573642048e-01;
            coef_C(3) = 6.222557631344432e-01; coef_C(4) = 9.582821306746903e-01;
          }
        
        // line to avoid warning
        VectReal_wp toto(1); toto(0) = alpha;

#ifdef MONTJOIE_WITH_MULTIPLE
        // refining coefficients with Newton's method
        param_scheme.Reallocate(1); param_scheme(0) = alpha; 
        
        VectReal_wp xsol(9);        
        FillGlobalVector(coef_alpha, coef_beta, xsol);
        
        /*         
          bool not_convergence = true;
          Matrix<Real_wp> a(5, 5); VectReal_wp b(5), c(5), fvec(9);
          while (not_convergence)
          {
          FillGlobalVector(coef_alpha, coef_beta, xsol);

          VectReal_wp fvec(xsol);
          Matrix<Real_wp> fjac(xsol.GetM(), xsol.GetM());
          IVect Control(10); VectReal_wp RControl(10);
          Control.Fill(0); RControl.Fill(Real_wp(0));
          Control(0) = 1000;
          Control(3) = 0;
          RControl(0) = 1e-12;
          RControl(1) = Real_wp(1);
            
          int info = SolveMinpack(*this, xsol, fvec, fjac, scale, Control, RControl);
        */  
            
        NewtonSolver<Real_wp> newton_solver;
	newton_solver.SetAbsoluteSchemeNorm(Real_wp(0.01));
	newton_solver.ForceReevaluationJacobian(false);
	newton_solver.EnableBacktracking(false);
	newton_solver.SetPrintLevel(0);
	newton_solver.SetStoppingCriterion(1e4*epsilon_machine, 1e8*epsilon_machine);
	newton_solver.SetMaxNumberOfIterations(40);
        
        newton_solver.Init(*this, xsol);
        newton_solver.Solve(*this, xsol);
        
        /* this->EvaluateFunction(xsol, fvec);
            
           FillVectorAB(xsol, coef_alpha, coef_beta);           
           
           GetCoefficientsC(coef_alpha, coef_beta, c);
           GetCoefficientsAB(coef_alpha, coef_beta, a, b);
           bool cond_not_satisfied = false;
           for (int i = 0; i < c.GetM(); i++)
           if ((c(i) < -1e-12) || (c(i) > 1.0))
           cond_not_satisfied = true;
           
           //for (int i = 0; i < c.GetM(); i++)
           //if (b(i) < -1e-12)
           //cond_not_satisfied = true;
            
           Real_wp residu = Norm2(fvec);
           if (cond_not_satisfied)
           residu = 1.0;
            
           if (residu < 1e-10)
           {
           DISP(Norm2(fvec));
           not_convergence = false;
           }
           else
           {
           cond_not_satisfied = true;
           while (cond_not_satisfied)
           {
           Real_wp lambda;
           for (int i = 0; i < coef_alpha.GetM(); i++)
           {
           lambda = Real_wp(rand())/RAND_MAX;
           coef_alpha(i) = -Real_wp(2)*lambda;

           lambda = Real_wp(rand())/RAND_MAX;
           coef_beta(i) = lambda;
           }
                    
           coef_alpha(0) = 0;
                    
           GetCoefficientsC(coef_alpha, coef_beta, c);
           GetCoefficientsAB(coef_alpha, coef_beta, a, b);
           cond_not_satisfied = false;
           for (int i = 0; i < c.GetM(); i++)
           if ((c(i) < -1e-12) || (c(i) > 1.0))
           cond_not_satisfied = true;
                    
           //for (int i = 0; i < c.GetM(); i++)
           //if (b(i) < -1e-12)
           //  cond_not_satisfied = true;
           }
           }
           }
        */
        
        FillVectorAB(xsol, coef_alpha, coef_beta);
        GetCoefficientsC(coef_alpha, coef_beta, coef_C);
        //DISP(coef_alpha); DISP(coef_beta); DISP(coef_C);
#endif
        
      }
    else
      {
	cout << "Order not implemented" << endl;
	abort();
      }
  }


#ifdef MONTJOIE_WITH_MULTIPLE
  //! Fills vector X with values of A and B
  template<class T> template<class Vector1, class Vector2>
  void LowStorageRK_Iterator<T>::FillGlobalVector(const Vector1& A, const Vector1& B, Vector2& X)
  {
    int N = A.GetM()-1;
    for (int i = 1; i < A.GetM(); i++)
      X(i-1) = A(i);

    for (int i = 0; i < B.GetM(); i++)
      X(N+i) = B(i);        
  }  
  
  
  //! fills values of A and B from the solution X
  template<class T> template<class Vector1, class Vector2>
  void LowStorageRK_Iterator<T>::FillVectorAB(const Vector1& X, Vector2& A, Vector2& B)
  {
    int N = A.GetM()-1;
    A(0) = 0;
    for (int i = 1; i < A.GetM(); i++)
      A(i) = X(i-1);

    for (int i = 0; i < B.GetM(); i++)
      B(i) = X(N+i);
  }

  
  //! computes Runge-Kutta coefficients c_i from values A and B
  //! (denoted alpha, beta in documentation)
  template<class T> template<class Vector1, class Vector2>
  void LowStorageRK_Iterator<T>::GetCoefficientsC(const Vector1& A, const Vector1& B, Vector2& C)
  {
    Real_wp one(1), coef = one+A(1);
    C(0) = 0;
    C(1) = B(0);
    C(2) = C(1) + B(1)*coef;
    for (int i = 3; i < C.GetM(); i++)
      {
        coef = A(i-1)*coef + one;
        C(i) = C(i-1) + B(i-1)*coef;
      }
  }

  
  //! computes Runge-Kutta coefficients a_{i, j}, b_i from values A and B
  template<class T>
  void LowStorageRK_Iterator<T>
  ::GetCoefficientsAB(const VectReal_wp& A, const VectReal_wp& B,
                      Matrix<Real_wp>& a, VectReal_wp& b)
  {
    int N = A.GetM();
    a.Fill(0);
    
    for (int i = 1; i < N; i++)
      {
        a(i, i-1) = B(i-1);
        for (int j = i-2; j >= 0; j--)
          a(i, j) = A(j+1)*a(i, j+1) + B(j);
      }
    
    b(N-1) = B(N-1);
    for (int i = N-2; i >= 0; i--)
      b(i) = A(i+1)*b(i+1) + B(i);
  }

  
  //! computes the non-linear system to solve in order to satisfy fourth-order scheme
  template<class T>
  void LowStorageRK_Iterator<T>::ComputeScheme(const VectReal_wp& X, VectReal_wp& F)
  {
    int N = 5;
    VectReal_wp A(N), B(N), c(N);
    Matrix<Real_wp> a(N, N); VectReal_wp b(N);
    FillVectorAB(X, A, B);
    GetCoefficientsC(A, B, c);
    GetCoefficientsAB(A, B, a, b);
    
    Real_wp one (1);
    // 9 equations : \sum b_i = 1, \sum b_i c_i = 1/2, \sum b_i c_i^2 = 1/3, \sum b_i c_i^3 = 1/4
    // \sum b_i a_i,j c_j = 1/6, \sum b_i c_i a_i,j c_j = 1/8, \sum b_i a_i,j c_j = 1/12
    // \sum b_i a_i,j a_j,k c_k = 1/24, \sum b_i a_i,j a_j,k a_k,l c_l = alpha 
    F(0) = one; F(1) = one/2; F(2) = one/3; F(3) = one/4; F(4) = one/6;
    F(5) = one/8; F(6) = one/12; F(7) = one/24; F(8) = param_scheme(0);
    for (int i = 0; i < b.GetM(); i++)
      {
        F(0) -= b(i);
        F(1) -= b(i)*c(i);
        F(2) -= b(i)*square(c(i));
        F(3) -= b(i)*square(c(i))*c(i);
        for (int j = 0; j < i; j++)
          {
            F(4) -= b(i)*a(i, j)*c(j);
            F(5) -= b(i)*c(i)*a(i, j)*c(j);
            F(6) -= b(i)*a(i, j)*square(c(j));
            for (int k = 0; k < j; k++)
              {
                F(7) -= b(i)*a(i, j)*a(j, k)*c(k);
                for (int l = 0; l < k; l++)
                  F(8) -= b(i)*a(i, j)*a(j, k)*a(k, l)*c(l);
              }
          }
      }
  }
  
  
  //! returns the 2-norm of scheme
  template<class T>
  Real_wp LowStorageRK_Iterator<T>::GetNorm2Vector(const VectReal_wp& scheme)
  {
    return Norm2(scheme);
  }
  
  
  //! computes and factorises the jacobian matrix associated the non-linear system to solve
  template<class T>
  void LowStorageRK_Iterator<T>::ComputeAndFactoriseDiff(const VectReal_wp& X, const VectReal_wp& scale)
  {
    EvaluateJacobian(X, mat_jacobian_lu);  
    GetLU(mat_jacobian_lu, pivot_jacobian);
  }

  
  //! computes the non-linear system to solve in order to satisfy fourth-order conditions
  template<class T>
  void LowStorageRK_Iterator<T>::EvaluateFunction(const VectReal_wp& X, VectReal_wp& F)
  {
    ComputeScheme(X, F);
  }

  
  //! computes the jacobian matrix associated the non-linear system to solve
  template<class T>
  void LowStorageRK_Iterator<T>::EvaluateJacobian(const VectReal_wp& X, Matrix<Real_wp>& DF)
  {
    Real_wp one(1), two(2), three(3);
    Real_wp h = pow(epsilon_machine, one/three);
    int N = X.GetM();
    VectReal_wp f_plus(N), f_minus(N), Xtmp(X);
    DF.Reallocate(N, N);
    for (int j = 0; j < N; j++)
      {
        Xtmp(j) += h;
        ComputeScheme(Xtmp, f_plus);
        
        Xtmp(j) -= two*h;
        ComputeScheme(Xtmp, f_minus);
        
        for (int i = 0; i < N; i++)
          DF(i, j) = (f_plus(i) - f_minus(i))/(two*h);
        
        Xtmp(j) += h;
      }
  }


  //! solves the jacobian matrix associated the non-linear system to solve (X  = DF^{-1} F )
  template<class T>  
  void LowStorageRK_Iterator<T>::SolveDifferential(const VectReal_wp& F, VectReal_wp& X)
  {
    Copy(F, X);
    SolveLU(mat_jacobian_lu, pivot_jacobian, X);
  }
#endif  
  
  
  //! initialisation of u^0
  template<class T>
  void LowStorageRK_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = Y.GetM();
    Y = u0; u0.Clear();
    sys.SetDirichletCondition(t, 0, Y);
    Y_tilde.Reallocate(N); Fy.Reallocate(N);
    Y_tilde.Fill(0); Fy.Fill(0);
  }
  

  //! releases memory used by the scheme
  template<class T>
  void LowStorageRK_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! releases memory used by the scheme
  template<class T>
  void LowStorageRK_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_tilde.Clear(); Fy.Clear();
    Yfinal = Y; Y.Clear();
  }
  

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> LowStorageRK_Iterator<T>::GetStabilityFunction() const
  {
    UnivariatePolynomial<Real_wp> Pol, PolTilde, z;
    z.SetOrder(1); z(0) = 0.0; z(1) = 1.0;
    Pol.SetOrder(0); Pol(0) = 1.0;
    PolTilde.SetOrder(0); PolTilde(0) = 0.0;
    for (int i = 0; i < coef_C.GetM(); i++)
      {
        PolTilde = z*Pol + coef_alpha(i)*PolTilde;
        Add(coef_beta(i), PolTilde, Pol);
      }
    
    return Pol;
  }

  
  //! returns the CFL of the scheme
  template<class T>
  Real_wp LowStorageRK_Iterator<T>::GetCFL(bool real_p) const
  {
    if (real_p)
      return Real_wp(3.34);
    
    return Real_wp(4.65);
  }
  
  
  //! Low-storage Runge-Kutta algorithm
  template<class T>
  void LowStorageRK_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {    
    int nb_unknowns = Fy.GetM();
    if (Y.GetM() > nb_unknowns)
      {
	// new unknowns, resizing vectors
	Fy.Resize(Y.GetM());
	Y_tilde.Resize(Y.GetM());
	for (int i = nb_unknowns; i < Y.GetM(); i++)
	  {
	    Fy(i) = 0;
	    Y_tilde(i) = 0;
	  }
      }
	
    for (int i = 0; i < coef_C.GetM(); i++)
      {
	sys.EvaluateFunction(t + coef_C(i)*dt, Y, Fy);
        glob_chrono.Start(VirtualTimer::SCHEME);
        
        // simple loop (more optimal than calls to functions)
#pragma omp parallel for
        for (int j = 0; j < Y.GetM(); j++)
          {
            Y_tilde(j) = dt*Fy(j) + coef_alpha(i)*Y_tilde(j);
            Y(j) += coef_beta(i)*Y_tilde(j);
          }
        
        // Old calls to functions
	//Mlt(coef_alpha(i), Y_tilde); Add(dt, Fy, Y_tilde);
	//Add(coef_beta(i), Y_tilde, Y);
        
        glob_chrono.Stop(VirtualTimer::SCHEME);
      }
  }


  //! changes time step during the simulation
  template<class T>
  void LowStorageRK_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& LowStorageRK_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& LowStorageRK_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int LowStorageRK_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& LowStorageRK_Iterator<T>::GetIterate(int k)
  {
    return Y;
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t LowStorageRK_Iterator<T>::GetMemorySize() const
  {
    size_t taille = coef_alpha.GetMemorySize() + coef_beta.GetMemorySize() + coef_C.GetMemorySize();
    taille += Y.GetMemorySize() + Y_tilde.GetMemorySize() + Fy.GetMemorySize();
    return taille;
  }
  
  
  /*****************************
   * MultiStepButcher_Iterator *
   *****************************/

  
  //! default constructor
  template<class T>
  MultiStepButcher_Iterator<T>::MultiStepButcher_Iterator()
  {
    order = 1; type_algo = NO_CORRECTION;
    dt = 0;
  }
  
  
  //! computation of coefficients
  template<class T>
  void MultiStepButcher_Iterator<T>::SetOrder(int order_, int type)
  { 
    type_algo = type;
    Real_wp one(1);
    order = order_; 
    alpha.Reallocate(order); beta.Reallocate(order);
    gamma.Reallocate(order); delta.Reallocate(order);
    if (order==1)
      {
	alpha(0) = Real_wp(3)/2; 
	beta(0) = -one/2;
	gamma(0) = one/2;
        delta(0) = one/2;
      }
    else if (order==2)
      {
	alpha(0) = -one/2; alpha(1) = Real_wp(17)/12;
	beta(0) = Real_wp(3)/2; beta(1) = Real_wp(7)/12;
	gamma(0) = one/2; gamma(1) = -one/12;
	delta(0) = one/2; delta(1) = one/12;
      }
    else if (order==3)
      {
	alpha(0) = Real_wp(15)/2; alpha(1) = Real_wp(-31)/10; alpha(2) = Real_wp(37)/40;
	beta(0) = Real_wp(-13)/2; beta(1) = Real_wp(-29)/10; beta(2) = Real_wp(-49)/120;
	gamma(0) = one/2; gamma(1) = -one/10; gamma(2) = one/120;
	delta(0) = one/2; delta(1) = one/10; delta(2) = one/120;
      }
    else if (order==4)
      {
	alpha(0) = Real_wp(-65)/2; alpha(1) = Real_wp(515)/28;
        alpha(2) = Real_wp(-107)/28; alpha(3) = Real_wp(769)/1680;
	beta(0) = Real_wp(67)/2; beta(1) = Real_wp(437)/28;
        beta(2) = Real_wp(239)/84; beta(3) = Real_wp(117)/560;
	gamma(0) = one/2; gamma(1) = -Real_wp(3)/28; gamma(2) = one/84; gamma(3) = -one/1680;
	delta(0) = one/2; delta(1) = Real_wp(3)/28; delta(2) = one/84; delta(3) = one/1680;
      }
    else
      {
	cout << "Order not implemented" << endl;
	abort();
      }
  }
  
  
  //! initialisation of the scheme with initial condition U^0
  template<class T>
  void MultiStepButcher_Iterator<T>::
  SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
		      Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM();
    Y = u0; u0.Clear();
    sys.SetDirichletCondition(t, 0, Y);
    fnm1.Reallocate(order); fnm2.Reallocate(order);
    Y_star.Reallocate(N); Y_star.Fill(0);
    
    for (int i = 0; i < order; i++)
      {
	fnm1(i).Reallocate(N); fnm1(i).Fill(0);
	fnm2(i).Reallocate(N); fnm2(i).Fill(0);
      }
  }
  

  //! memory is released
  template<class T>
  void MultiStepButcher_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory is released
  template<class T>
  void MultiStepButcher_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_star.Clear(); fnm1.Clear(); fnm2.Clear();
    Yfinal = Y; Y.Clear(); 
  }
  

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> MultiStepButcher_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

    
  //! Butcher algorithm
  template<class T>
  void MultiStepButcher_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  { 
    if (n_time == 0)
      {
        // use of RK4 scheme for the first iterate
        Copy(Y, Y_star);
        sys.EvaluateFunction(t, Y, fnm1(0));
        Add(dt/6, fnm1(0), Y_star);
        
        Copy(Y, fnm2(0)); Add(dt/2, fnm1(0), fnm2(0));
        sys.EvaluateFunction(t+0.5*dt, fnm2(0), fnm1(0));
        Add(dt/3, fnm1(0), Y_star);
        
        Copy(Y, fnm2(0)); Add(dt/2, fnm1(0), fnm2(0));
        sys.EvaluateFunction(t+0.5*dt, fnm2(0), fnm1(0));
        Add(dt/3, fnm1(0), Y_star);
        
        Copy(Y, fnm2(0)); Add(dt, fnm1(0), fnm2(0));
        sys.EvaluateFunction(t+dt, fnm2(0), fnm1(0));
        Add(dt/6, fnm1(0), Y_star);

        // evaluating fnm1 and fnm2
        sys.EvaluateDerivativeFunction(t, 0, Y, fnm2(0));
	for (int i = 1; i < order; i++)
          sys.EvaluateDerivativeFunction(t, i, fnm2(i-1), fnm2(i));	

        sys.EvaluateDerivativeFunction(t+dt, 0, Y_star, fnm1(0));
	for (int i = 1; i < order; i++)
          sys.EvaluateDerivativeFunction(t+dt, i, fnm1(i-1), fnm1(i));
        
        Copy(Y_star, Y);
        
        return;
      }
    
    Vector<Real_wp> dtk(order); 
    dtk(0) = dt;
    for (int i = 1; i < order; i++)
      dtk(i) = dt*dtk(i-1);
    
    Copy(Y, Y_star);
    for (int i = 0; i < order; i++)
      {
	Add(alpha(i)*dtk(i), fnm1(i), Y_star);
	Add(beta(i)*dtk(i), fnm2(i), Y_star);
      }
    
    if (type_algo == PECE)
      {
	Copy(fnm1(0), fnm2(0));
	sys.EvaluateDerivativeFunction(t+dt, 0, Y_star, fnm1(0));
	for (int i = 1; i < order; i++)
	  {
	    Copy(fnm1(i), fnm2(i));
	    sys.EvaluateDerivativeFunction(t+dt, i, fnm1(i-1), fnm1(i));
	  }
	
	for (int i = 0; i < order; i++)
	  {
	    Add(gamma(i)*dtk(i), fnm1(i), Y);
	    Add(delta(i)*dtk(i), fnm2(i), Y);
	  }
	
	sys.EvaluateDerivativeFunction(t+dt, 0, Y, fnm1(0));
        for (int i = 1; i < order; i++)
          sys.EvaluateDerivativeFunction(t+dt, i, fnm1(i-1), fnm1(i));	
      }
    else
      {
	Copy(Y_star, Y); Copy(fnm1(0), fnm2(0));
	sys.EvaluateDerivativeFunction(t+dt, 0, Y, fnm1(0));
	for (int i = 1; i < order; i++)
	  {
	    Copy(fnm1(i), fnm2(i));
	    sys.EvaluateDerivativeFunction(t+dt, i, fnm1(i-1), fnm1(i));
	  }
      }
  }


  //! changes time step during the simulation
  template<class T>
  void MultiStepButcher_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented " << endl;
    abort();
  }
  
    
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& MultiStepButcher_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& MultiStepButcher_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
  */
  template<class T>
  int MultiStepButcher_Iterator<T>::GetNumberOfIterates() const
  {
    return 2*order+1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& MultiStepButcher_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else if (k <= order)
      return fnm1(k-1);
    else
      return fnm2(k-order-1);
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t MultiStepButcher_Iterator<T>::GetMemorySize() const
  {
    size_t taille = alpha.GetMemorySize() + beta.GetMemorySize() + delta.GetMemorySize() + gamma.GetMemorySize();
    taille += Y.GetMemorySize() + Y_star.GetMemorySize() + Seldon::GetMemorySize(fnm1) + Seldon::GetMemorySize(fnm2);
    return taille;
  }

  
  /************************
   * SplitScheme_Iterator *
   ************************/
  
  
  //! default constructor
  template<class T>
  SplitScheme_Iterator<T>::SplitScheme_Iterator()
  {
    order = 1;
    dt = 0;
  }
  
  
  //! sets the order of approximation
  template<class T>
  void SplitScheme_Iterator<T>::SetOrder(int r)
  {
    order = r;
    
    Real_wp zero(0), one(1), two(2);
    
    if ((order == 1) || (order == 2))
      {
        // classical schemes, no coefficients stored
      }
    else if (order == 4)
      {
        // Neri's scheme
        Real_wp two_third = pow(two, one/3);
        coef_ci.Reallocate(4); coef_di.Reallocate(4);
        coef_ci(0) = one/(two*(two - two_third));
        coef_ci(1) = (one-two_third)/(two*(two - two_third));
        coef_ci(2) = coef_ci(1); coef_ci(3) = coef_ci(0);
        coef_di(0) = one/(two - two_third);
        coef_di(1) = -two_third/(two - two_third);
        coef_di(2) = coef_di(0); coef_di(3) = zero;
      }
    else if (order == 6)
      {
        // Yoshida's scheme
        VectReal_wp w(4);
        w(1) = -1.17767998417887;
        w(2) = 0.235573213359357;
        w(3) = 0.784513610477560;
        w(0) = one - two*(w(1) + w(2) + w(3));
        
        coef_ci.Reallocate(8); coef_di.Reallocate(8);
        coef_ci(0) = 0.5*w(3); coef_ci(7) = coef_ci(0);
        coef_ci(1) = 0.5*(w(2)+w(3)); coef_ci(6) = coef_ci(1);
        coef_ci(2) = 0.5*(w(1)+w(2)); coef_ci(5) = coef_ci(2);
        coef_ci(3) = 0.5*(w(0)+w(1)); coef_ci(4) = coef_ci(3);
        
        coef_di(0) = w(3); coef_di(6) = coef_di(0); coef_di(7) = zero;
        coef_di(1) = w(2); coef_di(5) = coef_di(1);
        coef_di(2) = w(1); coef_di(4) = coef_di(2);
        coef_di(3) = w(0);
      }
    else if (order == 8)
      {
        // Yoshida's scheme
        VectReal_wp w(8);
        w(1) = 0.311790812418427;
        w(2) = -1.55946803821447;
        w(3) = -1.67896928259640;
        w(4) = 1.66335809963315;
        w(5) = -1.06458714789183;
        w(6) = 1.36934946416871;
        w(7) = 0.629030650210433;
        w(0) = one - two*(w(1) + w(2) + w(3) + w(4) + w(5) + w(6) + w(7));
        
        coef_ci.Reallocate(16); coef_di.Reallocate(16);
        coef_ci(0) = 0.5*w(7); coef_ci(15) = coef_ci(0);
        coef_ci(1) = 0.5*(w(6) + w(7)); coef_ci(14) = coef_ci(1);
        coef_ci(2) = 0.5*(w(5) + w(6)); coef_ci(13) = coef_ci(2);
        coef_ci(3) = 0.5*(w(4) + w(5)); coef_ci(12) = coef_ci(3);
        coef_ci(4) = 0.5*(w(3) + w(4)); coef_ci(11) = coef_ci(4);
        coef_ci(5) = 0.5*(w(2) + w(3)); coef_ci(10) = coef_ci(5);
        coef_ci(6) = 0.5*(w(1) + w(2)); coef_ci(9) = coef_ci(6);
        coef_ci(7) = 0.5*(w(0) + w(1)); coef_ci(8) = coef_ci(7);
        
        coef_di(0) = w(7); coef_di(14) = coef_di(0); coef_di(15) = zero;
        coef_di(1) = w(6); coef_di(13) = coef_di(1);
        coef_di(2) = w(5); coef_di(12) = coef_di(2);
        coef_di(3) = w(4); coef_di(11) = coef_di(3);
        coef_di(4) = w(3); coef_di(10) = coef_di(4);
        coef_di(5) = w(2); coef_di(9) = coef_di(5);
        coef_di(6) = w(1); coef_di(8) = coef_di(6);
        coef_di(7) = w(0);
      }
    else
      {
        cout << "Order not implemented" << endl;
        abort();
      }
  }
    
  
  //! provides the initial condition to the scheme
  template<class T>
  void SplitScheme_Iterator<T>::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
                                                    Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    Un = u0;
    Un_next = u0;
    u0.Clear();
  }
    

  //! memory used by the scheme is released
  template<class T>
  void SplitScheme_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory used by the scheme is released
  template<class T>
  void SplitScheme_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Un_next.Clear();
    Yfinal = Un; Un.Clear(); 
  }
    

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> SplitScheme_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! computes u^n+1 from u^n, u^n is overwritten by u^n+1
  template<class T>
  void SplitScheme_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    if (order == 1)
      {
        // basic splitting
        sys.SolveSplitOperator(t+0.5*dt, dt, Un, Un_next, 0);
        sys.SolveSplitOperator(t+0.5*dt, dt, Un_next, Un, 1);
      }
    else if (order == 2)
      {
        // Strang's splitting
        sys.SolveSplitOperator(t+0.25*dt, 0.5*dt, Un, Un_next, 0);
        sys.SolveSplitOperator(t+0.5*dt, dt, Un_next, Un, 1);
        sys.SolveSplitOperator(t+0.75*dt, 0.5*dt, Un, Un_next, 0);
        Copy(Un_next, Un);
      }
    else
      {
        // Yoshida's splitting
        int N = coef_ci.GetM()-1;
        for (int i = 0; i < N; i++)
          {
            sys.SolveSplitOperator(t, coef_ci(N-i)*dt, Un, Un_next, 0);
            sys.SolveSplitOperator(t, coef_di(N-i-1)*dt, Un_next, Un, 1);
          }
        
        sys.SolveSplitOperator(t, coef_ci(0)*dt, Un, Un_next, 0);
        Copy(Un_next, Un);
      }
  }
    
  
  //! changes time step during the simulation
  template<class T>
  void SplitScheme_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }
  
    
  //! returns the current solution u^n
  template<class T>
  Vector<T>& SplitScheme_Iterator<T>::GetIterate()
  {
    return Un;
  }
  
  
  //! returns the current solution u^n
  template<class T>
  const Vector<T>& SplitScheme_Iterator<T>::GetIterate() const
  {
    return Un;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int SplitScheme_Iterator<T>::GetNumberOfIterates() const
  {
    return 1;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& SplitScheme_Iterator<T>::GetIterate(int k)
  {
    return Un;
  }

  
  //! returns the size used to store the object in bytes
  template<class T>
  size_t SplitScheme_Iterator<T>::GetMemorySize() const
  {
    size_t taille = coef_ci.GetMemorySize() + coef_di.GetMemorySize();
    taille += Un.GetMemorySize() + Un_next.GetMemorySize();
    return taille;
  }
  
  
  /********************
   * Nystrom_Iterator *
   ********************/
  

  //! default constructor
  template<class T>
  Nystrom_Iterator<T>::Nystrom_Iterator()
  {
    dt = 0; type_algo = NYSTROM;
  }
  
  
  //! computation of coefficients
  /*!
    \param[in] order order of Nystrom coefficients
  */
  template<class T>
  void Nystrom_Iterator<T>::SetOrder(int order, int type)
  {
    type_algo = type;
    
    // evaluation of coefficients gamma and gamma_star
    VectReal_wp gamma, gamma_star;
    gamma.Reallocate(order);
    gamma_star.Reallocate(order);
    gamma(0) = 2.0; gamma_star(0) = 2.0;
    if (order > 1)
      {
        gamma(1) = 0;
        gamma_star(1) = -2.0;
      }
    
    for (int m = 2; m < order; m++)
      {
        gamma(m) = 1.0;
        for (int j = 0; j < m; j++)
          gamma(m) -= gamma(j)/(m+1-j);

        gamma_star(m) = 0.0;
        for (int j = 0; j < m; j++)
          gamma_star(m) -= gamma_star(j)/(m+1-j);
      }
    
    // computation of beta and beta_star from gamma, gamma_star
    beta.Reallocate(order); beta_star.Reallocate(order);
    beta(0) = gamma_star(order-1); beta_star(0) = gamma(order-1);
    VectReal_wp new_beta(order), new_beta_star(order);
    for (int m = order-2; m >= 0; m--)
      {
        new_beta(0) = gamma_star(m) + beta(0);
        for (int j = 1; j <= order-2-m; j++)
          new_beta(j) = -beta(j-1) + beta(j);
        
        new_beta(order-1-m) = -beta(order-2-m);
        beta = new_beta;

        new_beta_star(0) = gamma(m) + beta_star(0);
        for (int j = 1; j <= order-2-m; j++)
          new_beta_star(j) = -beta_star(j-1) + beta_star(j);
        
        new_beta_star(order-1-m) = -beta_star(order-2-m);
        beta_star = new_beta_star;
      }
  }
  
  
  //! initialisation of the scheme with u^0
  template<class T>
  void Nystrom_Iterator<T>
  ::SetInitialCondition(const Real_wp& t, const Real_wp& deltat,
			Vector<T>& u0, VirtualOdeSystem<T>& sys)
  {
    dt = deltat;
    int N = u0.GetM(); int order = beta.GetM();
    Y = u0; Y_star.Reallocate(N); Y_star.Fill(0);
    Ynm1.Reallocate(N); Ynm1.Fill(0);
    sys.SetDirichletCondition(t, 0, Y);
    F_star.Reallocate(N); F_star.Fill(0);
    fi.Reallocate(order);
    for (int i = 0; i < order; i++)
      {
	fi(i).Reallocate(N);
	fi(i).Fill(0);
      }    
    
    if (order >= 7)
      {
        rk_scheme.SetOrder(6);
        rk_scheme.SetInitialCondition(t, deltat, u0, sys);
      }

    u0.Clear();
  }
  

  //! memory is released
  template<class T>
  void Nystrom_Iterator<T>::Clear()
  {
    Vector<T> Yfinal; ClearFirst(Yfinal);
  }
  
  
  //! memory is released
  template<class T>
  void Nystrom_Iterator<T>::ClearFirst(Vector<T>& Yfinal)
  {
    Y_star.Clear();
    fi.Clear(); F_star.Clear();
    Yfinal = Y; Y.Clear(); 
    rk_scheme.Clear();
  }
  

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> Nystrom_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! Adams-Bashforth-Moulton algorithm
  template<class T>
  void Nystrom_Iterator<T>
  ::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    int order = beta.GetM();
    
    if (n_time < order)
      {
        if (order < 7)
          {
            // use of RK4 scheme for the first iterates
            Copy(Y, Y_star);
            sys.EvaluateFunction(t, Y, fi(0));
            Add(dt/6, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt/2, fi(0), F_star);
            sys.EvaluateFunction(t+0.5*dt, F_star, fi(0));
            Add(dt/3, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt/2, fi(0), F_star);
            sys.EvaluateFunction(t+0.5*dt, F_star, fi(0));
            Add(dt/3, fi(0), Y_star);
            
            Copy(Y, F_star); Add(dt, fi(0), F_star);
            sys.EvaluateFunction(t+dt, F_star, fi(0));
            Add(dt/6, fi(0), Y_star);
            
            Copy(Y, Ynm1);
            Copy(Y_star, Y);
          }
        else
          {
            rk_scheme.Advance(t, n_time, sys);

            Copy(Y, Ynm1);
            Y = rk_scheme.GetIterate();
          }
        
	// once Y_star = y^n+1 is computed
	// f(y^n+1) is stored in the vector fi
	sys.EvaluateFunction(t+dt, Y, fi(order-1-n_time));
	
	return;
      }
    
    if (type_algo == NYSTROM)
      {
        // Y_star = dt \sum beta_star_i  f_{n-i}
        Y_star.Fill(0);
	for (int i = 0; i < order; i++)
	  Add(dt*beta_star(i), fi(i), Y_star);
        
        // updating scheme
        T y_next;
        for (int i = 0; i < Y_star.GetM(); i++)
          {
            y_next = Ynm1(i) + Y_star(i);
            Ynm1(i) = Y(i);
            Y(i) = y_next;
          }
	
	// udpate of f(t_{n-i},y_{n-i})
	for (int i = (order-1); i > 0; i--)
	  Copy(fi(i-1), fi(i));
	
	sys.EvaluateFunction(t+dt, Y, fi(0));
      }
    else
      {
	// Nystrom step
	// y*_n = y_n + h \sum_{i=1}^k  beta_star(i) fi(i)
	Copy(Ynm1, Y_star);
	for (int i = 0; i < order; i++)
	  Add(dt*beta_star(i), fi(i), Y_star);
	
	// computation of f(t_n+1, y*_n)
	sys.EvaluateFunction(t+dt, Y_star, F_star);
	
	// Milne correction
	// y_{n+1} = y_{n-1} + h \beta_0 f(t_n+1, y*_n) + h \sum_{i=0}^k  beta(i) fi(i)
	Copy(Ynm1, Y_star);
        Add(dt*beta(0), F_star, Y_star);
	for (int i = 0; i < order-1; i++)
	  Add(dt*beta(i+1), fi(i), Y_star);
	
        // updating scheme
        Copy(Y, Ynm1);
        Copy(Y_star, Y);
        
	// udpate of f(t_{n-i},y_{n-i})
	for (int i = (order-1); i > 0; i--)
	  Copy(fi(i-1), fi(i));
	
	sys.EvaluateFunction(t+dt, Y, fi(0));
      }
  }
  

  //! changes time step during the simulation
  template<class T>
  void Nystrom_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    cout << "not implemented " << endl;
    abort();
  }

  
  //! returns the main unknown at the current time
  template<class T>
  Vector<T>& Nystrom_Iterator<T>::GetIterate()
  {
    return Y;
  }


  //! returns the main unknown at the current time
  template<class T>
  const Vector<T>& Nystrom_Iterator<T>::GetIterate() const
  {
    return Y;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int Nystrom_Iterator<T>::GetNumberOfIterates() const
  {
    return beta.GetM()+2;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& Nystrom_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else if (k == 1)
      return Ynm1;
    else
      return fi(k-2);
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t Nystrom_Iterator<T>::GetMemorySize() const
  {
    size_t taille = beta_star.GetMemorySize() + beta.GetMemorySize();
    taille += Y_star.GetMemorySize() + F_star.GetMemorySize() + Y.GetMemorySize() + Ynm1.GetMemorySize();
    taille += rk_scheme.GetMemorySize();
    return taille;
  }

  
  /******************************
   * RungeKuttaNystrom_Iterator *
   ******************************/
  
  
  //! default constructor
  template<class T>
  RungeKuttaNystrom_Iterator<T>::RungeKuttaNystrom_Iterator()
  {
    f_dependent_from_yprime = false;
    dt = 0;
    type_pair = OPTIMAL_PAIR;
  }
  
  
  //! sets order of the method
  template<class T>
  void RungeKuttaNystrom_Iterator<T>::SetOrder(int order, bool dependent, int type)
  {
    type_pair = type;
    Real_wp zero(0), one(1);
    f_dependent_from_yprime = dependent;
    order_rk = order;
    if (!dependent)
      {
	if (order == 2)
	  {
	    // dt <= 2 / sqrt(|| A ||)
	    // optimal CFL and conservative method
	    Abar.Reallocate(1, 1); Bbar.Reallocate(1);
	    B.Reallocate(1); C.Reallocate(1);
	    
	    Abar.Fill(zero);
	    B(0) = one; C(0) = one/2; Bbar(0) = one/2;	    	    
	  }
	else if (order == 3)
	  {
	    Abar.Reallocate(2, 2); Bbar.Reallocate(2);
	    B.Reallocate(2); C.Reallocate(2);
	    Abar.Fill(zero);
	    
	    // solution family with parameter alpha
	    // naive choice (Bbar(1) = 0)	    
	    // dt <= 2.1411 / sqrt(|| A ||)
	    //Real_wp alpha = one/3;
	    	    
	    // optimized choice (Chawla/Sharma)
	    // dt <= 2.498 / sqrt( || A || )
	    // optimal CFL, dissipative method
	    Real_wp alpha = (Real_wp(3) - sqrt(Real_wp(3)))/6;
	    
	    C(0) = alpha; C(1) = (2.0-3.0*alpha) / (3.0*(one-2*alpha));
	    B(0) = (one/2*C(1) - one/3)/(C(0)*(C(1) - C(0))); B(1) = one - B(0);
	    Bbar(0) = (one/2*C(1) - one/6)/(C(1) - C(0)); Bbar(1) = one/2 - Bbar(0);
	    Abar(1, 0) = one/(6*B(1));
	  }
	else if (order == 4)
          {
            Abar.Reallocate(3, 3); Bbar.Reallocate(3);
            B.Reallocate(3); C.Reallocate(3);            
            Abar.Fill(zero);

	    // classical fourth-order method
            // dt <= 2.58/sqrt(|| A ||)
	    // Real_wp alpha = 0;
            //Abar(1, 0) = one/8; Abar(2, 1) = one/2;
            
            //Bbar(0) = one/6; Bbar(1) = one/3; Bbar(2) = zero;
            //B(0) = one/6; B(1) = Real_wp(2)/3; B(2) = one/6;
            
            //C(0) = zero; C(1) = one/2; C(2) = one;
	    
	    // optimized choice (Chawla/Sharma)
	    // dt <= 3.94/sqrt(|| A ||)
	    // optimal CFL, but dissipative method
	    Real_wp alpha = one/(Real_wp(4)*(one + cos(pi_wp/9)));	    

	    // alpha = 1/6 provide conservative method, but lower CFL (3.46)
	    //Real_wp alpha = one/6;
	    
	    C(0) = alpha; C(1) = one/2; C(2) = one-alpha;
	    
	    Abar(1, 0) = (one - 4*alpha)*(one - 2*alpha) / (Real_wp(8)*(6*alpha*(alpha-one)+one));
	    Abar(2, 0) = 2*alpha*(one-2*alpha);
	    Abar(2, 1) = (one-2*alpha)*(one-4*alpha)/2;
	    
	    B(0) = one/(Real_wp(6)*square(one-2*alpha));
	    B(1) = Real_wp(4)*(6*alpha*(alpha-one) + one)/(Real_wp(6)*square(one-2*alpha));
	    B(2) = B(0);
	    
	    Bbar(0) = B(0)*(one-C(0));
	    Bbar(1) = B(1)*(one-C(1));
	    Bbar(2) = B(2)*(one-C(2));
          }
        else if (order == 5)
          {
            Abar.Reallocate(4, 4); Bbar.Reallocate(4);
            B.Reallocate(4); C.Reallocate(4);
            
            Abar.Fill(zero);
	    // method not A-stable
	    //Abar(1, 0) = one/50; Abar(2, 1) = Real_wp(7)/27; Abar(3, 2) = Real_wp(9)/35;
            //Abar(2, 0) = -one/27; Abar(3, 1) = Real_wp(-2)/35;
            //Abar(3, 0) = Real_wp(3)/10;
            
            //Bbar(0) = Real_wp(14)/336; Bbar(1) = Real_wp(100)/336;
            //Bbar(2) = Real_wp(54)/336; Bbar(3) = zero;
            //B(0) = Bbar(0); B(1) = Real_wp(125)/336;
            //B(2) = Real_wp(162)/336; B(3) = Real_wp(35)/336;
            
            //C(0) = zero; C(1) = one/5; C(2) = Real_wp(2)/3; C(3) = one;

	    // alternative family
	    // CFL is 2.7404
	    /* Real_wp alpha = 0.883415917004925;
	    C(0) = alpha;
	    C(1) = (Real_wp(4) - sqrt(Real_wp(6)))/10;
	    C(2) = (Real_wp(4) + sqrt(Real_wp(6)))/10;
	    C(3) = one;
	    
	    B(0) = zero;
	    B(1) = (Real_wp(16) - sqrt(Real_wp(6)))/36;
	    B(2) = (Real_wp(16) + sqrt(Real_wp(6)))/36;
	    B(3) = one/9;
	    
	    Bbar(0) = B(0)*(one-C(0));
	    Bbar(1) = B(1)*(one-C(1));
	    Bbar(2) = B(2)*(one-C(2));
	    Bbar(3) = B(3)*(one-C(3));
	    
	    Abar(1, 0) = C(1)*C(1)/2;
	    Abar(2, 1) = (one-5*C(0))/(120*B(2)*(C(1)-C(0))*(one-C(2)));
	    Abar(2, 0) = C(2)*C(2)/2 - Abar(2, 1);
	    Abar(3, 1) = Real_wp(3)/Real_wp(40)
	      *(20-25*C(1)-(one-5*C(0))/((C(1)-C(0))*(one-C(2))));
	    
	    Abar(3, 2) = 15*C(1)/8;
	    Abar(3, 0) = one/2 - Abar(3, 1) - Abar(3, 2);	    
	    
	    return; */

	    // optimized choice (Chawla/Sharma)
	    // dt <=  2.908 / sqrt(|| A ||)
	    // conservative method and optimal CFL
	    C(1) = Real_wp(4) / (Real_wp(11) - sqrt(16*sqrt(10) - Real_wp(39)));
	    C(3) = (165*C(1)*C(1) - 195*C(1) + Real_wp(50)
		    + sqrt(Real_wp(5)*(-4 + C(1)*(36 + C(1)*(-105 + C(1)*(90 + 45*C(1))))))) 
	      / (Real_wp(225)*C(1)*C(1) - 240*C(1) + Real_wp(60));
	    C(0) = zero;
	    C(2) = (12 - 15*(C(1)+C(3)) + 20*C(1)*C(3)) / (15-20*(C(1)+C(3))+30*C(1)*C(3));
	    B(1) = (3-4*(C(2)+C(3))+6*C(2)*C(3)) / (12*C(1)*(C(1)-C(2))*(C(1)-C(3)));
	    B(2) = (3-4*(C(1)+C(3))+6*C(1)*C(3)) / (12*C(2)*(C(2)-C(1))*(C(2)-C(3)));
	    B(3) = (3-4*(C(1)+C(2))+6*C(1)*C(2)) / (12*C(3)*(C(3)-C(1))*(C(3)-C(2)));
	    B(0) = one - B(1) - B(2) - B(3);
	    
	    Bbar(0) = B(0)*(one-C(0));
	    Bbar(1) = B(1)*(one-C(1));
	    Bbar(2) = B(2)*(one-C(2));
	    Bbar(3) = B(3)*(one-C(3));
	    
	    Abar(1, 0) = C(1)*C(1)/2;
	    Abar(2, 1) = (Real_wp(4) - 5*C(3)) / (Real_wp(120)*B(2)*C(1)*(C(2)-C(3)));
	    Abar(2, 0) = C(2)*C(2)/2 - Abar(2, 1);
	    Abar(3, 1) = (C(2)*(6 - 5*C(2)) + 5*C(1)*C(3) - 4*C(1) - 2*C(3))
	      / (Real_wp(120)*B(3)*C(1)*(C(3)-C(2))*(C(2)-C(1)));
	    Abar(3, 2) = (Real_wp(2) - 5*C(1)) / (Real_wp(120)*B(3)*C(2)*(C(2)-C(1)));
	    Abar(3, 0) = C(3)*C(3)/2 - Abar(3, 1) - Abar(3, 2);
          }
        else if (order == 6)
          {
	    if (type_pair == OPTIMAL_PAIR)
	      {
		Abar.Reallocate(5, 5); Bbar.Reallocate(5);
		B.Reallocate(5); C.Reallocate(5);
		Abar.Fill(zero);
		
		// Hairer's scheme
		// optimal choice dt <= 3.0895 / sqrt(||A||_2)
		Real_wp alpha = 0.229183261481182;
		C(0) = zero; C(1) = alpha; C(2) = one/2; C(3) = one-alpha; C(4) = one;
		B(0) = -(10*alpha*(alpha-one)+one) / (60*alpha*(one-alpha));
		B(1) = one / (60*alpha*(one-alpha)*square(one-2*alpha));
		B(2) = 8*(5*alpha*(alpha-one) + one) / (15*square(one-2*alpha));
		B(3) = B(1); B(4) = B(0);
		
		Bbar(0) = B(0)*(one-C(0));
		Bbar(1) = B(1)*(one-C(1));
		Bbar(2) = B(2)*(one-C(2));
		Bbar(3) = B(3)*(one-C(3));
		Bbar(4) = B(4)*(one-C(4));
		
		Abar(1, 0) = alpha*alpha/2;
		Abar(2, 0) = (5*alpha-one)*(6*alpha*(alpha-one)+one) / (48*alpha*(5*alpha*(alpha-one)+one));
		Abar(2, 1) = (one-3*alpha)*(one-2*alpha) / (48*alpha*(5*alpha*(alpha-one)+one));
		Abar(3, 0) = (one-alpha)/(6*alpha)*(one+alpha*(-9+alpha*(29-24*alpha)));
		Abar(3, 1) = -(one-alpha)/(6*alpha)*(one-2*alpha)*(one-6*alpha);
		Abar(3, 2) = 2*(one-alpha)*(one-2*alpha)*(one-3*alpha)/3;
		
		Abar(4, 0) = (one-3*alpha)*(one-4*alpha)/(6*alpha*(10*alpha*(alpha-one)+one));
		Abar(4, 1) = (one-alpha)*(-one+alpha*(9+alpha*(-29+24*alpha)))
		  / (6*alpha*(10*alpha*(alpha-one)+one)*square(one-2*alpha));
		Abar(4, 2) = 2*(one-alpha)*(one-5*alpha)*(6*alpha*(alpha-one)+one)
		  / (3*(10*alpha*(alpha-one)+one)*square(one-2*alpha));
		Abar(4, 3) = -alpha*alpha / (2*(10*alpha*(alpha-one)+one)*square(one-2*alpha));
		
		return;
	      }
	    
            // Dormand Prince pair
            // method not A-stable
            Abar.Reallocate(6, 6); Bbar.Reallocate(6);
            B.Reallocate(6); C.Reallocate(6);
            
            Real_wp R = sqrt(Real_wp(8581));
            Abar.Fill(zero);
            Abar(1, 0) = (Real_wp(26131)-209*R)/810000;
            Abar(2, 0) = (Real_wp(26131)-209*R)/607500;
            Abar(2, 1) = (Real_wp(26131)-209*R)/303750;
#ifdef MONTJOIE_WITH_MPFR
            Abar(3, 0) = (Real_wp("980403512254")
                          +Real_wp("7781688431")*R)/Real_wp("11694469921875");
            Abar(3, 1) = -(Real_wp("1262884486208")
                           +Real_wp("15385481287")*R)/Real_wp("11694469921875");
            Abar(3, 2) = (Real_wp("7166233891441")
                          +Real_wp("78694563299")*R)/Real_wp("46777879687500");
            Abar(4, 0) = -Real_wp(9)*(Real_wp(329260)+3181*R)/27040000;
            Abar(4, 1) = Real_wp(27)*(Real_wp(35129)+3331*R)/13520000;
	    Abar(4, 2) = -Real_wp(27)*(Real_wp(554358343)+31040327*R)/Real_wp("464060480000");
            Abar(4, 3) = Real_wp(153)*(Real_wp(8555257)-67973*R)/Real_wp("2745920000");
#else
	    Abar(3, 0) = 1.45474742801091785895935232316e-1;
	    Abar(3, 1) = -2.29860640522647473120261456297e-1;
	    Abar(3, 2) = 3.09034987202967536528726080729e-1;
	    Abar(4, 0) = -2.07668262950789954335146205252e-1;
	    Abar(4, 1) = 6.86366784292514312273571851042e-1;
	    Abar(4, 2) = -1.99549277872349252201326358888e-1;
	    Abar(4, 3) = 1.25850756530624894262900713098e-1;      
#endif
            Abar(5, 0) = Real_wp(329)/4212;
            Abar(5, 1) = zero;
            Abar(5, 2) = (Real_wp(84119543)+366727*R)/409622616;
            Abar(5, 3) = (Real_wp(84119543)-366727*R)/409622616;
            Abar(5, 4) = Real_wp(200)/17901;
            
            Bbar(0) = Abar(5, 0); Bbar(1) = zero; Bbar(2) = Abar(5, 2);
            Bbar(3) = Abar(5, 3); Bbar(4) = Abar(5, 4); Bbar(5) = zero;
            
            B(0) = Bbar(0); B(1) = zero;
#ifdef MONTJOIE_WITH_MPFR
	    B(2) = (Real_wp(389225579)+96856*R)/1024056540;
            B(3) = (Real_wp(389225579)-96856*R)/1024056540;
#else
	    B(2) = 3.88843478705982602715952208523e-1;
	    B(3) = 3.71320757928842267403035557523e-1;
#endif
	    B(4) = Real_wp(2000)/17901; B(5) = one/20;
            
            C(0) = zero; C(1) = (Real_wp(209)-R)/900; C(2) = (Real_wp(209)-R)/450;
            C(3) = (Real_wp(209)+R)/450; C(4) = Real_wp(9)/10; C(5) = one;

            // embedded formula of order 4
            //Bbar(0) = (Real_wp(2701)+23*R)/4563; Bbar(1) = -(Real_wp(9829)+131*R)/9126;
            //Bbar(2) = Real_wp(5)*(Real_wp(1798)+17*R)/9126;
            // Bbar(3) = zero; Bbar(4) = zero; Bbar(5) = zero;

            //B(0) = Real_wp(115)/2106; B(1) = zero; B(2) = (Real_wp(84119543)+366727*R)/256014135;
            //B(3) = (Real_wp(84119543)-366727*R)/256014135;
            // B(4) = Real_wp(6950)/17901; B(5) = -one/10;
          }
	else if (order == 7)
	  {
	    // optimal CFL for x = 0.110451398065702 0.173816271367107 0.459433163929695 0.652002232653235
	    // dt <= 7.0875 / sqrt{||A||_2}
	    
	    Abar.Reallocate(7, 7); Bbar.Reallocate(7);
	    B.Reallocate(7); C.Reallocate(7);
	    Abar.Fill(zero);

	    C(0) = zero;
	    C(1) = to_num<Real_wp>("1.10451398065702005579780120569921564e-01");
	    C(2) = to_num<Real_wp>("1.73816271367107005874075298379466403e-01");
	    C(3) = to_num<Real_wp>("4.59433163929695020932797433488303795e-01");
	    C(4) = to_num<Real_wp>("6.52002232653234981007983606104971841e-01");
	    C(5) = to_num<Real_wp>("8.70316039276850641731485221667251659e-01");
	    C(6) = one;
	    
	    Abar(1, 0) = to_num<Real_wp>("6.09975566733408036984889161815912258e-03");
	    Abar(2, 0) = to_num<Real_wp>("7.18196866340384833275861380082859722e-03");
	    Abar(2, 1) = to_num<Real_wp>("7.92407943257804305363202701761677378e-03");
	    Abar(3, 0) = to_num<Real_wp>("4.26084378823758050365778606252709387e-02");
	    Abar(3, 1) = to_num<Real_wp>("-8.24485401360933505966230685930712040e-02");
	    Abar(3, 2) = to_num<Real_wp>("1.45379518312942554584717959247231131e-01");
	    Abar(4, 0) = to_num<Real_wp>("-5.37088465036575760649625031419185680e-02");
	    Abar(4, 1) = to_num<Real_wp>("3.68661957798455919895784614844303236e-01");
	    Abar(4, 2) = to_num<Real_wp>("-1.83888557354620941168652331094260465e-01");
	    Abar(4, 3) = to_num<Real_wp>("8.14889017522241751888711105979766147e-02");
	    Abar(5, 0) = to_num<Real_wp>("1.70008045806008422772667069945839198e-01");
	    Abar(5, 1) = to_num<Real_wp>("-5.31910248775804906934345153218386581e-01");
	    Abar(5, 2) = to_num<Real_wp>("6.41668264620689834986545041313623124e-01");
	    Abar(5, 3) = to_num<Real_wp>("3.86015134460949377057623569522366645e-02");
	    Abar(5, 4) = to_num<Real_wp>("6.03574290142840259140481436287346275e-02");
	    Abar(6, 0) = to_num<Real_wp>("-3.29087728094379976485470177712491783e-01");
	    Abar(6, 1) = to_num<Real_wp>("1.64392502759664414591706705102985975e+00");
	    Abar(6, 2) = to_num<Real_wp>("-1.21266590861644488523436114606010461e+00");
	    Abar(6, 3) = to_num<Real_wp>("3.83564492404034718609694779559539695e-01");
	    Abar(6, 4) = to_num<Real_wp>("-3.31512602882803733557862215352776886e-02");
	    Abar(6, 5) = to_num<Real_wp>("4.74153769984263705488557147184746397e-02");
	    
	    Bbar(0) = to_num<Real_wp>("5.15872003859607865072663944920218165e-02");
	    Bbar(1) = zero;
	    Bbar(2) = to_num<Real_wp>("2.19783461816409289345320164170947548e-01");
	    Bbar(3) = to_num<Real_wp>("1.37308131717759403763775261611777162e-01");
	    Bbar(4) = to_num<Real_wp>("6.45744543518716639683756710428985512e-02");
	    Bbar(5) = to_num<Real_wp>("2.67467517279988564152625086823549215e-02");
	    Bbar(6) = zero;
	    
	    B(0) = to_num<Real_wp>("5.15872003859607865072663944920218165e-02");
	    B(1) = zero;
	    B(2) = to_num<Real_wp>("2.66022501048393332859171930797456505e-01");
	    B(3) = to_num<Real_wp>("2.54007687034469496467423582204510828e-01");
	    B(4) = to_num<Real_wp>("1.85559967364750240943299246803018726e-01");
	    B(5) = to_num<Real_wp>("2.06245641934842618966934601003454012e-01");
	    B(6) = to_num<Real_wp>("3.65770022315835242559042446995381115e-02");
	  }
	else if (order == 8)
	  {
            // x = 0.135294127286225 0.24015308384744 0.453046953126355 0.695039606659698
            // dt <= 7.8525 / sqrt{||A||_2}

	    Abar.Reallocate(8, 8); Bbar.Reallocate(8);
	    B.Reallocate(8); C.Reallocate(8);
	    Abar.Fill(zero);
	    
	    C(0) = zero;
	    C(1) = to_num<Real_wp>("6.76470636431125010856746371246117633e-02");
	    C(2) = to_num<Real_wp>("1.35294127286225002171349274249223527e-01");
	    C(3) = to_num<Real_wp>("2.40153083847439996478811963243060745e-01");
	    C(4) = to_num<Real_wp>("4.53046953126355012031467595079448074e-01");
	    C(5) = to_num<Real_wp>("6.95039606659697994928137632086873055e-01");
	    C(6) = to_num<Real_wp>("8.77487149679790856084927881563740034e-01");
	    C(7) = one;

	    Abar(1, 0) = to_num<Real_wp>("2.28806260976765658382704510140668146e-03");
	    Abar(2, 0) = to_num<Real_wp>("3.05075014635687544510272680187557528e-03");
	    Abar(2, 1) = to_num<Real_wp>("6.10150029271375089020545360375115056e-03");
	    Abar(3, 0) = to_num<Real_wp>("7.93634919636077359237354983187964128e-03");
	    Abar(3, 1) = to_num<Real_wp>("7.67646192116656802385029978026138770e-03");
	    Abar(3, 2) = to_num<Real_wp>("1.32239407231904310547492613539671037e-02");
	    Abar(4, 0) = to_num<Real_wp>("-5.31871144289635360146073427016748354e-02");
	    Abar(4, 1) = to_num<Real_wp>("2.36044633344624708849805558381396245e-01");
	    Abar(4, 2) = to_num<Real_wp>("-1.79271982092878879587216984208660673e-01");
	    Abar(4, 3) = to_num<Real_wp>("9.90402340457545644571362971306571663e-02");
	    Abar(5, 0) = to_num<Real_wp>("4.90355260776335867108556529640790833e-01");
	    Abar(5, 1) = to_num<Real_wp>("-1.49943742478680635969109255005908820e+00");
	    Abar(5, 2) = to_num<Real_wp>("1.59886032314332738026576595345173096e+00");
	    Abar(5, 3) = to_num<Real_wp>("-4.64289148325384786114948230250881684e-01");
	    Abar(5, 4) = to_num<Real_wp>("1.16051016605361751123161797592983842e-01");
	    Abar(6, 0) = to_num<Real_wp>("-6.92502765701010234294253963376878460e-01");
	    Abar(6, 1) = to_num<Real_wp>("2.41422212747769107963400089411525982e+00");
	    Abar(6, 2) = to_num<Real_wp>("-2.33465912691602993940635418262261007e+00");
	    Abar(6, 3) = to_num<Real_wp>("9.38656627817305460910770919304652028e-01");
	    Abar(6, 4) = to_num<Real_wp>("6.03592318815413418272570336462726437e-03");
	    Abar(6, 5) = to_num<Real_wp>("5.32390630604713399264011018872648789e-02");
	    Abar(7, 0) = to_num<Real_wp>("1.57644508199114216198898729030405752e+00");
	    Abar(7, 1) = to_num<Real_wp>("-5.17271376418954017841975722381190163e+00");
	    Abar(7, 2) = to_num<Real_wp>("5.58338164728802516538989832308607150e+00");
	    Abar(7, 3) = to_num<Real_wp>("-1.98207946139056331757622229921653891e+00");
	    Abar(7, 4) = to_num<Real_wp>("4.67372806575426333681573057076549864e-01");
	    Abar(7, 5) = to_num<Real_wp>("-6.76565728608477766538264111427169840e-03");
	    Abar(7, 6) = to_num<Real_wp>("3.43593470115946126009034936760333554e-02");
	    
	    Bbar(0) = to_num<Real_wp>("4.22648292488055179219567776938578283e-02");
	    Bbar(1) = zero;
	    Bbar(2) = to_num<Real_wp>("1.56999374433680387929725938638660018e-01");
	    Bbar(3) = to_num<Real_wp>("6.74417722662260878866967061900736795e-02");
	    Bbar(4) = to_num<Real_wp>("1.52240380082141876010559273409320573e-01");
	    Bbar(5) = to_num<Real_wp>("5.95587501397858668741222820182381299e-02");
	    Bbar(6) = to_num<Real_wp>("2.14948938293602633769390220498497713e-02");
	    Bbar(7) = zero;

	    B(0) = to_num<Real_wp>("4.22648292488055179219567776938578283e-02");
	    B(1) = zero;
	    B(2) = to_num<Real_wp>("1.81563904430250720489223626708463830e-01");
	    B(3) = to_num<Real_wp>("8.87570520226804628073088145877463872e-02");
	    B(4) = to_num<Real_wp>("2.78342685816159041208762209867990998e-01");
	    B(5) = to_num<Real_wp>("1.95299951863994684924106477115390700e-01");
	    B(6) = to_num<Real_wp>("1.75450116238252002598000576542664045e-01");
	    B(7) = to_num<Real_wp>("3.83214603798575700506415174838862108e-02");
	  }
	else if (order == 10)
	  {
	    if (type_pair == OPTIMAL_PAIR)
	      {
		Abar.Reallocate(11, 11); C.Reallocate(11);
		Bbar.Reallocate(11); B.Reallocate(11);
		B.Fill(zero); Bbar.Fill(zero);
		
		// Hairer's family with 11 steps and optimized coefficients
		// CFL <= 4.7527 / sqrt{ ||A||_2 }
		C(0) = zero;
		C(1) = to_num<Real_wp>("1.42258516529654787280495387961238034e-01");
		C(2) = to_num<Real_wp>("2.84517033059309574560990775922476068e-01");
		C(3) = to_num<Real_wp>("8.82527661964732346425501486979669075e-01");
		C(4) = to_num<Real_wp>("6.42615758240322548157075497020439536e-01");
		C(5) = to_num<Real_wp>("1.17472338035267653574498513020330925e-01");
		C(6) = to_num<Real_wp>("3.57384241759677451842924502979560464e-01");
		C(7) = to_num<Real_wp>("8.82527661964732346425501486979669075e-01");
		C(8) = to_num<Real_wp>("2.84517033059309574560990775922476068e-01");
		C(9) = zero;
		C(10) = one;
		
		Abar(1, 0) = to_num<Real_wp>("1.01187427626090321710860464807833977e-02");
		Abar(2, 0) = to_num<Real_wp>("1.34916570168120428947813953077111970e-02");
		Abar(2, 1) = to_num<Real_wp>("2.69833140336240857895627906154223940e-02");
		Abar(3, 0) = to_num<Real_wp>("4.30435109016211931870206095384142825e-01");
		Abar(3, 1) = to_num<Real_wp>("-8.87310981777171926495514760028347620e-01");
		Abar(3, 2) = to_num<Real_wp>("8.46303409827428436772602452820221004e-01");
		Abar(4, 0) = to_num<Real_wp>("5.75065858705501067406592793593364471e-02");
		Abar(4, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(4, 2) = to_num<Real_wp>("1.45887694153335099629854550557425308e-01");
		Abar(4, 3) = to_num<Real_wp>("3.08322634550713210681420724060967468e-03");
		Abar(5, 0) = to_num<Real_wp>("5.39472302083110558506596566968858728e-03");
		Abar(5, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(5, 2) = to_num<Real_wp>("2.13598787811262291344809091667854013e-03");
		Abar(5, 3) = to_num<Real_wp>("2.82777036684632509713401825351367559e-04");
		Abar(5, 4) = to_num<Real_wp>("-9.13612833892265286435157215371360828e-04");
		Abar(6, 0) = to_num<Real_wp>("-1.36841264551267974546489158297160215e-02");
		Abar(6, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(6, 2) = to_num<Real_wp>("-3.99811188980866174583397519725051276e-02");
		Abar(6, 3) = to_num<Real_wp>("-5.23469791631888034986528661855080587e-03");
		Abar(6, 4) = to_num<Real_wp>("1.74841752842128615728346388711275402e-02");
		Abar(6, 5) = to_num<Real_wp>("1.05277516114389224010271855686576308e-01");
		Abar(7, 0) = to_num<Real_wp>("-4.73456184955383390795846091005089416e-02");
		Abar(7, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(7, 2) = to_num<Real_wp>("-4.23743858788686507940334767865457753e-01");
		Abar(7, 3) = to_num<Real_wp>("1.13221953604904855706417888603329222e-03");
		Abar(7, 4) = to_num<Real_wp>("5.92561155805353959627250552153525800e-02");
		Abar(7, 5) = to_num<Real_wp>("3.74758283782496303359144274620312867e-01");
		Abar(7, 6) = to_num<Real_wp>("4.25370395451612541288279656420284165e-01");
		Abar(8, 0) = to_num<Real_wp>("5.32167041818089493457742338866506165e-02");
		Abar(8, 1) = to_num<Real_wp>("2.69833140336240857895627906154223940e-02");
		Abar(8, 2) = to_num<Real_wp>("7.48461026890723138580763349076303298e-02");
		Abar(8, 3) = to_num<Real_wp>("5.37020444352510481334404510138484566e-03");
		Abar(8, 4) = to_num<Real_wp>("-2.57972156168540331175749450046521932e-02");
		Abar(8, 5) = to_num<Real_wp>("-9.56496528139340803405735343817643824e-02");
		Abar(8, 6) = to_num<Real_wp>("-8.58589657784671761979758994351161829e-04");
		Abar(8, 7) = to_num<Real_wp>("2.36410379097846009771501979281314244e-03");
		Abar(9, 0) = to_num<Real_wp>("-1.19178005011498929545408849135041619e-02");
		Abar(9, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(9, 2) = to_num<Real_wp>("-2.91994044050538961776405711153706332e-01");
		Abar(9, 3) = to_num<Real_wp>("5.12764121092785313976895638281852086e-03");
		Abar(9, 4) = to_num<Real_wp>("-9.76278465378848754006192439027094184e-03");
		Abar(9, 5) = to_num<Real_wp>("4.26839128606462241119437787322013304e-02");
		Abar(9, 6) = to_num<Real_wp>("8.62701043355193366659768389974564646e-02");
		Abar(9, 7) = to_num<Real_wp>("-3.54057408752857202153576716776132327e-03");
		Abar(9, 8) = to_num<Real_wp>("1.83133544885912500374854713512766443e-01");
		Abar(10, 0) = to_num<Real_wp>("-2.33825771656339011636621168264749170e-01");
		Abar(10, 1) = to_num<Real_wp>("0.00000000000000000000000000000000000e+00");
		Abar(10, 2) = to_num<Real_wp>("8.11654972822852162237641553647185415e-01");
		Abar(10, 3) = to_num<Real_wp>("-1.13649356716289902062383924463188010e-02");
		Abar(10, 4) = to_num<Real_wp>("1.41931560837344334785123015547930666e-01");
		Abar(10, 5) = to_num<Real_wp>("-5.48667455428601908486148173282159517e-01");
		Abar(10, 6) = to_num<Real_wp>("-7.80088351619539552382175211323637163e-01");
		Abar(10, 7) = to_num<Real_wp>("3.56197116284457610676534456335243798e-02");
		Abar(10, 8) = to_num<Real_wp>("5.84740269087467204620764930488224190e-01");
		Abar(10, 9) = to_num<Real_wp>("5.00000000000000000000000000000000000e-01");
		
		Bbar(0) = zero;
		Bbar(1) = zero;
		Bbar(2) = to_num<Real_wp>("-7.15482966940690465156292405311234388e-02");
		Bbar(3) = zero;
		Bbar(4) = to_num<Real_wp>("9.91488201804162591693788280356558222e-02");
		Bbar(5) = to_num<Real_wp>("1.67007309146871573763622194916580283e-01");
		Bbar(6) = to_num<Real_wp>("1.78280368337326917338981434524998518e-01");
		Bbar(7) = to_num<Real_wp>("2.22301690020519163946842091894320437e-02");
		Bbar(8) = to_num<Real_wp>("7.15482966940690465156292405311234388e-02");
		Bbar(9) = to_num<Real_wp>("3.33333333333333333333333333333333333e-02");
		Bbar(10) = zero;
		
		B(0) = zero;
		B(1) = zero;
		B(2) = to_num<Real_wp>("-1.00000000000000005551115123125782702e-01");
		B(3) = zero;
		B(4) = to_num<Real_wp>("2.77429188517743176508360262560654340e-01");
		B(5) = to_num<Real_wp>("1.89237478148923490158306404106012326e-01");
		B(6) = to_num<Real_wp>("2.77429188517743176508360262560654340e-01");
		B(7) = to_num<Real_wp>("1.89237478148923490158306404106012326e-01");
		B(8) = to_num<Real_wp>("1.00000000000000005551115123125782702e-01");
		B(9) = to_num<Real_wp>("3.33333333333333333333333333333333333e-02");
		B(10) = to_num<Real_wp>("3.33333333333333333333333333333333333e-02");
		
		return;
	      }
	    
	    // Sharp Qureshi Grazier pair
	    // CFL <= 5.291 / sqrt{||A||_2}
	    Abar.Reallocate(13, 13); C.Reallocate(13);
	    Bbar.Reallocate(13); B.Reallocate(13);
	    B.Fill(zero); Bbar.Fill(zero);

	    C(0) = zero;
	    C(1) = to_num<Real_wp>("3.5369578561715839852580662156128003290245e-02");
	    C(2) = to_num<Real_wp>("7.0739157123431679705161324312256006580491e-02");
	    C(3) = to_num<Real_wp>("1.9397227470895648005755022968171580307125e-01");
	    C(4) = to_num<Real_wp>("2.7465849059990290000000000000000000000000e-01");
	    C(5) = to_num<Real_wp>("1.9939545825206940000000000000000000000000e-01");
	    C(6) = to_num<Real_wp>("5.2332097109723180000000000000000000000000e-02");
	    C(7) = to_num<Real_wp>("4.1285926718800690000000000000000000000000e-01");
	    C(8) = to_num<Real_wp>("5.9440567007045230000000000000000000000000e-01");
	    C(9) = to_num<Real_wp>("6.9648498893199050000000000000000000000000e-01");
	    C(10) = to_num<Real_wp>("8.5840043338316930000000000000000000000000e-01");
	    C(11) = to_num<Real_wp>("9.5922053070763063933434705195204984252787e-01");
	    C(12) = to_num<Real_wp>("1.0000000000000000000000000000000000000000e+00");
	    
	    Abar(1, 0) = to_num<Real_wp>("6.2550354381669436926370260206784041201048e-04");
	    Abar(2, 0) = to_num<Real_wp>("8.3400472508892582568493680275712054934730e-04");
	    Abar(2, 1) = to_num<Real_wp>("1.6680094501778516513698736055142410986946e-03");
	    Abar(3, 0) = to_num<Real_wp>("1.4377592173651130455265029764307572975032e-02");
	    Abar(3, 1) = to_num<Real_wp>("-2.5520389586889438033660194426507730127969e-02");
	    Abar(3, 2) = to_num<Real_wp>("2.9955419091121746433503478890580675956400e-02");
	    Abar(4, 0) = to_num<Real_wp>("5.6606010103291956408355747332339448084576e-03");
	    Abar(4, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(4, 2) = to_num<Real_wp>("2.2438147205568485440414320313187753497836e-02");
	    Abar(4, 3) = to_num<Real_wp>("9.6198950134107937593128496677833016937067e-03");
	    Abar(5, 0) = to_num<Real_wp>("4.2497841145700712790685880544080818730288e-03");
	    Abar(5, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(5, 2) = to_num<Real_wp>("1.3737389558120657344878988260209037684744e-02");
	    Abar(5, 3) = to_num<Real_wp>("2.1090823734166889521084117969028056034572e-03");
	    Abar(5, 4) = to_num<Real_wp>("-2.1698166033104208350129690333992516122959e-04");
	    Abar(6, 0) = to_num<Real_wp>("8.4469918165805880106087924734363680346449e-04");
	    Abar(6, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(6, 2) = to_num<Real_wp>("7.6286170426770036487508428629504122228445e-04");
	    Abar(6, 3) = to_num<Real_wp>("-4.5274889794326362037916138257410810321144e-03");
	    Abar(6, 4) = to_num<Real_wp>("-9.4704498758374659029865654611073300453711e-05");
	    Abar(6, 5) = to_num<Real_wp>("4.3839567862160003018135640613696763068192e-03");
	    Abar(7, 0) = to_num<Real_wp>("4.0257188676144292921381323534074832686999e-02");
	    Abar(7, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(7, 2) = to_num<Real_wp>("2.8132558570671430000000000000000000000000e-01");
	    Abar(7, 3) = to_num<Real_wp>("-9.9086973311055342850364525233745309401480e-02");
	    Abar(7, 4) = to_num<Real_wp>("3.1558678646031991227831775524256371039974e-02");
	    Abar(7, 5) = to_num<Real_wp>("7.6263009309052972315471957912384344318273e-02");
	    Abar(7, 6) = to_num<Real_wp>("-2.4509110177537917817170493451316523864377e-01");
	    Abar(8, 0) = to_num<Real_wp>("-5.1852206149158582964016598801731374371401e-01");
	    Abar(8, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(8, 2) = to_num<Real_wp>("-4.0220496748396490000000000000000000000000e+00");
	    Abar(8, 3) = to_num<Real_wp>("1.3339542387088720000000000000000000000000e+00");
	    Abar(8, 4) = to_num<Real_wp>("-3.6280137620198678475540384353108923505414e-01");
	    Abar(8, 5) = to_num<Real_wp>("-4.4662932597538443824830112336373381632573e-01");
	    Abar(8, 6) = to_num<Real_wp>("4.1093353899979717196794462089910810944472e+00");
	    Abar(8, 7) = to_num<Real_wp>("8.3371860107714029551447509208700700646694e-02");
	    Abar(9, 0) = to_num<Real_wp>("4.5526515039304022037674417630889739137541e-01");
	    Abar(9, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(9, 2) = to_num<Real_wp>("3.4410244296399106992261498105244921193431e+00");
	    Abar(9, 3) = to_num<Real_wp>("-1.1420304634021581044575027059094514599562e+00");
	    Abar(9, 4) = to_num<Real_wp>("3.3122504529344460041163434350588640458224e-01");
	    Abar(9, 5) = to_num<Real_wp>("5.2510814510725835145952044841619177157062e-01");
	    Abar(9, 6) = to_num<Real_wp>("-3.4007204966989148959577640759197726181960e+00");
	    Abar(9, 7) = to_num<Real_wp>("1.4959025547025800806021814497275304454960e-02");
	    Abar(9, 8) = to_num<Real_wp>("1.7714834024190792778113334621606086825758e-02");
	    Abar(10, 0) = to_num<Real_wp>("-5.0738859071291316865190644559724141522647e-02");
	    Abar(10, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(10, 2) = to_num<Real_wp>("-8.5258315802766612821674340523993103859370e-01");
	    Abar(10, 3) = to_num<Real_wp>("2.9256281989303017562028747765972648095241e-01");
	    Abar(10, 4) = to_num<Real_wp>("-4.2631304538837970000000000000000000000000e-01");
	    Abar(10, 5) = to_num<Real_wp>("3.0845126792446679274386181736916000785280e-01");
	    Abar(10, 6) = to_num<Real_wp>("8.2068063064681728553592877499937152955118e-01");
	    Abar(10, 7) = to_num<Real_wp>("2.6353200561785124318114669620069412596405e-01");
	    Abar(10, 8) = to_num<Real_wp>("-3.8002959536498107924589290701514760515073e-02");
	    Abar(10, 9) = to_num<Real_wp>("5.0836949957876193531014830503462796310983e-02");
	    Abar(11, 0) = to_num<Real_wp>("-8.5506341904459305448936612697062868735670e-01");
	    Abar(11, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(11, 2) = to_num<Real_wp>("-4.4172967833111218167905547111397401957735e+00");
	    Abar(11, 3) = to_num<Real_wp>("1.4620414738250776987869127049027576880060e+00");
	    Abar(11, 4) = to_num<Real_wp>("1.5800603670627463865955782231377536485336e+00");
	    Abar(11, 5) = to_num<Real_wp>("-2.1500473087668603894015427990085825826253e+00");
	    Abar(11, 6) = to_num<Real_wp>("5.2192882951843358632420497699580187954903e+00");
	    Abar(11, 7) = to_num<Real_wp>("-7.0122224598841812446074882425859149461800e-01");
	    Abar(11, 8) = to_num<Real_wp>("4.0674284722476557322266460477065027849355e-01");
	    Abar(11, 9) = to_num<Real_wp>("-1.0995016400771572122089410865399735394434e-01");
	    Abar(11, 10) = to_num<Real_wp>("2.5498951087297871672803054502111066906939e-02");
	    Abar(12, 0) = to_num<Real_wp>("3.6124205767243278950403294321213352485132e+00");
	    Abar(12, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
	    Abar(12, 2) = to_num<Real_wp>("1.9614376424164347318861614225759179537283e+01");
	    Abar(12, 3) = to_num<Real_wp>("-6.5540954527470471068838634201654381264077e+00");
	    Abar(12, 4) = to_num<Real_wp>("-5.3634775178894355119823336908885299369771e+00");
	    Abar(12, 5) = to_num<Real_wp>("8.9549200633414078297191926517465420210063e+00");
	    Abar(12, 6) = to_num<Real_wp>("-2.2111999575685298066957085805778892539418e+01");
	    Abar(12, 7) = to_num<Real_wp>("3.0666418328538943745146953495973014350614e+00");
	    Abar(12, 8) = to_num<Real_wp>("-1.2974380337600327021222956536808230517716e+00");
	    Abar(12, 9) = to_num<Real_wp>("5.9822684827188431105752038905995447885348e-01");
	    Abar(12, 10) = to_num<Real_wp>("-2.4107078892562108246647516780668823095046e-02");
	    Abar(12, 11) = to_num<Real_wp>("4.5319136185137669988740390100397569527517e-03");
	    
	    Bbar(0) = to_num<Real_wp>("1.1445045431083081161076675575316257764110e-02");
	    Bbar(5) = to_num<Real_wp>("1.5182228814165001267592100569234113490800e-01");
	    Bbar(6) = to_num<Real_wp>("9.3833383282371058262139365435233240765713e-02");
	    Bbar(7) = to_num<Real_wp>("1.3138871401731356660181720771852165705253e-01");
	    Bbar(8) = to_num<Real_wp>("4.2452793993460570894743314052986770767710e-02");
	    Bbar(9) = to_num<Real_wp>("4.6614359052634087726314462069090004222681e-02");
	    Bbar(10) = to_num<Real_wp>("1.9739340751760337610363200447178885100858e-02");
	    Bbar(11) = to_num<Real_wp>("2.7040753297272850676247690093320494184034e-03");
	    
	    B(0) = to_num<Real_wp>("1.1445045431083081161076675575316257764110e-02");
	    B(5) = to_num<Real_wp>("1.8963455766836141912201425856209518713321e-01");
	    B(6) = to_num<Real_wp>("9.9015048411147152930074891966312987547087e-02");
	    B(7) = to_num<Real_wp>("2.2377720821386995442668671882695678061627e-01");
	    B(8) = to_num<Real_wp>("1.0466811506175315699681616849938105031123e-01");
	    B(9) = to_num<Real_wp>("1.5358172529459859690396485461784888933400e-01");
	    B(10) = to_num<Real_wp>("1.3940255061073114260364669651015505753606e-01");
	    B(11) = to_num<Real_wp>("6.6309723413523447970758037743751771753948e-02");
	    B(12) = to_num<Real_wp>("1.2166025894932047884961697698182018004080e-02");
	    
	    // embedded eighth order
	    //bh( 1)   =     4.5091769947912279554753974287168407711151e-03
	    //bh( 6)   =     1.3527220707133125439791607750061158575207e-01
	    //bh( 7)   =     1.0821493793737774061748262678958123203921e-01
	    //bh( 8)   =     1.5696421153248798279514458286583792716591e-01
	    //bh( 9)   =    -2.3852476313780400695870370918234366321058e-03
	    //bh(10)   =     8.3696394474180869285163548968221903261701e-02
	    //bh(11)   =     8.1055737581001137792995468759747065645341e-03
	    //bh(12)   =     5.6227458631088512391052566628792410775657e-03
	    //bh(13)   =     0.0000000000000000000000000000000000000000e+00

	    //bph( 1)  =     4.5091769947912279554753974287168407711151e-03
	    //bph( 6)  =     1.6896257767411160899294871388888369845498e-01
	    //bph( 7)  =     1.1419078097647369221167431000919672432731e-01
	    //bph( 8)  =     2.6733660732536693139786660143108056772225e-01
	    //bph( 9)  =    -5.8808702572157774128880948775823679755328e-03
	    //bph(10)  =     2.7575701834209040931496844904292248176983e-01
	    //bph(11)  =     5.7242927727553335856011122467538259863711e-02
	    //bph(12)  =     1.3788178121682857168394350060924379506634e-01
	    //bph(13)  =    -2.0000000000000000000000000000000000000000e-02
	  }
        else if (order == 12)
          {            
	    
            Abar.Reallocate(17, 17); C.Reallocate(17);
            Bbar.Reallocate(17); B.Reallocate(17);
            Bbar.Fill(zero); B.Fill(zero);

	    if (type_pair == SHARP_QURESHI_GRAZIER)
	      {
		C(0) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		C(1) = to_num<Real_wp>("2.2987911522171001779713237615438538119610e-02");
		C(2) = to_num<Real_wp>("4.5975823044342003559426475230877076239221e-02");
		C(3) = to_num<Real_wp>("1.2385385838913258856443744159596778850622e-01");
		C(4) = to_num<Real_wp>("1.7390230229219979999999999999999960000000e-01");
		C(5) = to_num<Real_wp>("1.3449139239361099999999999999999960000000e-01");
		C(6) = to_num<Real_wp>("6.9942050473186690000000000000000090000000e-02");
		C(7) = to_num<Real_wp>("2.4040563275957759999999999999999940000000e-01");
		C(8) = to_num<Real_wp>("2.9925801585104509999999999999999930000000e-01");
		C(9) = to_num<Real_wp>("3.7207169132292630000000000000000050000000e-01");
		C(10) = to_num<Real_wp>("5.6269822875885759999999999999999910000000e-01");
		C(11) = to_num<Real_wp>("6.0511133917257530000000000000000310000000e-01");
		C(12) = to_num<Real_wp>("7.8928236288669310000000000000000130000000e-01");
		C(13) = to_num<Real_wp>("9.0575623053613730000000000000000300000000e-01");
		C(14) = to_num<Real_wp>("9.5911819987361900885877549244129356798106e-01");
		C(15) = to_num<Real_wp>("1.0000000000000000000000000000000000000000e+00");
		C(16) = to_num<Real_wp>("1.0000000000000000000000000000000000000000e+00");
		
		Abar(1, 0) = to_num<Real_wp>("2.6422203807558115202416706164019592284447e-04");
		Abar(2, 0) = to_num<Real_wp>("3.5229605076744153603222274885359456379263e-04");
		Abar(2, 1) = to_num<Real_wp>("7.0459210153488307206444549770718912758525e-04");
		Abar(3, 0) = to_num<Real_wp>("5.5616305161226259562843795370385134036589e-03");
		Abar(3, 1) = to_num<Real_wp>("-9.5580465816756542096284718844433432423919e-03");
		Abar(3, 2) = to_num<Real_wp>("1.1666305184490682696288436830144992260589e-02");
		Abar(4, 0) = to_num<Real_wp>("2.3634420556285561401011018821299321480668e-03");
		Abar(4, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(4, 2) = to_num<Real_wp>("9.0339664428499038047721918341677026225487e-03");
		Abar(4, 3) = to_num<Real_wp>("3.7235968727853599617566678437222956684636e-03");
		Abar(5, 0) = to_num<Real_wp>("1.8624547311177989395250472281344554557062e-03");
		Abar(5, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(5, 2) = to_num<Real_wp>("6.1448860345187159719066099755428932961800e-03");
		Abar(5, 3) = to_num<Real_wp>("1.1457435952158885319853129717780726715124e-03");
		Abar(5, 4) = to_num<Real_wp>("-1.0911704686628006944016051495547521995554e-04");
		Abar(6, 0) = to_num<Real_wp>("9.4180611081097833155341752966915992193089e-04");
		Abar(6, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(6, 2) = to_num<Real_wp>("1.8411090719361789751615510687526024889725e-03");
		Abar(6, 3) = to_num<Real_wp>("-2.3672320089043180505789780821582798360570e-03");
		Abar(6, 4) = to_num<Real_wp>("-1.8988771241998357880073132661059796671045e-04");
		Abar(6, 5) = to_num<Real_wp>("2.2201497507740415659320626069251716866486e-03");
		Abar(7, 0) = to_num<Real_wp>("4.4811515899536082744011946275281854845397e-03");
		Abar(7, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(7, 2) = to_num<Real_wp>("4.4256990629448735507509946140853006366110e-03");
		Abar(7, 3) = to_num<Real_wp>("7.2846553696940604574233642704455287613172e-03");
		Abar(7, 4) = to_num<Real_wp>("7.1782081765624091353293259037536679128887e-03");
		Abar(7, 5) = to_num<Real_wp>("-6.5830701673077061483826142049328170387362e-03");
		Abar(7, 6) = to_num<Real_wp>("1.2110790099419199999999999999999990000000e-02");
		Abar(8, 0) = to_num<Real_wp>("6.9593181374308079774723096132815131405848e-03");
		Abar(8, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(8, 2) = to_num<Real_wp>("-3.1177443259978198952083678705131279450776e-03");
		Abar(8, 3) = to_num<Real_wp>("6.1032418595952454428011365524770693348181e-03");
		Abar(8, 4) = to_num<Real_wp>("1.0927595990740454242922223941675992465677e-02");
		Abar(8, 5) = to_num<Real_wp>("-6.7930937679711936043719211199166804766137e-03");
		Abar(8, 6) = to_num<Real_wp>("2.6472848566893970000000000000000020000000e-02");
		Abar(8, 7) = to_num<Real_wp>("4.2255135648607160000000000000000040000000e-03");
		Abar(9, 0) = to_num<Real_wp>("1.0252050304966763158356291178314013065472e-02");
		Abar(9, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(9, 2) = to_num<Real_wp>("-2.4129973490628836365443719414500318140269e-02");
		Abar(9, 3) = to_num<Real_wp>("-3.8755472247109303714439971656256992895260e-03");
		Abar(9, 4) = to_num<Real_wp>("5.2689010910517151712228074699134632225656e-02");
		Abar(9, 5) = to_num<Real_wp>("-3.7270093015784036846590444612012734704310e-02");
		Abar(9, 6) = to_num<Real_wp>("6.8690833039890067962362284424903866300045e-02");
		Abar(9, 7) = to_num<Real_wp>("-6.1673197410092857949677927520958092021650e-03");
		Abar(9, 8) = to_num<Real_wp>("9.0297109587105816070133012577272357809429e-03");
		Abar(10, 0) = to_num<Real_wp>("4.0839000487334870513260623735454672935169e-02");
		Abar(10, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(10, 2) = to_num<Real_wp>("4.9034416855323231577295962694814848879583e-02");
		Abar(10, 3) = to_num<Real_wp>("-1.7875287829096277825736823040412256255493e-03");
		Abar(10, 4) = to_num<Real_wp>("-2.2982955467709093819887468226913629673054e+00");
		Abar(10, 5) = to_num<Real_wp>("1.7000107425214392401643092768774703860059e+00");
		Abar(10, 6) = to_num<Real_wp>("-3.5456418905719625605713482050341828334589e-01");
		Abar(10, 7) = to_num<Real_wp>("1.6311241524436521595729735363150438472236e+00");
		Abar(10, 8) = to_num<Real_wp>("-7.9116857407379441688711884558508178519577e-01");
		Abar(10, 9) = to_num<Real_wp>("1.8312217470123800000000000000000000000000e-01");
		Abar(11, 0) = to_num<Real_wp>("7.3117582239584326764160210617739254388489e-03");
		Abar(11, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(11, 2) = to_num<Real_wp>("8.6489774963948022669952679177430863730906e-02");
		Abar(11, 3) = to_num<Real_wp>("1.0896419840881710915435825193663319835910e-02");
		Abar(11, 4) = to_num<Real_wp>("-5.4714415238133902124989300656545193564159e-01");
		Abar(11, 5) = to_num<Real_wp>("4.1076518650536265488246959230745741261579e-01");
		Abar(11, 6) = to_num<Real_wp>("-1.1915400866541802235615771541381182674910e-01");
		Abar(11, 7) = to_num<Real_wp>("5.0039262158596160435546622500371348310559e-01");
		Abar(11, 8) = to_num<Real_wp>("-2.8332989352730655341977384803927893730143e-01");
		Abar(11, 9) = to_num<Real_wp>("1.1236602367914300023588119112647321122841e-01");
		Abar(11, 10) = to_num<Real_wp>("4.4861361724219026664208532330773595818204e-03");
		Abar(12, 0) = to_num<Real_wp>("-4.8756368479806916005172964618355930802728e-02");
		Abar(12, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(12, 2) = to_num<Real_wp>("-1.5750049493581786657025299815149867173596e-01");
		Abar(12, 3) = to_num<Real_wp>("-9.2182483482169766903099366449178091083924e-03");
		Abar(12, 4) = to_num<Real_wp>("3.2416897610706500704927852859604584695406e+00");
		Abar(12, 5) = to_num<Real_wp>("-2.9298195952570410950093226062174009023345e+00");
		Abar(12, 6) = to_num<Real_wp>("1.0033057738136136627434303487783021181928e+00");
		Abar(12, 7) = to_num<Real_wp>("-8.2129770978549559091597509462618005882374e-01");
		Abar(12, 8) = to_num<Real_wp>("-4.6351380226962875016681813934827268715430e-01");
		Abar(12, 9) = to_num<Real_wp>("5.1781981712073242097803801494279404536411e-01");
		Abar(12, 10) = to_num<Real_wp>("-1.6009204686453711212618738328112314707084e-01");
		Abar(12, 11) = to_num<Real_wp>("1.3886623811754890000000000000000060000000e-01");
		Abar(13, 0) = to_num<Real_wp>("8.8519019091169314709300861131478052506914e-01");
		Abar(13, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(13, 2) = to_num<Real_wp>("-1.1586538612853209129448116708129654509206e-01");
		Abar(13, 3) = to_num<Real_wp>("-4.7036425888763808510571548962688126862829e-02");
		Abar(13, 4) = to_num<Real_wp>("-1.8362063545324825922941116074179573989096e+01");
		Abar(13, 5) = to_num<Real_wp>("2.1473609214157552885588996395019839611979e+01");
		Abar(13, 6) = to_num<Real_wp>("-6.3447854281975789699695674852994056905042e+00");
		Abar(13, 7) = to_num<Real_wp>("-5.4332043898859568106860458482491681423241e+00");
		Abar(13, 8) = to_num<Real_wp>("1.4562805773075549551724840032731143179332e+01");
		Abar(13, 9) = to_num<Real_wp>("-7.1877097267025992231588009492930964407492e+00");
		Abar(13, 10) = to_num<Real_wp>("2.2259573055575777734126249765480463848056e+00");
		Abar(13, 11) = to_num<Real_wp>("-1.2916906700097048206719895088388832841780e+00");
		Abar(13, 12) = to_num<Real_wp>("4.4990263013104439167203168515950234889611e-02");
		Abar(14, 0) = to_num<Real_wp>("-1.6868346747844423333984017651778594381518e+00");
		Abar(14, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(14, 2) = to_num<Real_wp>("9.4675926633898683010662692429980593442487e-01");
		Abar(14, 3) = to_num<Real_wp>("1.2626744623164166748157170500957771530878e-01");
		Abar(14, 4) = to_num<Real_wp>("2.5108789550496441711136851177723470164796e+01");
		Abar(14, 5) = to_num<Real_wp>("-3.4930811892262862771968512146948082044057e+01");
		Abar(14, 6) = to_num<Real_wp>("1.0704667740708699235509969145651899823899e+01");
		Abar(14, 7) = to_num<Real_wp>("2.0665965404332405335241737335252534324165e+01");
		Abar(14, 8) = to_num<Real_wp>("-3.4854037785719070454211867593326577973023e+01");
		Abar(14, 9) = to_num<Real_wp>("1.6123228230410273670449719393869120295288e+01");
		Abar(14, 10) = to_num<Real_wp>("-4.5375658696847941217969199818363325533831e+00");
		Abar(14, 11) = to_num<Real_wp>("2.8152698664874812842280569734715615441373e+00");
		Abar(14, 12) = to_num<Real_wp>("-2.7094518385659462182588238454815128725755e-02");
		Abar(14, 13) = to_num<Real_wp>("5.3510964953051006741559524410235616442448e-03");
		Abar(15, 0) = to_num<Real_wp>("-1.1899121194130261301945389206176085613039e+00");
		Abar(15, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(15, 2) = to_num<Real_wp>("-5.1418227389890362188340858018430934696529e-01");
		Abar(15, 3) = to_num<Real_wp>("9.1267925296460964589863388314990027522180e-03");
		Abar(15, 4) = to_num<Real_wp>("2.9030883819232492308114409848665340760069e+01");
		Abar(15, 5) = to_num<Real_wp>("-3.3569598167824625327922689651452037204801e+01");
		Abar(15, 6) = to_num<Real_wp>("1.0765669845391384149613078204678650121616e+01");
		Abar(15, 7) = to_num<Real_wp>("8.6355662506415796611718976577301923660219e+00");
		Abar(15, 8) = to_num<Real_wp>("-2.3395963065200365516704305270609713997495e+01");
		Abar(15, 9) = to_num<Real_wp>("1.2151002776771489532848226382586024154885e+01");
		Abar(15, 10) = to_num<Real_wp>("-4.0177980459424927702182618024059579123895e+00");
		Abar(15, 11) = to_num<Real_wp>("2.6203747622601793145100602536608201757836e+00");
		Abar(15, 12) = to_num<Real_wp>("-4.0832925209436904252496980431263952317502e-02");
		Abar(15, 13) = to_num<Real_wp>("1.4424863031084061722555832338595182868311e-02");
		Abar(15, 14) = to_num<Real_wp>("1.2374876309951467364866872097692112758611e-03");
		Abar(16, 0) = to_num<Real_wp>("1.0335383824149195801153977034753233284874e+00");
		Abar(16, 1) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");
		Abar(16, 2) = to_num<Real_wp>("-1.4312314563199904195770210335082868280635e+00");
		Abar(16, 3) = to_num<Real_wp>("-1.2156128454760159934902435859245851269225e-01");
		Abar(16, 4) = to_num<Real_wp>("-3.8726117099919094881739101219309007559289e+00");
		Abar(16, 5) = to_num<Real_wp>("1.2959385668956505697956499297659173683569e+01");
		Abar(16, 6) = to_num<Real_wp>("-3.7810015800145339728936039525971235578005e+00");
		Abar(16, 7) = to_num<Real_wp>("-1.9996342845852889481527024712922749822547e+01");
		Abar(16, 8) = to_num<Real_wp>("2.4408966656424686666529832053542102363298e+01");
		Abar(16, 9) = to_num<Real_wp>("-9.8742991528915224876172458844402012564335e+00");
		Abar(16, 10) = to_num<Real_wp>("2.4650914319163148310301713920391543869935e+00");
		Abar(16, 11) = to_num<Real_wp>("-1.3509150380790346553797321167282400407704e+00");
		Abar(16, 12) = to_num<Real_wp>("5.0662525147977496945068862292087246659488e-02");
		Abar(16, 13) = to_num<Real_wp>("7.9579559461211008765965338037242021138025e-03");
		Abar(16, 14) = to_num<Real_wp>("2.3604468909567310639963379083955631142852e-03");
		Abar(16, 15) = to_num<Real_wp>("0.0000000000000000000000000000000000000000e+00");

		Bbar(0) = to_num<Real_wp>("1.9380954768361851166532694062993208190969e-02");
		Bbar(6) = to_num<Real_wp>("1.0762918742385763743772663992457421533648e-01");
		Bbar(7) = to_num<Real_wp>("2.4491264594176250841899888042955020242794e-01");
		Bbar(8) = to_num<Real_wp>("-2.0374922115399364806432751137961934609925e-01");
		Bbar(9) = to_num<Real_wp>("2.0740727320952827962107356647182199493656e-01");
		Bbar(10) = to_num<Real_wp>("2.7037653506820505081757399074068434420151e-03");
		Bbar(11) = to_num<Real_wp>("7.9679738352602164691660988746774951275038e-02");
		Bbar(12) = to_num<Real_wp>("3.3167284463365381427652957493559180018599e-02");
		Bbar(13) = to_num<Real_wp>("6.6900219258499383396602439405513013549996e-03");
		Bbar(14) = to_num<Real_wp>("2.1783497179838364528458004023874491166387e-03");

		B(0) = to_num<Real_wp>("1.9380954768361851166532694062993208190969e-02");
		B(6) = to_num<Real_wp>("1.1572309819900605895692766053163426860441e-01");
		B(7) = to_num<Real_wp>("3.2242556883553637429514273903799695148412e-01");
		B(8) = to_num<Real_wp>("-2.9076211467683832448784819269422304773424e-01");
		B(9) = to_num<Real_wp>("3.3030406551744133799449147500918406992040e-01");
		B(10) = to_num<Real_wp>("6.1828364952839545598258700474356962806252e-03");
		B(11) = to_num<Real_wp>("2.0177773194511659405537873497173112951798e-01");
		B(12) = to_num<Real_wp>("1.5740155839698742820716683845706162883867e-01");
		B(13) = to_num<Real_wp>("7.0986357654287094725349594264370868642753e-02");
		B(14) = to_num<Real_wp>("5.3284094908975136275468492612906036065700e-02");
		B(15) = to_num<Real_wp>("-1.1704152044157505748435906301090809811391e-02");
		B(16) = to_num<Real_wp>("2.5000000000000000000000000000000000000000e-02");

		// embedded tenth order
		//bh( 1)   =     2.0049392762716821389909832000985406153996e-02
		//bh( 7)   =     1.0512024642076356191992432055317523624786e-01
		//bh( 8)   =     2.7473594011357684776282528527103878081291e-01
		//bh( 9)   =    -2.6662621561251224553826210851278294490370e-01
		//bh(10)   =     2.5409586567177420162813815081470954934107e-01
		//bh(11)   =    -4.7112041111540887585242653806025146012547e-02
		//bh(12)   =     1.2142856389149294758476535487915928917733e-01
		//bh(13)   =     2.7611825649944075458158180164908918900778e-02
		//bh(14)   =     9.4856265764403256941333566886450972124144e-03
		//bh(15)   =     1.2107956373443516856502819461858130698927e-03
		
		//bph( 1)  =     2.0049392762716821389909832000985406153996e-02
		//bph( 7)  =     1.1302548026630567999853968647855055214417e-01
		//bph( 8)  =     3.6168770064960081875142012129331276995984e-01
		//bph( 9)  =    -3.8049128158965312507692344324021843504142e-01
		//bph(10)  =     4.0465744601817073585423548616024469848232e-01
		//bph(11)  =    -1.0773347882362401395663403883833990008068e-01
		//bph(12)  =     3.0750076144769318595812869074465546295734e-01
		//bph(13)  =     1.3103708843838576826000971033691401446771e-01
		//bph(14)  =     1.0064990641187735553289060926986931455985e-01
		//bph(15)  =     2.9616984418526773288423345794026116396883e-02
		//bph(16)  =     2.0000000000000000000000000000000000000000e-02
		//bph(17)  =     0.0000000000000000000000000000000000000000e+00
		
		return;
	      }
            
            // Dormand Prince pair
            // dt <= 8.32/sqrt(|| A ||) for a linear system and stiff problem
            C(0) = zero; C(1) = one/50; C(2) = one/25; C(3) = one/10;
            C(4) = Real_wp(4)/30; C(5) = Real_wp(4)/25; 
            C(6) = one/20; C(7) = one/5; C(8) = one/4; C(9) = one/3;
            C(10) = one/2; C(11) = Real_wp(5)/9;
            C(12) = Real_wp(3)/4; C(13) = Real_wp(6)/7; C(14) = Real_wp(8437)/8926;
            C(15) = one; C(16) = one;

            Abar.Fill(zero);
            Abar(1, 0) = one/5000; Abar(2, 0) = one/3750; Abar(2, 1) = one/1875;
            Abar(3, 0) = Real_wp(7)/2400; Abar(3, 1) = -one/240; Abar(3, 2) = one/160;
            Abar(4, 0) = Real_wp(2)/1215; Abar(4, 2) = Real_wp(4)/729;
            Abar(4, 3) = Real_wp(32)/18225;
            Abar(5, 0) = Real_wp(152)/78125; Abar(5, 2) = Real_wp(1408)/196875;
            Abar(5, 3) = Real_wp(2048)/703125; Abar(5, 4) = Real_wp(432)/546875;
            Abar(6, 0) = Real_wp(29)/51200; Abar(6, 2) = Real_wp(341)/387072;
            Abar(6, 3) = -Real_wp(151)/345600; Abar(6, 4) = Real_wp(243)/716800;
            Abar(6, 5) = -Real_wp(11)/110592;
            Abar(7, 0) = Real_wp(37)/12000; Abar(7, 3) = Real_wp(2)/1125;
            Abar(7, 4) = Real_wp(27)/10000;
            Abar(7, 5) = Real_wp(5)/3168; Abar(7, 6) = Real_wp(224)/20625;

            Abar(8, 0) = to_num<Real_wp>("3.65183937480112971375119150338e-3");
            Abar(8, 2) = to_num<Real_wp>("3.96517171407234306617557289807E-3");
            Abar(8, 3) = to_num<Real_wp>("3.19725826293062822350093426091E-3");
            Abar(8, 4) = to_num<Real_wp>("8.22146730685543536968701883401E-3");
            Abar(8, 5) = to_num<Real_wp>("-1.31309269595723798362013884863E-3");
            Abar(8, 6) = to_num<Real_wp>("9.77158696806486781562609494147E-3");
            Abar(8, 7) = to_num<Real_wp>("3.75576906923283379487932641079E-3");
            
            Abar(9, 0) = to_num<Real_wp>("3.70724106871850081019565530521E-3");
            Abar(9, 2) = to_num<Real_wp>("5.08204585455528598076108163479E-3");
            Abar(9, 3) = to_num<Real_wp>("1.17470800217541204473569104943E-3");
            Abar(9, 4) = to_num<Real_wp>("-2.11476299151269914996229766362E-2");
            Abar(9, 5) = to_num<Real_wp>("6.01046369810788081222573525136E-2");
            Abar(9, 6) = to_num<Real_wp>("2.01057347685061881846748708777E-2");
            Abar(9, 7) = to_num<Real_wp>("-2.83507501229335808430366774368E-2");
            Abar(9, 8) = to_num<Real_wp>("1.48795689185819327555905582479E-2");
            
            Abar(10, 0) = to_num<Real_wp>("3.51253765607334415311308293052E-2");
            Abar(10, 2) = to_num<Real_wp>("-8.61574919513847910340576078545E-3");
            Abar(10, 3) = to_num<Real_wp>("-5.79144805100791652167632252471E-3");
            Abar(10, 4) = to_num<Real_wp>("1.94555482378261584239438810411E0");
            Abar(10, 5) = to_num<Real_wp>("-3.43512386745651359636787167574E0");
            Abar(10, 6) = to_num<Real_wp>("-1.09307011074752217583892572001E-1");
            Abar(10, 7) = to_num<Real_wp>("2.3496383118995166394320161088E0");
            Abar(10, 8) = to_num<Real_wp>("-7.56009408687022978027190729778E-1");
            Abar(10, 9) = to_num<Real_wp>("1.09528972221569264246502018618E-1");

            Abar(11, 0) = to_num<Real_wp>("2.05277925374824966509720571672E-2");
            Abar(11, 2) = to_num<Real_wp>("-7.28644676448017991778247943149E-3");
            Abar(11, 3) = to_num<Real_wp>("-2.11535560796184024069259562549E-3");
            Abar(11, 4) = to_num<Real_wp>("9.27580796872352224256768033235E-1");
            Abar(11, 5) = to_num<Real_wp>("-1.65228248442573667907302673325E0");
            Abar(11, 6) = to_num<Real_wp>("-2.10795630056865698191914366913E-2");
            Abar(11, 7) = to_num<Real_wp>("1.20653643262078715447708832536E0");
            Abar(11, 8) = to_num<Real_wp>("-4.13714477001066141324662463645E-1");
            Abar(11, 9) = to_num<Real_wp>("9.07987398280965375956795739516E-2");
            Abar(11, 10) = to_num<Real_wp>("5.35555260053398504916870658215E-3");
            
            Abar(12, 0) = to_num<Real_wp>("-1.43240788755455150458921091632E-1");
            Abar(12, 2) = to_num<Real_wp>("1.25287037730918172778464480231E-2");
            Abar(12, 3) = to_num<Real_wp>("6.82601916396982712868112411737E-3");
            Abar(12, 4) = to_num<Real_wp>("-4.79955539557438726550216254291E0");
            Abar(12, 5) = to_num<Real_wp>("5.69862504395194143379169794156E0");
            Abar(12, 6) = to_num<Real_wp>("7.55343036952364522249444028716E-1");
            Abar(12, 7) = to_num<Real_wp>("-1.27554878582810837175400796542E-1");
            Abar(12, 8) = to_num<Real_wp>("-1.96059260511173843289133255423E0");
            Abar(12, 9) = to_num<Real_wp>("9.18560905663526240976234285341E-1");
            Abar(12, 10) = to_num<Real_wp>("-2.38800855052844310534827013402E-1");
            Abar(12, 11) = to_num<Real_wp>("1.59110813572342155138740170963E-1");
            
            Abar(13, 0) = to_num<Real_wp>("8.04501920552048948697230778134E-1");
            Abar(13, 2) = to_num<Real_wp>("-1.66585270670112451778516268261E-2");
            Abar(13, 3) = to_num<Real_wp>("-2.1415834042629734811731437191E-2");
            Abar(13, 4) = to_num<Real_wp>("1.68272359289624658702009353564E1");
            Abar(13, 5) = to_num<Real_wp>("-1.11728353571760979267882984241E1");
            Abar(13, 6) = to_num<Real_wp>("-3.37715929722632374148856475521E0");
            Abar(13, 7) = to_num<Real_wp>("-1.52433266553608456461817682939E1");
            Abar(13, 8) = to_num<Real_wp>("1.71798357382154165620247684026E1");
            Abar(13, 9) = to_num<Real_wp>("-5.43771923982399464535413738556E0");
            Abar(13, 10) = to_num<Real_wp>("1.38786716183646557551256778839E0");
            Abar(13, 11) = to_num<Real_wp>("-5.92582773265281165347677029181E-1");
            Abar(13, 12) = to_num<Real_wp>("2.96038731712973527961592794552E-2");
            
            Abar(14, 0) = to_num<Real_wp>("-9.13296766697358082096250482648E-1");
            Abar(14, 2) = to_num<Real_wp>("2.41127257578051783924489946102E-3");
            Abar(14, 3) = to_num<Real_wp>("1.76581226938617419820698839226E-2");
            Abar(14, 4) = to_num<Real_wp>("-1.48516497797203838246128557088E1");
            Abar(14, 5) = to_num<Real_wp>("2.15897086700457560030782161561E0");
            Abar(14, 6) = to_num<Real_wp>("3.99791558311787990115282754337E0");
            Abar(14, 7) = to_num<Real_wp>("2.84341518002322318984542514988E1");
            Abar(14, 8) = to_num<Real_wp>("-2.52593643549415984378843352235E1");
            Abar(14, 9) = to_num<Real_wp>("7.7338785423622373655340014114E0");
            Abar(14, 10) = to_num<Real_wp>("-1.8913028948478674610382580129E0");
            Abar(14, 11) = to_num<Real_wp>("1.00148450702247178036685959248E0");
            Abar(14, 12) = to_num<Real_wp>("4.64119959910905190510518247052E-3");
            Abar(14, 13) = to_num<Real_wp>("1.12187550221489570339750499063E-2");
            
            Abar(15, 0) = to_num<Real_wp>("-2.75196297205593938206065227039E-1");
            Abar(15, 2) = to_num<Real_wp>("3.66118887791549201342293285553E-2");
            Abar(15, 3) = to_num<Real_wp>("9.7895196882315626246509967162E-3");
            Abar(15, 4) = to_num<Real_wp>("-1.2293062345886210304214726509E1");
            Abar(15, 5) = to_num<Real_wp>("1.42072264539379026942929665966E1");
            Abar(15, 6) = to_num<Real_wp>("1.58664769067895368322481964272E0");
            Abar(15, 7) = to_num<Real_wp>("2.45777353275959454390324346975E0");
            Abar(15, 8) = to_num<Real_wp>("-8.93519369440327190552259086374E0");
            Abar(15, 9) = to_num<Real_wp>("4.37367273161340694839327077512E0");
            Abar(15, 10) = to_num<Real_wp>("-1.83471817654494916304344410264E0");
            Abar(15, 11) = to_num<Real_wp>("1.15920852890614912078083198373E0");
            Abar(15, 12) = to_num<Real_wp>("-1.72902531653839221518003422953E-2");
            Abar(15, 13) = to_num<Real_wp>("1.93259779044607666727649875324E-2");
            Abar(15, 14) = to_num<Real_wp>("5.20444293755499311184926401526E-3");
            
            Abar(16, 0) = to_num<Real_wp>("1.30763918474040575879994562983E0");
            Abar(16, 2) = to_num<Real_wp>("1.73641091897458418670879991296E-2");
            Abar(16, 3) = to_num<Real_wp>("-1.8544456454265795024362115588E-2");
            Abar(16, 4) = to_num<Real_wp>("1.48115220328677268968478356223E1");
            Abar(16, 5) = to_num<Real_wp>("9.38317630848247090787922177126E0");
            Abar(16, 6) = to_num<Real_wp>("-5.2284261999445422541474024553E0");
            Abar(16, 7) = to_num<Real_wp>("-4.89512805258476508040093482743E1");
            Abar(16, 8) = to_num<Real_wp>("3.82970960343379225625583875836E1");
            Abar(16, 9) = to_num<Real_wp>("-1.05873813369759797091619037505E1");
            Abar(16, 10) = to_num<Real_wp>("2.43323043762262763585119618787E0");
            Abar(16, 11) = to_num<Real_wp>("-1.04534060425754442848652456513E0");
            Abar(16, 12) = to_num<Real_wp>("7.17732095086725945198184857508E-2");
            Abar(16, 13) = to_num<Real_wp>("2.16221097080827826905505320027E-3");
            Abar(16, 14) = to_num<Real_wp>("7.00959575960251423699282781988E-3");
            
            Bbar(0) = to_num<Real_wp>("1.21278685171854149768890395495E-2");
            Bbar(6) = to_num<Real_wp>("8.62974625156887444363792274411E-2");
            Bbar(7) = to_num<Real_wp>("2.52546958118714719432343449316E-1");
            Bbar(8) = to_num<Real_wp>("-1.97418679932682303358307954886E-1");
            Bbar(9) = to_num<Real_wp>("2.03186919078972590809261561009E-1");
            Bbar(10) = to_num<Real_wp>("-2.07758080777149166121933554691E-2");
            Bbar(11) = to_num<Real_wp>("1.09678048745020136250111237823E-1");
            Bbar(12) = to_num<Real_wp>("3.80651325264665057344878719105E-2");
            Bbar(13) = to_num<Real_wp>("1.16340688043242296440927709215E-2");
            Bbar(14) = to_num<Real_wp>("4.65802970402487868693615238455E-3");
            
            B(0) = to_num<Real_wp>("1.21278685171854149768890395495E-2");
            B(6) = to_num<Real_wp>("9.08394342270407836172412920433E-2");
            B(7) = to_num<Real_wp>("3.15683697648393399290429311645E-1");
            B(8) = to_num<Real_wp>("-2.63224906576909737811077273181E-1");
            B(9) = to_num<Real_wp>("3.04780378618458886213892341513E-1");
            B(10) = to_num<Real_wp>("-4.15516161554298332243867109382E-2");
            B(11) = to_num<Real_wp>("2.46775609676295306562750285101E-1");
            B(12) = to_num<Real_wp>("1.52260530105866022937951487642E-1");
            B(13) = to_num<Real_wp>("8.14384816302696075086493964505E-2");
            B(14) = to_num<Real_wp>("8.50257119389081128008018326881E-2");
            B(15) = to_num<Real_wp>("-9.15518963007796287314100251351E-3");
            B(16) = to_num<Real_wp>("2.5E-2");
            
            // embedded 10th order
            //Bbar(0) = Real_wp("1.70087019070069917527544646189E-2");
            //Bbar(6) = Real_wp("7.22593359308314069488600038463E-2");
            //Bbar(7) = Real_wp("3.72026177326753045388210502067E-1");
            //Bbar(8) = Real_wp("-4.01821145009303521439340233863E-1");
            //Bbar(9) = Real_wp("3.35455068301351666696584034896E-1");
            //Bbar(10) = Real_wp("-1.31306501075331808430281840783E-1");
            //Bbar(11) = Real_wp("1.89431906616048652722659836455E-1");
            //Bbar(12) = Real_wp("2.68408020400290479053691655806E-2");
            //Bbar(13) = Real_wp("1.63056656059179238935180933102E-2");
            //Bbar(14) = Real_wp("3.79998835669659456166597387323E-3");
            
            //B(0) = Real_wp("1.70087019070069917527544646189E-2");
            //B(6) = Real_wp("7.60624588745593757356421093119E-2");
            //B(7) = Real_wp("4.65032721658441306735263127583E-1");
            //B(8) = Real_wp("-5.35761526679071361919120311817E-1");
            //B(9) = Real_wp("5.03182602452027500044876052344E-1");
            //B(10) = Real_wp("-2.62613002150663616860563681567E-1");
            //B(11) = Real_wp("4.26221789886109468625984632024E-1");
            //B(12) = Real_wp("1.07363208160116191621476662322E-1");
            //B(13) = Real_wp("1.14139659241425467254626653171E-1");
            //B(14) = Real_wp("6.93633866500486770090602920091E-2");
            //B(15) = Real_wp("2.0E-2");
            //B(16) = zero;
          }
	else
	  {
	    cout << "Order not implemented" << endl;
	    abort();
	  }
      }
    else
      {
        // f does depend on y', we use formula 14.4, p283 of Hairer's book (vol I)
        RungeKutta_Iterator<Real_wp> rk;
        rk.SetOrder(order);
        rk.GetCoefficients(A, B, C);
        
        // Abar = A^2
        // Bbar = A^T b
        int N = A.GetM();
        Abar.Reallocate(N, N); Bbar.Reallocate(N);
        Mlt(A, A, Abar);
        Mlt(SeldonTrans, A, B, Bbar);
      }
  }
  
  
  //! fills arrays A, \bar{A}, b, \bar{b} and C
  template<class T>
  void RungeKuttaNystrom_Iterator<T>
  ::GetCoefficients(Matrix<Real_wp>& coefA, Matrix<Real_wp>& coefAbar,
		    Vector<Real_wp>& coefB, Vector<Real_wp>& coefBbar,
		    Vector<Real_wp>& coefC) const
  {
    coefA = A; coefAbar = Abar;
    coefB = B; coefBbar = Bbar; coefC = C;    
  }
  
  
  //! inits scheme (needed before calling Advance)
  template<class T>
  void RungeKuttaNystrom_Iterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                         Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys)
  {
    T zero; SetComplexZero(zero);
    int N = u0.GetM();
    Y = u0; Yprime = du0_dt;
    u0.Clear(); du0_dt.Clear();
    ki.Reallocate(Abar.GetM());
    for (int i = 0; i < ki.GetM(); i++)
      {
        ki(i).Reallocate(N);
        ki(i).Fill(zero);
      }

    Ytilde.Reallocate(N); Ytilde.Fill(zero);
    if (f_dependent_from_yprime)
      {
        Yprime_tilde.Reallocate(N);
        Yprime_tilde.Fill(zero);
      }
    
    dt = deltat;
  }
    

  //! releases memory used by the scheme
  template<class T>
  void RungeKuttaNystrom_Iterator<T>::Clear()
  {
    Y.Clear(); Ytilde.Clear(); Yprime.Clear(); Yprime_tilde.Clear();
    ki.Clear();
  }
  
  
  //! releases memory used by the scheme
  template<class T>
  void RungeKuttaNystrom_Iterator<T>::ClearSecond(Vector<T>& Uf, Vector<T>& dUf)
  {
    Ytilde.Clear(); Yprime_tilde.Clear(); ki.Clear();
    Uf = Y; dUf = Yprime; Y.Clear(); Yprime.Clear();
  }
    

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> RungeKuttaNystrom_Iterator<T>::GetStabilityFunction() const
  {
    UnivariatePolynomial<Real_wp> Pol;
    Pol.SetOrder(0);
    if (!f_dependent_from_yprime)
      {
        switch (order_rk)
          {
          case 2 : Pol(0) = 2.0; return Pol;
          case 3 : Pol(0) = 2.498; return Pol;
          case 4 : Pol(0) = 3.94; return Pol;
          case 5 : Pol(0) = 2.908; return Pol;
          case 6 : Pol(0) = 3.0895; return Pol;
          case 7 : Pol(0) = 7.0875; return Pol;
          case 8 : Pol(0) = 7.8525; return Pol;
          case 10 : Pol(0) = 4.7527; return Pol;
          }
      }    
    
    cout << "Case not implemented" << endl;
    abort(); 
    return Pol;
  }
  
  
  //! computes y^{n+1} and yprime^{n+1} from y^n and yprime^n, values are overwritten
  template<class T>
  void RungeKuttaNystrom_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    // evaluating k_i = f(t^n + c_i dt, y^n + c_i dt yprime^n 
    //        + h^2 \sum abar_{i,j} k_j, yprime^n + dt \sum a_{i,j} k_j)
    for (int i = 0; i < B.GetM(); i++)
      {
        Copy(Y, Ytilde);
        Add(dt*C(i), Yprime, Ytilde);
        for (int j = 0; j < i; j++)
          if (Abar(i, j) != Real_wp(0))
            Add(dt*dt*Abar(i, j), ki(j), Ytilde);
        
        if (f_dependent_from_yprime)
          {
            Copy(Yprime, Yprime_tilde);
            for (int j = 0; j < i; j++)
              if (A(i, j) != Real_wp(0))
                Add(dt*A(i, j), ki(j), Yprime_tilde);
          }
        
        sys.EvaluateFunctionS(t + C(i)*dt, Ytilde, Yprime_tilde, ki(i));
      }
    
    // y^{n+1} = y^n + dt yprime^n + dt^2 \sum bbar_i k_i
    Add(dt, Yprime, Y);
    for (int i = 0; i < Bbar.GetM(); i++)
      Add(dt*dt*Bbar(i), ki(i), Y);
    
    // yprime^{n+1} = yprime^n + dt \sum b_i k_i
    for (int i = 0; i < B.GetM(); i++)
      Add(dt*B(i), ki(i), Yprime);

    // setting Dirichlet condition
    sys.SetDirichletCondition(t+dt, 0, Y);
    sys.SetDirichletCondition(t+dt, 1, Yprime);
  }
  
  
  //! changes time step during the simulation
  template<class T>
  void RungeKuttaNystrom_Iterator<T>::ChangeTimeStep(const Real_wp& dt_)
  {
    dt = dt_;
  }

    
  //! returns the current iterate y^n
  template<class T> 
  Vector<T>& RungeKuttaNystrom_Iterator<T>::GetIterate()
  {
    return Y;
  }
  
  
  //! returns the current iterate y^n
  template<class T>
  const Vector<T>& RungeKuttaNystrom_Iterator<T>::GetIterate() const
  {
    return Y;
  }

  
  //! returns the current derivative dy/dt(t^n)
  template<class T>
  Vector<T>& RungeKuttaNystrom_Iterator<T>::GetDeriveIterate()
  {
    return Yprime;
  }
  
  
  //! returns the current derivative dy/dt(t^n)
  template<class T>
  const Vector<T>& RungeKuttaNystrom_Iterator<T>::GetDeriveIterate() const
  {
    return Yprime;
  }


  //! returns the number of vectors needed to advance to time t^{n+1}
  /*!
    Intermediary vectors do not count
   */
  template<class T>
  int RungeKuttaNystrom_Iterator<T>::GetNumberOfIterates() const
  {
    return 2;
  }
  
  
  //! returns the k-th vector needed to advance to time t^{n+1}
  template<class T>
  Vector<T>& RungeKuttaNystrom_Iterator<T>::GetIterate(int k)
  {
    if (k == 0)
      return Y;
    else
      return Yprime;
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t RungeKuttaNystrom_Iterator<T>::GetMemorySize() const
  {
    size_t taille = A.GetMemorySize() + Abar.GetMemorySize() + B.GetMemorySize() + Bbar.GetMemorySize()
      + C.GetMemorySize();
    
    taille += Y.GetMemorySize() + Yprime.GetMemorySize() + Ytilde.GetMemorySize()
      + Yprime_tilde.GetMemorySize();

    taille += Seldon::GetMemorySize(ki);
    return taille;
  }

  
  /*******************************
   * SymmetricMultistep_Iterator *
   *******************************/
  
  
  //! default constructor
  template<class T>
  SymmetricMultistep_Iterator<T>::SymmetricMultistep_Iterator()
  {
    dt = 0.0;
  }
    
  
  //! Retrieves cofficients alpha and beta stored in the class
  template<class T>
  void SymmetricMultistep_Iterator<T>::GetAlphaBeta(VectReal_wp& A, VectReal_wp& B)
  {
    A = alpha;
    B = beta;
  }


  //! finds coefficients alpha and beta from parameters theta
  template<class T>
  void SymmetricMultistep_Iterator<T>::FindCoefAlphaBeta(const VectReal_wp& theta)
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

    // symmetric coefficients : 
    // alpha = (alpha_N, alpha_N-1, ..., alpha_0, alpha_1, ..., alpha_N)
    // beta = (beta_N, beta_N-1, ..., beta_0, beta_1, ..., beta_N)
    // taille du systeme : on resout que la moitie des inconnues beta_0, beta_N-1
    int N = r/2;

    // on explicite le systeme lineaire que doit resoudre beta_0, ..., beta_N-1
    Matrix<Real_wp> sys(N, N);
    Vector<Real_wp> powI(N+1), rhs(N);
    powI.Zero();
    rhs.Zero();
    sys.Zero();
    
    // premiere equation : \alpha_0 + 2 \sum_{i=1}^N alpha_i = 0
    // deja verifiee par 1 est racine double de rho

    // deuxieme equation : 2 \sum_{i=1}^N i^2/2 alpha_i = \beta_0 + 2 \sum_{i=1}^N \beta_i
    sys(0, 0) = one;  
    
    // terme N^2 alpha_N
    rhs(0) = N*N*rho(0); powI(N) = N*N;
    for (int i = 1; i < N; i++)
      {
        // terme 2 beta_i
        sys(0, i) = Real_wp(2);

        // terme i^2 alpha_i
        rhs(0) += i*i*rho(N-i);

        // on initialise powI(i) = i^2 pour la boucle suivante
        powI(i) = i*i;
      }
    
    // equation pour obtenir un ordre s
    // 1/[s(s-1)] \sum_{i=1}^N i^s alpha_i = \sum_{i=1}^N i^{s-2} \beta_i 
    int num = 1;
    for (int s = 4; s <= r; s += 2)
      {
        // coef contient 1 / (s(s-1))
        Real_wp coef = one / (s*(s-1));
        
	// powI(i) contient i^{s-2}
        
        // terme N^s alpha_N / (s (s-1))
        rhs(num) = powI(N) * N * N * rho(0) * coef;
        for (int i = 1; i < N; i++)
          {
	    // terme beta_i i^{s-2} 
            sys(num, i) = powI(i);
	    
	    // terme alpha_i i^s / [s(s-1)]
            rhs(num) += rho(N-i) * powI(i) * i * i * coef;
          }
	
	// on multiplie row(i) par i*i pour le s+2 suivant
        for (int i = 1; i <= N; i++)
          powI(i) *= i*i;        
        
        num++;
      }
    
    // on resout le systeme lineaire
    IVect pivot(N);
    GetLU(sys, pivot);
    
    SolveLU(sys, pivot, rhs);
    
    // on remplit alpha et beta
    for (int i = 0; i < alpha.GetM(); i++)
      alpha(i) = rho(i);    
    
    beta(0) = 0; beta(2*N) = 0;
    beta(N) = rhs(0);
    for (int i = 1; i < N; i++)
      {
        beta(N+i) = rhs(i);
        beta(N-i) = rhs(i);
      }
  }

  
  //! inits the time scheme with an order and parameters theta
  template<class T>
  void SymmetricMultistep_Iterator<T>::SetOrder(int order, const VectReal_wp& theta)
  {
    if (order == 2)
      {
        alpha.Reallocate(3);
        beta.Reallocate(3);
        alpha(0) = 1.0; alpha(1) = -2.0; alpha(2) = 1.0;
        beta(0) = 0.0; beta(1) = 1.0; beta(2) = 0.0;
      }
    else if (order%2 == 0)
      FindCoefAlphaBeta(theta);
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
  void SymmetricMultistep_Iterator<T>
  ::SetInitialConditionS(const Real_wp& t, const Real_wp& deltat,
                         Vector<T>& u0, Vector<T>& du0_dt, VirtualOdeSystem<T>& sys)
  {
    int N = u0.GetM();
    if (beta(0) != 0)
      {
        cout << "Only explicit scheme is implemented" << endl;
        abort();
      }
    
    int order = alpha.GetM();
    dt = deltat;
    yk.Reallocate(order-1);
    fk.Reallocate(order-2);
    for (int i = 0; i < yk.GetM(); i++)
      {
        yk(i).Reallocate(N);
        yk(i).Fill(0);
      }
    
    yk(0) = u0;
    
    for (int i = 0; i < fk.GetM(); i++)
      {
        fk(i).Reallocate(N);
        fk(i).Fill(0);
      }
        
    Y_star.Reallocate(N);
    F_star.Reallocate(N);
    Y_star.Fill(0);
    F_star.Fill(0);

    // initializing RKN scheme
    rk_nystrom.SetInitialConditionS(t, deltat, u0, du0_dt, sys);
    u0.Clear(); du0_dt.Clear();
  }
  
  
  template<class T>
  void SymmetricMultistep_Iterator<T>::Clear()
  {
    yk.Clear();
    fk.Clear();
  }
  
  
  template<class T>
  void SymmetricMultistep_Iterator<T>::ClearSecond(Vector<T>& y0, Vector<T>& dy0)
  {
    y0 = yk(0);
    dy0.Clear();
    
    yk.Clear();
    fk.Clear();
  }
    

  //! returns the stability function of scheme
  template<class T>
  UnivariatePolynomial<Real_wp> SymmetricMultistep_Iterator<T>::GetStabilityFunction() const
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  template<class T>
  void SymmetricMultistep_Iterator<T>::Advance(const Real_wp& t, int n_time, VirtualOdeSystem<T>& sys)
  {
    int order = beta.GetM();
    
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
        // on evalue Y_star = \sum \beta_k fk
        Y_star.Zero();
        for (int k = 0; k < order-2; k++)
          Add(dt*dt*beta(k+1), fk(k), Y_star);
        
        // on retranche la partie -\sum \alpha_k yk
        for (int k = 0; k < order-1; k++)
          Add(-alpha(k+1), yk(k), Y_star);        
      }

    sys.SetDirichletCondition(t+dt, 0, Y_star);
    
    // on decale les y_{n-k} stockes
    for (int k = order-3; k >= 0; k--)
      Copy(yk(k), yk(k+1));
    
    // on decale les f_{n-k} stockes
    for (int k = order-4; k >= 0; k--)
      Copy(fk(k), fk(k+1));
    
    // on calcule le nouveau f_{n+1}
    Copy(Y_star, yk(0));
    Vector<T> yp_null;
    sys.EvaluateFunctionS(t+dt, yk(0), yp_null, fk(0));
  }
    

  template<class T>
  void SymmetricMultistep_Iterator<T>::ChangeTimeStep(const Real_wp&)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  template<class T>
  Vector<T>& SymmetricMultistep_Iterator<T>::GetIterate()
  {
    return yk(0);
  }
  
  
  template<class T>
  const Vector<T>& SymmetricMultistep_Iterator<T>::GetIterate() const
  {
    return yk(0);
  }
  
    
  template<class T>
  int SymmetricMultistep_Iterator<T>::GetNumberOfIterates() const
  {
    return alpha.GetM();
  }
  
  
  template<class T>
  Vector<T>& SymmetricMultistep_Iterator<T>::GetIterate(int k)
  {
    return yk(k);
  }


  //! returns the size used to store the object in bytes
  template<class T>
  size_t SymmetricMultistep_Iterator<T>::GetMemorySize() const
  {
    size_t taille = alpha.GetMemorySize() + beta.GetMemorySize() + rk_nystrom.GetMemorySize();
    taille += Y_star.GetMemorySize() + F_star.GetMemorySize();
    taille += Seldon::GetMemorySize(yk) + Seldon::GetMemorySize(fk);
    return taille;
  }
  
}
  
#define MONTJOIE_FILE_EXPLICIT_TIME_SCHEMES_CXX
#endif
