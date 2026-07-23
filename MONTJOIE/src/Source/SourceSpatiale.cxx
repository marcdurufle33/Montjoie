#ifndef MONTJOIE_FILE_SOURCE_SPATIALE_CXX

namespace Montjoie
{


  template<class T, class Dimension>
  void VirtualSourceField<T, Dimension>::EvaluateGradient(const R_N& x, Vector<T>& df) const
  {
    cout << "Not implemented for this source" << endl;
    abort();
  }
  
    
  /******************
   * GaussianSource *
   ******************/
  
  
  //! default constructor
  template<class Dimension>
  GaussianSource<Dimension>::GaussianSource()
  {
    this->origin.Fill(Real_wp(0));
    RadiusSource = 1.0;
    RadiusSource2 = 1.0;
    RadiusSourceCutOff = 1.0;
    alpha = -log(1.e-6)/RadiusSource2;
    beta = sqrt(alpha/pi_wp);
  }
  
  
  //! constructor with the origin and the radius
  template<class Dimension>
  GaussianSource<Dimension>::GaussianSource(const R_N& pt, const Real_wp& radius)
  {
    origin = pt;
    RadiusSource = radius;
    RadiusSourceCutOff = radius;
    RadiusSource2 = radius*radius;
    alpha = -log(1.e-6)/RadiusSource2;
    beta = sqrt(alpha/pi_wp);
  }

  
  //! returns the radius of the gaussian
  template<class Dimension>
  Real_wp GaussianSource<Dimension>::GetRadius() const
  {
    return RadiusSource;
  }


  //! returns the radius of the gaussian
  template<class Dimension>
  Real_wp GaussianSource<Dimension>::GetCutOffRadius() const
  {
    return RadiusSourceCutOff;
  }


  //! returns the center of the gaussian
  template<class Dimension>
  const typename Dimension::R_N& GaussianSource<Dimension>::GetCenter() const
  {
    return origin;
  }

  
  //! initializing origin and radius
  template<class Dimension>
  void GaussianSource<Dimension>::Init(const R_N& pt, const Real_wp& radius,
				       const Real_wp& radius_cut_off)
  {
    origin = pt;
    RadiusSource = radius;
    RadiusSourceCutOff = radius_cut_off;
    RadiusSource2 = radius*radius;
    alpha = -log(1.e-6)/RadiusSource2;
    beta = sqrt(alpha/pi_wp);
  }

  
  //! returns only gaussian amplitude
  template<class Dimension>
  Real_wp GaussianSource<Dimension>::GetAmplitude(const R_N & x) const
  {
    Real_wp r = Distance(x, this->origin);
    Real_wp aj(0);
    if (r < RadiusSourceCutOff)
      {
	Real_wp t2 = r*r;
	aj = beta*exp(-alpha*t2);
      }
    
    return aj;
  }

  
  //! evaluation of gradient of the amplitude
  template<class Dimension>
  void GaussianSource<Dimension>
  ::GetGradAmplitude(const R_N & x, Real_wp& aj, R_N& grad_aj) const
  {
    Real_wp r = Distance(x, this->origin);
    aj = Real_wp(0);
    grad_aj.Fill(0);
    if (r < RadiusSourceCutOff)
      {
	Real_wp t2 = r*r;
        aj = beta*exp(-alpha*t2);
	Real_wp coef_aj = -2.0*alpha*aj;
        for (int j = 0; j < Dimension::dim_N; j++)
          grad_aj(j) = coef_aj*(x(j)-this->origin(j));
      }
  }


  //! evaluation of gradient and hessian of the amplitude
  template<class Dimension>
  void GaussianSource<Dimension>
  ::GetHessianAmplitude(const R_N & x, Real_wp& aj, R_N& grad_aj,
                        TinyMatrix<Real_wp, Symmetric, dim_N, dim_N>& hess_aj) const
  {
    Real_wp r = Distance(x, this->origin);
    aj = Real_wp(0);
    grad_aj.Fill(0);
    hess_aj.Fill(0);
    if (r < RadiusSourceCutOff)
      {
	Real_wp t2 = r*r;
        aj = beta*exp(-alpha*t2);
	Real_wp coef_aj = -2.0*alpha*aj;
        for (int j = 0; j < Dimension::dim_N; j++)
          {
            grad_aj(j) = coef_aj*(x(j)-this->origin(j));
            hess_aj(j, j) = coef_aj;
            for (int k = j; k < Dimension::dim_N; k++)
              hess_aj(j, k) -= 2.0*alpha*grad_aj(j)*(x(k)-this->origin(k));
          }
      }
  }
  

  //! sets the polarization
  template<class T, class Dimension>
  void GaussianSourceField<T, Dimension>
  ::SetPolarization(const Vector<T>& P)
  {
    polar = P;
  }
  

  //! evaluation of gaussian
  template<class T, class Dimension>
  void GaussianSourceField<T, Dimension>
  ::EvaluateFunction(const R_N & x, Vector<T>& fvec) const
  {
    Real_wp f = this->GetAmplitude(x);
    for (int i = 0; i < polar.GetM(); i++)
      fvec(i) = f*polar(i);
  }


  //! evaluation of gradient of the gaussian
  template<class T, class Dimension>
  void GaussianSourceField<T, Dimension>
  ::EvaluateGradient(const R_N & x, Vector<T>& fvec) const
  {
    Real_wp f; R_N df;
    TinyVector<T, Dimension::dim_N> vec_u;
    this->GetGradAmplitude(x, f, df);
    for (int i = 0; i < polar.GetM(); i++)
      {
	vec_u = df;
	vec_u *= polar(i);
	CopyVector(vec_u, i, fvec);
      }
  }

  
  /**********************
   * UniformSourceField *
   **********************/
  

  //! sets the polarization
  template<class T, class Dimension>
  void UniformSourceField<T, Dimension>
  ::SetPolarization(const Vector<T>& P)
  {
    polar = P;
  }
  

  //! evaluation of uniform function
  template<class T, class Dimension>
  void UniformSourceField<T, Dimension>
  ::EvaluateFunction(const R_N & x, Vector<T>& fvec) const
  {
    for (int i = 0; i < polar.GetM(); i++)
      fvec(i) = polar(i);
  }


  //! evaluation of gradient of the uniform field
  template<class T, class Dimension>
  void UniformSourceField<T, Dimension>
  ::EvaluateGradient(const R_N & x, Vector<T>& fvec) const
  {
    fvec.Fill(0);
  }


  /**************************
   * PlaneWaveIncidentField *
   **************************/
  
  
  //! default constructor
  template<class Dimension>
  PlaneWaveIncidentField<Complex_wp, Dimension>::PlaneWaveIncidentField(const R_N& pt, const R_N& u)
  {
    origin = pt;
    k_wave = u;
  }
  

  //! kind of constructor
  template<class Dimension>
  void PlaneWaveIncidentField<Complex_wp, Dimension>::Init(const R_N& pt, const R_N& u)
  {
    origin = pt;
    k_wave = u;
  }
  
  
  //! evaluation of a plane wave
  template<class Dimension>
  void PlaneWaveIncidentField<Complex_wp, Dimension>
  ::EvaluateFunction(const R_N& pt, Complex_wp& f) const
  {
    R_N diff; diff = pt - this->origin;
    Real_wp arg = DotProd(k_wave, diff);
    f = exp(Iwp*arg);
  }
  
  
  //! gradient of plane wave
  template<class Dimension>
  void PlaneWaveIncidentField<Complex_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& pt, Complex_wp& uinc,
			     TinyVector<Complex_wp, Dimension::dim_N>& grad_uinc) const
  {
    R_N diff; diff = pt - this->origin;
    Real_wp arg = DotProd(k_wave, diff);
    uinc = exp(Iwp*arg);
    for (int j = 0; j < Dimension::dim_N; j++)
      grad_uinc(j) = Iwp*k_wave(j)*uinc;
  }


  /**************************
   * PlaneWaveIncidentFieldComplex *
   **************************/
  

  //! default constructor
  template<class Dimension>
  PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>
  ::PlaneWaveIncidentFieldComplex(const R_N& pt, const TinyVector<Complex_wp, Dimension::dim_N>& u)
  {
    origin = pt;
    k_wave = u;
  }
  

  //! kind of constructor
  template<class Dimension>
  void PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>::Init(const R_N& pt, const TinyVector<Complex_wp, Dimension::dim_N>& u)
  {
    origin = pt;
    k_wave = u;
  }
  
  
  //! evaluation of a plane wave
  template<class Dimension>
  void PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>
  ::EvaluateFunction(const R_N& pt, Complex_wp& f) const
  {
    R_N diff; diff = pt - this->origin;
    Complex_wp arg = DotProd(k_wave, diff);
    f = exp(Iwp*arg);
  }
  
  
  //! gradient of plane wave
  template<class Dimension>
  void PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& pt, Complex_wp& uinc,
			     TinyVector<Complex_wp, Dimension::dim_N>& grad_uinc) const
  {
    R_N diff; diff = pt - this->origin;
    Complex_wp arg = DotProd(k_wave, diff);
    uinc = exp(Iwp*arg);
    for (int j = 0; j < Dimension::dim_N; j++)
      grad_uinc(j) = Iwp*k_wave(j)*uinc;
  }



  /*********************************
   * LayeredPlaneWaveIncidentField *
   *********************************/


  //! plane wave incident field in layered medium
  /*!
    \param[in] kwave wave vector (kx, ky)
    \param[in] omega pulsation
    \param[in] a_infty coefficient a at infinity
    \param[in] b_infty coefficient b at infinity
    \param[in] parameters definition of d, a, b
   */
  template<class Dimension>
  LayeredPlaneWaveIncidentField<Complex_wp, Dimension>
  ::LayeredPlaneWaveIncidentField(const R_N& ptO, const R_N& kwave, const Real_wp& omega,
                                  const Real_wp& a_infty, const Real_wp& b_infty, 
                                  const VectString& parameters)
  {
    num_layer = -1;
    Init(ptO, kwave, omega, a_infty, b_infty, parameters);
  }

  
  //! this method is called before EvaluateFunction or EvaluateFunctionGradient
  template<class Dimension>
  void LayeredPlaneWaveIncidentField<Complex_wp, Dimension>
  ::InitElement(int num_elem, const typename Dimension::VectR_N& s)
  {
    // center of the element
    R_N center;
    for (int k = 0; k < s.GetM(); k++)
      center += s(k);
    
    center /= Real_wp(s.GetM());

    // then computing the layer where the element i is
    int i = 0;
    for (i = 0; i < d.GetM(); i++)
      if( center(1) <= d(i) ) 
	break;
    
    num_layer = i;
  }
  

  //! plane wave incident field in layered medium
  /*!
    \param[in] kwave wave vector (kx, ky)
    \param[in] omega pulsation
    \param[in] a_infty coefficient a at infinity
    \param[in] a_infty coefficient b at infinity
    \param[in] parameters definition of d, a, b
   */
  template<class Dimension>
  void LayeredPlaneWaveIncidentField<Complex_wp, Dimension>
  ::Init(const R_N& ptO, const R_N& kwave, const Real_wp& omega,
         const Real_wp& a_infty, const Real_wp& b_infty, 
         const VectString& parameters)
  {  
    origin = ptO;
    kx = kwave(0);
    
    int N = to_num<int>(parameters(2));
    int nb = 3;
    a.Reallocate(N+2);
    b.Reallocate(N+2);
    d.Reallocate(N+1);
    b(0) = to_num<Complex_wp>(parameters(nb++));
    a(0) = to_num<Complex_wp>(parameters(nb++));        
    for (int i = 0; i < N; i++)
      {
	d(i) = to_num<Real_wp>(parameters(nb++));
	b(i+1) = to_num<Complex_wp>(parameters(nb++));
	a(i+1) = to_num<Complex_wp>(parameters(nb++));
      }
    
    d(N) = to_num<Real_wp>(parameters(nb++));
    a(N+1) = a_infty;
    b(N+1) = b_infty;
    
    //DISP(a); DISP(b); DISP(d); DISP(kwave);
    Matrix<Complex_wp, General, BandedCol> M;
    M.Reallocate(2*N+4, 2*N+4, 2, 2);
    M.Zero();

    ky.Reallocate(N+2);
    Real_wp signe = sign(kwave(1)); // signe des ky
    for (int i=0; i<N+2; i++)
      {
	Complex_wp epsilon = b(i)/a(i);
	ky(i) = signe * sqrt(omega*omega*epsilon - kx*kx);
      }
    
    //DISP(ky);
    
    // the unknowns are ordered as follows
    // A : 2*i 
    // B : 2*i+1
    // A(0) and B(0) correspond to the lower infinite medium 
    // A(i) and B(i) correspond to the layer between d(i-1) and d(i)
    // A(N+1) and B(N+1) correspond to the upper infinite medium
    for (int i = 0; i <= N; i++)
      {
        // continuity of u at interface
	M.Get(2*i+1,2*i)       = exp( Iwp * ky(i) * (d(i) - origin(1)) ); // *A(i)
	M.Get(2*i+1,2*i+1)     = Real_wp(1) / M(2*i+1, 2*i); // *B(i)

	M.Get(2*i+1,2*i+2)     = - exp(Iwp * ky(i+1) * (d(i) - origin(1)) ); // *A(i+1)
	M.Get(2*i+1, 2*i+3)    = Real_wp(1) / M(2*i+1,2*i+2);	// *B(i+1)
	
        // continuity of (a du/dy) at interface
	M.Get(2*i+2, 2*i)    = Iwp * ky(i) * a(i) * M(2*i+1, 2*i) ; 
	M.Get(2*i+2, 2*i+1)  = -Iwp * ky(i) * a(i) * M(2*i+1, 2*i+1) ; 

	M.Get(2*i+2, 2*i+2)  = Iwp * ky(i+1)  * a(i+1) * M(2*i+1,2*i+2) ;
	M.Get(2*i+2, 2*i+3)  = -Iwp * ky(i+1) * a(i+1) * M(2*i+1,2*i+3) ;
      } 
    
    // the plane wave is supposed to come from above, we set 
    // A(N+1) = 1
    // and no reflected wave from low infty 
    // B(0) = 0
    M.Get(0, 0) = 0.0; 
    M.Get(0, 1) = 1.0; // * B(0)

    M.Get(2*N+3, 2*N+2) = 1.0; // *A(N+1)
    
    Vector<int> pivot; 
    M.Factorize(pivot);
    
    VectComplex_wp rhs(2*N+4);
    rhs.Zero();
    rhs(2*N+3) = 1.0;

    M.Solve(pivot,rhs);
    
    A.Reallocate(N+2);
    B.Reallocate(N+2);

    for (int i = 0; i < N+2; i++)
      {
	A(i) = rhs(2*i);
	B(i) = rhs(2*i+1);
      }
    
    //DISP(A); DISP(B);
  }

  template<class Dimension>
  void LayeredPlaneWaveIncidentField<Complex_wp, Dimension>
  ::EvaluateFunction(const R_N& ptM, Complex_wp& f) const
  {
    int i = num_layer;
    R_N pt(ptM - origin);
    
    // use coefs A(i) and B(i) , and ky(i)
    f = A(i)*exp (Iwp*(kx*pt(0)+ky(i)*pt(1))) + B(i)*exp (Iwp*(kx*pt(0)-ky(i)*pt(1)) );
    
  }

  template<class Dimension>
  void LayeredPlaneWaveIncidentField<Complex_wp, Dimension>    
  ::EvaluateFunctionGradient(const R_N& ptM, Complex_wp& f,
                             TinyVector<Complex_wp, Dimension::dim_N>& df) const
  {
    int i = num_layer;
    R_N pt(ptM - origin);
      
    // use coefs A(i) and B(i) , and ky(i)
    Complex_wp tmp1 = exp (Iwp*(kx*pt(0)+ky(i)*pt(1)));
    Complex_wp tmp2 = exp (Iwp*(kx*pt(0)-ky(i)*pt(1)));
    f = A(i)*tmp1 + B(i)*tmp2;
    df(0) = Iwp * kx * A(i) * tmp1 + Iwp * kx * B(i)*tmp2;
    df(1) = Iwp * ky(i) * A(i) * tmp1 - Iwp * ky(i) * B(i) * tmp2;
  }
  
  
  template<class Dimension>
  void LayeredPlaneWaveIncidentField<Complex_wp, Dimension>    
  ::UpdateCoefAB(Complex_wp& a_, Complex_wp& b_) const
  { 
    int i = num_layer;
    
    a_ = a(i);
    b_ = b(i);
  }


  /***********************
   * HankelIncidentField *
   ***********************/
  
  
  //! default constructor
  template<class Dimension>
  HankelIncidentField<Complex_wp, Dimension>::
  HankelIncidentField(const R_N& pt, const R_N& u)
  {
    Init(pt, u);
  }
  
  
  //! kind of constructor
  template<class Dimension>
  void HankelIncidentField<Complex_wp, Dimension>::
  Init(const R_N& pt, const R_N& u)
  {
    origin = pt;
    k_wave = u;
    norm_kwave = Norm2(k_wave);
  }

  
  //! evaluation of Hankel function
  template<>
  void HankelIncidentField<Complex_wp, Dimension2>
  ::EvaluateFunction(const R_N& pt, Complex_wp& f) const
  {
    Real_wp ka = norm_kwave*pt.Distance(this->origin);
    f = jn(0, ka) + Iwp* yn(0, ka);
  }
  
  
  //! gradient of Hankel function
  template<>
  void HankelIncidentField<Complex_wp, Dimension2>
  ::EvaluateFunctionGradient(const R_N& pt, Complex_wp& f,
			     TinyVector<Complex_wp, 2>& df) const
  {
    Real_wp radius = pt.Distance(this->origin);
    Real_wp ka = norm_kwave*radius;
    Real_wp coef = norm_kwave/radius;
    f = jn(0, ka) + Iwp* yn(0, ka);
    Complex_wp dH0 = -(jn(1, ka) + Iwp* yn(1, ka)) * coef;
    df(0) =  dH0 * (pt(0) - this->origin(0));    
    df(1) =  dH0 * (pt(1) - this->origin(1));    
  }


  //! evaluation of spherical Hankel function
  template<>
  void HankelIncidentField<Complex_wp, Dimension3>
  ::EvaluateFunction(const R_N& pt, Complex_wp& f) const
  {
    Real_wp ka = norm_kwave*pt.Distance(this->origin);
    f = Complex_wp(sin(ka), -cos(ka)) / ka; 
  }
  
  
  //! gradient of spherical Hankel function
  template<>
  void HankelIncidentField<Complex_wp, Dimension3>
  ::EvaluateFunctionGradient(const R_N& pt, Complex_wp& f,
			     TinyVector<Complex_wp, 3>& df) const
  {
    Real_wp radius = pt.Distance(this->origin);
    Real_wp ka = norm_kwave*radius;
    Real_wp sin_ka = sin(ka), cos_ka = cos(ka);
    Real_wp coef = norm_kwave / radius;
    f = Complex_wp(sin_ka, -cos_ka) / ka;
    Complex_wp dH0 = (-f / ka + Complex_wp(cos_ka, sin_ka) / ka) * coef;
    df(0) = dH0 * (pt(0) - this->origin(0));
    df(1) = dH0 * (pt(1) - this->origin(1));
    df(2) = dH0 * (pt(2) - this->origin(2));
  }
  
  
  /*****************************
   * GaussianBeamIncidentField *
   *****************************/
  
  
#ifdef MONTJOIE_WITH_TWO_DIM  
  //! constructor with parameter
  GaussianBeamIncidentField<Complex_wp, Dimension2>::
  GaussianBeamIncidentField(const R2& pt, const R2& u, const Real_wp& w)
  {
    Init(pt, u, w);
  }


  //! kind of constructor
  void GaussianBeamIncidentField<Complex_wp, Dimension2>
  ::Init(const R2& pt, const R2& u, const Real_wp& w)
  {
    origin = pt;
    k_wave = u;
    omega = Norm2(k_wave);
    waist = w;
  }
  
  
  //! evaluation of gaussian beam
  void GaussianBeamIncidentField<Complex_wp, Dimension2>
  ::EvaluateFunction(const R2& point, Complex_wp& f) const
  {
    // on change de coordonnees
    // on s'aligne dans la direction du faisceau gaussien
    R2 pt = point - this->origin;
    Real_wp z = DotProd(k_wave, pt)/omega;
    Real_wp x = (k_wave(0)*pt(1)-k_wave(1)*pt(0))/omega;
    
    // z -> direction de propagation de l'onde gaussienne
    // x -> direction d'atténuation
    // u_incident = exp(-0.5*log(1+iz/z0) + ikz) * exp(-0.5 kx²/(z0+iz))
    // on en deduit le laplacien + k² U
    Real_wp z0 = 0.5*waist*waist*omega;
    Complex_wp facteur = exp(Complex_wp(-0.5)*log(Complex_wp(1)+Iwp*z/z0)+Iwp*omega*z);
    f = facteur*exp(Complex_wp(-0.5*omega*x*x)/(Complex_wp(z0)+Iwp*z));
  }
  
  
  //! gradient of gaussian beam
  void GaussianBeamIncidentField<Complex_wp, Dimension2>
  ::EvaluateFunctionGradient(const R2& point, Complex_wp& f, TinyVector<Complex_wp, 2>& df) const
  {
    R2 pt = point - this->origin;
    
    // on change de coordonnees
    // on s'aligne dans la direction du faisceau gaussien
    Real_wp z = DotProd(k_wave, pt)/omega;
    Real_wp x = (k_wave(0)*pt(1)-k_wave(1)*pt(0))/omega;
    // z -> direction de propagation de l'onde gaussienne
    // x -> direction d'attenuation

    // changement de variables  z = cos(\teta) \hat{x} + sin(teta) \hat{y}
    //                          x = -sin(\teta) \hat{x} + cos(\teta) \hat{y}
    // rotation de teta des coordonnees ou teta est l'angle d'incidence

    // u_incident = exp(-0.5*log(1+iz/z0) + ikz) * exp(-0.5 kx/(z0+iz))

    Real_wp z0 = 0.5*waist*waist*omega;
    Complex_wp facteur = exp(Complex_wp(-0.5)*log(Complex_wp(1.0)+Iwp*z/z0)
			     + Iwp*omega*z);

    f = facteur*exp(Complex_wp(-0.5*omega*x*x) / (Complex_wp(z0) + Iwp*z));

    // on calcule ensuite le gradient
    // \dvp{U}{x} = -\frac{k x}{z_0+\ii z} U
    // \dvp{U}{z} = ( -\frac{\ii }{2(z_0+\ii z)} + \ii k
    ///                   + \ii k \frac{x^2}{2(z_0+\ii z)^2} ) U
    Complex_wp z0Plusiz = Complex_wp(z0, z);
    Complex_wp du_dx = Complex_wp(-omega*x)/z0Plusiz * f;
    Complex_wp du_dz = Iwp * ( Complex_wp(-0.5)/z0Plusiz + omega 
			       + Complex_wp(0.5*omega*x*x)
			       /(z0Plusiz*z0Plusiz)) * f;

    // on fait ensuite la rotation de teta
    // \dvp{U}{\hat{x}} = cos(\teta) \dvp{U}{z} - sin(teta) \dvp{U}{x}
    // \dvp{U}{\hat{x}} = sin(\teta) \dvp{U}{z} + cos(teta) \dvp{U}{x}
    df(0) = (Complex_wp(-k_wave(1)) * du_dx + k_wave(0) * du_dz)/omega;
    df(1) = (k_wave(0) * du_dx + k_wave(1) * du_dz)/omega;
  }
#endif


#ifdef MONTJOIE_WITH_THREE_DIM
  //! constructor with parameter
  GaussianBeamIncidentField<Complex_wp, Dimension3>::
  GaussianBeamIncidentField(const R3& pt, const R3& u, const Real_wp& w)
  {
    Init(pt, u, w);
  }


  //! kind of constructor
  void GaussianBeamIncidentField<Complex_wp, Dimension3>
  ::Init(const R3& pt, const R3& u, const Real_wp& w)
  {
    origin = pt;
    k_wave = u;
    omega = Norm2(k_wave);
    waist = w;
    k_unit = k_wave; Mlt(1.0/omega, k_unit);
  }
  
  
  //! evaluation of gaussian beam
  void GaussianBeamIncidentField<Complex_wp, Dimension3>
  ::EvaluateFunction(const R3& point, Complex_wp& f) const
  {
    // on change de coordonnees
    // on s'aligne dans la direction du faisceau gaussien
    R3 pt, vec_u; pt = point - this->origin;
    Real_wp z = DotProd(k_unit, pt);
    TimesProd(k_unit, pt, vec_u);
    Real_wp r2 = vec_u(0)*vec_u(0) + vec_u(1)*vec_u(1) + vec_u(2)*vec_u(2);
    
    // z -> direction de propagation du faisceau gaussien
    // x,y -> directions d'atténuation, \f$ r^2 = x^2 + y^2 \f$
    // \f$ u = exp(-log(1+\frac{iz}{z0} ) + ikz) \cdot exp(-0.5 \frac{kr^2}{(z0+iz)})
    Real_wp z0 = 0.5*waist*waist*omega;
    Complex_wp facteur = exp(-log(Complex_wp(1)+Iwp*z/z0)+Iwp*omega*z);
    f = facteur*exp(Complex_wp(-0.5*omega*r2)/(z0+Iwp*z));
  }
  
  
  //! gradient of gaussian beam
  void GaussianBeamIncidentField<Complex_wp, Dimension3>
  ::EvaluateFunctionGradient(const R3& point, Complex_wp& u_inc, TinyVector<Complex_wp, 3>& f) const
  {
    // on change de coordonnees
    // on s'aligne dans la direction du faisceau gaussien
    R3 pt, vec_u, dr, dz; pt = point - this->origin;
    Real_wp z = DotProd(k_unit, pt);
    TimesProd(k_unit, pt, vec_u);
    Real_wp r2 = vec_u(0)*vec_u(0) + vec_u(1)*vec_u(1) + vec_u(2)*vec_u(2);
    Real_wp r = sqrt(r2);
    
    // z -> direction de propagation du faisceau gaussien
    // x,y -> directions d'atténuation, \f$ r^2 = x^2 + y^2 \f$
    // \f$ u = exp(-log(1+\frac{iz}{z0} + ikz)) \cdot exp(-0.5 \frac{kr^2}{(z0+iz)})
    Real_wp z0 = 0.5*waist*waist*omega;
    Complex_wp z0Plusiz = Complex_wp(z0, z);
    Complex_wp facteur = exp(-log(z0Plusiz/z0)+Iwp*omega*z);
    u_inc = facteur*exp(Complex_wp(-0.5*omega*r2)/z0Plusiz);
    Complex_wp du_dz =(Iwp/z0Plusiz*(Complex_wp(-1.0) + Complex_wp(omega*r2)
				       /(Complex_wp(2.0)*z0Plusiz)) + Iwp*omega)*u_inc;
    
    Complex_wp du_dr = Complex_wp(-omega*r)/z0Plusiz*u_inc;
    // on change de coordonnes 
    // du/dx' = du/dr dr/dx' + du/dz dz/dx'
    // avec r = || k \times x' || / omega, z = k \cdot x' / omega
    TimesProd(vec_u, k_unit, dr); 
    if (abs(r) <= R3::threshold)
      dr.Zero();
    else
      Mlt(1.0/r, dr);
    
    dz = k_unit;
    f(0) = du_dz*dz(0) + du_dr*dr(0);
    f(1) = du_dz*dz(1) + du_dr*dr(1);
    f(2) = du_dz*dz(2) + du_dr*dr(2);
  }
#endif
  

    //! Applies FFT only for a single mode
  template<>
  void VirtualProjectorFEM_Base<Real_wp>::ApplyFFT_Point(int nx, int ny, int nz,
							 Vector<Vector<Vector<VectReal_wp> > >& feval)
  {
  }

  
  //! Applies FFT only for a single mode
  template<>
  void VirtualProjectorFEM_Base<Complex_wp>::ApplyFFT_Point(int nx, int ny, int nz,
							    Vector<Vector<Vector<VectComplex_wp> > >& feval)
  {
    Vector<Complex_wp> feval_point(feval.GetM()); Complex_wp u_point;
    SetComplexZero(u_point);
    for (int k = 0; k < feval(0).GetM(); k++)
      for (int p = 0; p < feval(0)(k).GetM(); p++)
        for (int q = 0; q < feval(0)(k)(p).GetM(); q++)
          {
            for (int n = 0; n < feval.GetM(); n++)
              feval_point(n) = feval(n)(k)(p)(q);
            
            if (nz < 0)
              {
                if (ny < 0)
                  {
                    fft_interface.ApplyForwardPoint(nx, feval_point, u_point);
                    feval(0)(k)(p)(q) = u_point;
                  }
                else
                  {
                    fft_interface.ApplyForwardPoint(nx, ny, feval_point, u_point);
                    feval(0)(k)(p)(q) = u_point;
                  }
              }
            else
              {
                fft_interface.ApplyForwardPoint(nx, ny, nz, feval_point, u_point);
                feval(0)(k)(p)(q) = u_point;
              }
          }
  }


  //! Applies standard FFT
  template<>
  void VirtualProjectorFEM_Base<Real_wp>::ApplyFFT_Global(Vector<Vector<Vector<VectReal_wp> > >& feval)
  {
  }

  
  //! Applies standard FFT
  template<>
  void VirtualProjectorFEM_Base<Complex_wp>::ApplyFFT_Global(Vector<Vector<Vector<VectComplex_wp> > >& feval)
  {
    Vector<Complex_wp> feval_point(feval.GetM());
    for (int k = 0; k < feval(0).GetM(); k++)
      for (int p = 0; p < feval(0)(k).GetM(); p++)
        for (int q = 0; q < feval(0)(k)(p).GetM(); q++)
          {
            for (int n = 0; n < feval.GetM(); n++)
              feval_point(n) = feval(n)(k)(p)(q);
	  
            fft_interface.ApplyForward(feval_point);
            
            for (int n = 0; n < feval.GetM(); n++)
              feval(n)(k)(p)(q) = feval_point(n);
          }
  }
  

  /***********************
   * VirtualProjectorFEM *
   ***********************/

  
  //! returns the number of unknowns
  template<class T, class Dimension>
  int VirtualProjectorFEM<T, Dimension>::GetNbUnknowns() const
  {
    return nb_unknowns;
  }
  
  
  template<class T, class Dimension>
  size_t VirtualProjectorFEM<T, Dimension>::GetMemorySize() const
  {
    size_t taille = this->fft_interface.GetMemorySize();
    return taille;
  }

  
  //! Applies FFT to recover the components of function on the current mode
  template<class T, class Dimension>
  void VirtualProjectorFEM<T, Dimension>
  ::ApplyFFT(Vector<Vector<Vector<Vector<T> > > >& feval)
  {
    if (var_boundary.ModesNotStored())
      {
        int nx, ny, nz; bool teta_sym;
        var_boundary.GetPeriodicModes(var_boundary.GetCurrentModeNumber(), nx, ny, nz, teta_sym);
	this->ApplyFFT_Point(nx, ny, nz, feval);
      }
    else
      {
	this->ApplyFFT_Global(feval);
      }
  }
  
    
  //! modification of points in order to compute the source for all the modes
  template<class T, class Dimension>
  void VirtualProjectorFEM<T, Dimension>
  ::ModifyPoints(int nsrc, VectR_N& s, SetPoints<Dimension>& Pts,
		 SetMatrices<Dimension>& Mat)
  {
    if (var_boundary.GetSymmetryType() == var_boundary.NO_SYMMETRY)
      {
        // nothing to do 
        return;
      }
    
    int nx, ny, nz; bool teta_sym;
    var_boundary.GetPeriodicNumberModes(nx, ny, nz, teta_sym);
    
    int ix, iy, iz;
    var_boundary.GetPeriodicModes(nsrc, ix, iy, iz, teta_sym);
    
    switch (var_boundary.GetSymmetryType())
      {
      case VarBoundaryCondition_Base::PERIODIC_THETA :
        {
          if (nsrc > 0)
            {
              Real_wp teta = nsrc*var_problem.mesh.GetPeriodicAlpha();
              Pts.RotatePoints(teta);
              Mat.RotateNormale(teta);
            }
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_ZTHETA :
        {
          Real_wp dz = var_problem.GetZmax() - var_problem.GetZmin();
          int posz = (ny + 1)/2 - ny + iy;
          R_N vec_u(0, 0, posz*dz);
          
          Pts.TranslatePoints(vec_u);
          
          if (ix > 0)
            {
              Real_wp teta = ix*var_problem.mesh.GetPeriodicAlpha();
              Pts.RotatePoints(teta);
              Mat.RotateNormale(teta);
            }
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_X :
        {
          Real_wp dx = var_problem.GetXmax() - var_problem.GetXmin();
          int posx = (nx + 1)/2 - nx + ix;
          R_N vec_u(posx*dx, 0, 0);
          
          Pts.TranslatePoints(vec_u);          
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_Y :
        {
          Real_wp dy = var_problem.GetYmax() - var_problem.GetYmin();
          int posy = (nx + 1)/2 - nx + ix;
          R_N vec_u(0, posy*dy, 0);
          
          Pts.TranslatePoints(vec_u);
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_Z :
        {
          Real_wp dz = var_problem.GetZmax() - var_problem.GetZmin();
          int posz = (nx + 1)/2 - nx + ix;
          R_N vec_u(0, 0, posz*dz);
          
          Pts.TranslatePoints(vec_u);          
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_XY :
        {
          Real_wp dx = var_problem.GetXmax() - var_problem.GetXmin();
          int posx = (nx + 1)/2 - nx + ix;
          Real_wp dy = var_problem.GetYmax() - var_problem.GetYmin();
          int posy = (ny + 1)/2 - ny + iy;
          R_N vec_u(posx*dx, posy*dy, 0);
          
          Pts.TranslatePoints(vec_u);          
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_XZ :
        {
          Real_wp dx = var_problem.GetXmax() - var_problem.GetXmin();
          int posx = (nx + 1)/2 - nx + ix;
          Real_wp dz = var_problem.GetZmax() - var_problem.GetZmin();
          int posz = (ny + 1)/2 - ny + iy;
          R_N vec_u(posx*dx, 0, posz*dz);
          
          Pts.TranslatePoints(vec_u);
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_YZ :
        {
          Real_wp dy = var_problem.GetYmax() - var_problem.GetYmin();
          int posy = (nx + 1)/2 - nx + ix;
          Real_wp dz = var_problem.GetZmax() - var_problem.GetZmin();
          int posz = (ny + 1)/2 - ny + iy;
          R_N vec_u(0, posy*dy, posz*dz);
          
          Pts.TranslatePoints(vec_u);
        }
        break;
      case VarBoundaryCondition_Base::PERIODIC_XYZ :
        {
          Real_wp dx = var_problem.GetXmax() - var_problem.GetXmin();
          int posx = (nx + 1)/2 - nx + ix;
          Real_wp dy = var_problem.GetYmax() - var_problem.GetYmin();
          int posy = (ny + 1)/2 - ny + iy;
          Real_wp dz = var_problem.GetZmax() - var_problem.GetZmin();
          int posz = (nz + 1)/2 - nz + iz;
          R_N vec_u(posx*dx, posy*dy, posz*dz);
          
          Pts.TranslatePoints(vec_u);
        }
        break;
      }
  }

  
  //! computes components of the function on each mode
  template<class T, class Dimension>
  void VirtualProjectorFEM<T, Dimension>
  ::ModifyEvaluationProjection(Vector<Vector<Vector<Vector<T> > > >& feval,
                               Vector<bool>& is_f_vec)
  {
    if (var_boundary.GetSymmetryType() == var_boundary.NO_SYMMETRY)
      {
        // nothing to do 
        return;
      }

    int nx, ny, nz; bool teta_sym;
    var_boundary.GetPeriodicNumberModes(nx, ny, nz, teta_sym);

    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    if (teta_sym) 
      {
        for (int ix = 0; ix < nx; ix++)
          for (int iy = 0; iy < ny; iy++)
            {
              int n = ix*ny + iy;
              Real_wp cos_nalpha(1), sin_nalpha(0);
              this->fft_interface.GetCosSinAlpha(ix, cos_nalpha, sin_nalpha);
              
              for (int k = 0; k < feval(n).GetM(); k++)
                if (is_f_vec(k))
                  for (int m = 0; m < feval(n)(k).GetM(); m++)
                    {
                      int nb_vec = feval(n)(k)(m).GetM() / Dimension::dim_N;
                      for (int p = 0; p < nb_vec; p++)
                        {
                          CopyVector(feval(n)(k)(m), p, vec_u);
                          vec_v(0) = cos_nalpha * vec_u(0) + sin_nalpha * vec_u(1);
                          vec_v(1) = -sin_nalpha * vec_u(0) + cos_nalpha * vec_u(1);
                          if (Dimension::dim_N == 3)
                            vec_v(2) = vec_u(2);
                          
                          CopyVector(vec_v, p, feval(n)(k)(m));
                        }
                    }
            }
      }
    
    this->ApplyFFT(feval);    
  }


  /***************************
   * IncidentWaveSourceField *
   ***************************/

  
  template<class T, class Dimension>
  void IncidentWaveSourceField<T, Dimension>::EvaluateFunction(const R_N& x, Vector<T>& f) const
  {
    T u_inc;
    incident_wave->EvaluateFunction(x, u_inc);
    f = polar; Mlt(u_inc, f);
  }
  

  template<class T, class Dimension>
  void IncidentWaveSourceField<T, Dimension>::EvaluateGradient(const R_N& x, Vector<T>& fvec) const
  {
    T f;
    TinyVector<T, Dimension::dim_N> vec_u, df;
    incident_wave->EvaluateFunctionGradient(x, f, df);
    for (int i = 0; i < polar.GetM(); i++)
      {
	vec_u = df;
	vec_u *= polar(i);
	CopyVector(vec_u, i, fvec);
      }
  }
  
  
  /*************************
   * IncidentWaveProjector *
   *************************/
  

  //! initialization of the projector with the initial time
  template<class T, class Dimension>
  void IncidentWaveProjector<T, Dimension>::Init(const Real_wp& t0)
  {
    this->t = t0;
  }

  
  //! method called for each element i
  template<class T, class Dimension>
  void IncidentWaveProjector<T, Dimension>::InitElement(int i, const typename Dimension::VectR_N& s)
  {
    incident_wave.InitElement(i, s);
  }
  
  
  //! evaluates incident field
  template<class T, class Dimension>
  void IncidentWaveProjector<T, Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f)
  {
    T u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
    f = polarization; Mlt(u_inc, f);
  }
  
  
  /********************
   * VirtualSourceFEM *
   ********************/

  
  //! overloaded for time-dependent sources
  template<class T, class Dimension>
  Real_wp VirtualSourceFEM<T, Dimension>
  ::Init(const Real_wp& t, const Real_wp& dt, int print_level,
	 int nb_deriv, bool scalar_eq)
  {
    cout << "This method should be overloaded" << endl;
    abort();
  }


  //! initialization for surfacic integration
  /*!
    \param[in] i referenced edge number or face number
    \param[in] num_face global edge number or global face number
    \param[in] num_elem element number, element adjacent to the edge or face
    \param[in] num_loc local edge number in the element
  */
  template<class T, class Dimension>
  void VirtualSourceFEM<T, Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    num_bound_ref_ = i;
    num_boundary_ = num_face;
    num_elem_ = num_elem;
    ref_boundary_ = var_problem.mesh.BoundaryRef(i).GetReference();
    num_loc_ = num_loc;
    // DISP(i); DISP(this->vars.mesh.BoundaryRef(i).GetReference());
  }
  
  
  //! returns true is the boundary owns to an inhomogeneous dirichlet boundary
  /*!
    \param[in] cond_ref boundary condition
  */
  template<class T, class Dimension>
  bool VirtualSourceFEM<T, Dimension>
  ::IsNonNullDirichletSource(int ref)
  {
    int cond_ref = this->var_problem.mesh.GetBoundaryCondition(ref);
    if ((cond_ref == BoundaryConditionEnum::LINE_DIRICHLET)
        || (cond_ref == BoundaryConditionEnum::LINE_SUPPORTED))
      {
        if (this->var_problem.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          return true;
        
        if (!var_problem.FirstOrderFormulationDG())
          return true;
      }
    
    return false;
  }


  //! returns true if a Dirichlet condition may be present
  template<class T, class Dimension>
  bool VirtualSourceFEM<T, Dimension>::PresenceDirichlet() const
  {
    return true;
  }
  
  
  //! projects volumetric source on modes
  template<class T, class Dimension>
  void VirtualSourceFEM<T, Dimension>
  ::ModifyEvaluationVolume(Vector<bool>& int_phi, Vector<bool>& int_grad,
			   Vector<Vector<Vector<Vector<T> > > >& feval,
			   Vector<Vector<Vector<Vector<T> > > >& feval_grad,
                           Vector<bool>& is_f_vec, Vector<bool>& is_df_vec)
  {
    if (var_boundary.GetSymmetryType() == var_boundary.NO_SYMMETRY)
      {
        // nothing to do 
        return;
      }

    int nx, ny, nz; bool teta_sym;
    var_boundary.GetPeriodicNumberModes(nx, ny, nz, teta_sym);

    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    if (teta_sym) 
      {
        for (int ix = 0; ix < nx; ix++)
          for (int iy = 0; iy < ny; iy++)
            {
              int n = ix*ny + iy;
              Real_wp cos_nalpha(1), sin_nalpha(0);
              this->fft_interface.GetCosSinAlpha(ix, cos_nalpha, sin_nalpha);
              
              if (int_phi(n))
                for (int k = 0; k < feval(n).GetM(); k++)
                  if (is_f_vec(k))
                    for (int m = 0; m < feval(n)(k).GetM(); m++)
                      {
                        int nb_vec = feval(n)(k)(m).GetM() / Dimension::dim_N;
                        for (int p = 0; p < nb_vec; p++)
                          {
                            CopyVector(feval(n)(k)(m), p, vec_u);
                            vec_v(0) = cos_nalpha * vec_u(0) + sin_nalpha * vec_u(1);
                            vec_v(1) = -sin_nalpha * vec_u(0) + cos_nalpha * vec_u(1);
                            if (Dimension::dim_N == 3)
                              vec_v(2) = vec_u(2);
                            
                            CopyVector(vec_v, p, feval(n)(k)(m));
                          }
                      }
              
              if (int_grad(n))
                for (int k = 0; k < feval_grad(n).GetM(); k++)
                  if (is_df_vec(k))
                    for (int m = 0; m < feval_grad(n)(k).GetM(); m++)
                      {
                        int nb_vec = feval_grad(n)(k)(m).GetM() / Dimension::dim_N;
                        for (int p = 0; p < nb_vec; p++)
                          {
                            CopyVector(feval_grad(n)(k)(m), p, vec_u);
                            vec_v(0) = cos_nalpha * vec_u(0) + sin_nalpha * vec_u(1);
                            vec_v(1) = -sin_nalpha * vec_u(0) + cos_nalpha * vec_u(1);
                            if (Dimension::dim_N == 3)
                              vec_v(2) = vec_u(2);
                            
                            CopyVector(vec_v, p, feval_grad(n)(k)(m));
                          }
                      }
            }
      }
    
    if (int_phi(0))
      this->ApplyFFT(feval);

    if (int_grad(0))
      this->ApplyFFT(feval_grad);
  }
    
  
  //! projects surfacic source on modes  
  template<class T, class Dimension>
  void VirtualSourceFEM<T, Dimension>
  ::ModifyEvaluationSurface(bool& int0_phi, bool& int0_grad,
			    Vector<Vector<Vector<Vector<T> > > >& feval,
                            Vector<Vector<Vector<Vector<T> > > >& feval_diff,
                            Vector<bool>& is_f_vec, Vector<bool>& is_df_vec)
  {
    Vector<bool> int_phi(feval.GetM()), int_grad(feval_diff.GetM());
    int_phi.Fill(int0_phi);
    int_grad.Fill(int0_grad);
    ModifyEvaluationVolume(int_phi, int_grad, feval, feval_diff, is_f_vec, is_df_vec);
  }
  
  
  /********************
   * VolumetricSource *
   ********************/
  

  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::InitDefaultValues()
  {
    f_vol.Reallocate(var_problem.GetNbPhysicalIndices()+1);
    g_vol.Reallocate(var_problem.GetNbPhysicalIndices()+1);
    f_surf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    g_surf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    VirtualSourceField<T, Dimension>* null_ptr = NULL; 
    f_vol.Fill(null_ptr);
    g_vol.Fill(null_ptr);
    f_surf.Fill(null_ptr);
    g_surf.Fill(null_ptr);
    
    variable_fvol.Reallocate(f_vol.GetM());
    variable_gvol.Reallocate(g_vol.GetM());
    variable_fsurf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    variable_fsurf.Fill(false);
    variable_fvol.Fill(false);
    variable_gvol.Fill(false);
    
    evalSurf.Reallocate(var_problem.GetNbComponentsAll());
    for (int k = 0; k < evalSurf.GetM(); k++)
      evalSurf(k).Reallocate(var_problem.mesh.GetNbBoundaryRef());
    
  }
  
  
  //! destructor
  template<class T, class Dimension>
  VolumetricSource_Base<T, Dimension>::~VolumetricSource_Base()
  {
    RemoveDuplicate(f_vol);
    RemoveDuplicate(g_vol);
    RemoveDuplicate(f_surf);
    
    for (int i = 0; i < f_vol.GetM(); i++)
      if (f_vol(i) != NULL)
	delete f_vol(i);

    for (int i = 0; i < g_vol.GetM(); i++)
      if (g_vol(i) != NULL)
	delete g_vol(i);
    
    for (int i = 0; i < f_surf.GetM(); i++)
      if (f_surf(i) != NULL)
	delete f_surf(i);

    for (int i = 0; i < g_surf.GetM(); i++)
      if (g_surf(i) != NULL)
	delete g_surf(i);
  }
  
  
  //! returns the size used by the object in bytes
  template<class T, class Dimension>
  size_t VolumetricSource_Base<T, Dimension>::GetMemorySize() const
  {
    size_t taille = VirtualSourceFEM<T, Dimension>::GetMemorySize();
    taille += Seldon::GetMemorySize(evalS) + Seldon::GetMemorySize(evalG) + Seldon::GetMemorySize(evalSurf);
    taille += sizeof(*this);
    return taille;
  }
  

  //! sets the volume source for references ref
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::SetVolumeSource(const IVect& ref, VirtualSourceField<T, Dimension>* f)
  {
    for (int i = 0; i < ref.GetM(); i++)
      {
	variable_fvol(ref(i)) = false;
	f_vol(ref(i)) = f;
      }
  }


  //! sets the gradient volume source for references ref
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::SetVolumeSourceGrad(const IVect& ref, VirtualSourceField<T, Dimension>* f)
  {
    for (int i = 0; i < ref.GetM(); i++)
      {
	variable_gvol(ref(i)) = false;
	g_vol(ref(i)) = f;
      }
  }


  //! sets the volume source
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::SetVolumeSourceFunction(VirtualSourceField<T, Dimension>& f)
  {
    RemoveDuplicate(f_vol);
    
    for (int i = 0; i < f_vol.GetM(); i++)
      if (f_vol(i) != NULL)
	delete f_vol(i);

    f_vol.Reallocate(var_problem.GetNbPhysicalIndices()+1);
    f_vol.Fill(&f);
  }
  
  
  //! nullifies the volume source
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::NullifyVolumeSourceFunction()
  {
    VirtualSourceField<T, Dimension>* null_ptr = NULL; 
    f_vol.Fill(null_ptr);
    g_vol.Fill(null_ptr);
  }
  

  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::SetModalSource(const IVect& ref)
  {
    for (int i = 0; i < ref.GetM(); i++)
      variable_fsurf(ref(i)) = true;
  }

  
  //! sets the surface source for references ref
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::SetSurfaceSource(const IVect& ref, VirtualSourceField<T, Dimension>* f)
  {
    for (int i = 0; i < ref.GetM(); i++)
      {
	variable_fsurf(ref(i)) = false;
	f_surf(ref(i)) = f;
      }
  }


  //! sets the surface source for references ref
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::SetSurfaceSourceGrad(const IVect& ref, VirtualSourceField<T, Dimension>* f)
  {
    for (int i = 0; i < ref.GetM(); i++)
      g_surf(ref(i)) = f;
  }


  //! sets a variable volume source
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::SetVariableSource(const IVect& ref)
  {
    for (int i = 0; i < ref.GetM(); i++)
      variable_fvol(ref(i)) = true;
  }


  //! sets a variable volume gradient source
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::SetVariableGradientSource(const IVect& ref)
  {
    for (int i = 0; i < ref.GetM(); i++)
      variable_gvol(ref(i)) = true;
  }


  //! sets a variable surface source
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>::SetVariableSurfaceSource(const IVect& ref)
  {
    for (int i = 0; i < ref.GetM(); i++)
      variable_fsurf(ref(i)) = true;
  }

  
  //! volumetric source
  template<class T, class Dimension>
  bool VolumetricSource_Base<T, Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (f_vol(ref) != NULL)
      return true;
    
    if (variable_fvol(ref))
      return true;
    
    return false;
  }
  
  
  //! evaluation of f in \f$ \int f \varphi \f$
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f)
  {
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (variable_fvol(ref))
      {
	for (int k = 0; k < f.GetM(); k++)
	  f(k) = evalS(k)(i)(j);
      }
    else
      {
	f_vol(ref)->EvaluateFunction(x, f);
        var_source.ModifyVolumetricSource(i, j, x, *f_vol(ref), f);
	Mlt(coef_vol, f);
      }    
  }
  
  
  //! volumetric with gradient of basis functions ?
  template<class T, class Dimension>
  bool VolumetricSource_Base<T, Dimension>::IsNonNullGradientSource(const VectR_N& s)
  {
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (g_vol(ref) != NULL)
      return true;
    
    if (variable_gvol(ref))
      return true;
    
    return false;
  }
    
  
  //! evaluation of f in \f$ \int f \nabla \varphi \f$
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f)
  {
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (variable_gvol(ref))
      {
	for (int k = 0; k < f.GetM(); k++)
	  f(k) = evalG(k)(i)(j);
      }
    else
      {
	g_vol(ref)->EvaluateFunction(x, f);
      }
  }
  

  //! returns true for an hetereogeneous Dirichlet function
  template<class T, class Dimension>
  bool VolumetricSource_Base<T, Dimension>::IsNonNullDirichletSource(int ref)
  {
    if (VirtualSourceFEM<T, Dimension>::IsNonNullDirichletSource(ref))
      {
	if (ref <= 0)
	  return false;
	
	if (variable_fsurf(ref))
	  return true;
	
	if (f_surf(ref) != NULL)
	  return true;
      }
    
    return false;
  }
  
  
  //! evaluation of function f to project on dofs (for Dirichlet function)
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f)
  {
    if (variable_fsurf(this->ref_boundary_))
      {
	for (int k = 0; k < f.GetM(); k++)
	  f(k) = evalSurf(k)(this->num_bound_ref_)(j);
      }
    else
      f_surf(this->ref_boundary_)->EvaluateFunction(x, f);
  }
  

  //! initializes surface source for referenced edge/face i
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<T, Dimension>::InitSurface(i, num_face, num_elem, num_loc);

    if (variable_fsurf(this->ref_boundary_))
      return;
    
    //if (f_surf(this->ref_boundary_))
    //f_surf.Init(i);
  }
  

  //! returns true if there is a source term \int_\Gamma f \varphi 
  template<class T, class Dimension>
  bool VolumetricSource_Base<T, Dimension>
  ::IsNonNullSurfacicSource(int ref)
  {
    if (this->IsNonNullDirichletSource(ref))
      return false;
    
    if (ref <= 0)
      return false;
    
    if (variable_fsurf(ref))
      return true;
    
    if (f_surf(ref) != NULL)
      return true;
    
    return false;
  }
  
  
  //! evaluates f in source term \int_\Gamma f \varphi 
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
  {
    R_N x = PointsElem.GetPointQuadratureBoundary(k);
    if (variable_fsurf(this->ref_boundary_))
      {
	for (int n = 0; n < f.GetM(); n++)
	  f(n) = evalSurf(n)(this->num_bound_ref_)(k);
      }
    else
      f_surf(this->ref_boundary_)->EvaluateFunction(x, f);

  }


  //! returns true if there is a source term \int_\Gamma f \varphi 
  template<class T, class Dimension>
  bool VolumetricSource_Base<T, Dimension>
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    if (this->IsNonNullDirichletSource(ref))
      return false;
    
    if (ref <= 0)
      return false;
    
    if (g_surf(ref) != NULL)
      return true;
    
    return false;
  }
  
  
  //! evaluates f in source term \int_\Gamma f \nabla \varphi 
  template<class T, class Dimension>
  void VolumetricSource_Base<T, Dimension>
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
                                   const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
  {
    R_N x = PointsElem.GetPointQuadratureBoundary(k);
    g_surf(this->ref_boundary_)->EvaluateFunction(x, f);
  }
  
} // end namespace

#define MONTJOIE_FILE_SOURCE_SPATIALE_CXX
#endif
