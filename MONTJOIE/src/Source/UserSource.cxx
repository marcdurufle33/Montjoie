#ifndef MONTJOIE_FILE_USER_SOURCE_CXX

namespace Montjoie
{

  ////////////////////
  // SPATIAL SOURCE //


  //! inhomogeneous dirichlet condition
  /*!
    \param[in] i number of the element
    \param[in] j number of the quadrature point
    \param[in] x point where the function needs to be evaluated
    \param[out] f evaluation of the function
    */
  template<class T, class Dimension>
    void UserDefinedSource<T, Dimension>
    ::EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f)
    {
      if (this->num_loc_ >= 0)
      {
        Real_wp amplitude = to_num<Real_wp>(source_space_param(1));
        f(2) = amplitude;
      }
    }


  //! true if volumetric source
  template<class T, class Dimension>
    bool UserDefinedSource<T, Dimension>::
    IsNonNullVolumetricSource(const VectR_N& s)
    {
      //int ref = to_num<int>(source_space_param(1));
      //int ref_elem = this->vars.mesh.Element(this->num_elem_).GetReference();
      //if (ref_elem == ref)

      return true;

      //return false;
    }


  //! used for volumic integrals \f$ \int f \varphi \f$
  /*!
    \param[in] i number of the element
    \param[in] j number of the quadrature point
    \param[in] x point where the function f needs to be evaluated
    \param[out] f evaluation of function f
    */
  template<class T, class Dimension>
    void UserDefinedSource<T, Dimension>
    ::EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f)
    {
      R_N center;
      center(0) = 1.5;
      center(1) = 1.05;
      Real_wp r = 0.3;
      Real_wp gauss = exp(-(x(0)-center(0))*(x(0)-center(0))/(r*r)-(x(1)-center(1))*(x(1)-center(1))/(r*r));
      f(0) = -200.0*((x(0)-center(0))/(r*r))*gauss;
      f(1) = -200.0*((x(1)-center(1))/(r*r))*gauss;

      //T miomega,omega;
      //var_problem.GetMiomega(miomega);
      //var.GetOmega(omega);

      //f(0) = miomega*cos(x(0))-x(1)*sin(x(0))+1j*sin(x(1));
      //f(1) = omega*sin(x(1))+1j*x(0)*cos(x(1))+cos(x(0));
    }


  //! true if volumetric source with \f$ \int f \nabla \varphi \f$
  template<class T, class Dimension>
    bool UserDefinedSource<T, Dimension>::IsNonNullGradientSource(const VectR_N& s)
    {
      return false;
    }


  //! used for volumic integrals \f$ \int f \nabla \varphi \f$
  /*!
    \param[in] i number of the element
    \param[in] j number of the quadrature point
    \param[in] x point where the function f needs to be evaluated
    \param[out] f evaluation of function f
    For H(curl) elements \nabla \varphi is \nabla \times \varphi
    For H(div) elements \nabla \varphi is \nabla \cdot \varphi
    */
  template<class T, class Dimension>
    void UserDefinedSource<T, Dimension>
    ::EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f)
    {
    }


  //! true if surfacic source
  template<class T, class Dimension>
    bool UserDefinedSource<T, Dimension>::IsNonNullSurfacicSource(int ref)
    {
      // DISP(this->ref_boundary_);
      //if (ref == 3)
      //return true;

      return false;
    }


  //! surfacic integral \f$ \int_\Gamma f \varphi \f$
  /*!
    \param[in] k number of the quadrature point on the edge/face
    \param[in] PointsElem transformation of references points
    \param[in] MatricesElem jacobian matrices
    \param[out] f evaluation of f on quadrature points
    */
  template<class T, class Dimension>
    void UserDefinedSource<T, Dimension>
    ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
        const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
    {
    }


  //! returns true if a term \int_\Gamma f \cdot \nabla \varphi ds is present
  /*!
    \param[int] ref reference number of the face of Gamma
    */
  template<class T, class Dimension>
    bool UserDefinedSource<T, Dimension>::IsNonNullSurfacicSourceGradient(int ref)
    {
      return false;
    }


  //! surface term \int_\Gamma f \cdot \nabla \varphi ds is present
  /*!
    \param[in] k number of the quadrature point on the edge/face
    \param[in] PointsElem transformation of references points
    \param[in] MatricesElem jacobian matrices
    \param[out] f evaluation of f on quadrature points
    For H(curl) elements \nabla \varphi is \nabla \times \varphi
    For H(div) elements \nabla \varphi is \nabla \cdot \varphi
    */
  template<class T, class Dimension>
    void UserDefinedSource<T, Dimension>
    ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
        const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
    {
    }


  // SPATIAL SOURCE //
  ////////////////////


  ////////////////////////
  // INITIAL CONDITIONS //


  template<class T, class Dimension>
    void InitialUserFunction<T, Dimension>::InitElement(int, const typename Dimension::VectR_N& s)
    {
    }


  //! evaluation of initial vector u0, involved in a condition like u(t = 0) = u0
  /*!
    \param[in] i element number
    \param[in] j dof point number
    \param[in] x dof point
    \param[out] f u0(x)
    */
  template<class T, class Dimension>
    void InitialUserFunction<T, Dimension>
    ::EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f)
    {
      R_N center(to_num<Real_wp>(param_initial_condition(0)),
          to_num<Real_wp>(param_initial_condition(1)));  

      Real_wp r0 = to_num<Real_wp>(param_initial_condition(2));    
      //Real_wp ampl = param_initial_condition(3);

      Real_wp r = center.Distance(x);
      Real_wp alpha = square(r/r0);
      f(0) = exp(-7.0*alpha)/(r0*r0);

      // expressions of Kronbichler
      /*
         Real_wp alpha = to_num<Real_wp>(param_initial_condition(0));
         Real_wp l = to_num<Real_wp>(param_initial_condition(1));
         Real_wp xsi = to_num<Real_wp>(param_initial_condition(2));
         Real_wp rho = to_num<Real_wp>(param_initial_condition(3));

         Real_wp gaussian = exp(-alpha*(x(0)*x(0) + square(x(1)-l)));
         f(0) = xsi*alpha*l*x(0)*gaussian; //DISP(alpha); DISP(l); DISP(xsi); DISP(rho); DISP(gaussian);
         f(1) = -alpha/rho*gaussian*x(0);
         f(2) = (l-x(1))*alpha/rho*gaussian; */

      // random
      //f(0) = Real_wp(rand())/RAND_MAX;
      //f(1) = Real_wp(rand())/RAND_MAX;
      //f(2) = Real_wp(rand())/RAND_MAX;

      // oscillatory
      /*Real_wp L = 0.0012;
        Real_wp kx = 2.0*pi_wp/L;
        f(0) = sin(kx*x(0));*/
    }


  // INITIAL CONDITIONS //
  ////////////////////////


  /////////////////////
  // TEMPORAL SOURCE //


  //! virtual destructor
  template<class T>
    VirtualTimeSource<T>::~VirtualTimeSource()
    {
    }


  //! Returns the derivative of f
  /*!
    The user should not modify this functions
    */
  template<class T>
    T VirtualTimeSource<T>::EvaluateDerivative(const Real_wp& f)
    {
      cout << "EvaluateDerivative not implemented for this source" << endl;
      abort();
      T zero; SetComplexZero(zero);
      return zero;
    }


  //! initializes the time source
  /*!
    freq is the frequency as given in the data file
    param is the list of parameters given in the data file (keyword TemporalSource)
    t_final can be modified if you know when the source stops
    */
  TimeUserSource
    ::TimeUserSource(const Real_wp& freq, Real_wp& t_final, const VectString& param)
    {

    }


  //! evaluation of the source at time t
  /*!
    \param[in] t fixed time where the source has to be computed
    \returns value of the source
    */  
  Real_wp TimeUserSource::Evaluate(const Real_wp&)
  {
    // here you can place your expression of the source
    return Real_wp(0);
  }


  // TEMPORAL SOURCE //
  /////////////////////


  ////////////////////
  // VARIABLE MEDIA //


  void ComputeIndexGalbrun(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 3>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 3, 3>& hess_coef)
  {
    abort();
  }

  void ComputeIndexGalbrun(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, TinyVector<Real_wp, 3>& grad_coef,
      TinyMatrix<Real_wp, Symmetric, 3, 3>& hess_coef)
  {
    abort();
  }

  void ComputeIndexGalbrun(const Real_wp& x, const Real_wp& z, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 2>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 2, 2>& hess_coef)
  {
    Real_wp r = sqrt(x*x + z*z);
    Real_wp rho = exp(-r*r/16);
    Real_wp c = exp(-0.6*r)*(0.8*r+1.0+0.1*r*r);

    switch (num_index)
    {
      case 0:
        {
          // coefficient 1/(rho c^2)
          coef = Real_wp(1) / (rho*c*c);  
        }
        break;
      case 1:
        {
          // coefficient 1/rho
          coef = Real_wp(1) / rho;
        }
        break;
      case 2:
        {
          // damping
          coef = Real_wp(1) / (100*rho*c*c);
        }
        break;
      case 3:
        {
          Real_wp Rt = 14.475500686554541;
          Real_wp f = jn(10, r*Rt/4.0);

          if (r == Real_wp(0))
            SetComplexZero(coef);
          else
          {	 
            Real_wp sin_teta = z/r, cos_teta = x/r;
            if (num_component == 0)
              coef = -f*sin_teta;
            else
              coef = f*cos_teta;
          }

          //Real_wp df = 0.5*(jn(9, r*Rt/4) - jn(11, r*Rt/4))*Rt/4.0;
          //Real_wp u_s = f/10, v_s = (f + r/2*df)/10;

          // flow M
          /*if (num_component == 0)
            {
          // coef Mx / (rho c^2)
          Real_wp mr = u_s * x*z/(r*r) - v_s*x*z/(r*r);	      
          coef = mr / (rho*rho*c*c);
          }
          else if (num_component == 1)
          {
          // coef Mteta / (rho c^2)
          Real_wp mteta = r*r/20;
          coef = mteta / (rho*rho*c*c);
          }
          else
          {
          // coef Mz / (rho c^2)
          Real_wp mz = u_s*z*z/(r*r) + v_s*x*x/(r*r);
          coef = mz / (rho*rho*c*c);
          }*/

          if (r == Real_wp(0))
            SetComplexZero(coef);
        }
        break;
      case 4:
        {
          // coefficient beta = rho c^2
          coef = rho*c*c;
        }
    }    
  }

  void ComputeIndexGoldstein(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 3>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 3, 3>& hess_coef)
  {
    abort();
  }

  void ComputeIndexGoldstein(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, TinyVector<Real_wp, 3>& grad_coef,
      TinyMatrix<Real_wp, Symmetric, 3, 3>& hess_coef)
  {
    abort();
  }
  void ComputeIndexGoldstein(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, R2& grad_coef, Matrix2_2sym& hess_coef)
  {
    switch (num_index)
    {
      case 0:
        {
          // flow M
          if (num_component ==0)
          {
            coef = cosh(x);
            grad_coef(0) = sinh(x);
            grad_coef(1) = 0;
            hess_coef.Zero();
          }
          else
          {
            coef = -y*sinh(x);
            grad_coef(0) = -y*cosh(x);
            grad_coef(1) = -sinh(x);
            hess_coef.Zero();
          } 
        }
      break;
      case 1:
        {
          //damping
          coef = 0.1;
          grad_coef.Zero();
          hess_coef.Zero();
        }
      break;
      case 2:
      {
        //density
        coef = 1;
        grad_coef.Zero();
        hess_coef.Zero();
      }
      break;
      case 3:
      {
        //sound speed
        coef = sqrt(0.25*cosh(2*x)-0.5*y*y+10);
        grad_coef.Zero();
        hess_coef.Zero();
      }
      break;
      case 4:
      {
        //pressure
        coef = 0.25*cosh(2*x)-0.5*y*y+10;
        grad_coef(0) = 0.5*sinh(2*x);
        grad_coef(1) = -y;
        hess_coef.Zero();
      }
    }
  }


  void ComputeIndexGoldstein(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 2>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 2, 2>& hess_coef)
  {
    switch (num_index)
    {
      case 0:
        {
          // flow M
          if (num_component ==0)
          {
            coef = cosh(x);
            grad_coef(0) = sinh(x);
            grad_coef(1) = 0;
            hess_coef.Zero();
          }
          else
          {
            coef = -y*sinh(x);
            grad_coef(0) = -y*cosh(x);
            grad_coef(1) = -sinh(x);
            hess_coef.Zero();
          } 
        }
      break;
      case 1:
        {
          //damping
          coef = 0.1;
          grad_coef.Zero();
          hess_coef.Zero();
        }
      break;
      case 2:
      {
        //density
        coef = 1;
        grad_coef.Zero();
        hess_coef.Zero();
      }
      break;
      case 3:
      {
        //sound speed
        coef = sqrt(0.25*cosh(2*x)-0.5*y*y+10);
        grad_coef.Zero();
        hess_coef.Zero();
      }
      break;
      case 4:
      {
        //pressure
        coef = 0.25*cosh(2*x)-0.5*y*y+10;
        grad_coef(0) = 0.5*sinh(2*x);
        grad_coef(1) = -y;
        hess_coef.Zero();
      }
    }
  }

  void ComputeIndexGalbrun_old(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, R2& grad_coef, Matrix2_2sym& hess_coef)
  {
    switch (num_index)
    {
      case 0:
        {
          // flow M
          Real_wp a,b;

          a = 50.0;
          b = 10.0;
          if (num_component == 0)
          {
            // Mx

            coef = 1.0+(a*a)/(x*x+(y+b)*(y+b))-(2*a*a*x*x)/((x*x+(y+b)*(y+b))*(x*x+(y+b)*(y+b)));
            grad_coef(0) = 0.0;
            grad_coef(1) = 0.0;
            hess_coef(0, 0) = 0.0;
            hess_coef(0,1) = 0.0;
            hess_coef(1,0) = 0.0;
          }
          else
          {
            // My
            coef = -2*a*a*x*(y+b)/((x*x+(y+b)*(y+b))*(x*x+(y+b)*(y+b)));
            grad_coef(0) = 0.0;
            grad_coef(1) = 0.0;
            hess_coef(0, 0) = 0.0;
            hess_coef(0,1) = 0.0;
            hess_coef(1,0) = 0.0;
          }
        }
        break;
      case 1:
        {
          // damping sigma
          coef = 0.1;
          grad_coef.Zero();
          hess_coef.Zero();
        }
        break;
      case 2:
        {
          // density rho0
          coef = 1.5; 
          grad_coef.Zero();
          hess_coef.Zero();
        }
        break;
      case 3:
        {
          // sound speed c0
          coef = 12.3;
          grad_coef.Zero();
          hess_coef.Zero();
        }
        break;
      case 4:
        {
          // pressure p0
          coef = 7.6;
          grad_coef.Zero();
          hess_coef.Zero();
        }
        break;
    }
  }

  void ComputeIndexGalbrunAxi(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 2>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 2, 2>& hess_coef)
  {
    abort();
  }


  void ComputeIndexGalbrunAxi(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, R2& grad_coef, Matrix2_2sym& hess_coef)
  {
    switch (num_index)
    {
      case 3:
        {
          // flow M
          Real_wp rho = 1.6 + 0.3*cos(pi_wp/4*x)*sin(pi_wp/4*y);
          Real_wp drho_dr = -0.3*pi_wp/4*sin(pi_wp/4*x)*sin(pi_wp/4*y);
          Real_wp drho_dz = 0.3*pi_wp/4*cos(pi_wp/4*x)*cos(pi_wp/4*y);

          if (num_component == 0)
          {
            // Mr
            coef = 0.0;
            grad_coef(0) = 0.0;
            grad_coef(1) = 0.0;
          }
          else if (num_component == 1)
          {
            // Mphi
            coef = 0.0;
            grad_coef(0) = 0.0;
            grad_coef(1) = 0.0;
          }
          else
          {
            // Mz
            coef = cos(pi_wp/4*x)/rho;
            grad_coef(0) = (-pi_wp/4*sin(pi_wp/4*x)/rho - cos(pi_wp/4*x)*drho_dr/(rho*rho));
            grad_coef(1) = (-cos(pi_wp/4*x)*drho_dz/(rho*rho));
          }

          hess_coef.Fill(0);
        }
        break;
      case 4:
        {
          // damping sigma
          coef = 0.1;
          grad_coef.Fill(0);
          hess_coef.Fill(0);
        }
        break;
      case 0:
        {
          // density rho0
          coef = 1.6 + 0.3*cos(pi_wp/4*x)*sin(pi_wp/4*y);
          grad_coef(0) = -0.3*pi_wp/4*sin(pi_wp/4*x)*sin(pi_wp/4*y);
          grad_coef(1) = 0.3*pi_wp/4*cos(pi_wp/4*x)*cos(pi_wp/4*y);
          hess_coef.Fill(0);
        }
        break;
      case 1:
        {
          // sound speed c0
          coef = 1.2 + 0.2*(x*x/10 + sin(pi_wp*x/4)*cos(pi_wp*y/4));
          grad_coef(0) = 0.4*x/10 + 0.2*pi_wp/4*cos(pi_wp*x/4)*cos(pi_wp*y/4);
          grad_coef(1) = -0.2*pi_wp/4*sin(pi_wp*x/4)*sin(pi_wp*y/4);
          hess_coef.Fill(0);
        }
        break;
      case 2:
        {
          // pressure p0
          coef = 2.2 + 0.4*sin(pi_wp/4*x)*cos(pi_wp/4*y);
          grad_coef(0) = 0.4*pi_wp/4*cos(pi_wp/4*x)*cos(pi_wp/4*y);
          grad_coef(1) = -0.4*pi_wp/4*sin(pi_wp/4*x)*sin(pi_wp/4*y);
          hess_coef.Fill(0);
        }
        break;
    }
  }


  void ComputeIndexGalbrunAxi(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, R3& grad_coef, Matrix3_3sym& hess_coef)
  {
    abort();
  }


  void ComputeIndexGalbrunAxi(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 3>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 3, 3>& hess_coef)
  {
    abort();
  }


  void ComputeIndexAcous(const Real_wp& x, const Real_wp& y, const Real_wp& z, int num_index,
      int num_component, Complex_wp& coef, TinyVector<Complex_wp, 2>& grad_coef,
      TinyMatrix<Complex_wp, Symmetric, 2, 2>& hess_coef)
  {
    Real_wp r = sqrt(x*x+y*y);
    if (r < 0.7)
      coef = 1.0-0.5*r/0.7;
    else
      coef = sqrt(Real_wp(1.0005-r))/sqrt(Real_wp(1.0005-0.7))*0.5;

    coef = Real_wp(1)/square(coef);
  }

  void ComputeIndexAcous(const Real_wp& x, const Real_wp& y, const Real_wp& z, int num_index,
      int num_component, Real_wp& coef, R2& grad_coef, Matrix2_2sym& hess_coef)
  {
    Real_wp r = sqrt(x*x+y*y);
    if (r < 0.7)
      coef = 1.0-0.5*r/0.7;
    else
      coef = sqrt(Real_wp(1.0005-r))/sqrt(Real_wp(1.0005-0.7))*0.5;

    coef = Real_wp(1)/square(coef);
  }


  //Computing bensalah's indexes
  //void ComputeIndexGalbrunRing(const Real_wp& x, const Real_wp& y, int num_index, int num_component, Complex_wp& coef)
  void ComputeIndexGalbrun(const Real_wp& x, const Real_wp& y, int num_index,
      int num_component, Real_wp& coef, R2& grad_coef, Matrix2_2sym& hess_coef)
  {

    Real_wp r = sqrt(x*x+y*y);
    Real_wp gamma = 1.4;
    Real_wp mu = 0.3;
    Real_wp ri = 0.5;

    Real_wp v0 = 1.75;
    Real_wp rho = 10.2+v0*v0*(r-ri)*(gamma-1.0)/(mu*gamma);
    //Real_wp ep = 0.001;
    //Real_wp int_vit = v0*v0*(1.0/(2*ri*ri)-1.0/(2*r*r))+2*ep*v0*(1.0/ri-1.0/r)+ep*ep*log(r/ri);
    //Real_wp rho = 10.2+int_vit*(gamma-1.0)/(mu*gamma);
    //


    //rotating cylinder
    Real_wp theta = atan2(y,x);
    //Real_wp Gs2pir = 2.0/(2.0*M_PI*r);
    Real_wp R = 0.5;
    //Real_wp Uinf = 0.2;
    Real_wp U = 1.0;
    Real_wp Gamma = 4.0*M_PI*R*U+1.5;
    //Real_wp Gamma = 4.0*M_PI*R*U+0.00001;
    //Real_wp Gamma = 0.0;

    switch(num_index)
    {

      case 0: // flow M
        {
          Real_wp Vr = (U/r)*(r-(R*R)/r)*cos(theta);
          Real_wp Vt = -U*(1+(R*R)/(r*r))*sin(theta)-Gamma/(2*M_PI*r);
          if(num_component==0)
          {
            //orthoradial
            //coef = -v0*r*y;
            //grad_coef(0) = v0*y*x/r;
            //grad_coef(1) = -v0*r-v0*r*y*y/r;

            //radial
            //coef = v0*r*x;
            //grad_coef(0) = v0*r+v0*x*x/r;
            //grad_coef(1) = v0*x*y/r;


            //rotating cylinder
            coef = (x*Vr-y*Vt)/r;
            // d(atan2)/dx = -y/(r*r)
            Real_wp r3 = (x*x+y*y)*(x*x+y*y)*(x*x+y*y);
            grad_coef(0) = (x*(2.0*M_PI*R*R*U*(x*x-3.0*y*y)-Gamma*y*(x*x+y*y)))/(M_PI*r3);
            grad_coef(1) = (Gamma*(x*x*x*x-y*y*y*y)-4.0*M_PI*R*R*U*y*(y*y-3.0*x*x))/(2.0*M_PI*r3);
          }
          else
          {
            //orthoradial
            //coef = v0*r*x;
            //grad_coef(0) = v0*r+v0*x*x/r;
            //grad_coef(1) = v0*x*y/r;

            //radial
            //coef = v0*r*y;
            //grad_coef(0) = v0*y*x/r;
            //grad_coef(1) = +v0*r+v0*r*y*y/r;

            //test for transport
            coef = (y*Vr+x*Vt)/r;
            Real_wp r3 = (x*x+y*y)*(x*x+y*y)*(x*x+y*y);
            grad_coef(0) = (x*(2.0*M_PI*R*R*U*(x*x-3.0*y*y)-Gamma*y*(x*x+y*y)))/(M_PI*r3);
            grad_coef(1) = (Gamma*(x*x*x*x-y*y*y*y)-4.0*M_PI*R*R*U*y*(y*y-3.0*x*x))/(2.0*M_PI*r3);
          }
        }
        break;
      case 1: //damping sigma
        {
          coef = 0.1;
        }
        break;
      case 2: //density rho
        {
          coef = pow(rho,1.0/(gamma-1.0));
        }
        break;

      case 3: //sound spped c
        {
          coef = gamma*mu*rho;
        }
        break;
      case 4: //pressure p
        {
          coef = mu*pow(pow(rho,1.0/(gamma-1.0)),gamma);
        }
        break;
    }

  }

  //! in this function, you can specify your own variable index (rho, mu, epsilon, etc)
  /*!
    \param[in] var considered problem
    \param[in] ElementRho intermediary array used to access to rho, grad_rho
    \param[in] ref reference of the domain where the variable media is defined
    \param[in] num_index number of the index (depends on the equation)
    for example 0 -> rho, 1 -> mu for Helmholtz equation
    \param[in] num_component number of the component of the index
    (for a vectorial index, tensorial, etc)
    \param[in] PointsQuadrature quadrature points where you have to compute the value of index
    \param[in] compute_grad true if the gradient needs to be provided
    (e.g. for aeroacoustic equation)
    \param[inout] rho evaluation of the index on all quadrature points
    \param[inout] grad_rho evaluation of the gradient of index on all quadrature points
    (if compute_grad = true)
    */
  template<class Dimension, class T>
  void ComputeVariableUserIndex(const VarGeometryProblem<Dimension>& var,
                                const IVect& ElementRho, int ref, int num_index, int num_component,
                                const Vector<typename Dimension::VectR_N>& PointsQuadrature,
                                bool compute_grad, bool compute_hess,
                                Vector<Vector<T> >& rho,
                                Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
                                Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N,
                                Dimension::dim_N> > >& hess_rho,
                                const T& offset, const T& amplitude)
  {
    typedef typename Dimension::R_N R_N;
    
    // rho and grad_rho are already allocated
    // loop over elements
    T coef(0); TinyVector<T, Dimension::dim_N> grad_coef;
    TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> hess_coef;
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      if (var.mesh.Element(i).GetReference() == ref)
        {
          // element number used in arrays rho, rho_grad
          int i1 = ElementRho(i);
          
          // you can check which component, index you are modifying
          //if ((num_index == 0) && (num_component == 0))
          {
            // then you loop over quadrature points of the element
            for (int j = 0; j < PointsQuadrature(i).GetM(); j++)
              {
                // retrieving the point
                R_N point = PointsQuadrature(i)(j);
                Real_wp x = point(0), y = point(1);
                
                if (num_component == 0)
                  {
                    coef = y;
                    grad_coef.Init(0.0, 1.0);
                  }
                else
                  {
                    coef = -x;
                    grad_coef.Init(-1.0, 0.0);
                  }
                
                //DISP(offset); DISP(amplitude); DISP(x); DISP(y); DISP(coef);
                // computing value of index
                //ComputeIndexGoldstein(x, y, num_index, num_component,
                //                      coef, grad_coef, hess_coef);
                //ComputeIndexGalbrun(x, y, num_index, num_component,
                //    coef, grad_coef, hess_coef);
                //ComputeIndexGalbrunRing(x,y,num_index,num_component,coef);
                //              DISP(typeid(coef).name());
                //              double aa;
                //              DISP(typeid(&aa).name());
                //              std::cout << typeid(coef).name() <<std::endl;
                //              abort();
                
                /* rho(i1)(j) = 1.0 + exp(-(x*x + y*y));
                // and gradient if necessary
                if (compute_grad)
                {
                grad_rho(i1)(j).Init(-2.0*x, -2.0*y);
                grad_rho(i1)(j) *= exp(-(x*x + y*y));
                }*/
                
                // final modification of the index
                coef = offset + amplitude*coef;
                grad_coef *= amplitude;
                hess_coef *= amplitude;
                
                rho(i1)(j) = coef;
                if (compute_grad)
                  grad_rho(i1)(j) = grad_coef;
                
                if (compute_hess)
                  hess_rho(i1)(j) = hess_coef;
                
              }
          }
        }
  }
  
  
  // VARIABLE MEDIA //
  ////////////////////

}

#define MONTJOIE_FILE_USER_SOURCE_CXX
#endif
