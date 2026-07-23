#ifndef MONTJOIE_FILE_ADVECTION_PROBLEM_CXX

namespace Montjoie
{
  
  /******************************
   * Generic Advection Equation *
   ******************************/
  
  
  //! default constructor
  AdvectionEquation::AdvectionEquation() : GenericWaterWaveEquation()
  {
  }
    
  
  //! computation of c and c' for the collocation/quadrature points of the method
  template<class GenericPb>
  void AdvectionEquation::ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax)
  {      
    // we approximate bottom with higher degree in order to have
    // a correct approximation of Cx
    InitOrder(var.mesh.GetOrder() + 3, var);
    int Nlob = lob.GetNbPointsQuad();
    VectReal_wp Xhigh(Nlob*var.mesh.GetNbElt());
    int nb = 0;
    
    // Xhigh : position of nodal points used to interpolate c
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	Real_wp h = var.mesh.Vertex(i+1) - var.mesh.Vertex(i);
	for (int j = 0; j < Nlob; j++)
	  Xhigh(nb++) = var.mesh.Vertex(i) + h*lob.Points(j);
      }
      
    // computation of the bottom
    VectReal_wp Bhigh, Chigh;
    ComputeBottom(xmin, xmax, Xhigh, Bhigh, Chigh);
    
    // computation of Cx
    VectReal_wp CxDof;
    GetGradient(var, Chigh, CxDof);
    
    GetUquadrature(var, Chigh, Cquad);
    GetUquadrature(var, CxDof, CxQuad);
  }

  
  /**************************
   * LDG Advection Equation *
   **************************/

  
  //! default constructor
  AdvectionStiffnessOperatorDG::AdvectionStiffnessOperatorDG(AdvectionEquationDG& var_kdv)
    : var(var_kdv)
  {
  }
    
    
  //! numerical fluxes for advection equation
  template<class Vector1, class Vector2>
  void AdvectionStiffnessOperatorDG
  ::EvaluateNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up, Vector2& flux) const
  {
    // centered fluxes : 1/2 ( c^+ u^+  - c^- u^-)
    flux.Zero();
    Real_wp flux_m = var.Cquad(i)*Um(0);
    Real_wp flux_p = var.Cquad(i)*Up(0);
    flux(0) = 0.5*(flux_p - flux_m);    
    
    // penalisation term -1/2 alpha (u^-  -  u^+)
    if (k == 0)
      {
	flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
    else
      {
	flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
  }
  
  
  //! numerical fluxes and derivatives for advection equation
  template<class Vector1, class Vector2, class Matrix1>
  void AdvectionStiffnessOperatorDG
  ::EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um, const Vector1& Up,
                                    Vector2& flux,
				    Matrix1& df_dum, Matrix1& df_dup) const
  {
    flux.Zero(); df_dum.Zero(); df_dup.Zero();
    Real_wp flux_m = var.Cquad(i)*Um(0);
    Real_wp flux_p = var.Cquad(i)*Up(0);
    flux(0) = 0.5*(flux_p - flux_m);    
    if (k == 0)
      {
	flux(0) -= 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
    else
      {
	flux(0) += 0.5*var.alpha_lax*(Um(0) - Up(0));
      }
          
    df_dum(0,0) = -0.5*var.Cquad(i);
    df_dup(0,0) = 0.5*var.Cquad(i);
          
    if (k == 0)
      {
	df_dum(0,0) -= 0.5*var.alpha_lax;
	df_dup(0,0) += 0.5*var.alpha_lax;	
      }
    else
      {
	df_dum(0,0) += 0.5*var.alpha_lax;
	df_dup(0,0) -= 0.5*var.alpha_lax;
      }
  }
  
  
  //! operator involved in volume integrals for advection equation
  template<class Vector1, class Vector2>
  void AdvectionStiffnessOperatorDG
  ::EvaluateFunction(int i, const Vector1& u, const Vector1& du,
		     Vector2& g, Vector2& h) const
  {
    // the volume integral is equal to
    // \int (c du/dx + 1/2 c' u) \varphi dx 
    h(0) = du(0)*var.Cquad(i) + 0.5*var.CxQuad(i)*u(0);
  }
  
  
  //! operator involved in volume integrals for advection equation and its derivatives
  template<class Vector1, class Vector2, class Matrix1>
  void AdvectionStiffnessOperatorDG
  ::EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
		       Vector2& g, Vector2& h, Matrix1& dg_du,
		       Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    h(0) = 0.5*var.CxQuad(i)*u(0) + du(0)*var.Cquad(i);
              
    dh_du(0, 0) = 0.5*var.CxQuad(i);
    dh_ddu(0, 0) = var.Cquad(i);
  }

  
  //! default constructor
  AdvectionEquationDG::AdvectionEquationDG() : var_stiff(*this)
  {
    //alpha_lax = 1.0;
    alpha_lax = 0.0;
  }
  
  
  //! copy constructor
  AdvectionEquationDG::AdvectionEquationDG(const AdvectionEquationDG& var)
    : AdvectionEquation(var), var_stiff(*this)
  {
    alpha_lax = var.alpha_lax;
  }
  

  /*************************************
   * Finite Element Advection Equation *
   *************************************/
  
  
  //! default constructor
  AdvectionStiffnessOperatorFem::
  AdvectionStiffnessOperatorFem(AdvectionEquationFem& var_eq) : var(var_eq)
  {
  }
  
  
  //! operators involved in volume integrals for advection equation
  template<class Vector1, class Vector2>
  void AdvectionStiffnessOperatorFem::
  EvaluateFunction(int i, const Vector1& u, const Vector1& du,
                   Vector2& g, Vector2& h) const
  {
    h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0);
  }
   
  
  //! operators involved in volume integrals for advection equation
  template<class Vector1, class Vector2, class Matrix1>
  void AdvectionStiffnessOperatorFem::
  EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                     Vector2& g, Vector2& h, Matrix1& dg_du,
                     Matrix1& dh_du, Matrix1& dg_ddu, Matrix1& dh_ddu) const
  {
    h(0) = var.Cquad(i)*du(0) + 0.5*var.CxQuad(i)*u(0);
    
    dh_du(0, 0) = 0.5*var.CxQuad(i);
    dh_ddu(0, 0) = var.Cquad(i);
  }

  
  //! default constructor
  AdvectionEquationFem::AdvectionEquationFem() : var_stiff(*this)
  {
  }
  

  AdvectionEquationFem::AdvectionEquationFem(const AdvectionEquationFem& A)
    : AdvectionEquation(A), var_stiff(*this)
  {
  }

  
  /****************************************
   * Finite Difference Advection Equation *
   ****************************************/

  
  //! default constructor
  FiniteDifferenceScheme1D<AdvectionEquation>::FiniteDifferenceScheme1D()
    : FiniteDifferenceScheme1D_Base<AdvectionEquation>()
  {
  }
  
    
  /*******************************
   * Spectral Advection Equation *
   *******************************/

  
  //! default constructor
  AdvectionEquationSpectral::AdvectionEquationSpectral()
  {
    constant_advec = false;
  }
  
  
  //! computation of c and cx 
  void AdvectionEquationSpectral
  ::ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
			    double xmin, double xmax, int N)
  {      
    VectReal_wp Xdof(N);
    Real_wp dx = (xmax-xmin)/N;
    for (int i = 0; i < N; i++)
      Xdof(i) = xmin + i*dx;
    
    VectReal_wp Bdof;
    ComputeBottom(xmin, xmax, Xdof, Bdof, Cquad);
    
    // computation of Cx
    VectComplex_wp Chat(N/2+1);
    fft.ApplyForward(Cquad, Chat);
    
    for (int i = 0; i < nu.GetM(); i++)
      Chat(i) *= Iwp*nu(i);
    
    CxQuad.Reallocate(N);
    fft.ApplyInverse(Chat, CxQuad);
    
    ProdReal.Reallocate(N);
    Xreal.Reallocate(N);
    dXreal.Reallocate(N);
    dXn.Reallocate(N/2+1);
    XnHat.Reallocate(N/2+1);
  }
  
  
  //! evaluation of the operator that evolves in time
  void AdvectionEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectComplex_wp& Xn,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectComplex_wp& ProdXn)
  {
    if (constant_advec)
      for (int i = 0; i < Xn.GetM(); i++)
	ProdXn(i) = -Iwp*nu(i)*Xn(i);
    else
      {
	for (int i = 0; i < Xn.GetM(); i++)
	  dXn(i) = Iwp*nu(i)*Xn(i);
	
	fft.ApplyInverse(dXn, dXreal);
	fft.ApplyInverse(Xn, Xreal);
	
	for (int i = 0; i < Xreal.GetM(); i++)
	  ProdReal(i) = -(Cquad(i)*dXreal(i) + 0.5*CxQuad(i)*Xreal(i));
	
	fft.ApplyForward(ProdReal, ProdXn);
      }
  }


  //! evaluation of the operator that evolves in time
  void AdvectionEquationSpectral
  ::EvaluateFunction(const Real_wp& t, const VectReal_wp& Xn,
		     FftRealInterface& fft, VectReal_wp& nu,
		     VectReal_wp& ProdXn)
  {
    fft.ApplyForward(Xn, XnHat);
    
    for (int i = 0; i < XnHat.GetM(); i++)
      dXn(i) = Iwp*nu(i)*XnHat(i);
	
    fft.ApplyInverse(dXn, dXreal);
    
    for (int i = 0; i < Xreal.GetM(); i++)
      ProdXn(i) = -(Cquad(i)*dXreal(i) + 0.5*CxQuad(i)*Xn(i));
    
  }
  
}

#define MONTJOIE_FILE_ADVECTION_PROBLEM_CXX
#endif
