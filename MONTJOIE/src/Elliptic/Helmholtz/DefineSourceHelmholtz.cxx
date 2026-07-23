#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_CXX

#include "DefineSourceHelmholtzInline.cxx"

namespace Montjoie
{
  
  /*****************************
   * DiffractedWaveSource_Helm *
   *****************************/
  
  
  //! Sets attributes to default values
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>::InitDefaultValues()
  {
    coef_vol = 1.0; rho0 = var_helm.rho0; rho = 1.0;
    
    // variables for gradient source
    coef_grad.SetIdentity(); mu0 = var_helm.mu0;
    mu.SetIdentity();
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    // variables for high conductivity condition
    order_condition = var_boundary.GetHighConductivityOrder();
    epsilon = 1e-5;
    alpha = ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >::alpha;
#endif
    omega2 = var_problem.GetSquareOmega();
    coef_mu0 = 0.0; coef_mu = 1.0;
    
    varying_rho = false; varying_mu = false;
  }
  
  
  //! initialization for each element
  /*!
    \param[in] num_elem element number
    if physical properties are constant inside the element
    coefficients \f$ -\omega^2(\rho - \rho_0) \f$ and (mu-mu0) are computed
   */
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
    incident_wave.InitElement(num_elem, s);
    
    // if the incident wave is defined for a non constant medium we update the value
    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);

    int ref = var_problem.mesh.Element(num_elem).GetReference();
    rho = Complex_wp(var_helm.ref_rho(ref))
      + Iwp*Complex_wp(var_helm.ref_sigma(ref))/var_problem.GetOmega();
    
    mu = var_helm.ref_mu(ref).GetConstant();

    Complex_wp sigma_p;
    var_helm.ref_drude(ref)
      .ModifyCoefficientHelmholtz(var_problem.GetOmega(), rho, sigma_p, mu);
    
    // coef_vol = omega2 (rho-rho_0) with rho complex !
    if (rho != rho0)
      {
        if (var_problem.FirstOrderFormulation())
          coef_vol = Iwp*var_problem.GetOmega()*(rho-rho0);
        else
          coef_vol = var_problem.GetSquareOmega()*(rho-rho0);
      }
    else
      coef_vol = Complex_wp(0);
    
    // coef_grad = mu_0 - mu with mu complex
    Complex_wp zero; SetComplexZero(zero);
    coef_grad.SetDiagonal(mu0);
    coef_grad -= mu;
    
    varying_rho = (var_helm.ref_rho(ref).IsVarying() || var_helm.ref_sigma(ref).IsVarying());
    varying_mu = var_helm.ref_mu(ref).IsVarying();
    this->num_elem_ = num_elem;
    this->num_loc_ = -1;
    
    // for pml element, no source
    if (var_problem.InsidePML(num_elem))
      {
        varying_rho = false;
        varying_mu = false;
        coef_vol = 0;
        coef_grad.Fill(zero);
      }

    //DISP(coef_vol); DISP(coef_grad); DISP(rho0); DISP(mu0); DISP(rho); DISP(mu);
  }

  
  //! volumetric source if rho <> rho0
  template<class Dimension>
  bool DiffractedWaveSource_Helm<Dimension>::IsNonNullVolumetricSource(const VectR_N& s)
  {
    // if rho different from rho 0, we have a volumic source \int f \phi
    if ((coef_vol != Complex_wp(0))||varying_rho)
      return true;
    
    return false;
  }
  
  
  //! Dirichlet condition
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    // incident wave 
    incident_wave.EvaluateFunction(x, f(0));
    
    // for Dirichlet condition u = -u_inc
    if (this->num_loc_ >= 0)
      {
        f(0) = -f(0);
        // for PML, no source
        if (var_problem.InsidePML(i))
          SetComplexZero(f(0));
      }
  }
  
  
  //! evaluation of volumetric source
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    // volumetric source -> omega2 \int (rho-rho_0) u^inc \phi
    incident_wave.EvaluateFunction(x, f(0));
    if (varying_rho)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
	rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, i, j);
        rho += Iwp*var_helm.ref_sigma(ref).GetCoefficient(var_problem, i, j)
          /var_problem.GetOmega();
        
	coef_vol = var_problem.GetSquareOmega()*(rho - rho0);
      }
    
    f(0) *= coef_vol;
  }
  
  
  //! volumetric source for \f$ \nabla phi \f$ if mu <> mu0
  template<class Dimension>
  bool DiffractedWaveSource_Helm<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    // gradient source if mu different from mu_0
    if ((coef_grad.IsZero())&&(!varying_mu))
      return false;
    
    return true;
  }

  
  //! evaluation of volumetric source g in \f$ \int g \nabla \varphi \f$
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    // gradient source -> \int (mu_0 - mu) \nabla u^{inc} \cdot \nabla \phi exp(ikx)
    R_N_Complex_wp grad_u_inc, fgrad; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_u_inc);
    
    // coef_grad = mu_0 - mu with mu complex
    Complex_wp zero; SetComplexZero(zero);
    if (varying_mu)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
	mu = var_helm.ref_mu(ref).GetCoefficient(var_problem, i, j);
	coef_grad.SetDiagonal(mu0);
        coef_grad -= mu;
      }
    
    Mlt(coef_grad, grad_u_inc, fgrad);
    CopyVector(fgrad, 0, f);
  }
  
  
  //! initialization before evaluation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);

    VectR_N s; 
    var_problem.mesh.GetVerticesElement(num_elem, s);
    incident_wave.InitElement(num_elem, s);
    
    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);

    // for neumann condition and high conductivity condition
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    mu = var_helm.ref_mu(ref).GetConstant();
    coef_mu0 = -mu0;
    coef_mu = mu(0, 0);
    if (var_problem.InsidePML(num_elem))
      {
        coef_mu = 0;
        coef_mu0 = 0;
      }
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    ref = var_problem.mesh.BoundaryRef(i).GetReference();
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if (cond_ref == BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY)
      epsilon = var_boundary.GetParamCondition(ref, 0);
#endif
  }

  
  //! surfacic for Neumann and High-Conductivity boundary conditions
  template<class Dimension>
  bool DiffractedWaveSource_Helm<Dimension>
  ::IsNonNullSurfacicSource(int ref)
  {    
    // for neumann condition and high conductivity condition
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref;
    if ( (cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
	 || (cond_ref== BoundaryConditionEnum::LINE_IMPEDANCE)
	 || (cond_ref == BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY))
      return true;
    
    return false;
  }
  
  
  //! evaluation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    R_N x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    R_N_Complex_wp grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    Complex_wp du_inc_dn = DotProd(grad_uinc, MatricesElem.GetNormaleQuadratureBoundary(k));
    SetComplexZero(f(0));
    if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
      {
	// neumann condition -\int \mu_0 du_inc/dn \phi 	
	f(0) = coef_mu0*du_inc_dn;
      }
    else if (reference_condition == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
        // retrieving value of Z
	int ref_domain = var_problem.mesh.Element(this->num_elem_).GetReference();
	Complex_wp imped = fct_imped.GetCoefficient(this->num_bound_ref_, this->num_elem_, this->num_loc_,
						    k, ref_domain, this->ref_boundary_, PointsElem, MatricesElem);
	
	// neumann condition -\int \mu_0 du_inc/dn \phi 	
	f(0) = coef_mu0*du_inc_dn - imped*u_inc;
      }
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    else if (reference_condition == BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY)
      {
	// high conductivity boundary 
	// -\int \mu_0 du_inc/dn \phi - \int mu alpha/epsilon u_inc \phi   for order 1
	// -\int \mu_0 du_inc/dn \phi - \int mu (alpha/epsilon + h) u_inc \phi   for order 2
	// -\int \mu_0 du_inc/dn \phi 
        // - \int mu [alpha/epsilon + h -epsilon/(2 alpha) (h^2 -g +omega2)] u_inc \phi
	//  for order 3
	// h is the mean curvature, g the gaussian curvature
	f(0) = coef_mu0*du_inc_dn;
	Real_wp k1 = MatricesElem.GetK1QuadratureBoundary(k),
	  k2 = MatricesElem.GetK2QuadratureBoundary(k);
	Real_wp h = 0.5*(k1+k2), g = k1*k2;
	
	if (order_condition == 1)
	  f(0) -= coef_mu*alpha/epsilon*u_inc;
	else if (order_condition == 2)
	  f(0) -= coef_mu*(alpha/epsilon+h)*u_inc;
	else if (order_condition == 3)
	  f(0) -= coef_mu*(alpha/epsilon+h-Real_wp(0.5)*epsilon/alpha*(h*h-g+omega2))*u_inc;

      }
#endif
  }

  
  //! surfacic source \f$ \int f \frac{d\varphi}{ds} \f$  for high-conductivity boundary condition
  template<class Dimension>
  bool DiffractedWaveSource_Helm<Dimension>
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    // for high conductivy boundary if order >= 3
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if ( (cond_ref == BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY)&&(order_condition>=3) )
      return true;
#endif
    
    return false;
  }
  
  //! evaluation of surfacic source f in \f$ \int f \frac{d \varphi}{ds} \f$ 
  template<class Dimension>
  void DiffractedWaveSource_Helm<Dimension>
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
				   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    f.Fill(Real_wp(0));
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    if (order_condition == 3)
      {
        // for r=3 -\int epsilon/alpha mu du_inc/ds dphi/ds
        R_N_Complex_wp grad_uinc, fvec; Complex_wp u_inc;
        R_N x = PointsElem.GetPointQuadratureBoundary(k);
        incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
        
        Complex_wp du_ds; MatrixN_Nsym proj;
        GetTangentialProjector(MatricesElem.GetNormaleQuadratureBoundary(k), proj);
        Mlt(proj, grad_uinc, fvec);
        Mlt(-Real_wp(0.5)*coef_mu*epsilon/alpha, fvec);
	CopyVector(fvec, 0, f);
      }
#endif
  }
  

  /************************
   * TotalWaveSource_Helm *
   ************************/


  //! Is the element in the PML ?
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
    num_quad = num_elem;
    element_insidePML = var_problem.InsidePML(num_elem);    
    incident_wave.InitElement(num_elem, s);

    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);
  }

  
  //! for Dirichlet condition
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    if (element_insidePML)
      {
	// incident plane wave 
	incident_wave.EvaluateFunction(x, f(0));
      }
    else
      {
	f(0) = 0.0;
      }
  }
  
  
  //! volumetric source if PML element
  template<class Dimension>
  bool TotalWaveSource_Helm<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return element_insidePML;
  }

  
  //! evaluation of volumetric source
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    int ref = var_problem.mesh.Element(num_quad).GetReference();
    rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, i, j);
    mu = var_helm.ref_mu(ref).GetCoefficient(var_problem, i, j);
    Complex_wp sigma;
    TinyVector<Complex_wp, Dimension::dim_N> v0;
    int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    var_helm.ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, i, i1, j);
    
    incident_wave.EvaluateFunction(x, f(0));
    f(0) *= -var_problem.GetSquareOmega()*(rho-rho0);
  }

  
  //! volumetric source if PML element
  template<class Dimension>
  bool TotalWaveSource_Helm<Dimension>::
  IsNonNullGradientSource(const VectR_N& s)
  {
    return element_insidePML;
  }

  
  //! evaluation of volumetric source f in \f$ \int f \nabla \varphi \f$
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    int ref = var_problem.mesh.Element(num_quad).GetReference();
    Complex_wp sigma;
    TinyVector<Complex_wp, Dimension::dim_N> v0;
    mu = var_helm.ref_mu(ref).GetCoefficient(var_problem, i, j);
    int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    var_helm.ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, i, i1, j);
        
    // gradient source -> \int (mu_0 - mu) \nabla u^{inc} \cdot \nabla \phi exp(ikx)
    typename Dimension::R_N_Complex_wp fgrad, fvec; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, fgrad);
    for (int k = 0; k < Dimension::dim_N; k++)
      mu(k, k) -= mu0;
    
    Mlt(mu, fgrad, fvec);
    CopyVector(fvec, 0, f);
  }

  
  //! initialization before computation of surfacic source
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension>
      ::InitSurface(i, num_face, num_elem, num_loc);

    VectR_N s; 
    var_problem.mesh.GetVerticesElement(num_elem, s);
    incident_wave.InitElement(num_elem, s);

    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);
    
    Complex_wp rho; typename Dimension::MatrixN_Nsym_Complex_wp mu;
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    rho = var_helm.ref_rho(ref);
    mu = var_helm.ref_mu(ref).GetConstant();
    coef_u = omega*sqrt(abs(rho/mu(0,0)));
    num_quad = num_elem;
    element_insidePML = var_problem.InsidePML(num_elem);
  }

  
  //! non null surfacic source for absorbing condition
  template<class Dimension>
  bool TotalWaveSource_Helm<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref; 
    if (cond_ref== BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    if (cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
      return true;
    
    return false;
  }

  
  //! evaluation of surfacic source
  template<class Dimension>
  void TotalWaveSource_Helm<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    f.Fill(0); 
    if (!element_insidePML)
      {
        if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
          return;
      }
    
    int ref_element = var_problem.mesh.Element(this->num_elem_).GetReference();
    Complex_wp imped = 
      fct_imped_abc.GetCoefficient(this->num_bound_ref_, this->num_elem_, this->num_loc_,
				   k, ref_element, this->ref_boundary_, PointsElem, MatricesElem);
    
    coef_u = imped;
   
    // mu_0 \int_\Sigma (du_inc - i k u_inc) 
    typename Dimension::R_N x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc; 
    typename Dimension::R_N_Complex_wp grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
	
    Complex_wp du_inc_dn = DotProd(grad_uinc, MatricesElem.GetNormaleQuadratureBoundary(k));
    if (element_insidePML)
      {
	if (reference_condition == BoundaryConditionEnum::LINE_ABSORBING)
	  f(0) = coef*du_inc_dn + coef_u*u_inc;
	else
          f(0) = coef*du_inc_dn;
      }
    else if (reference_condition == BoundaryConditionEnum::LINE_ABSORBING)
      {		
	f(0) = coef*du_inc_dn + coef_u*u_inc;
      }
  }


  /*******************************
   * DiffractedWaveSource_HelmDG *
   *******************************/

  
  //! Kind of constructor
  template<class Dimension>
  void DiffractedWaveSource_HelmDG<Dimension>::InitDefaultValues()
  {
    coef_rho = 1.0; rho0 = var_helm.rho0; rho = 1.0;
    
    // variables for gradient source
    coef_mu.SetIdentity(); mu0 = var_helm.mu0;
    mu.SetIdentity();
    
    omega = var_problem.GetOmega();
    varying_rho = false; varying_mu = false;
  }
  
  
  //! initialization before computation of volumetric source
  template<class Dimension>
  void DiffractedWaveSource_HelmDG<Dimension>
  ::InitElement(int num_elem, const VectR_N& s)
  {
    incident_wave.InitElement(num_elem, s);
    
    // if the incident wave is defined for a non constant medium we update the value
    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);

    int ref = var_problem.mesh.Element(num_elem).GetReference();
    rho = var_helm.ref_rho(ref).GetConstant()
      + Iwp*var_helm.ref_sigma(ref).GetConstant()/var_problem.GetOmega();
    
    mu = var_helm.ref_invMu(ref).GetConstant();
    
    // coef_rho = omega (rho-rho_0) with rho complex !
    if (rho != rho0)
      coef_rho = var_problem.GetOmega()*(rho-rho0);
    else
      coef_rho = Complex_wp(0);
    
    // coef_grad = mu_0 - mu with mu complex
    Complex_wp zero; SetComplexZero(zero);

    coef_mu.SetDiagonal(Real_wp(1)/mu0);
    coef_mu -= mu;
    coef_mu *= var_problem.GetOmega();

    if (var_problem.FirstOrderFormulation())
      {
	coef_rho *= Iwp;
	coef_mu *= -Real_wp(1);
      }
    
    varying_rho = (var_helm.ref_rho(ref).IsVarying() || var_helm.ref_sigma(ref).IsVarying());
    varying_mu = var_helm.ref_mu(ref).IsVarying();
    
    if (var_problem.InsidePML(num_elem))
      {
        varying_rho = false;
        varying_mu = false;
        coef_rho = 0;
        coef_mu.Fill(zero);
      }
  }
  
  
  //! volumetric source if rho != rho0 or mu != mu0
  template<class Dimension>
  bool DiffractedWaveSource_HelmDG<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    // if rho different from rho 0, we have a volumic source \int f \phi
    if ((coef_rho!=Complex_wp(0)) || (!coef_mu.IsZero()) || varying_rho || varying_mu)
      return true;
    
    return false;
  }
  
  
  //! evaluation of volumetric source
  template<class Dimension>
  void DiffractedWaveSource_HelmDG<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    // volumetric source -> omega2 \int (rho-rho_0) u^inc \phi
    Complex_wp u_inc;
    TinyVector<Complex_wp, Dimension::dim_N> grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    R_N_Complex_wp fvec, fvol;
    fvec = -grad_uinc/var_problem.GetOmega();
    
    if (varying_rho)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
	rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, i, j);
        rho += Iwp*var_helm.ref_sigma(ref).GetCoefficient(var_problem, i, j)
          /var_problem.GetOmega();
        
        coef_rho = var_problem.GetOmega()*(rho-rho0);
	if (var_problem.FirstOrderFormulation())
	  coef_rho *= Iwp;
      }
    
    if (varying_mu)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
	mu = var_helm.ref_invMu(ref).GetCoefficient(var_problem, i, j);
	coef_mu.SetDiagonal(mu0);
        coef_mu -= mu;
        Mlt(var_problem.GetOmega(), coef_mu);
	if (var_problem.FirstOrderFormulation())
	  coef_mu *= -Real_wp(1);
      }
    
    Mlt(coef_mu, fvec, fvol);
    
    f(0) = u_inc*coef_rho;
    for (int k = 0; k < x.GetM(); k++)
      f(k+1) = fvol(k);
    
  }

  
  //! initialization before computation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_HelmDG<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);

    VectR_N s; 
    var_problem.mesh.GetVerticesElement(num_elem, s);
    incident_wave.InitElement(num_elem, s);

    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);
  }
  
  
  //! surfacic source for Neumann or Dirichlet conditions
  template<class Dimension>
  bool DiffractedWaveSource_HelmDG<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    // for neumann condition and high conductivity condition
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref;
    if ((cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
	||(cond_ref == BoundaryConditionEnum::LINE_DIRICHLET)
	||(cond_ref == BoundaryConditionEnum::LINE_IMPEDANCE))
      return true;
    
    return false;
  }

  
  //! evaluation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_HelmDG<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    R_N x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    R_N_Complex_wp grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);

    R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp du_inc_dn = DotProd(grad_uinc, normale);
    f.Zero();

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int num_point = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
	num_point = var_problem.GetNbPointsQuadratureInside(this->num_elem_);
	for (int pos = 0; pos < this->num_loc_; pos++)
	  num_point += mesh_num.
	    GetNbPointsQuadratureBoundary(var_problem.mesh.Element(this->num_elem_).numBoundary(pos));
      }
    
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
      {
	if (var_problem.FirstOrderFormulation())
	  {
	    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	      f(0) = -Iwp*mu0*du_inc_dn/omega;
	    else
	      f(0) = Iwp*mu0*du_inc_dn/omega;
	  }
	else
	  f(0) = -mu0*du_inc_dn/omega;
      }
    else if (reference_condition == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
	Complex_wp imped; SetComplexOne(imped);
	var_boundary.MltParamCondition(this->ref_boundary_, 0, imped);
	if (var_problem.FirstOrderFormulation())
	  f(0) = Iwp*mu0*(du_inc_dn + imped*u_inc)/omega;
	else
	  f(0) = -mu0*(du_inc_dn + imped*u_inc)/omega;
	
	if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	  {
	    Complex_wp invMu = var_helm.ref_invMu(ref).GetCoefficient(var_problem, this->num_elem_, num_point+k)(0, 0);
		
	    if (var_problem.FirstOrderFormulation())
	      f(0) *= -Real_wp(0.5);
	    else
	      f(0) *= Real_wp(0.5);
	    
	    Complex_wp f0 = f(0)*invMu/imped*omega;
	    if (var_problem.FirstOrderFormulation())
	      f0 = Iwp*f0;
	    
	    for (int i = 0; i < x.GetM(); i++)
	      f(i+1) = -f0*normale(i);
	  }
      }
    else
      {
	if (var_problem.FirstOrderFormulation())
	  for (int i = 0; i < x.GetM(); i++)
	    f(i+1) = -u_inc*normale(i);
	else
	  for (int i = 0; i < x.GetM(); i++)
	    f(i+1) = u_inc*normale(i);
      }
    
    if (var_problem.InsidePML(this->num_elem_))
      f.Zero();
  }
  
  
  /**************************
   * TotalWaveSource_HelmDG *
   **************************/

  
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>::InitDefaultValues()
  {
    omega = var_problem.GetOmega();
    rho0 = var_helm.rho0; mu0 = var_helm.mu0;
    element_insidePML = false;
    k_wave = var_problem.GetWaveVector();
  }
  

  //! Is the element in the PML ?
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
    num_quad = num_elem;
    element_insidePML = var_problem.InsidePML(num_elem);

    incident_wave.InitElement(num_elem, s);
    
    // if the incident wave is defined for a non constant medium we update the value
    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);
  }
  
  
  //! nothing special
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);

    VectR_N s; 
    var_problem.mesh.GetVerticesElement(num_elem, s);
    incident_wave.InitElement(num_elem, s);
    
    rho0 = var_helm.rho0;
    mu0 = var_helm.mu0;
    incident_wave.UpdateCoefAB(mu0, rho0);
  }

  
  //! is element inside PML ?
  template<class Dimension>
  bool TotalWaveSource_HelmDG<Dimension>::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return element_insidePML;
  }
  
  
  //! evaluation of volumetric source
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    int ref = var_problem.mesh.Element(num_quad).GetReference();
    rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, i, j);
    invMu = var_helm.ref_invMu(ref).GetCoefficient(var_problem, i, j);
    Complex_wp sigma;
    TinyVector<Complex_wp, Dimension::dim_N> v0;
    int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    var_helm.ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, i, i1, j);
    
    Complex_wp u_inc;
    typename Dimension::R_N_Complex_wp fvec, fgrad;
    incident_wave.EvaluateFunctionGradient(x, u_inc, fvec);
    for (int k = 0; k < Dimension::dim_N; k++)
      invMu(k, k) -= 1.0/mu0;

    f(0) = -var_problem.GetOmega()*(rho-rho0)*u_inc;

    Mlt(invMu, fvec, fgrad);
    for (int i = 0; i < x.GetM(); i++)
      f(i+1) = -fgrad(i);
  }
  

  //! volumetric source for \f$ \nabla phi \f$ if mu <> mu0
  template<class Dimension>
  bool TotalWaveSource_HelmDG<Dimension>::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
  

  //! evaluation of volumetric source g in \f$ \int g \nabla \varphi \f$
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    f.Fill(0);
  }
  
  
  //! surfacic source for absorbing boundary condition
  template<class Dimension>
  bool TotalWaveSource_HelmDG<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref;
    if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    if ((cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
	|| (cond_ref == BoundaryConditionEnum::LINE_DIRICHLET))
      return true;
    
    return false;
  }

  
  //! evaluation of surfacic source
  template<class Dimension>
  void TotalWaveSource_HelmDG<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    f.Fill(0);
    TinyVector<Complex_wp, Dimension::dim_N> fac_tau;
    fac_tau.Fill(1.0);
    if (!var_problem.InsidePML(this->num_elem_))
      {
        if (reference_condition != BoundaryConditionEnum::LINE_ABSORBING)
          return;
      }

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    typename Dimension::R_N x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    typename Dimension::R_N_Complex_wp grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp du_inc_dn = DotProd(grad_uinc, normale);    

    if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
      {
	if (var_problem.FirstOrderFormulation())
	  f(0) = -Iwp*mu0*du_inc_dn/omega;
	else
	  f(0) = -mu0*du_inc_dn/omega;
      }
    else if (reference_condition == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        if (var_problem.FirstOrderFormulation())
	  for (int i = 0; i < x.GetM(); i++)
	    f(i+1) = u_inc*normale(i);
	else
	  for (int i = 0; i < x.GetM(); i++)
	    f(i+1) = -u_inc*normale(i);
      }
    else
      {
	int num_point = var_problem.GetNbPointsQuadratureInside(this->num_elem_);
	for (int pos = 0; pos < this->num_loc_; pos++)
	  num_point += mesh_num.
	    GetNbPointsQuadratureBoundary(var_problem.mesh.Element(this->num_elem_).numBoundary(pos));

	int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
	Complex_wp rho = var_helm.GetRhoTilde(ref, this->num_elem_, num_point+k);
	Complex_wp invMu = var_helm.ref_invMu(ref).GetCoefficient(var_problem, this->num_elem_, num_point+k)(0, 0);

        Complex_wp vloc = (du_inc_dn/omega - Iwp*sqrt(rho*invMu)*u_inc)/invMu;
        if (var_problem.FirstOrderFormulation())
	  {
	    vloc *= Iwp;
	    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
	      vloc = -vloc;
	  }

	f(0) = vloc;
	vloc *= sqrt(invMu/rho);
	
	if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	  {
	    if (var_problem.FirstOrderFormulation())
	      for (int i = 0; i < x.GetM(); i++)
		f(i+1) = vloc*normale(i);
	    else
	      for (int i = 0; i < x.GetM(); i++)
		f(i+1) = -Iwp*vloc*normale(i);

	    Mlt(Real_wp(0.5), f);
	  }
      }    
  }
  
  
  //! incident wave for DG
  template<class Dimension>
  void IncidentWaveProjector_HelmDG<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    Complex_wp u_inc;
    TinyVector<Complex_wp, Dimension::dim_N> grad_uinc;
    this->incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    f(0) = u_inc;
    Complex_wp coef_v = -var_helm.mu0/var_problem.GetOmega();
    for (int n = 0; n < x.GetM(); n++)
      f(n+1) = coef_v*grad_uinc(n);
  }


  /***************************
   * VolumetricSource_HelmDG *
   ***************************/
      
  template<class T, class Dimension>
  bool VolumetricSource_HelmDG<T, Dimension>::IsNonNullDirichletSource(int cond_ref)
  {
    return false;
  }

  
  template<class T, class Dimension>
  void VolumetricSource_HelmDG<T, Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
  {
    VolumetricSource_Base<T, Dimension>::EvaluateSurfacicSource(k, PointsElem, MatricesElem, f);
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
 
    // only first component is kept
    T f0 = f(0);
    f.Zero();

    if (var_problem.mesh.GetBoundaryCondition(this->ref_boundary_)
        == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// case of an inhomogeneous Dirichlet
	typename Dimension::R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
	if (var_problem.FirstOrderFormulation())
	  for (int p = 1; p <= Dimension::dim_N; p++)
	    f(p) = f0*normale(p-1);
	else
	  for (int p = 1; p <= Dimension::dim_N; p++)
	    f(p) = -f0*normale(p-1);
      }
    else
      {
	T coef_v(1);
	if (sizeof(T) != sizeof(Real_wp))
	  {
	    coef_v = Real_wp(1) / var_problem.GetOmega();
	    if (var_problem.FirstOrderFormulation())
	      to_complex(-Iwp*coef_v, coef_v);
	  }
	
	// case of an inhomogeneous Neumann or Robin condition
	if (var_problem.FormulationDG() == ElementReference_Base::HDG)
	  {
	    if ((sizeof(T) == sizeof(Real_wp)) && (var_problem.FirstOrderFormulation()))
	      coef_v = -coef_v;
	    
	    f(0) = coef_v*f0;
	  }
	else
	  {
	    if ((sizeof(T) != sizeof(Real_wp)) && var_problem.FirstOrderFormulation())
	      coef_v = -coef_v;	    
	    
	    // for LDG, the source depends on the boundary condition
	    int cond = var_problem.mesh.GetBoundaryCondition(this->ref_boundary_);
	    if (cond == BoundaryConditionEnum::LINE_NEUMANN)
	      f(0) = coef_v*f0;
	    else
	      {
		int num_point = var_problem.GetNbPointsQuadratureInside(this->num_elem_);
		for (int pos = 0; pos < this->num_loc_; pos++)
		  num_point += mesh_num.
		    GetNbPointsQuadratureBoundary(var_problem.mesh.Element(this->num_elem_).numBoundary(pos));
		
		int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
		T invMu = var_helm.ref_invMu(ref).GetCoefficient(var_problem, this->num_elem_, num_point)(0, 0);
		typename Dimension::R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
		
		f(0) = Real_wp(0.5)*coef_v*f0;
		
		if (cond == BoundaryConditionEnum::LINE_IMPEDANCE)
		  {
		    T Z; SetComplexOne(Z);
		    var_boundary.MltParamCondition(this->ref_boundary_, 0, Z);
		    f0 *= Real_wp(0.5)*invMu/Z;
		    if (var_problem.FirstOrderFormulation())
		      f0 = -f0;
		    
		    for (int p = 1; p <= Dimension::dim_N; p++)
		      f(p) = -f0*normale(p-1);
		  }
		else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
		  {
		    T rho = var_helm.GetRhoTilde(ref, this->num_elem_, num_point), Z;
		    if (sizeof(T) != sizeof(Real_wp))
		      to_complex(Complex_wp(0, 0.5)/ var_problem.GetOmega(), Z);
		    else
		      SetComplexReal(Real_wp(0.5), Z);

		    f0 *= Z*sqrt(invMu/rho);
		    if (var_problem.FirstOrderFormulation())
		      f0 = -f0;
		    
		    for (int p = 1; p <= Dimension::dim_N; p++)
		      f(p) = -f0*normale(p-1);		
		  }
	      }
	  }
      }
  }
  
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_CXX
#endif
