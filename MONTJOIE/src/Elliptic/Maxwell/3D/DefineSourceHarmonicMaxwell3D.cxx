#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_CXX

namespace Montjoie
{
  //! constructor with given problem
  void DiffractedWaveSource_Maxwell3D
  ::InitDefaultValues()
  {
    // other initializations
    coef_vol.SetIdentity(); epsilon.SetIdentity();
    
    var_problem.GetPolarization(E_0);    
    coef_curl.SetIdentity(); invMu.SetIdentity(); 
    
    varying_eps = false;
    varying_mu = false;
  }
  
  
  //! computation of coefficients for volumetric source
  void DiffractedWaveSource_Maxwell3D
  ::InitElement(int num_elem, const VectR3& s)
  {
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    TinyMatrix<Complex_wp, Symmetric, 3, 3> sigma_p, mu_p;

    epsilon = var_maxwell.ref_epsilon(ref).GetConstant();
    mu_p = var_maxwell.ref_mu(ref).GetConstant();

    var_maxwell.ref_drude(ref)
      .ModifyCoefficientMaxwell(var_problem.GetOmega(), epsilon, sigma_p, mu_p);

    invMu = mu_p; GetInverse(invMu);
    
    coef_vol.SetDiagonal(-var_maxwell.epsilon0);
    coef_vol += epsilon;

    if (var_problem.FirstOrderFormulation())
      {
        Mlt(Iwp*var_problem.GetOmega(), coef_vol);
        coef_vol += var_maxwell.ref_sigma(ref).GetConstant();
      }
    else
      {
        Mlt(var_problem.GetSquareOmega(), coef_vol);
        Add(Iwp*var_problem.GetOmega(), var_maxwell.ref_sigma(ref).GetConstant(), coef_vol);
      }

    coef_curl.SetDiagonal(-var_maxwell.invMu0);
    coef_curl += invMu;
    
    varying_eps = var_maxwell.ref_epsilon(ref).IsVarying() || var_maxwell.ref_sigma(ref).IsVarying();
    varying_mu = var_maxwell.ref_mu(ref).IsVarying();
    
    this->num_elem_ = num_elem; this->num_loc_ = -1;
  }
  
  
  //! computation of incident electric field
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] fvec resulting electric field
   */
  void DiffractedWaveSource_Maxwell3D
  ::EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& fvec)
  {
    // inhomogeneous Dirichlet condition
    // E = -E^{inc}  (on tangential dofs)
    fvec.Fill(0);
    if (!var_problem.InsidePML(this->num_elem_))
      {
	Complex_wp u_inc;
	incident_wave.EvaluateFunction(x, u_inc);
	
	fvec(0) = -u_inc*E_0(0);
	fvec(1) = -u_inc*E_0(1);
	fvec(2) = -u_inc*E_0(2);
      }
  }
  
  
  //! volumetric source for dielectric material
  /*!
    \param[in] s vertices
   */
  bool DiffractedWaveSource_Maxwell3D
  ::IsNonNullVolumetricSource(const VectR3& s)
  {
    // there is a volumetric integral if epsilon is different from epsilon0
    // (or sigma different from 0) 
    if ((coef_vol.IsZero())&&(!varying_eps))
      return false;
    
    return true;
  }
  
  
  //! evaluation of volumetric source
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] fvec source vector
   */
  void DiffractedWaveSource_Maxwell3D
  ::EvaluateVolumetricSource(int i, int j, const R3& x,
			     VectComplex_wp& fvec)
  {
    R3_Complex_wp E_inc = E_0, f;
    
    Complex_wp u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
    E_inc *= u_inc;
    
    // source is ( \omega^2 (eps-eps0) + i omega sigma) E^inc
    if (varying_eps)
      {
        // case where epsilon or sigma depends on the quadrature point
	int ref = var_problem.mesh.Element(i).GetReference();
        epsilon.SetDiagonal(-var_maxwell.epsilon0);
	epsilon += var_maxwell.ref_epsilon(ref).GetCoefficient(var_problem, i, j);
	
	Mlt(var_problem.GetSquareOmega(), epsilon);
        Add(Iwp*var_problem.GetOmega(), var_maxwell.ref_sigma(ref)
            .GetCoefficient(var_problem, i, j), epsilon);
        
	Mlt(epsilon, E_inc, f);
      }
    else
      Mlt(coef_vol, E_inc, f);
    
    CopyVector(f, 0, fvec);
  }
  
  
  //! volumetric source \f$ \int f curl \varphi \f$  if  mu <> mu0
  /*!
    \param[in] s vertices
   */
  bool DiffractedWaveSource_Maxwell3D
  ::IsNonNullGradientSource(const VectR3& s)
  {
    // non-null term if mu is different from mu0
    if ((coef_curl.IsZero())&&(!varying_mu))
      return false;
    
    return true;
  }
  
  
  //! volumetric source \f$ \int f curl(\varphi) \f$
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] fvec source vector
   */
  void DiffractedWaveSource_Maxwell3D
  ::EvaluateGradientSource(int i, int j, const R3& x, VectComplex_wp& fvec)
  {
    R3_Complex_wp H_inc, f;    
    R3_Complex_wp grad_uinc; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    TimesProd(E_0, grad_uinc, H_inc);
    
    // source -(1/mu - 1/mu0) curl(E^inc)
    if (varying_mu)
      {
        // case where mu depends on the quadrature point
	int ref = var_problem.mesh.Element(i).GetReference();
        invMu.SetDiagonal(-var_maxwell.invMu0);
	invMu += var_maxwell.ref_invMu(ref).GetCoefficient(var_problem, i, j);
	
	Mlt(invMu, H_inc, f);
      }
    else
      Mlt(coef_curl, H_inc, f);
    
    CopyVector(f, 0, fvec);
  }
  
  
  //! surfacic source for Neumann or impedance boundary condition
  bool DiffractedWaveSource_Maxwell3D
  ::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if ((cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
	||(cond_ref == BoundaryConditionEnum::LINE_IMPEDANCE))
      return true;
    
    return false;
  }
  
  
  //! empty method
  void DiffractedWaveSource_Maxwell3D
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension3>::
      InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! evaluation of surfacic source
  void DiffractedWaveSource_Maxwell3D
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
			   const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f)
  {
    // for Neumann condition :
    // \int_\Gamma 1/mu0 n \times curl(E_inc) \cdot \varphi
    R3_Complex_wp rotE_inc;
    R3_Complex_wp E_inc = E_0, Esrc;
    R3 x = PointsElem.GetPointQuadratureBoundary(k);
    R3 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    
    Complex_wp u_inc;
    R3_Complex_wp grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    E_inc *= u_inc;
    
    TimesProd(grad_uinc, E_0, rotE_inc);
    int cond_ref = var_problem.mesh.GetBoundaryCondition(this->ref_boundary_);
    if (cond_ref == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
        // for impedance condition :
        // \int_\Gamma (1/mu0 n \times curl(E_inc) - Z (n \times E^inc) \times n) \cdot \varphi
        // Z is the impedance and retrieved by calling EvaluateImpedancePhi
	
        // Esrc = E^{inc} \times n
        TimesProd(E_inc, normale, Esrc);
        
        // retrieving value of Z
        int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
	Complex_wp imped = fct_imped.GetCoefficient(this->num_bound_ref_, this->num_elem_, this->num_loc_,
						    k, ref, this->ref_boundary_, PointsElem, MatricesElem);
	
        // Z_Etangent = Z (E^{inc} \times n)
        R3_Complex_wp Z_Etangent;
        Z_Etangent = imped*Esrc; 
        Mlt(var_maxwell.invMu0, rotE_inc);
	
        // rotE_inc = 1/mu0 curl(E^{inc}) - Z (E^{inc} \times n) 
        rotE_inc -= Z_Etangent;
      }
    else 
      Mlt(var_maxwell.invMu0, rotE_inc);
    
    // f = n \times rotE_inc
    TimesProd(normale, rotE_inc, Esrc);
    CopyVector(Esrc, 0, f);
  }
  
  
  //! Sets values of attributes to default values
  void TotalWaveSource_Maxwell3D::
  InitDefaultValues()
  {
    var_problem.GetPolarization(E_0);    
    Real_wp c0 = var_maxwell.GetVelocityOfInfinity();
    kwave = var_problem.GetOmega() / c0;    
  }
  
  
  //! Is the element in the PML ?
  void TotalWaveSource_Maxwell3D
  ::InitElement(int num_elem, const VectR3& s)
  {
    this->num_elem_ = num_elem; this->num_loc_ = -1;
    element_insidePML = var_problem.InsidePML(num_elem);
  }
  
  
  //! for Dirichlet condition
  void TotalWaveSource_Maxwell3D
  ::EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f)
  {
    f.Zero();
    if (element_insidePML)
      {
	// incident plane wave, we impose E = E^inc on the boundary of the PML
	Complex_wp u_inc;
	incident_wave.EvaluateFunction(x, u_inc);

	R3_Complex_wp fvec;
	fvec = E_0;
	
	Mlt(u_inc, fvec);
	CopyVector(fvec, 0, f);
      }
  }
  
  
  //! volumetric source if inside a PML layer
  bool TotalWaveSource_Maxwell3D
  ::IsNonNullVolumetricSource(const VectR3& s)
  {
    return element_insidePML;
  }
  
  
  //! evaluation of volumetric source
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] f vector source
  */
  void TotalWaveSource_Maxwell3D
  ::EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f)
  {
    // source \omega^2 (eps-eps0) E^{inc}
    int ref = var_problem.mesh.Element(i).GetReference();
    epsilon = var_maxwell.ref_epsilon(ref).GetCoefficient(var_problem, i, j);
    mu = var_maxwell.ref_mu(ref).GetCoefficient(var_problem, i, j);
    sigma = var_maxwell.ref_sigma(ref).GetCoefficient(var_problem, i, j);

    int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    var_maxwell.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, j);
    
    Complex_wp u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
    epsilon(0, 0) -= var_maxwell.epsilon0;
    epsilon(1, 1) -= var_maxwell.epsilon0;
    epsilon(2, 2) -= var_maxwell.epsilon0;
    Complex_wp coef = Iwp/var_problem.GetOmega();
    epsilon += coef*sigma;
    
    R3_Complex_wp fvec;
    Mlt(epsilon, E_0, fvec);
    Mlt(var_problem.GetSquareOmega()*u_inc, fvec);

    CopyVector(fvec, 0, f);
  }
  
  
  //! volumetric source if inside a PML layer
  bool TotalWaveSource_Maxwell3D
  ::IsNonNullGradientSource(const VectR3& s)
  {
    return element_insidePML;
  }
  
  
  //! evaluation of volumetric source involving a curl
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] f vector source
  */
  void TotalWaveSource_Maxwell3D
  ::EvaluateGradientSource(int i, int j, const R3& x, VectComplex_wp& f)
  {
    int ref = var_problem.mesh.Element(i).GetReference();
    epsilon = var_maxwell.ref_epsilon(ref).GetCoefficient(var_problem, i, j);
    mu = var_maxwell.ref_mu(ref).GetCoefficient(var_problem, i, j);
    sigma = var_maxwell.ref_sigma(ref).GetCoefficient(var_problem, i, j);

    int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    var_maxwell.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, j);
    GetInverse(mu, invMu);

    R3_Complex_wp rotE_inc;
    // source (1/mu - 1/mu0) curl(E^inc)
     
    R3_Complex_wp grad_uinc; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    TimesProd(grad_uinc, E_0, rotE_inc);
    
    // gradient source -> \int (mu_0 - mu) H^{inc} \cdot \nabla \times \phi
    invMu(0, 0) -= var_maxwell.invMu0;
    invMu(1, 1) -= var_maxwell.invMu0;
    invMu(2, 2) -= var_maxwell.invMu0;
    
    R3_Complex_wp fvec;
    Mlt(invMu, rotE_inc, fvec);
    CopyVector(fvec, 0, f);
  }
  
  
  //! initialization before evaluation of surfacic source
  void TotalWaveSource_Maxwell3D
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension3>::InitSurface(i, num_face, num_elem, num_loc);
    
    // int ref = this->vars.mesh.elements(num_elem).GetReference();
    // epsilon = this->vars.ref_epsilon(ref).GetConstant();
    // mu = this->vars.ref_mu(ref).GetConstant();
    // coef_u = omega*sqrt(abs(epsilon(0,0)*mu(0,0)));
    element_insidePML = var_problem.InsidePML(num_elem);
  }
  
  
  //! surfacic source for absorbing boundary condition
  bool TotalWaveSource_Maxwell3D
  ::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref; 
    //DISP(cond_ref); DISP(BoundaryConditionEnum::LINE_ABSORBING);
    if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    if (element_insidePML && (cond_ref== BoundaryConditionEnum::LINE_NEUMANN))
      return true;
    
    return false;
  }
  
  
  //! evaluation of surfacic source
  /*!
    \param[in] k quadrature point number in the edge
    \param[in] PointsElem transformation of points with respect to Fi
    \param[in] MatricesElem jacobian matrices DF_i
    \param[out] f source vector
   */
  void TotalWaveSource_Maxwell3D
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
			   const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f)
  {
    R3_Complex_wp fvec;
    if (reference_condition == BoundaryConditionEnum::LINE_ABSORBING)
      {		
        R3 x = PointsElem.GetPointQuadratureBoundary(k);
        R3 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
        
        R3_Complex_wp grad_uinc; Complex_wp u_inc;
        incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
	
        // E^inc  = normale ^ ( E^0 u^inc) ^ normale
        R3_Complex_wp E_src = E_0; Mlt(u_inc, E_src);
        Matrix3_3sym ProjTangent; R3_Complex_wp E_inc;
        GetTangentialProjector(normale, ProjTangent);
        Mlt(ProjTangent, E_src, E_inc); //DISP(normale); DISP(E_src); DISP(E_inc);
        
        // rot E^inc = grad_uinc \times E0
        R3_Complex_wp rotE_inc;                     
        TimesProd(grad_uinc, E_0, rotE_inc); //DISP(rotE_inc);
        
        //int ref = this->vars.mesh.Element(this->num_elem_).GetReference();
        Complex_wp coef = Iwp*kwave;
        //DISP(x); DISP(u_inc); DISP(coef); DISP(this->vars.coefficient_impedance_absorbing(ref));
        // -\int  ( n \times rot E^inc + coef (n \times E^inc) \times n ) \cdot phi
        TimesProd(normale, rotE_inc, fvec); // DISP(f);
        Add(coef, E_inc, fvec); //DISP(E_inc);
        Mlt(-var_maxwell.invMu0, fvec);  //DISP(f);      
        // DISP(this->ref_boundary_); DISP(x); DISP(E_inc); DISP(rotE_inc); DISP(f);
      } 
    
    CopyVector(fvec, 0, f);
  }
  
  
  //! default constructor
  template<class T>
  MaxwellGaussianSource<T>::MaxwellGaussianSource()
    : GaussianSourceField<T, Dimension3>()
  {
    presence_curl = false;
  }
  

  template<class T>
  void MaxwellGaussianSource<T>
  ::Init(const R3& x, const Real_wp& r, const Real_wp& rmax, bool curl, Vector<T>& P)
  {
    GaussianSource<Dimension3>::Init(x, r, rmax);
    this->polar = P;
    presence_curl = curl;
  }
  
  
  //! evaluation of the volume source
  template<class T>
  void MaxwellGaussianSource<T>::EvaluateFunction(const R3 & x, Vector<T>& fvec) const
  {
    if (!presence_curl)
      {
	GaussianSourceField<T, Dimension3>::EvaluateFunction(x, fvec);
	return;
      }
    
    Real_wp f; R3 grad_f;
    this->GetGradAmplitude(x, f, grad_f);

    SetComplexZero(fvec(0)); SetComplexZero(fvec(1)); SetComplexZero(fvec(2));
    if (f != Real_wp(0))
      {
	TinyVector<T, 3> g, polarization;
	CopyVector(this->polar, 0, polarization);
	TimesProd(grad_f, polarization, g);
	CopyVector(g, 0, fvec);
      }
  }
  
  
  //! constructor with given problem
  void DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::InitDefaultValues()
  {
    coef_vol.SetIdentity(); epsilon.SetIdentity();
    
    var_problem.GetPolarization(E_0);    
    
    TimesProd(var_problem.GetWaveVector(), E_0, H_0);
    coef_curl.SetIdentity(); invMu.SetIdentity(); 
    
    varying_eps = false; varying_mu = false;
  }
  
  
  //! computation of coefficients before evaluation of volumetric sources
  void DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::InitElement(int num_elem, const VectR3& s)
  {
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    epsilon = var_maxwell.ref_epsilon(ref).GetConstant();
    mu = var_maxwell.ref_mu(ref).GetConstant();
    
    coef_vol.SetDiagonal(-var_maxwell.epsilon0);
    coef_vol += epsilon;
    
    Mlt(var_problem.GetOmega(), coef_vol);
    
    coef_curl.SetDiagonal(-var_maxwell.mu0);
    coef_curl += mu;
    
    Mlt(var_problem.GetOmega(), coef_curl);
	
    varying_eps = var_maxwell.ref_epsilon(ref).IsVarying() || var_maxwell.ref_sigma(ref).IsVarying();
    varying_mu = var_maxwell.ref_mu(ref).IsVarying();
    
    this->num_elem_ = num_elem; this->num_loc_ = -1;
  }
  
  
  //! volumetric source if dielectric
  bool DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::
  IsNonNullVolumetricSource(const VectR3& s)
  {
    if ((coef_vol.IsZero())&&(!varying_eps)&&(!varying_mu))
      return false;
    
    return true;
  }
  
  
  //! evaluation of volumetric source
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] fvec vector source
  */
  void DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::
  EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& fvec)
  {
    R3_Complex_wp E_inc = E_0, E_src, H_inc = H_0, H_src;
    
    Complex_wp u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
    
    Mlt(u_inc, E_inc);  Mlt(u_inc, H_inc);
    if (varying_eps)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
        epsilon.SetDiagonal(-var_maxwell.epsilon0);
	epsilon += var_maxwell.ref_epsilon(ref).GetCoefficient(var_problem, i, j);
	
	Mlt(var_problem.GetOmega(), epsilon);
	Mlt(epsilon, E_inc, E_src);
      }
    else
      Mlt(coef_vol, E_inc, E_src);
    
    if (varying_mu)
      {
	int ref = var_problem.mesh.Element(i).GetReference();
        mu.SetDiagonal(-var_maxwell.mu0);
	mu += var_maxwell.ref_mu(ref).GetCoefficient(var_problem, i, j);
	
	Mlt(var_problem.GetOmega(), mu);
	Mlt(mu, H_inc, H_src);
      }
    else
      Mlt(coef_curl, H_inc, H_src);
    
    CopyVector(E_src, 0, fvec);
    CopyVector(H_src, 1, fvec);
  }
  
  //! surfacic source for Dirichlet or Neumann condition
  bool DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::
  IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if ((cond_ref==BoundaryConditionEnum::LINE_NEUMANN)||(cond_ref==BoundaryConditionEnum::LINE_DIRICHLET))
      return true;
    
    return false;
  }
  
  //! empty method
  void DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension3>::InitSurface(i, num_face, num_elem, num_loc);
  }
  
  //! to be done
  void DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>::
  EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
			  const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f)
  {
  }
  
  //! constructor with given problem
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>::InitDefaultValues()
  {
    var_problem.GetPolarization(E_0);    
    
    TimesProd(var_problem.GetWaveVector(), E_0, H_0);
    Mlt(1.0/var_problem.GetOmega(), H_0);
    element_insidePML = false;
  }
  
  
  // Is the element in the PML ?
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::InitElement(int num_elem, const VectR3& s)
  {
    this->num_elem_ = num_elem; this->num_loc_ = -1;
    element_insidePML = var_problem.InsidePML(num_elem);
  }
  
  //! to be done
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f)
  {
    
  }
  
  //! to be done
  bool TotalWaveSource<HarmonicMaxwellEquation_3D_DG>::
  IsNonNullVolumetricSource(const VectR3& s)
  {
    return element_insidePML;
  }
  
  //! to be done
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f)
  {
  }
  
  //! to be done
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension3>::InitSurface(i, num_face, num_elem, num_loc);
    
    // int ref = this->vars.mesh.elements(num_elem).GetReference();
    // epsilon = this->vars.ref_epsilon(ref).GetConstant();
    // mu = this->vars.ref_mu(ref).GetConstant();
    // coef_u = omega*sqrt(abs(epsilon(0,0)*mu(0,0)));
    
    element_insidePML = var_problem.InsidePML(num_elem);
  }
  
  //! to be done
  bool TotalWaveSource<HarmonicMaxwellEquation_3D_DG>::
  IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref; 
    if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    if (element_insidePML&&(cond_ref==BoundaryConditionEnum::LINE_NEUMANN))
      return true;
    
    return false;
  }
  
  //! to be done
  void TotalWaveSource<HarmonicMaxwellEquation_3D_DG>::
  EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
			  const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f)
  {
    R3 x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
        
    R3_Complex_wp H_inc = H_0, E_inc = E_0, Einc_N, Hinc_N, H_src;
    Mlt(u_inc, H_inc); Mlt(u_inc, E_inc);
    TimesProd(E_inc, MatricesElem.GetNormaleQuadratureBoundary(k), Einc_N);
    TimesProd(MatricesElem.GetNormaleQuadratureBoundary(k), Einc_N, E_inc);
    TimesProd(H_inc, MatricesElem.GetNormaleQuadratureBoundary(k), Hinc_N);
    TimesProd(MatricesElem.GetNormaleQuadratureBoundary(k), Hinc_N, H_inc);
    
    f(0) = -Iwp*(E_inc(0) + Hinc_N(0));
    f(1) = -Iwp*(E_inc(1) + Hinc_N(1));
    f(2) = -Iwp*(E_inc(2) + Hinc_N(2));
    
    f(3) = Einc_N(0) - H_inc(0);
    f(4) = Einc_N(1) - H_inc(1);
    f(5) = Einc_N(2) - H_inc(2);
  }
  

  //! evaluation of incident wave
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] fvec vector source
  */
  void IncidentWaveProjector_Maxwell3D_DG
  ::EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& fvec)
  {
    R3_Complex_wp E_inc = E_0;
    
    Complex_wp u_inc;
    this->incident_wave.EvaluateFunction(x, u_inc);
    
    Mlt(u_inc, E_inc);
    
    R3_Complex_wp H_inc = H_0;
    Mlt(u_inc, H_inc);
    
    CopyVector(E_inc, 0, fvec);
    CopyVector(H_inc, 1, fvec);
  }


  template<class Complexe>
  CircularCoaxialModeField<Complexe>
  ::CircularCoaxialModeField(const IVect& ref_surf,
			     const DistributedProblem<Dimension3>& var_problem,
			     const HarmonicMaxwell_3D<Complexe>& var_maxwell,
			     const VectString& param, int& nb)
  {
    bool search_axis = false;
    if (param(nb) == "AUTO")
      {
	search_axis = true;
	nb++;
      }
    else
      {
	pt_axe_1(0)  = to_num<Real_wp>(param(nb++));
	pt_axe_1(1)  = to_num<Real_wp>(param(nb++));
	pt_axe_1(2)  = to_num<Real_wp>(param(nb++));
	pt_axe_2(0)  = to_num<Real_wp>(param(nb++));
	pt_axe_2(1)  = to_num<Real_wp>(param(nb++));
	pt_axe_2(2)  = to_num<Real_wp>(param(nb++));

	a = to_num<Real_wp>(param(nb++));
	b = to_num<Real_wp>(param(nb++));
      }
    
    volt = to_num<Complexe>(param(nb++));    
    
    Vector<bool> is_ref_on_surf(var_problem.mesh.GetNbReferences()+1);
    is_ref_on_surf.Fill(false);
    for (int i = 0; i < ref_surf.GetM(); i++)
      is_ref_on_surf(ref_surf(i)) = true;
    
    int ref = -1;
    Vector<int> other_ref;
    for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
      if (var_problem.mesh.BoundaryRef(i).GetReference() > 0)
	if (is_ref_on_surf(var_problem.mesh.BoundaryRef(i).GetReference()))
	  {
	    int num_elem = var_problem.mesh.BoundaryRef(i).numElement(0);
	    ref = var_problem.mesh.Element(num_elem).GetReference();
	    
	    if (search_axis)
	      {
		for (int j = 0; j < var_problem.mesh.Element(num_elem).GetNbBoundary(); j++)
		  {
		    int nf = var_problem.mesh.Element(num_elem).numBoundary(j);
		    int ref2 = var_problem.mesh.Boundary(nf).GetReference();
		    if ((ref2 > 0) && (!is_ref_on_surf(ref2)))
		      other_ref.PushBack(ref2);
		  }
	      }
	  }
    
    // wave number for the selected media
    kwave = var_problem.GetOmega()*var_maxwell.coefficient_impedance_absorbing(ref);
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
    if (search_axis)
      {
	// here we perform manipulations such that in each processor
	// other_ref will be the same as if the simulation was sequential
	if (rank_proc == 0)
	  {
	    // in processor 0, we retrieve all the other references found
	    MPI_Status status;
	    int nb_ref; Vector<int> all_ref(other_ref);
	    for (int p = 1; p < nb_proc; p++)
	      {
		MPI_Recv(&nb_ref, 1, MPI_INTEGER, p, 122, var_problem.comm_group_mode, &status);
		other_ref.Reallocate(nb_ref);
		MPI_Recv(other_ref.GetData(), nb_ref, MPI_INTEGER, p, 123, var_problem.comm_group_mode, &status);
		all_ref.PushBack(other_ref);
	      }
	    
	    RemoveDuplicate(all_ref);
	    other_ref = all_ref;
	    
	    // sending back the array other_ref to all processors
	    nb_ref = other_ref.GetM();
	    for (int p = 1; p < nb_proc; p++)
	      {
		MPI_Send(&nb_ref, 1, MPI_INTEGER, p, 124, var_problem.comm_group_mode);
		MPI_Send(other_ref.GetData(), nb_ref, MPI_INTEGER, p, 125, var_problem.comm_group_mode);
	      }
	  }
	else
	  {
	    // sending references of the current processor to the master
	    int nb_ref = other_ref.GetM();
	    MPI_Send(&nb_ref, 1, MPI_INTEGER, 0, 122, var_problem.comm_group_mode);
	    MPI_Send(other_ref.GetData(), nb_ref, MPI_INTEGER, 0, 123, var_problem.comm_group_mode);
	    
	    // retrieving the final references (common to everyone)
	    MPI_Status status;
	    MPI_Recv(&nb_ref, 1, MPI_INTEGER, 0, 124, var_problem.comm_group_mode, &status);
	    other_ref.Reallocate(nb_ref);
	    MPI_Recv(other_ref.GetData(), nb_ref, MPI_INTEGER, 0, 125, var_problem.comm_group_mode, &status);
	  }
      }
#else
    int rank_proc(0);
#endif    
        
    // now we retrieve the curve parameters of other_ref
    if (search_axis)
      {
	a = 1e300; b = 0;
	R3 vec_axis, pt_axis;
	for (int i = 0; i < other_ref.GetM(); i++)
	  {
	    ref = other_ref(i);
	    if (var_problem.mesh.GetCurveType(ref) == Mesh<Dimension3>::CURVE_CYLINDER)
	      {
		const VectReal_wp& param = var_problem.mesh.GetCurveParameter(ref);
		R3 ptA(param(0), param(1), param(2));
		R3 vec_u(param(3), param(4), param(5));
		Real_wp radius = param(6);
		
		if (radius < a)
		  {
		    a = radius;
		    vec_axis = vec_u;
		    pt_axis = ptA;
		  }
		
		if (radius > b)
		  b = radius;
	      }
	  }
	
	if ((a > 1e200) || (b < a))
	  {
	    cout << "Parameters of the two cylinders have not been found" << endl;
	    abort();
	  }
	
	pt_axe_1 = pt_axis;
	pt_axe_2 = pt_axis + vec_axis;	
      }
    
    if ( pt_axe_1 == pt_axe_2)
      {
	cout<<"WARNING : Points on the cylinder axis are the same"<<endl;
	abort();
      }
        
    
    if ((rank_proc == 0) && (var_problem.print_level >= 2))
      {
        cout << "You have selected a coaxial source on a cylinder" << endl;
        cout << "Wave number k = " << kwave << endl;
        cout << "Points on the cylinder " << pt_axe_1 << "  " << pt_axe_2 << endl;
        cout << "Cylinder radii = " << a << ", " << b << endl;
        cout << "Voltage = " << volt << endl;
      }
  }
  
  
  template<class Complexe>
  void CircularCoaxialModeField<Complexe>
  ::ProjOrtho (const R3& pt, R3& pt_proj) const
  {
    Real_wp prod_scal_1 = 0;
    Real_wp prod_scal_2 = 0;
    R3 A,B;
    A = pt-pt_axe_1;
    B = pt_axe_2-pt_axe_1;
    prod_scal_1 = DotProd(A, B);
    prod_scal_2 = DotProd(B, B);
    pt_proj = pt_axe_1 + (prod_scal_1/prod_scal_2)*B;
  }
  
  
  
  template<class Complexe>
  void CircularCoaxialModeField<Complexe>
  ::CalculLaplacianMode(const R3& pt, Complexe& grad_phi_x, Complexe& grad_phi_y, Complexe& grad_phi_z) const
  {
    R3 proj, vec_proj;
    ProjOrtho(pt, proj);
    vec_proj = pt-proj;
    Real_wp r = Norm2(vec_proj);
    Complexe coef = volt/(r*r*log(b/a));
    grad_phi_x = vec_proj(0)*coef;
    grad_phi_y = vec_proj(1)*coef;
    grad_phi_z = vec_proj(2)*coef;
  }

  
  template<class Complexe>
  void CircularCoaxialModeField<Complexe>
  ::EvaluateFunction(const R3& x, Vector<Complexe>& f) const
  {
    Complexe grad_phi_x, grad_phi_y, grad_phi_z;
    CalculLaplacianMode(x, grad_phi_x, grad_phi_y,grad_phi_z); 
    
    // coef = -2 i k
    Complexe coef;
    to_complex(-2.0*Iwp*kwave, coef);
    
    // Ici, on se place dans un repere local tel que z=0
    f(0) = coef*grad_phi_x;
    f(1) = coef*grad_phi_y;
    f(2) = coef*grad_phi_z;    
    //DISP(x); DISP(f);
  }


  /*********************************
   * VolumetricSource_MaxwellHdg3D *
   *********************************/
      
  template<class T>
  bool VolumetricSource_MaxwellHdg3D<T>::IsNonNullDirichletSource(int cond_ref)
  {
    // Dirichlet condition is treated weakly
    return false;
  }

  
  //! evaluates surface source
  template<class T>
  void VolumetricSource_MaxwellHdg3D<T>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
			   const SetMatrices<Dimension3>& MatricesElem, Vector<T>& f)
  {
    VolumetricSource_Base<T, Dimension3>::EvaluateSurfacicSource(k, PointsElem, MatricesElem, f);

    // only three first components are kept
    TinyVector<T, 3> f0(f(0), f(1), f(2)), vec_u;
    f.Zero();

    R3 normale = MatricesElem.GetNormaleQuadratureBoundary(k);

    if (var_problem.mesh.GetBoundaryCondition(this->ref_boundary_)
	== BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// case of an inhomogeneous Dirichlet
        TimesProd(f0, normale, vec_u);
        f(3) = vec_u(0);
        f(4) = vec_u(1);
        f(5) = vec_u(2);
      }
    else
      {
        // we subtract normal component of f0 to have a tangential vector
        T f0_dot_n = DotProd(f0, normale);
        Add(-f0_dot_n, normale, f0);
        
	// case of an inhomogeneous Neumann or Robin condition
        T m_iomega, coef_v;
        var_problem.GetMiomega(m_iomega);
        coef_v = Real_wp(1) / m_iomega;
	
        f(0) = coef_v*f0(0);
        f(1) = coef_v*f0(1);
        f(2) = coef_v*f0(2);
      }
  }
  
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_CXX
#endif
