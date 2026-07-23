#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_2D_CXX

namespace Montjoie
{
  //! constructor with a given problem
  /*!
    The given problem is used to get the wave vector of the place wave
   */
  template<class TypeEquation>
  DiffractedWaveSource_Maxwell2D
  ::DiffractedWaveSource_Maxwell2D(const EllipticProblem<TypeEquation>& vars,
				   IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension2>(vars), incident_wave(u_inc),
      var_problem(vars), var_maxwell(vars)
  {
    coef_vol.SetIdentity(); epsilon0 = vars.epsilon0;
    epsilon.SetIdentity();
    coef_curl = 1.0; coef_neu = 1.0;
    mu0 = vars.mu0;
  }
  
  
  //! initialization before evaluating volumetric source
  void DiffractedWaveSource_Maxwell2D::InitElement(int num_elem, const VectR2& s)
  {
    Complex_wp mu; int ref = var_problem.mesh.Element(num_elem).GetReference();
    epsilon = var_maxwell.ref_epsilon(ref).GetConstant();
    mu = var_maxwell.ref_mu(ref);
    coef_vol(0, 0) = var_problem.GetSquareOmega()*(epsilon(0, 0)-epsilon0);
    coef_vol(1, 1) = var_problem.GetSquareOmega()*(epsilon(1, 1)-epsilon0);
    coef_vol(0, 1) = var_problem.GetSquareOmega()*epsilon(1, 0);
    
    coef_curl = -1.0/mu + 1.0/mu0;
    coef_curl *= -mu0;
  }
  
  
  //! evaluation of incident wave
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x point where the incident wave has to be evaluated
    \param[out] f incident electric field
    E^inc = rot ( exp( i k x))
   */
  void DiffractedWaveSource_Maxwell2D
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    TinyVector<Complex_wp, 2> grad_uinc; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    f(0) = Iwp*grad_uinc(1)/var_problem.GetOmega();
    f(1) = -Iwp*grad_uinc(0)/var_problem.GetOmega();
    
    // for Dirichlet condition -E_inc
    if (this->num_loc_ >= 0)
      {
	f(0) = -f(0);
	f(1) = -f(1);
      }
  }
  
  
  //! volumetric source if dielectric material, epsilon <> eps0
  /*!
    \param[in] s vertices
   */
  bool DiffractedWaveSource_Maxwell2D
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    if (coef_vol.IsZero())
      return false;
    
    return true;
  }
  
  
  //! evaluation of volumetric source
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] f source vector 
   */
  void DiffractedWaveSource_Maxwell2D
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    TinyVector<Complex_wp, 2> Einc, grad_inc; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_inc);
    Einc(0) = Iwp*grad_inc(1)/var_problem.GetOmega();
    Einc(1) = -Iwp*grad_inc(0)/var_problem.GetOmega();
    f(0) = Einc(0)*coef_vol(0,0) + Einc(1)*coef_vol(0,1);
    f(1) = Einc(0)*coef_vol(0,1) + Einc(1)*coef_vol(1,1);
  }
  
  
  //! volumetric source if dielectric material, mu <> mu0
  /*!
    \param[in] s vertices
   */
  bool DiffractedWaveSource_Maxwell2D::IsNonNullGradientSource(const VectR2& s)
  {
    if (coef_curl==Real_wp(0))
      return false;
    
    return true;
  }
  
  
  //! evaluation of volumetric source
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x quadrature point
    \param[out] f source vector 
   */
  void DiffractedWaveSource_Maxwell2D
  ::EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    incident_wave.EvaluateFunction(x, f(0));
    f(0) *= Iwp*var_problem.GetOmega()*coef_curl;
  }
  
  
  //! surfacic source for Neumann condition
  bool DiffractedWaveSource_Maxwell2D
  ::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if (cond_ref== BoundaryConditionEnum::LINE_NEUMANN)
      return true;
    
    return false;
  }
  
  //! surfacic source 
  void DiffractedWaveSource_Maxwell2D
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp uinc;
    incident_wave.EvaluateFunction(x, uinc);
    uinc *= coef_neu*Iwp*var_problem.GetOmega();
    f(0) = uinc*normale(1);
    f(1) = -uinc*normale(0);
  }


  template<class TypeEquation>
  DiffractedWaveSource<HarmonicMaxwellEquation_2D>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : DiffractedWaveSource_Maxwell2D(var, u_inc)
  {
  }
  
  
  //! constructor with given problem
  template<class TypeEquation>
  TotalWaveSource_Maxwell2D
  ::TotalWaveSource_Maxwell2D(const EllipticProblem<TypeEquation>& vars,
			      IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension2>(vars), incident_wave(u_inc),
      var_problem(vars)
  {
  }
  
  
  //! surfacic source for edges with an absorbing boundary condition
  bool TotalWaveSource_Maxwell2D::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }
  
  
  //! surfacic source
  void TotalWaveSource_Maxwell2D
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp uinc;
    TinyVector<Complex_wp, 2> Einc, grad_uinc;
    incident_wave.EvaluateFunctionGradient(x, uinc, grad_uinc);
    Einc(0) = Iwp*grad_uinc(1)/var_problem.GetOmega();
    Einc(1) = -Iwp*grad_uinc(0)/var_problem.GetOmega();
    Complex_wp rotEinc = uinc*Iwp*var_problem.GetOmega();
    
    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp Einc_n = DotProd(Einc, normale);
    Einc(0) -= Einc_n*normale(0);
    Einc(1) -= Einc_n*normale(1);
    
    f(0) = -rotEinc*normale(1) - Iwp*var_problem.GetOmega()*Einc(0);
    f(1) = rotEinc*normale(0) - Iwp*var_problem.GetOmega()*Einc(1);
  }
  

  template<class TypeEquation>
  TotalWaveSource<HarmonicMaxwellEquation_2D>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& var,
		    IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : TotalWaveSource_Maxwell2D(var, u_inc)
  {
  }


  template<class TypeEquation>
  IncidentWaveProjector_Maxwell2D
  ::IncidentWaveProjector_Maxwell2D(const EllipticProblem<TypeEquation>& var,
				    IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension2>(var, u_inc), var_problem(var)
  {
  }
  
  
  void IncidentWaveProjector_Maxwell2D
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    TinyVector<Complex_wp, 2> grad_uinc; Complex_wp u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    f(0) = Iwp*grad_uinc(1)/var_problem.GetOmega();
    f(1) = -Iwp*grad_uinc(0)/var_problem.GetOmega();
  }
  
}

#define MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_2D_CXX
#endif
