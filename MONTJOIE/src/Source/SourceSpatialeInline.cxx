#ifndef MONTJOIE_FILE_SOURCE_SPATIALE_INLINE_CXX

namespace Montjoie
{

  template<class T, class Dimension>
  inline VirtualSourceField<T, Dimension>::~VirtualSourceField() 
  {
  }


  template<class T, class Dimension>
  inline IncidentWaveField<T, Dimension>::IncidentWaveField() 
  {
    time_pulse = NULL;
  }

  
  template<class T, class Dimension>
  inline IncidentWaveField<T, Dimension>::~IncidentWaveField() 
  {
  }


  template<class T, class Dimension>
  inline VirtualTimeSource<T>* IncidentWaveField<T, Dimension>::GetTimeSource()
  {
    return time_pulse;
  }


  template<class T, class Dimension>
  inline void IncidentWaveField<T, Dimension>::SetTimeSource(VirtualTimeSource<T>* fsrc)
  {
    time_pulse = fsrc;
  }


  template<class T, class Dimension>
  inline void IncidentWaveField<T, Dimension>::InitElement(int i, const typename Dimension::VectR_N& s)
  {
  }
  

  template<class T, class Dimension>
  inline void IncidentWaveField<T, Dimension>::UpdateCoefAB(T& a, T& b) const
  {
  }
  

  template<class Dimension>
  inline PlaneWaveIncidentField<Real_wp, Dimension>::PlaneWaveIncidentField(const R_N&, const R_N&) 
  {
  }

  
  template<class Dimension>
  inline void PlaneWaveIncidentField<Real_wp, Dimension>
  ::EvaluateFunction(const R_N& x, Real_wp& f) const
  {
    f = 0;
    abort();
  }

  
  template<class Dimension>
  inline void PlaneWaveIncidentField<Real_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& x, Real_wp& f,
			     TinyVector<Real_wp, Dimension::dim_N>& df) const
  {
    f = 0;
    df.Fill(0);
    abort();
  }


  template<class Dimension>
  inline PlaneWaveIncidentFieldComplex<Real_wp, Dimension>
  ::PlaneWaveIncidentFieldComplex(const R_N&, const TinyVector<Complex_wp, Dimension::dim_N>&) 
  {
  }

  
  template<class Dimension>
  inline void PlaneWaveIncidentFieldComplex<Real_wp, Dimension>
  ::EvaluateFunction(const R_N& x, Real_wp& f) const
  {
    f = 0;
    abort();
  }

  
  template<class Dimension>
  inline void PlaneWaveIncidentFieldComplex<Real_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& x, Real_wp& f,
			     TinyVector<Real_wp, Dimension::dim_N>& df) const
  {
    f = 0;
    df.Fill(0);
    abort();
  }


  template<class Dimension>
  inline HankelIncidentField<Real_wp, Dimension>::HankelIncidentField(const R_N&, const R_N&) 
  {
  }

  
  template<class Dimension>
  inline void HankelIncidentField<Real_wp, Dimension>
  ::EvaluateFunction(const R_N& x, Real_wp& f) const
  {
    f = 0;
    abort();
  }

  
  template<class Dimension>
  inline void HankelIncidentField<Real_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& x, Real_wp& f,
			     TinyVector<Real_wp, Dimension::dim_N>& df) const
  {
    f = 0;
    df.Fill(0);
    abort();
  }


  template<class Dimension>
  inline GaussianBeamIncidentField<Real_wp, Dimension>
  ::GaussianBeamIncidentField(const R_N&, const R_N&, const Real_wp&) 
  {
  }

  
  template<class Dimension>
  inline void GaussianBeamIncidentField<Real_wp, Dimension>
  ::EvaluateFunction(const R_N& x, Real_wp& f) const
  {
    f = 0;
    abort();
  }

  
  template<class Dimension>
  inline void GaussianBeamIncidentField<Real_wp, Dimension>
  ::EvaluateFunctionGradient(const R_N& x, Real_wp& f,
			     TinyVector<Real_wp, Dimension::dim_N>& df) const
  {
    f = 0;
    df.Fill(0);
    abort();
  }

  
  template<class T, class Dimension>
  inline IncidentWaveSourceField<T, Dimension>
  ::IncidentWaveSourceField(IncidentWaveField<T, Dimension>* u_inc)
  {
    incident_wave = u_inc;
  }


  template<class T, class Dimension>
  inline IncidentWaveSourceField<T, Dimension>::~IncidentWaveSourceField()
  {
    delete incident_wave;
  }
  
  
  template<class T, class Dimension>
  inline void IncidentWaveSourceField<T, Dimension>::SetPolarization(const Vector<T>& P)
  {
    polar = P;
  }
  

  template<class T, class Dimension>
  inline IncidentWaveProjector<T, Dimension>::~IncidentWaveProjector() 
  {
  }


  template<class T, class Dimension>
  inline const IncidentWaveField<T, Dimension>& IncidentWaveProjector<T, Dimension>
  ::GetIncidentWaveField() const 
  {
    return incident_wave;
  }

  
  template<class T, class Dimension>
  inline void IncidentWaveProjector<T, Dimension>::SetTime(const Real_wp& tn, int n) 
  {
    t = tn;
  }

  
  template<class T, class Dimension>
  inline void IncidentWaveProjector<T, Dimension>
  ::EvaluateFunctionGradient(int i, int j, const R_N& x, Vector<T>& f, Vector<T>& df)
  { 
    abort(); 
  }


  //! constructor with EllipticProblem
  template<class T, class Dimension> template<class TypeEquation>
  inline VirtualProjectorFEM<T, Dimension>
  ::VirtualProjectorFEM(const EllipticProblem<TypeEquation>& var)
    : var_boundary(var), var_problem(var)
  {
    nb_unknowns = TypeEquation::nb_unknowns;
    var.InitFftComputation(this->fft_interface);
  }

  
  //! destructor
  template<class T, class Dimension>
  inline VirtualProjectorFEM<T, Dimension>::~VirtualProjectorFEM()
  {
  }
  
  
  //! constructor with a given problem and an incident field
  template<class T, class Dimension> template<class TypeEquation>
  inline IncidentWaveProjector<T, Dimension>
  ::IncidentWaveProjector(const EllipticProblem<TypeEquation>& var,
			  IncidentWaveField<T, Dimension>& fsrc_)
    : VirtualProjectorFEM<T, Dimension>(var), incident_wave(fsrc_)
  {
    var.GetPolarization(polarization);
  }
  
  
  //! sets the polarization of the incident wave
  template<class T, class Dimension>
  inline void IncidentWaveProjector<T, Dimension>::SetPolarization(const Vector<T>& p)
  {
    polarization = p;
  }
  
  
  //! constructor to overload in derived classes
  template<class T, class Dimension> template<class TypeEquation>
  inline VirtualSourceFEM<T, Dimension>::VirtualSourceFEM(const EllipticProblem<TypeEquation>& var)
    : VirtualProjectorFEM<T, Dimension>(var), var_problem(var), var_boundary(var)
  {
    this->grad_dirac = false;
    num_bound_ref_ = -1; num_boundary_ = -1; num_loc_ = -1;
    num_elem_ = -1; ref_boundary_ = -1;
    var.GetPolarization(polarization);
    var.GetPolarizationGrad(polarization_grad);
  }

  
  //! default destructor
  template<class T, class Dimension>
  inline VirtualSourceFEM<T, Dimension>::~VirtualSourceFEM()
  {
  }


  //! returns polarization associated with the source (for Dirac)
  template<class T, class Dimension>
  inline const Vector<T>& VirtualSourceFEM<T, Dimension>::GetPolarization() const
  {
    return polarization;
  }


  //! returns polarization associated with the gradient source (for Dirac)
  template<class T, class Dimension>
  inline const Vector<T>& VirtualSourceFEM<T, Dimension>::GetPolarizationGrad() const
  {
    return polarization_grad;
  }


  //! returns true if it is a dirac source with gradients
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>::IsGradientDirac() const
  {
    return grad_dirac;
  }

  
  //! returns origin associated with the source
  template<class T, class Dimension>
  inline const typename Dimension::R_N& VirtualSourceFEM<T, Dimension>::GetOrigin() const
  {
    return var_problem.GetPhaseOrigin();
  }

  
  //! sets polarization associated with the source
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>::SetPolarization(const Vector<T>& polar)
  {
    polarization = polar;
  }


  //! sets polarization associated with the source
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>::SetPolarizationGrad(const Vector<T>& polar)
  {
    polarization_grad = polar;
  }

 
  //! initialization of object for each volumic element of integration
  /*!
    \param[in] num_elem element number
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
    num_elem_ = num_elem;
  }
  
  
  // dirichlet condition
  //! used for the evaluation of incident wave and for inhomogeneous Dirichlet condition
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x point where the function needs to be evaluated
    \param[out] f evaluation of the function
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f)
  {
  }
  
  
  // volumetric integral \int f \varphi
  //! used for volumic integrals \f$ \int f \varphi \f$
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x point where the function f needs to be evaluated
    \param[out] f evaluation of function f
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f) 
  {
  }
  
  
  //! returns true if the source is non-null
  /*!
    \param[in] s list of vertices of the element
  */
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return false;
  }
  
  
  // volumetric integral \int f grad(\varphi)
  //! used for volumetric integrals \f$ \int f \nabla \varphi \f$
  /*!
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] x where the function f needs to be evaluated
    \param[out] f evaluation of function f
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f)
  {
  }
  
  
  //! returns true if the source is non-null
  /*!
    \param[in] s list of vertices of the element
  */
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
    
  
  //! surfacic integral \f$ \int_\Gamma f \varphi \f$
  /*!
    \param[in] k quadrature point number on the edge/face
    \param[in] PointsElem transformation of references points
    \param[in] MatricesElem jacobian matrices
    \param[out] f evaluation of f on quadrature points
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
			   const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
  {
  }
  
  
  //! returns true if a surfacic integral needs to be computed
  /*!
    \param[in] cond_ref boundary condition
  */
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>
  ::IsNonNullSurfacicSource(int ref)
  {
    return false;
  }
  
  
  //! surfacic integral \f$ \int_\Gamma f \cdot \nabla_\Gamma \varphi \f$
  /*!
    \param[in] k quadrature point number on the edge/face
    \param[in] PointsElem transformation of references points
    \param[in] MatricesElem jacobian matrices
    \param[out] f evaluation of f on quadrature points
  */
  template<class T, class Dimension>
  inline void VirtualSourceFEM<T, Dimension>
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
				   const SetMatrices<Dimension>& MatricesElem, Vector<T>& f)
  {
  }
  
  
  //! returns true if a surfacic gradient integral needs to be computed
  /*!
    \param[in] cond_ref boundary condition
  */
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    return false;
  }
  
  
  //! return true if there is a dirac delta to add to the source
  template<class T, class Dimension>
  inline bool VirtualSourceFEM<T, Dimension>::IsDiracSource() const
  {
    return false;
  }


  template<class T, class Dimension>
  inline T VirtualSourceFEM<T, Dimension>::GetCoefficientVolume() const
  {
    T one; SetComplexOne(one);
    return one;
  }
  
  
  /*************************
   * VolumetricSource_Base *
   *************************/


  //! constructor with given problem
  template<class T, class Dimension> template<class TypeEquation>
  inline VolumetricSource_Base<T, Dimension>
  ::VolumetricSource_Base(const EllipticProblem<TypeEquation>& var,
			  const Vector<VectString>& param)  
    : VirtualSourceFEM<T, Dimension>(var), var_source(var), var_problem(var)
  {
    InitDefaultValues();
    
    var.ConstructVolumeSourceFunctions(param, *this);
    coef_vol = var.GetCoefficientVolumeSource();
    var.ConstructSurfaceSourceFunctions(param, *this);
  }
  

  //! constructor with EllipticProblem
  template<class TypeEquation>
  inline VolumetricSource<TypeEquation>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_Base<typename TypeEquation::Complexe,
			    typename TypeEquation::Dimension>(var, param)
  {
  }


  template<class T, class Dimension>
  inline T VolumetricSource_Base<T, Dimension>::GetCoefficientVolume() const
  {
    return coef_vol;
  }
  
  
  /****************
   * Dirac Source *
   ****************/
  
  
  //! constructor with given problem
  template<class TypeEquation>
  inline DiracSource<TypeEquation>::
  DiracSource(const EllipticProblem<TypeEquation>& var,
	      const Vector<VectString>& param)
    : VirtualSourceFEM<typename TypeEquation::Complexe,
		       typename TypeEquation::Dimension>(var)
  {
    coef_vol = var.GetCoefficientVolumeSource();
    if (param(0).GetM() > 1)
      for (int k = 0; k < Dimension::dim_N; k++)
	origin(k) = to_num<Real_wp>(param(0)(k+1));
    else
      for (int k = 0; k < Dimension::dim_N; k++)
	origin(k) = var.GetPhaseOrigin()(k);

    for (int k = 0; k < param(0).GetM(); k++)
      if (param(0)(k) == "Polarization")
	{
	  this->polarization.Zero();
	  int pmax = min(this->polarization.GetM()+1, param(0).GetM()-k);
	  for (int p = 1; p < pmax; p++)
	    this->polarization(p-1) = to_num<Complexe>(param(0)(k+p));
	  
	  break;
	}
    
    for (int k = 0; k < param(0).GetM(); k++)
      if (param(0)(k) == "PolarizationGrad")
	{
          this->grad_dirac = true;
	  this->polarization_grad.Zero();
	  int pmax = min(this->polarization_grad.GetM()+1, param(0).GetM()-k);
	  for (int p = 1; p < pmax; p++)
	    this->polarization_grad(p-1) = to_num<Complexe>(param(0)(k+p));
	  
	  break;
	}
  }
  
  
  //! dirac source only
  template<class TypeEquation>
  inline bool DiracSource<TypeEquation>::IsDiracSource() const
  {
    return true;
  }


  //! returns the origin point of the Dirac
  template<class TypeEquation>
  inline const typename TypeEquation::Dimension::R_N& DiracSource<TypeEquation>::GetOrigin() const
  {
    return origin;
  }

  
  template<class TypeEquation>
  inline typename TypeEquation::Complexe DiracSource<TypeEquation>::GetCoefficientVolume() const
  {
    return coef_vol;
  }
  

}

#define MONTJOIE_FILE_SOURCE_SPATIALE_INLINE_CXX
#endif
