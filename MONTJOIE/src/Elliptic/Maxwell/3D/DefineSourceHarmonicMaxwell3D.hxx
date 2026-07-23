#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_HXX

namespace Montjoie
{
  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed.
    For Maxwell's equations, the diffracted field solves :
    - omega^2 \varepsilon E + curl (1/mu curl E)
             = omega^2 \varepsilon E^inc - curl(1/mu curl E^inc)   in \Omega
   
    E x n = -E^{inc} x n   on Dirichlet boundaries

    n/mu x curl E = - n/mu x curl E^{inc}   on Neumann boundaries

    n/mu x curl E - Z (n x E) x n = n/mu x curl E^{inc} - Z (n x E^{inc}) x n 
          on Impedance boundaries
    
    There are no source terms on other boundaries, E^{inc} is the incident field
    \varepsilon may be complex (and include the damping sigma)
    In order to have a null source when epsilon = epsilon0, mu = mu0, we subtract
    the equation satisfied by the incident field :

    -\omega^2 \varpepsilon_0 E + curl(1/mu_0 curl E) = 0
    
    For H(curl) usual formulation, by performing integration by parts,
    and keeping only the source terms, we obtain :

    \omega^2 \int_\Omega (\varepsilon - \varepsilon_0) E^{inc} \cdot \varphi \, dx
    - \int_\Omega (1/mu - 1/mu_{0}) curl E^{inc} \cdot curl \varphi \, dx
    + \int_{Gamma_Neumann} n/mu_0 x curl E^{inc} \cdot \varphi \, dx
    + \int_{Gamma_Impedance} ( n/mu_0 x curl E^{inc} -Z (n x E^{inc}) x n ) \cdot \varphi \, dx
    
    And an inhomogeneous Dirichlet condition E x n = -E^{inc} x n
  */
  class DiffractedWaveSource_Maxwell3D : public VirtualSourceFEM<Complex_wp, Dimension3>
  {
  protected :
    // coefficients for physical properties
    Matrix3_3sym_Complex_wp mu; //!< magnetic permeability
    Matrix3_3sym_Complex_wp epsilon; //!< dielectric permittivity
    Matrix3_3sym_Complex_wp coef_vol; //!< coefficient for volumetric integral
    Matrix3_3sym_Complex_wp coef_curl; //!< coefficient for integral with curl
    Matrix3_3sym_Complex_wp invMu; //!< inverse of mu
    //! polarization of electric field
    R3_Complex_wp E_0;
    R3_Complex_wp  H_0; //!< polarization of magnetic field
    //! boundary condtion
    int reference_condition;
    //! varying physical properties ?
    bool varying_eps;
    bool varying_mu; //!< varying mu ?

  private:
    IncidentWaveField<Complex_wp, Dimension3>& incident_wave;
    const DistributedProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    //const VarBoundaryCondition_Dim<Complex_wp, Dimension3>& var_boundary;
    VirtualImpedanceFunction<Complex_wp, Dimension3>& fct_imped;
     
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_Maxwell3D(const EllipticProblem<TypeEquation>& vars_helm,
				   IncidentWaveField<Complex_wp, Dimension3>& u_inc);
    
    void InitElement(int num_elem, const VectR3& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR3& s);
    
    // volumic integral \int f grad(\phi)
    void EvaluateGradientSource(int i, int j, const R3& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR3& s);
    
    // initialization for surfacic integration
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  template<>
  class DiffractedWaveSource<HarmonicMaxwellEquation_3D>
    : public DiffractedWaveSource_Maxwell3D
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension3>& fsrc);
    
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  class TotalWaveSource_Maxwell3D
    : public VirtualSourceFEM<Complex_wp, Dimension3>
  {
  protected :
    //! polarization of electric field
    R3_Complex_wp E_0;
    R3_Complex_wp H_0; //!< polarization of magnetic field
    Matrix3_3sym_Complex_wp mu; //!< magnetic permeability
    Matrix3_3sym_Complex_wp epsilon, sigma; //!< dielectric permittivity
    Matrix3_3sym_Complex_wp invMu; //!< inverse of magnetic permeability
    //! element belonging to a PML layer ?
    bool element_insidePML;
    //! boundary condition
    int reference_condition;
    Real_wp kwave;

  private:
    IncidentWaveField<Complex_wp, Dimension3>& incident_wave;
    const DistributedProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension3>& var_boundary;
    //VirtualImpedanceFunction<Complex_wp, Dimension3>& fct_imped_abc;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    TotalWaveSource_Maxwell3D(const EllipticProblem<TypeEquation>& vars_helm,
			      IncidentWaveField<Complex_wp, Dimension3>& u_inc);
    
    void EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f);
    
    // volumic source in PML
    void InitElement(int num_elem, const VectR3& s);
    bool IsNonNullVolumetricSource(const VectR3& s);
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f);

    bool IsNonNullGradientSource(const VectR3& s);
    void EvaluateGradientSource(int i, int j, const R3& x, VectComplex_wp& f);
    
    // surfacic source on absorbing boundary condition
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				 const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };

  //! curl of a gaussian source around a point
  template<class Complexe>
  class MaxwellGaussianSource : public GaussianSourceField<Complexe, Dimension3>
  {
  protected :
    bool presence_curl;    
    
  public :
    MaxwellGaussianSource();

    void Init(const R3&, const Real_wp& r, const Real_wp& rmax, bool curl, Vector<Complexe>&);
    
    void EvaluateFunction(const R3& x, Vector<Complexe>& f) const;
    
  };
  
  
  //! scattering by an incident wave
  template<>
  class TotalWaveSource<HarmonicMaxwellEquation_3D>
    : public TotalWaveSource_Maxwell3D
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension3>& u_inc);
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  template<>
  class DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>
    : public VirtualSourceFEM<Complex_wp, Dimension3>
  {
  protected :
    // physical properties
    Matrix3_3sym_Complex_wp mu; //!< magnetic permeability
    Matrix3_3sym_Complex_wp epsilon; //!< dielectric permittivity
    Matrix3_3sym_Complex_wp coef_vol; //!< coefficient for volumetric integral
    Matrix3_3sym_Complex_wp coef_curl; //!< coefficient for integral in curl
    Matrix3_3sym_Complex_wp invMu; //!< inverse of mu
    //! polarization of electric field
    R3_Complex_wp E_0;
    R3_Complex_wp H_0;     //!< polarization of magnetic field
    //! boundary condition
    int reference_condition;
    //! varying epsilon  ?
    bool varying_eps;
    bool varying_mu; //!< varying my ?

  private:
    IncidentWaveField<Complex_wp, Dimension3>& incident_wave;
    const DistributedProblem<Dimension3>& var_problem;
    const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    //const VarBoundaryCondition_Dim<Complex_wp, Dimension3>& var_boundary;
    //ImpedanceFunction_Base<Complex_wp, Dimension3, 1>& fct_imped;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension3>& u_inc);
    
    void InitElement(int num_elem, const VectR3& s);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR3& s);
    
    // initialization for surfacic integration
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				 const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f);

    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  template<>
  class TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
    : public VirtualSourceFEM<Complex_wp, Dimension3>
  {
  protected :
    //! polarization of electric field
    R3_Complex_wp E_0;
    R3_Complex_wp H_0; //!< polarization of magnetic field

    Matrix3_3sym_Complex_wp mu; //!< magnetic permeability
    Matrix3_3sym_Complex_wp epsilon; //!< dielectric permittivity
    Matrix3_3sym_Complex_wp invMu; //!< inverse of mu
    //! element belonging to a PML layer ?
    bool element_insidePML;
    //! boundary condition
    int reference_condition;

  private:
    IncidentWaveField<Complex_wp, Dimension3>& incident_wave;
    const DistributedProblem<Dimension3>& var_problem;
    //const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    //const VarBoundaryCondition_Dim<Complex_wp, Dimension3>& var_boundary;
    //ImpedanceFunction_Base<Complex_wp, Dimension3, 1>& fct_imped_abc;

    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm, 
		    IncidentWaveField<Complex_wp, Dimension3>& u_inc);
    
    void EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f);
    
    // volumic source in PML
    void InitElement(int num_elem, const VectR3& s);

    bool IsNonNullVolumetricSource(const VectR3& s);
    void EvaluateVolumetricSource(int i, int j, const R3& x, VectComplex_wp& f);
    
    // surfacic source on absorbing boundary condition
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				 const SetMatrices<Dimension3>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };


  //! projector for incident wave and DG
  class IncidentWaveProjector_Maxwell3D_DG : public IncidentWaveProjector<Complex_wp, Dimension3>
  {
  private:
    //const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    const DistributedProblem<Dimension3>& var_problem;
    
    R3_Complex_wp E_0, H_0;
    
  public :
    template<class TypeEquation>
    IncidentWaveProjector_Maxwell3D_DG(const EllipticProblem<TypeEquation>& var,
				       IncidentWaveField<Complex_wp, Dimension3>& u_inc);
    
    void EvaluateFunction(int i, int j, const R3& x, VectComplex_wp& f);
    
  };
  

  //! class to specify a mode in a coaxial cable (circular section)
  template<class Complexe>
  class CircularCoaxialModeField : public VirtualSourceField<Complexe, Dimension3>
  {
  private:
    R3 pt_axe_1, pt_axe_2;
    Real_wp a, b;
    Complexe kwave;
    Complexe volt;
    
  public:
    CircularCoaxialModeField(const IVect& ref_surf,
			     const DistributedProblem<Dimension3>& mesh,
			     const HarmonicMaxwell_3D<Complexe>& var_maxwell,
			     const VectString& param, int& nb);
    
    void ProjOrtho (const R3& pt, R3& pt_proj) const;
    void CalculLaplacianMode(const R3& pt, Complexe& grad_phi_x, Complexe& grad_phi_y, Complexe& grad_phi_z) const;
    
    void EvaluateFunction(const R3& x, Vector<Complexe>& f) const;
    
  };


  //! volumetric source for DG
  template<class T>
  class VolumetricSource_MaxwellHdg3D : public VolumetricSource_Base<T, Dimension3>
  {
  private :
    const HarmonicMaxwell_3D<T>& var_maxwell;
    const DistributedProblem<Dimension3>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
  public:
    template<class TypeEquation>
    VolumetricSource_MaxwellHdg3D(const EllipticProblem<TypeEquation>& var,
                                  const Vector<VectString>& param);
    
    bool IsNonNullDirichletSource(int cond_ref);        
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
				const SetMatrices<Dimension3>& MatricesElem, Vector<T>& f);
        
  };

  
  template<>
  class VolumetricSource<HarmonicMaxwellEquationHdg_3D>
    : public VolumetricSource_MaxwellHdg3D<Complex_wp>
  {
  public :    
    VolumetricSource(const EllipticProblem<HarmonicMaxwellEquationHdg_3D>&,
		     const Vector<VectString>& param);
    
  };
  
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_HXX
#endif
