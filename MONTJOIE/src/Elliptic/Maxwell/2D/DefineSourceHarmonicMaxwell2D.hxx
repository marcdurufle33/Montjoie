#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_2D_HXX

namespace Montjoie
{
    
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  class DiffractedWaveSource_Maxwell2D : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  protected :
    Complex_wp coef_curl; //!< coefficient for curl integral
    Complex_wp coef_neu; //!< coefficient for Neumann condition
    //! dielectric permittivity in empty space
    Complex_wp epsilon0;
    //! magnetic permeability in empty space
    Complex_wp mu0;
    //! relative magnetic permeability
    Complex_wp mu;
    //! relative dielectric permittivity
    Matrix2_2sym_Complex_wp epsilon;
    //! coefficient for volumetric integral
    Matrix2_2sym_Complex_wp coef_vol;
    
  private:
    IncidentWaveField<Complex_wp, Dimension2>& incident_wave;
    const DistributedProblem<Dimension2>& var_problem;
    const HarmonicMaxwell_2D<Complex_wp>& var_maxwell;
    //const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;
    
  public:
    template<class TypeEquation>
    DiffractedWaveSource_Maxwell2D(const EllipticProblem<TypeEquation>& vars_helm,
				   IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void InitElement(int num_elem, const VectR2& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR2& s);
    
    // volumic integral \int f grad(\phi)
    void EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR2& s);
    
    // surfacic integral \int_\Gamma f \phi
    bool IsNonNullSurfacicSource(int cond_ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
				const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
  };


  template<>
  class DiffractedWaveSource<HarmonicMaxwellEquation_2D>
    : public DiffractedWaveSource_Maxwell2D
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension2>& fsrc);
    
  };

  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
  */
  class TotalWaveSource_Maxwell2D : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  protected :
    //! coefficient
    Complex_wp coef;

  private:
    IncidentWaveField<Complex_wp, Dimension2>& incident_wave;
    const DistributedProblem<Dimension2>& var_problem;
    //const HarmonicMaxwell_2D<Complex_wp>& var_maxwell;
    //const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;
    
  public :
    template<class TypeEquation>
    TotalWaveSource_Maxwell2D(const EllipticProblem<TypeEquation>& vars_helm,
			      IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
                                const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  

  //! scattering by an incident wave
  template<>
  class TotalWaveSource<HarmonicMaxwellEquation_2D>
    : public TotalWaveSource_Maxwell2D
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension2>& u_inc);
  };


  //! projector for incident wave and DG
  class IncidentWaveProjector_Maxwell2D : public IncidentWaveProjector<Complex_wp, Dimension2>
  {
  private:
    const DistributedProblem<Dimension2>& var_problem;
    
  public :
    template<class TypeEquation>
    IncidentWaveProjector_Maxwell2D(const EllipticProblem<TypeEquation>& var,
				    IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
  };
  
} // end namespace

#define MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_2D_HXX 
#endif
