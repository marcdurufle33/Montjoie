#ifndef MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_HXX

namespace Montjoie
{
  
  template<class T>
  class VirtualMaxwellAxiSourceFEM : public VirtualSourceFEM<T, Dimension2>
  {
  public:
    template<class TypeEquation>
    VirtualMaxwellAxiSourceFEM(const EllipticProblem<TypeEquation>&);
    
    virtual bool IsNonNull_SourceH(const VectR2& s);
    virtual void EvaluateSource_H(int i, int j, const R2& x, Vector<T>& f);
    
  };


  class PlaneWaveMaxwellAxi
  {
  protected:
    Real_wp k_bot; //!< norm of transverse component of wave vector (ie sqrt(kx^2+ky^2))
    Real_wp kwave_z; //!< kz
    Real_wp teta0; //!< angle between kx and ky
    Real_wp threshold; //!< threshold to determine if kx = ky = 0
    bool incidence_axial; //!< true if kx = ky = 0    

    Complex_wp E0_x; //!< x-component of electric field
    Complex_wp E0_y; //!< y-component of electric field
    Complex_wp E0_z; //!< z-component of electric field
    Complex_wp H0_x; //!< x-component of magnetic field
    Complex_wp H0_y; //!< y-component of magnetic field
    Complex_wp H0_z; //!< z-component of magnetic field

    Complex_wp mu0; //!< magnetic permeability
    Complex_wp omega; //!< pulsation
    Complex_wp epsilon0; //!< dielectric permittivity

    bool modified_form;
    int number_mode; //!< mode number
    
  private:
    const EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;
    
    void InitDefaultValues();
    
  public:
    template<class TypeEquation>
    PlaneWaveMaxwellAxi(const EllipticProblem<TypeEquation>& var);
    
    void Evaluate_E_H(const R2& x, R3_Complex_wp& E, R3_Complex_wp& H, int nb);    
    void Evaluate_E(const R2& x, R3_Complex_wp& f);
    void Evaluate_H(const R2& x, R3_Complex_wp& f);

    void EvaluateIncidentField(int i, int j, const R2& x, R3_Complex_wp& f);
    
  };


  class IncidentWaveProjector_MaxwellAxi
    : public IncidentWaveProjector<Complex_wp, Dimension2>, public PlaneWaveMaxwellAxi
  {
  private:
    //const EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
				     IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
  };

  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  class DiffractedWaveSource_MaxwellAxi
    : public VirtualMaxwellAxiSourceFEM<Complex_wp>, public PlaneWaveMaxwellAxi
  {
  public :
    //! coef_vol -> coefficient for volumetric integrals  \f$ \int f \varphi \f$
    TinyMatrix<Complex_wp, Symmetric, 3, 3> coef_E; 
    TinyMatrix<Complex_wp, Symmetric, 3, 3> coef_H; //!< coefficients for volumetric integral
    int ref_domain;
    Vector<TinyVector<Complex_wp, 3> > EvalF, EvalH;    
    bool pml_element;

  private:
    const EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_MaxwellAxi(const EllipticProblem<TypeEquation>& vars_helm,
				    IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void InitElement(int num_elem, const VectR2& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR2& s);

    bool IsNonNullGradientSource(const VectR2& s);    
    void EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f);
    
    void EvaluateSource_H(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNull_SourceH(const VectR2& s);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
                                const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);

    void EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension2>& PointsElem,
                                        const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    
  };


  //! source for the scattering of an incident wave  
  template<>
  class DiffractedWaveSource<HarmonicMaxwellEquation_HcurlAxi>
    : public DiffractedWaveSource_MaxwellAxi
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
  };


  //! source for the scattering of an incident wave
  /*!
    the total field is computed
  */
  class TotalWaveSource_MaxwellAxi
    : public VirtualMaxwellAxiSourceFEM<Complex_wp>, public PlaneWaveMaxwellAxi
  {
  protected:
    int ref_domain;
    Vector<TinyVector<Complex_wp, 3> > EvalF, EvalH;    
    bool pml_element;

  private:
    const EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;

  public :
    template<class TypeEquation>
    TotalWaveSource_MaxwellAxi(const EllipticProblem<TypeEquation>& vars_helm,
			       IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
    void EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f);    
    bool IsNonNullVolumetricSource(const VectR2& s);

    void EvaluateSource_H(int i, int j, const R2& x, VectComplex_wp& f);    
    bool IsNonNull_SourceH(const VectR2& s);

    bool IsNonNullGradientSource(const VectR2& s);
    void EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
                                const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);

    void EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension2>& PointsElem,
                                        const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);

    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    
  };


  //! source for the scattering of an incident wave  
  template<>
  class TotalWaveSource<HarmonicMaxwellEquation_HcurlAxi> : public TotalWaveSource_MaxwellAxi
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
  };
  
  
  //! volumetric source
  template<class T>
  class VolumetricSource_MaxwellAxi : public VirtualMaxwellAxiSourceFEM<T>
  {
  private:
    const VarProblem<Dimension2>& var_problem;
    const VarSourceProblem<T, Dimension2>& var_source;
    const VarAxisymProblem& var_axi;
    const VarBoundaryCondition_Base& var_boundary;
    
    void InitDefaultValues(bool);
    
  public :
    T coef_vol; //!< coefficient
    int number_mode;
    Globatto<Real_wp> gauss_teta; //!< quadrature formula over theta    
    Vector<TinyVector<T, 3> > EvalF;
    bool modified_form;
    VirtualSourceField<T, Dimension3>* fsrc;

    VectR2 left_corner, right_corner;
    Vector<TinyVector<T, 2> > liste_current;
    int reference_source, ref_domain;
    bool uniform_source;

    template<class TypeEquation>
    VolumetricSource_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
				const Vector<VectString>&);
    
    ~VolumetricSource_MaxwellAxi();

    void EvaluateUniformSource(int i, int j, const R2& x, TinyVector<T, 3>& f);
    
    bool IsNonNullVolumetricSource(const VectR2& s);
    void EvaluateVolumetricSource(int i, int j, const R2& x, Vector<T>& f);
   
    bool IsNonNullGradientSource(const VectR2& s);
    void EvaluateGradientSource(int i, int j, const R2& x, Vector<T>& f);

    void InitElement(int i, const VectR2& s);
    
  };
  
  
  //! volumetric source
  template<>
  class VolumetricSource<HarmonicMaxwellEquationAxiDG>
    : public VolumetricSource_MaxwellAxi<Complex_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>&);
    
  };

  
  //! volumetric source
  template<>
  class VolumetricSource<HarmonicMaxwellEquation_HcurlAxi>
    : public VolumetricSource_MaxwellAxi<Complex_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>&);
    
  };


  //! volumetric source
  template<class T>
  class DiracSource_MaxwellAxi : public VirtualMaxwellAxiSourceFEM<T>
  {
  private:
    
  public:
    template<class TypeEquation>
    inline DiracSource_MaxwellAxi(const EllipticProblem<TypeEquation>& var)
      : VirtualMaxwellAxiSourceFEM<T>(var) {}

    inline bool IsDiracSource() const { return true; }
    
  };


  template<>
  class DiracSource<HarmonicMaxwellEquation_HcurlAxi>
    : public DiracSource_MaxwellAxi<Complex_wp>
  {
  public :
    template<class TypeEquation>
    inline DiracSource(const EllipticProblem<TypeEquation>& var,
                       const Vector<VectString>&) : DiracSource_MaxwellAxi<Complex_wp>(var)
    {}
    
    
  };
  
} 

#define MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_HXX
#endif
