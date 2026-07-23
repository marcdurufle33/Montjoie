#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_HXX

namespace Montjoie
{
  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  template<class Dimension>
  class DiffractedWaveSource_Helm : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::R_N_Complex_wp R_N_Complex_wp; //!< "complex" point
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!< 2x2 or 3x3 complex symmetric matrix
    
    //! coef_vol -> coefficient for volumic integrals  \f$ \int f \varphi \f$
    Complex_wp coef_vol;
    // physical indices : rho, mu, rho0, mu0, sigma_rho, sigma_mu
    Complex_wp rho0; //!< rho in empty media
    Complex_wp rho; //!< rho in element
    Complex_wp mu0; //!< mu in empty media
    Complex_wp omega2; //!< omega^2
    MatrixN_Nsym_Complex_wp mu; //!< mu in element
    
    //! coefficient for volumic integrals coef_grad \f$\int f \cdot grad(phi) \f$
    MatrixN_Nsym_Complex_wp coef_grad;
    
    // in order to treat high conductivity boundary conditions
    int reference_condition; //!< boundary condition
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    int order_condition; //!< order for high conductivity condition
    Complex_wp epsilon; //!< thickness of high conducting object
    Complex_wp alpha; //!< parameter = sqrt(i) 
#endif
    Complex_wp coef_mu0; //!< coefficient
    Complex_wp coef_mu; //!< coefficient
    bool varying_rho; //!< variable rho ?
    bool varying_mu; //!< variable mu ?
    
  private:
    IncidentWaveField<Complex_wp, Dimension>& incident_wave;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension>& var_boundary;
    ImpedanceFunction_Base<Complex_wp, Dimension>& fct_imped;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_Helm(const EllipticProblem<TypeEquation>& vars_helm,
			      IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void InitElement(int num_elem, const VectR_N& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    // volumic integral \int f grad(\phi)
    void EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR_N& s);
    
    // initialization for surfacic integration
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    bool IsNonNullSurfacicSource(int cond_ref);
    
    // surfacic integral \int_\Gamma f \cdot \nabla_\Gamma \phi
    void EvaluateSurfacicSourceGradient(int k,const SetPoints<Dimension>& PointsElem,
					 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    
  };

  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class DiffractedWaveSource<HelmholtzEquation<Dimension> >
    : public DiffractedWaveSource_Helm<Dimension>
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension>& f_inc);
    
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  template<class Dimension>
  class TotalWaveSource_Helm : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!< 2x2 or 3x3 complex symmetric matrix
    
    R_N k_wave; //!< wave vector
    Complex_wp coef; //!< coefficient
    Complex_wp coef_u; //!< coefficient
    Complex_wp rho; //!< coefficient
    Complex_wp rho0; //!< coefficient
    Real_wp omega; //!< pulsation
    MatrixN_Nsym_Complex_wp mu, invMu; //!< physical index mu
    Complex_wp mu0; //!< physical index mu in vacuum
    bool element_insidePML; //!< element of a PML layer ?
    int num_quad; //!< element number
    int reference_condition; //!<  boundary condition

  private:
    IncidentWaveField<Complex_wp, Dimension>& incident_wave;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension>& var_boundary;
    ImpedanceFunction_Base<Complex_wp, Dimension>& fct_imped_abc;
    
  public :
    template<class TypeEquation>
    TotalWaveSource_Helm(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f);
    
    // volumic source in PML
    void InitElement(int num_elem, const VectR_N& s);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR_N& s);
    void EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f);
    
    // surfacic source on absorbing boundary condition
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };


  //! source for the scattering of an incident wave
  template<class Dimension>
  class TotalWaveSource<HelmholtzEquation<Dimension> >
    : public TotalWaveSource_Helm<Dimension>
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension>& f_inc);
    
  };
  
  
  //! source for the scattering of an incident wave for DG discretization of Helmholtz equation
  /*!
    the diffracted field is computed
   */
  template<class Dimension>
  class DiffractedWaveSource_HelmDG : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::R_N_Complex_wp R_N_Complex_wp; //!< complex point
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!< 2x2 or 3x3 complex symmetric matrix
    
    int reference_condition; //!< boundary condition
    //! coef_vol -> coefficient for volumic integrals  \f$ \int f \varphi \f$
    Complex_wp coef_rho;
    // physical indices : rho, mu, rho0, mu0, sigma_rho, sigma_mu
    Complex_wp rho0; //!< rho in empty media
    Complex_wp rho; //!< rho of current element
    Complex_wp mu0; //!< mu in empty media
    Complex_wp omega; //!< pulsation
    MatrixN_Nsym_Complex_wp mu; //!< mu in current element
    
    //! coefficient for volumic integrals coef_grad \f$ \int f \cdot \varphi \f$
    MatrixN_Nsym_Complex_wp coef_mu;
    bool varying_rho; //!< variable physical index ?
    bool varying_mu; //!< variable physical index ?
   
  private:
    IncidentWaveField<Complex_wp, Dimension>& incident_wave;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension>& var_boundary;

    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_HelmDG(const EllipticProblem<TypeEquation>& vars_helm,
				IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void InitElement(int num_elem, const VectR_N& s);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    // initialization for surfacic integration
    void InitSurface(int i,int num_face,int num_elem,int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  

  //! source for the scattering of an incident wave
  template<class Dimension>
  class DiffractedWaveSource<HelmholtzEquationDG<Dimension> >
    : public DiffractedWaveSource_HelmDG<Dimension>
  {
  public :
    
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension>& f_inc);
  };


  //! source for the scattering of an incident wave for DG methods on Helmholtz equation
  /*!
    the total field is computed
   */
  template<class Dimension>
  class TotalWaveSource_HelmDG : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!< 2x2 or 3x3 complex symmetric matrix

    int reference_condition; //!< boundary condition    
    R_N k_wave; //!< wave vector
    Complex_wp rho; //!< rho in current element
    Complex_wp rho0; //!< rho in empty media
    int num_quad; //!< element number
    Real_wp omega; //!< pulsation
    bool element_insidePML; //!< element inside PML ?
    MatrixN_Nsym_Complex_wp mu; //!< mu in current element
    Complex_wp mu0; //!< mu in empty media
    MatrixN_Nsym_Complex_wp invMu; //!< inverse of mu
    
  private:
    IncidentWaveField<Complex_wp, Dimension>& incident_wave;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension>& var_boundary;

    void InitDefaultValues();

  public :
    template<class TypeEquation>
    TotalWaveSource_HelmDG(const EllipticProblem<TypeEquation>& vars_helm,
			   IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void InitElement(int num_elem, const VectR_N& s);
    
    bool IsNonNullVolumetricSource(const VectR_N& s);
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f);

    void EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR_N& s);

    // surfacic source on absorbing boundary condition
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class TotalWaveSource<HelmholtzEquationDG<Dimension> >
    : public TotalWaveSource_HelmDG<Dimension>
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension>& f_inc);
    
  };

  
  //! projector for incident wave and DG
  template<class Dimension>
  class IncidentWaveProjector_HelmDG : public IncidentWaveProjector<Complex_wp, Dimension>
  {
  private:
    typedef typename Dimension::R_N R_N;

    const VarHelmholtz_Base<Complex_wp, Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_HelmDG(const EllipticProblem<TypeEquation>& var,
				 IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f);
    
  };


  //! volumetric source for DG
  template<class T, class Dimension>
  class VolumetricSource_HelmDG : public VolumetricSource_Base<T, Dimension>
  {
  private :
    const VarHelmholtz_Base<T, Dimension>& var_helm;
    const DistributedProblem<Dimension>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
  public:
    template<class TypeEquation>
    VolumetricSource_HelmDG(const EllipticProblem<TypeEquation>& var,
			    const Vector<VectString>& param);
    
    bool IsNonNullDirichletSource(int cond_ref);        
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);
        
  };


  template<class Dimension>
  class VolumetricSource<HelmholtzEquationDG<Dimension> >
    : public VolumetricSource_HelmDG<Complex_wp, Dimension>
  {
  public :    
    VolumetricSource(const EllipticProblem<HelmholtzEquationDG<Dimension> >&,
		     const Vector<VectString>& param);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_HXX
#endif
