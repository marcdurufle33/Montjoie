#ifndef MONTJOIE_FILE_AXISYM_HELMHOLTZ_HXX

namespace Montjoie
{
  //! class to solve Helmholtz equation on axisymmetric domains
  template<class Complexe>
  class VarHelmholtz_Axi : public VarAxisymProblem, public VarHelmholtzIndex_Base<Complexe, Dimension2>
  {
  protected:

    //! envelope is computed ?
    bool calcul_enveloppe;
    bool formulation_r3;
    
  public:
    //! physical coefficients : rho and sigma
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_rho, ref_sigma, ref_beta;
    //! tensor coefficient mu
    /*!
      mu_teta = ref_mu(2, 2)
      mu_zz = ref_mu(1, 1)
      mu_rr = ref_mu(0, 0)
     */
    Vector<TensorPhysicalIndice<Dimension2, 3, Complexe> > ref_mu;
    //! eventual flow
    Vector<VectorPhysicalIndice<Dimension2, 3, Complexe> > ref_v0;
    //! if true the flow term -2i omega M \cdot \nabla u is included
    bool include_flow_term;
    //! if true we add the term M \cdot \nabla (M \cdot \nabla)
    bool include_double_gradM_flow;

    //! diagonal of mass matrix
    Vector<Complexe> mass_matrix_diag;

    //! coefficient of impedance for each physical domain (used for absorbing boundary condition)
    Vector<Complexe> coefficient_impedance_absorbing;

    //! arrays used to compute elementary matrix
    Vector<Vector<Complexe> > Glob_rtilde, Glob_radius, Glob_DrPML, Glob_DzPML;
    
    //! arrays used to perform a fast matrix vector product
    Vector<Vector<Complexe> > Glob_mass_omega2, Glob_mass_m2, Glob_mass_iomega,
      Glob_mass_m_omega;
    Vector<Vector<TinyVector<Complexe, 2> > > Glob_grad_iomega, Glob_grad_im;
    Vector<Vector<TinyMatrix<Complexe, Symmetric, 2, 2> > > Glob_stiff;
    bool use_iterative_solver, same_pattern_LU;
    
  private:
    VarProblem<Dimension2>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary;
    VarOutputProblem<Dimension2>& var_output;
    VarSourceProblem_Base& var_source;
    VarGeneralizedImpedance_Base& var_gibc_base;
    
    void InitDefaultValues();
    
  public:
    template<class TypeEquation>
    VarHelmholtz_Axi(EllipticProblem<TypeEquation>& var);

    // Inline methods
    TinyMatrix<Complexe, Symmetric, 2, 2>
    GetMu_QuadraturePoint(int ref, int i, int num_point) const;

    int GetNbPhysicalIndices() const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    Complexe GetRhoTilde(int ref, int num_elem, int npoint) const;
    Complexe GetMuNormale(int ref, int num_elem, int npoint,
			  const R2& normale) const;
    
    bool UseFormulationR3() const;
    virtual All_LinearSolver* GetNewLinearSolver() = 0;
    
    // other methods
    void SetInputData(const string& description_field, const VectString& parameters);

    bool UseAdditionalUnknownForABC() const;        
    bool IsSymmetricProblem(bool eigen = false) const;
    
    void ConstructFiniteElement(const string& name_elt);
    
    void ConstructAll(const string& name_file, const string& name_elt,
		      Vector<string>& lines_data_file);
    
    void ConstructAll(const Vector<string>& lines_data_file, const string& input_file,
		      const string&, const string&, All_LinearSolver*&,
		      bool compute_rho = true, bool delete_points = true, int num = -1);
    
    // methods related to indexes
    void WriteAllIndices();
    void PerformAdimensionalization(const Real_wp& t0_adim, const Real_wp& z0_adim);
    
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex,
			   IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
        
    bool IsVaryingMedia(int i) const;

    void SetPhysicalIndexAtInfinity(const Vector<bool>&);
    
    void UpdateDirichlet(int n);    
    void AllocateMassMatrices();
    
    void ComputeLocalMassMatrix(int i);

    void ModifyCoefficientPML(bool form_r3, TinyMatrix<Complexe, Symmetric, 2, 2>& muPML,
			      Complexe& rtilde, Complexe& invRtilde, Complexe& coef,
			      Complexe& coef_cte, Complexe& coef_grad, TinyVector<Complexe, 3>& v0,
			      int iquad, int i1, int j);
    
    template<class T>
    void ComputeElementaryMatrixHelmAxi(int iquad, IVect& num_dof, VirtualMatrix<T>& mat,
					const GlobalGenericMatrix<T>& nat_mat,
					const ElementReference<Dimension2, 1>& Fb);

    Complexe ComputeImpedanceCoefficient(const IVect& ref, Vector<Complexe>& U);
    
  };

  template<class TypeEquation>
  class VarHelmholtzAxi_Eq : public VarHarmonic<TypeEquation>,
    public VarHelmholtz_Axi<typename TypeEquation::Complexe>
  {
    typedef typename TypeEquation::Complexe Complexe;

  public:
    VarHelmholtzAxi_Eq();
    
    // Inline methods
    bool UseNumericalIntegration(int i) const;
    const R3& GetPhaseOrigin() const;
    int GetNbModesSource() const;
    
    IncidentWaveField<Real_wp, Dimension2>* GetNewIncidentField(int n, const Vector<VectString>& param,
								const Real_wp& val) const;
    
    IncidentWaveField<Complex_wp, Dimension2>* GetNewIncidentField(int n, const Vector<VectString>& param,
								   const Complex_wp& val) const;

    void SetInputData(const string& description_field, const VectString& parameters);
    void CheckInputMesh();

    int GetNbPhysicalIndices() const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void ConstructFiniteElement(const string& name_elt);

    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    void SetPhysicalIndexAtInfinity(const Vector<bool>& RefUsed);

    void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp& coef) const;
    void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp& coef) const;

    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;

    void AllocateMassMatrices();
    
    void ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
				SetPoints<Dimension2>& PointsElem,
                                SetMatrices<Dimension2>& MatricesElem,
                                IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb);
    
    // other methods
    void PerformOtherInitializations();
    
    int GetBoundaryConditionId(const IVect&, int pos, const VectString&, bool& per);

    void AddDiracSource(const Real_wp& alpha_, Vector<Vector<Real_wp> > & b_source,
			Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const;
    
    void AddDiracSource(const Complex_wp& alpha_, Vector<Vector<Complex_wp> > & b_source,
			Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const;
    
    void ComputeDiracSource(Vector<Complexe>& b_source, const VectR3& point_source3D,
                            const Vector<VectString>& param_src) const;
    
    void ComputeRightHandSide(Vector<Real_wp>& b_rhs, bool assemble = false);
    void ComputeRightHandSide(Vector<Complex_wp>& b_rhs, bool assemble = false);

    void ComputeRightHandSide(Vector<Vector<Complexe> >& b_rhs, bool assemble = false);
    
    void ModifyOutputUnknown(Vector<Complex_wp>& val_u,
			     Vector<Complex_wp>& grad_u,
			     int i, const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(Vector<Real_wp>& val_u,
			     Vector<Real_wp>& grad_u,
			     int i, const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;
      
    void ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
			     int, bool, bool) const;
    
    void ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
			     int, bool, bool) const;

    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			       Vector<VectReal_wp>& grad_quadrature,
			       int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			       bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;
    
    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			       Vector<VectComplex_wp>& grad_quadrature,
			       int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			       bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;
    
    void ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const;
    
    void ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const;
    
    bool IsSymmetricProblem(bool eigen = false) const;

    void ConstructAll(const string& name_file, const string& name_element,
		      Vector<string>& lines_data_file);
    
    void ConstructAll(const Vector<string>& lines_data_file, const string& input_file,
		      const string&, const string&, All_LinearSolver*&,
		      bool compute_rho = true, bool delete_points = true, int num = -1);

    All_LinearSolver* GetNewLinearSolver();    
    
  };
  

  template<class TypeEquation>
  class VarHelmholtzAxi_Cplx : public VarHelmholtzAxi_Eq<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe;
    
    // Inline methods
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    IncidentWaveProjector<Complex_wp, Dimension2>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const;
    
    Complex_wp ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U);

    void ComputeRightHandSide(Vector<Complex_wp>& b_rhs, bool assemble = false);    
  };


  template<class TypeEquation>
  class VarHelmholtzAxi_Real : public VarHelmholtzAxi_Eq<TypeEquation>
  {
  public:
    typedef typename TypeEquation::Complexe Complexe;
    
    // Inline methods
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
        
  };
  

  class HelmholtzEquationAxi : public GenericEquation<Complex_wp>
  {
  public :
    typedef Dimension2 Dimension;

    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 3};

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    static void ComputeMassMatrix(EllipticProblem<HelmholtzEquationAxi>& var,
				  int i, const ElementReference_Dim<Dimension2>& Fb);
    
    template<class TypeEquation, class T, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass);

    static TinyMatrix<Complex_wp, Symmetric, 2, 2>
    GetStiffCoefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j, 
			const GlobalGenericMatrix<Complex_wp>&, int ref);
    
    static Complex_wp GetMassCoefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
				const GlobalGenericMatrix<Complex_wp>&, int ref);
    
    static TinyVector<Complex_wp, 2>
    GetNablaU_Coefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref);

    static TinyVector<Complex_wp, 2>
    GetNablaPhi_Coefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref);
    
  };


  class LaplaceEquationAxi : public GenericEquation<Real_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 3};

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();
    
    template<class TypeEquation, class T, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass);
    
  };


  //! class for absorbing boundary condition for Helmholtz equation in axisymmetric domain
  template<>
  class ImpedanceABC<Complex_wp, HelmholtzEquationAxi> : public ImpedanceABC_Helm<Dimension2>
  {    
  private:
    const VarHelmholtz_Axi<Complex_wp>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;
    
  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
    //! evaluation of impedance coefficient
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                              const SetPoints<Dimension2>& Pts,
                              const SetMatrices<Dimension2>& Mat);
  };


  //! class for absorbing boundary condition for Laplace equation in axisymmetric domain
  template<class T>
  class ImpedanceABC<T, LaplaceEquationAxi>
    : public ImpedanceFunction_Base<T, Dimension2>
  {
  private:
    const VarHelmholtz_Axi<Real_wp>& var_helm;
    
  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
    //! evaluation of impedance coefficient
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<T>& nat_mat, int ref,
                              const SetPoints<Dimension2>& Pts,
                              const SetMatrices<Dimension2>& Mat);
    
  };


  //! base class for impedance boundary condition
  template<class Complexe>
  class ImpedanceGenericAxiHelm : public ImpedanceFunction_Base<Complexe, Dimension2>
  {
  private:
    const VarHelmholtz_Axi<Complexe>& var_helm;
    const VarBoundaryCondition_Base& var_boundary;
    
  public :
    template<class TypeEquation>
    ImpedanceGenericAxiHelm(const EllipticProblem<TypeEquation>& var);
    
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complexe>& nat_mat, int ref_d,
                              const SetPoints<Dimension2>& Pts,
                              const SetMatrices<Dimension2>& Mat);
    
  };
  
  
  //! class for impedance boundary conditions and axisymmetric Helmholtz
  template<>
  class ImpedanceGeneric<Complex_wp, HelmholtzEquationAxi> : public ImpedanceGenericAxiHelm<Complex_wp>
  {
  public:
    template<class TypeEquation>
    ImpedanceGeneric(const EllipticProblem<TypeEquation>& var);

  };


  //! class for impedance boundary conditions and axisymmetric Helmholtz
  template<>
  class ImpedanceGeneric<Real_wp, LaplaceEquationAxi> : public ImpedanceGenericAxiHelm<Real_wp>
  {
  public:
    template<class TypeEquation>
    ImpedanceGeneric(const EllipticProblem<TypeEquation>& var);
    
  };
  
  
  //! plane wave incident field
  /*!
    the function f is equal to \f$ f = \mbox{exp}( i k \cdot (x-x_0) ) \f$
    where k is the wave vector, and x_0 the origine of the phase     
  */
  class PlaneWaveAxiIncidentField
    : public IncidentWaveField<Complex_wp, Dimension2>
  {
  protected :
    Real_wp k_perp, kz;
    Complex_wp phase;
    int number_mode;
    bool incidence_axial;
    
  public :
    PlaneWaveAxiIncidentField(const R3&, const Complex_wp&, const Complex_wp&, const Complex_wp&,
			      const Real_wp&, int m);
    
    void EvaluateFunction(const R2& pt, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R2& pt, Complex_wp& f,
				  TinyVector<Complex_wp, 2>& df) const;
    
  };


  class IncidentWaveProjector_HelmAxi
    : public IncidentWaveProjector<Complex_wp, Dimension2>
  {
  private:
    const VarHelmholtz_Axi<Complex_wp>& var_helm;
    Real_wp omega; int number_mode;
    Real_wp k_perp, kz; bool incidence_axial; Complex_wp phase;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    IncidentWaveProjector_HelmAxi(const EllipticProblem<TypeEquation>& var,
				  IncidentWaveField<Complex_wp, Dimension2>& u_inc);
    
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
  };


  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  class DiffractedWaveSource_HelmAxi
    : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  protected :
    //! components of the wave vector : k_\perp and k_z
    Real_wp k_perp, kz;
    Real_wp number_mode;
    //! true if the wave vector is oriented in e_z
    bool incidence_axial;    
    //! coefficients for the incident field
    Complex_wp rho0, mu0;
    //! coefficients for volume integrals
    Complex_wp coef_vol, coef_grad, coef_m2, phase;
    int reference_condition;
    bool formulation_r3;
    
  private:
    IncidentWaveField<Complex_wp, Dimension2>& incident_wave;
    const VarProblem<Dimension2>& var_problem;
    const VarHelmholtz_Axi<Complex_wp>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;
    ImpedanceGenericAxiHelm<Complex_wp>& fct_imped;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEq>
    DiffractedWaveSource_HelmAxi(const EllipticProblem<TypeEq>& vars_helm,
				 IncidentWaveField<Complex_wp, Dimension2>& fsrc_);
    
    void InitElement(int num_elem, const VectR2& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR2& s);
    
    // volumic integral \int f grad(\phi)
    void EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR2& s);
    
    // initialization for surfacic integration
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
				 const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
    // surfacic integral \int_\Gamma f \cdot \nabla_\Gamma \phi
    void EvaluateSurfacicSourceGradient(int k,const SetPoints<Dimension2>& PointsElem,
					 const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    
  };
  
  //! source for the scattering of an incident wave
  template<>
  class DiffractedWaveSource<HelmholtzEquationAxi>
    : public DiffractedWaveSource_HelmAxi
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension2>& fsrc_);
    
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  class TotalWaveSource_HelmAxi : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  protected :
    //! components of the wave vector : k_\perp and k_z
    Real_wp k_perp, kz;
    Real_wp number_mode;
    //! coefficients for the incident field
    Complex_wp rho0, mu0, k0;
    //! true if the wave vector is oriented in e_z
    bool incidence_axial;    
    Complex_wp phase;
    bool formulation_r3;
    
  private:
    IncidentWaveField<Complex_wp, Dimension2>& incident_wave;
    const VarProblem<Dimension2>& var_problem;
    const VarHelmholtz_Axi<Complex_wp>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    TotalWaveSource_HelmAxi(const EllipticProblem<TypeEquation>& vars_helm,
			    IncidentWaveField<Complex_wp, Dimension2>&);

    void EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f);
    
    // volumic source in PML
    void InitElement(int num_elem, const VectR2& s);
    bool IsNonNullVolumetricSource(const VectR2& s);
    
    void EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR2& s);
    
    void EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f);
    
    // surfacic source on absorbing boundary condition
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
				 const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<>
  class TotalWaveSource<HelmholtzEquationAxi> : public TotalWaveSource_HelmAxi
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension2>& fsrc_);
    
  };
  

  //! class for GIBC and Helmholtz equation
  template<>
  class VarGeneralizedImpedance<HelmholtzEquationAxi> 
    : public VarGeneralizedImpedance_Helm<Dimension2>
  {    
  public:
    template<class TypeEquation>
    VarGeneralizedImpedance(EllipticProblem<TypeEquation>& var);
    
  };
  

  //! volumetric source
  template<class T>
  class VolumetricSource_AxiHelm
    : public VirtualSourceFEM<T, Dimension2>
  {
  private:
    const VarBoundaryCondition_Base& var_boundary;
    const VarHelmholtz_Axi<T>& var_helm;
    const DistributedProblem<Dimension2>& var_problem;
    
  protected:
    VirtualSourceField<T, Dimension3>* fsrc;
    bool surface_source; int ref_surf;
    bool external_source;
    Vector<bool> unif_source_vol, unif_source_surf;
    Vector<T> coef_source_vol, coef_source_surf;
    
  public :
    T coef_vol; Real_wp rmax; R2 center_source;
    Globatto<Real_wp> gauss_teta; //!< quadrature formula over theta
    
    template<class TypeEquation>
    VolumetricSource_AxiHelm(const EllipticProblem<TypeEquation>& var,
			     const Vector<VectString>& param);
    
    ~VolumetricSource_AxiHelm();


    void SetVolumeSourceFunction(VirtualSourceField<T, Dimension3>* f);
    void SetSurfaceSource(int ref, VirtualSourceField<T, Dimension3>* f);
    
    bool IsNonNullVolumetricSource(const VectR2& s);
    void EvaluateSourceMode(const R2& x, Vector<T>& f);
    void EvaluateVolumetricSource(int i, int j, const R2& x, Vector<T>& f);
    
    bool IsNonNullSurfacicSource(int ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
				const SetMatrices<Dimension2>& MatricesElem, Vector<T>& f);

    void EvaluateFunction(int i, int j, const R2& x, Vector<T>& f);    
    
  };
  
  
  //! volumetric source
  template<>
  class VolumetricSource<HelmholtzEquationAxi> : public VolumetricSource_AxiHelm<Complex_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
    
  };


  //! volumetric source
  template<>
  class VolumetricSource<LaplaceEquationAxi> : public VolumetricSource_AxiHelm<Real_wp>
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
    
  };
  
  
  //! matrix-free class
  template<>
  class FemMatrixFreeClass<Complex_wp, HelmholtzEquationAxi>
    : public FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquationAxi>
  {
  public:
    FemMatrixFreeClass(const EllipticProblem<HelmholtzEquationAxi>& var_);

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };
  

  //! class to solve Helmholtz equation in axisymmetric domain 
  template<>
  class EllipticProblem<HelmholtzEquationAxi>
    : public VarHelmholtzAxi_Cplx<HelmholtzEquationAxi>
  {
  };


  //! class to solve Laplace equation in axisymmetric domain 
  template<>
  class EllipticProblem<LaplaceEquationAxi>
    : public VarHelmholtzAxi_Real<LaplaceEquationAxi>
  {
  };

}

#define MONTJOIE_FILE_AXISYM_HELMHOLTZ_HXX
#endif
