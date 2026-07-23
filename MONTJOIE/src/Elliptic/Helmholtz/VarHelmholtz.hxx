#ifndef MONTJOIE_FILE_VAR_HELMHOLTZ_HXX

namespace Montjoie
{

  template<class Complexe, class Dimension>
  class VarHelmholtz_Base;

  
  template<class Complexe, class Dimension>
  class VarHelmholtzIndex_Base
  {
    typedef typename Dimension::R_N R_N;
  public:
    //! physical index in empty media
    Complexe rho0, mu0, sigma0;

    //! Drude material
    Vector<DrudeMaterial> ref_drude;
    int nb_dof_drude_vec, nb_dof_drude_vec_all;
    Vector<int> OffsetDofDrudeV;
    bool use_symm_drude; bool linearize_drude;
    
  public:
    VarHelmholtzIndex_Base();
    
    virtual Complexe GetRhoTilde(int ref, int num_elem, int npoint) const = 0;
    virtual Complexe GetMuNormale(int ref, int num_elem, int npoint, const R_N& normale) const = 0;
    
    virtual bool UseAdditionalUnknownForABC() const = 0;
    
    inline int GetNbVectorialDofDrude() const { return nb_dof_drude_vec; }
    inline int GetNbVectorialDofDrudeAll() const { return nb_dof_drude_vec_all; }
        
  };


  //! class to specify the resolution of Helmholtz equation
  /*!
    The Helmholtz equation solved is
    \f$ -\omega^2 \rho \, u - \mbox{div}( \nabla u) \, = \, f \f$
  */
  template<class T, class Dim>
  class HelmholtzEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    static bool store_dfjm1;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg = 1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = Dimension::dim_N};
    
    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();
    
    static bool ComputeDFjm1();

    static void ComputeMassMatrixGen(VarProblem<Dimension>& var_problem,
				     VarBoundaryCondition_Dim<T, Dimension>& var_boundary,
				     VarHelmholtz_Base<T, Dimension>& var_helm,
				     int i, const ElementReference<Dimension, 1>&);
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension>&);
    
    template<class T0>
    static TinyMatrix<T0, Symmetric,
		     Dimension::dim_N, Dimension::dim_N>
    GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j, 
			const GlobalGenericMatrix<T0>&, int ref);
    
    template<class T0>
    static T0 GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
				const GlobalGenericMatrix<T0>&, int ref);
    
    template<class T0>
    static TinyVector<T0, Dimension::dim_N>
    GetNablaU_Coefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
			  const GlobalGenericMatrix<T0>& nat_mat, int ref);

    template<class T0>
    static TinyVector<T0, Dimension::dim_N>
    GetNablaPhi_Coefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
			    const GlobalGenericMatrix<T0>& nat_mat, int ref);
    
    template<class T0, class MatStiff>
    static void GetMassPhiDFiTensor(const VarHelmholtz_Base<T, Dimension>& vars,
				    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
				    MatStiff& Amass, bool variable, bool affine,
				    const ElementReference_Dim<Dimension>& Fb);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>&,
                                  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class T0, class MatStiff>
    static void GetGradGradDFiTensor(const VarHelmholtz_Base<T, Dimension>&,
                                     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, MatStiff& Cgrad_grad,
                                     bool variable, bool affine,
				     const ElementReference_Dim<Dimension>& Fb);
    
    template<class T0, class MatStiff>
    static void GetGradPhiDFiTensor(const VarHelmholtz_Base<T, Dimension>&,
                                    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                    MatStiff& Dgrad_phi, MatStiff& Ephi_grad,
                                    bool variable, bool affine,
				    const ElementReference_Dim<Dimension>& Fb);

    template<class Matrix1, class TypeEquation, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
                           int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
			   int ref, int ref2, const EllipticProblem<TypeEquation>& vars,
                           const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class TypeEquation, class T0>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
			   int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);

  };
  
  template<class T, class Dimension>
  bool HelmholtzEquation_Base<T, Dimension>::store_dfjm1(false);
  
  //! class specifying Helmholtz equation
  /*!
    The Helmholtz equation solved is
    \f$ -\omega^2 \rho \, u - i \omega \sigma \, u - \mbox{div}( \mu \nabla u) \, = \, f \f$
  */
  template<class Dimension>
  class HelmholtzEquation : public HelmholtzEquation_Base<Complex_wp, Dimension>
  {
  public :
  };

  
  //! base class for the solution of Helmholtz equation
  /*!
    this class implements specificity of the Helmholtz Equation,
    physical indices and boundary conditions
  */
  template<class Complexe, class Dimension>
  class VarHelmholtz_Base : public VarHelmholtzIndex_Base<Complexe, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!<  d x d complex symmetric matrix
        
    typedef Vector<Complexe> VectComplexe; //!< "complex" vector
    typedef TinyMatrix<Complexe, Symmetric,
                       Dimension::dim_N, Dimension::dim_N> MatrixN_Nsym_Complexe;
    //!< d x d "complex" symmetric matrix
    typedef TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N> MatrixN_N_Complexe;
    //!< d x d "complex" matrix
    typedef TinyVector<Complexe, Dimension::dim_N> R_N_Complexe; //!< "complex" point in R2 or R3 
    
    //! physical index
    Vector<ScalarPhysicalIndice<Dimension, Complexe> > ref_rho, ref_invRho, ref_sigma, ref_beta;
    Vector<TensorPhysicalIndice<Dimension, Dimension::dim_N, Complexe> > ref_mu, ref_invMu;
    //! eventual flow
    Vector<VectorPhysicalIndice<Dimension, Dimension::dim_N, Complexe> > ref_v0;
    //! if true the flow term -2i omega M \cdot \nabla u is included
    bool include_flow_term;
    //! if true we add the term M \cdot \nabla (M \cdot \nabla)
    bool include_double_gradM_flow;
    //! geometric coefficients \omega^2 \rho \omega_k J(xi_k)  and -i omega sigma \omega_k J(xi_k)
    Vector<Vector<Complexe> > Glob_matMass_Dh, Glob_matMass_DhSigma;
    //! coefficients for PML (first-order formulation)
    Vector<Vector<Complexe> > Glob_matMass_DhSigmaDiff,
      Glob_matMass_DhSigmaDiff2, Glob_matMass_DhSigmaDiff3;
    
    Vector<Vector<R_N_Complexe> > Glob_matDamp_Mh;
    //! geometric coefficients \mu J(\xi_k) DF^-1(\xi_k) DF^*-1(\xi_k)
    Vector<Vector<MatrixN_Nsym_Complexe> > Glob_matMass_Bh;
    //! coefficient of impedance for each physical domain (used for absorbing boundary condition)
    Vector<Complexe> coefficient_impedance_absorbing;
    //! coefficient to use in the penalty term
    Vector<Real_wp> coefficient_mu_penalty;

    template<class T>
    void ComputeElementaryMatrixHelm(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				     const GlobalGenericMatrix<T>& nat_mat,
				     const VarProblem<Dimension>& vars,
				     const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
				     const ElementReference<Dimension, 1>& Fb);
    
  private :
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary;
    
    void InitDefaultValues();

  public :
    template<class TypeEquation>
    VarHelmholtz_Base(EllipticProblem<TypeEquation>&);
    
    inline const VarProblem<Dimension>& GetVarProblem() const { return var_problem; }
    
    TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N>
    GetMu_QuadraturePoint(int ref, int i, int num_point) const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    void CopyInputData(const VarHelmholtz_Base<Complexe, Dimension>& var_base);
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complexe>* >& rho_complex,
			   IVect& num_ref, IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);    

    void FinalizeComputationVaryingIndices();    
    bool IsVaryingMedia(int i) const;
    void SetPhysicalIndexAtInfinity(const Vector<bool>&);
    
    Complexe GetRhoTilde(int ref, int num_elem, int npoint) const;
    Complexe GetMuNormale(int ref, int num_elem, int npoint, const R_N& normale) const;

    bool UseAdditionalUnknownForABC() const;    
    bool IsSymmetricProblem(bool eigen = false) const;
    bool IsSymmetricMassMatrix() const;
    
    // PML only with complex coefficients
    void ModifyPMLCoefficient(Complexe& rho, Complexe& sigma, MatrixN_Nsym_Complexe& mu,
                              MatrixN_Nsym_Complexe& invMu, R_N_Complexe& v0,
                              int n, int i1, int i) const;
    
    void AllocateMassMatrices();

    void ComputeDrudeDofs();
    void PutOtherGlobalDofs();

    void FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const;
    
    template<class T>
    void GetPenalizationTauIPP(const VarProblem<Dimension>& vars, int num_face, T& coef_tau) const;

    template<class T>
    void GetCoefficientIPP(const VarProblem<Dimension>& vars,
			   Complexe& m_iomega, Real_wp& coef_ipp, T& coef_u, T& coef_v, T& coef_u_damp) const;
    
    template<class T>
    void ModifyCoefficientBC_HDG(const VarProblem<Dimension>& vars,
				 const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
				 const Mesh<Dimension>& mesh, int ref_boundary, int k,
				 const T& s, const Real_wp&, const T&,
				 int ref_domain, int iquad, int num_loc, T& coef_abc) const;
    
    void ComputeTauCoefficient();

    Complexe ComputeImpedanceCoefficient(const IVect& ref, const VectComplexe& U0);
    
    template<class T>
    void ComputeElementaryMatrixHelmDG(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				       CondensationBlockSolver_Base<T>&,
				       const GlobalGenericMatrix<T>& nat_mat,
				       const VarProblem<Dimension>&,
				       const VarBoundaryCondition_Dim<Complexe, Dimension>&,
				       const ElementReference<Dimension, 1>& Fb);
    
    template<class T, class TypeEquation>
    void AddElementaryFluxesHelmDG(VirtualMatrix<T>& mat_sp, const GlobalGenericMatrix<T>& nat_mat,
				   const EllipticProblem<TypeEquation>& vars,
				   int offset_row, int offset_col);
    
#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
    // implementation in Model_ThinSlot_Helmholtz.cxx
    //! slots present in the domain
    Vector<SlotModelParameters<Dimension2> > list_slots;
    
    void AddContributions1D_Slot(VirtualMatrix<Complex_wp> & mat2D);
    void TreatThinSlotCondition(VirtualMatrix<Complex_wp>& mat_sp);
    
    void LocalizeSlots();
    Real_wp FirstModeSlot2D(const Real_wp& t0, const Real_wp& t1, const Real_wp& ds, int k,
                            const ElementReference<Dimension, 1>& Fb) const;
#endif
    
  };
  

  template<class Dimension>
  class HelmholtzElementaryMatrixClass
  {
  public:
    template<class Complexe, class T>
    static void ComputeMatrixV(const VarProblem<Dimension>& vars,
			       const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
			       const VarHelmholtz_Base<Complexe, Dimension>& var_helm,
			       int iquad, int i, int nb_points_quad, bool affine,
			       int offset_uloc, int offset_uloc2, int offset_vloc,
			       const GlobalGenericMatrix<T>& nat_mat, bool stiff, bool variable,
			       const VectReal_wp& val_phi, const Vector<TinyVector<Real_wp, Dimension::dim_N> >& grad_phi,
			       Vector<T>& feval, Vector<T>& contrib, VirtualMatrix<T>& mat_interac);
    
    template<class Complexe, class T>
    static void ComputeMatrixVdiag(const VarProblem<Dimension>& vars,
                                   const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
                                   const VarHelmholtz_Base<Complexe, Dimension>& var_helm,
                                   int iquad, int nb_points_quad, bool affine, bool variable,
                                   int offset_uloc, int offset_uloc2,
                                   int offset_vloc, const GlobalGenericMatrix<T>& nat_mat,
                                   VirtualMatrix<T>& mat_interac);
    
  };
  
  
  //! base class for the resolution of Helmholtz equation
  /*!
    this class adds specificity of the Helmholtz Equation,
    physical indices and boundary conditions
  */
  template<class TypeEquation>
  class VarHelmholtz_Eq : public VarHelmholtz_Base<typename TypeEquation::Complexe,
						   typename TypeEquation::Dimension>,
			  public VarHarmonic<TypeEquation>
  {
  private :
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename Dimension::R_N R_N;

  public :
    VarHelmholtz_Eq();
    
    // Inline methods;
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void SetTypeEquation(const string& type_equation);
    void CopyInputData(const VarProblem_Base& var_base);

    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    void SetPhysicalIndexAtInfinity(const Vector<bool>&);
    
    inline void GetCoefAB_Infinity( Real_wp& a_infty, Real_wp& b_infty) {
      a_infty = realpart(this->mu0); 
      b_infty = realpart(this->rho0);
    }

    int GetMassMatrixType(Vector<bool>&) const;
    
    void FinalizeComputationVaryingIndices();
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;

    void AllocateMassMatrices();
    void PutOtherGlobalDofs();

    void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp& coef) const;
    void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp& coef) const;

    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;

    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;

    bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;

    bool IsSymmetricProblem(bool eigen = false) const;
    bool IsSymmetricMassMatrix() const;

    void SetFirstOrderFormulation(bool mix_);

    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    Real_wp GetCoefficientPenaltyStiffness(int ref) const;

    ModalSourceBoundary_Dim<Complexe, Dimension>* GetNewModalSourceEquation() const;
    
    // other methods
    void ComputeNumberOfDofs();
        
    void TreatDirichletCondition();
    void PerformOtherInitializations();    
    void ComputeTauCoefficient();

    int GetBoundaryConditionId(const IVect&, int pos, const VectString&, bool& per);
    
    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			       Vector<VectReal_wp>& grad_quadrature,
			       int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			       bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;

    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			       Vector<VectComplex_wp>& grad_quadrature,
			       int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			       bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;

    void ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const;
    
    void ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const;

    void FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const;
    
  };
  
  
  template<class TypeEquation>
  class VarHelmholtz_Cplx : public VarHelmholtz_Eq<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename TypeEquation::Complexe Complexe;

    // Inline methods
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);

    void ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const;

    void ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const;

    void ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;
    
    void ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;

    Complex_wp ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U);
    
  };

  
  //! class to specify the resolution of Helmholtz equation with dg method
  /*!
    The Helmholtz equation solved is
    \f$ -\omega \rho \, u - \mbox{div}( v) \, = \, f \f$
    \f$ -\frac{\omega }{\mu}  v + \mbox{grad}(v) \, = \, g \f$
  */
  template<class T, class Dim>
  class HelmholtzEquationDG_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;

    static bool store_dfjm1;
    static const bool FirstOrderFormulation = true;    
    
    enum{ nb_unknowns = 1+Dimension::dim_N, nb_unknowns_scal = 1, nb_unknowns_hdg = 1,
	  nb_unknowns_vec = Dimension::dim_N, nb_components_en = 1, nb_components_hn = 1};
    
    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int num_elem, const ElementReference_Dim<Dimension>& Fb);

    template<class T0, class Vector1>
    static void GetNeededDerivative(const VarHelmholtz_Base<T, Dimension>& vars,
				    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);    

    template<class T0, class MatStiff>
    static void GetGradPhiTensor(const VarHelmholtz_Base<T, Dimension>& vars,
				 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Dgrad_phi, MatStiff& Ephi_grad); 
    
    template<class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const VarHelmholtz_Base<T, Dimension>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Vn, Vector2& Un);
    
    template<class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const VarHelmholtz_Base<T, Dimension>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Un, Vector2& Vn);
    
  };

  template<class T, class Dimension>
  bool HelmholtzEquationDG_Base<T, Dimension>::store_dfjm1(false);


  //! class to specify the resolution of Helmholtz equation with dg method
  /*!
    The Helmholtz equation solved is
    \f$ -\omega \rho \, u - \mbox{div}( v) \, = \, f \f$
    \f$ -\frac{\omega }{\mu}  v + \mbox{grad}(v) \, = \, g \f$
  */
  template<class Dimension>
  class HelmholtzEquationDG : public HelmholtzEquationDG_Base<Complex_wp, Dimension>
  {
  public :
    template<class TypeEquation, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      int ref, MatMass& mass);
    
    template<class TypeEquation, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<Complex_wp>& nat_mat,
				int ref, Vector1& Un, Vector1& Vn);

    template<class Matrix1, class TypeEquation>
    static void GetNabc(Matrix1& Nabc, typename Dimension::R_N& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<Complex_wp>& nat_mat,
			int ref_d, const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension, 1>& Fb);

    template<class TypeEquation, class Matrix1>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
			   int iquad, int k, int num_face,
			   const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);

    template<class Vector1, class Vector2, class TypeEquation>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int num_face,
			   const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation>
    static void MltNabc(typename Dimension::R_N& normale, int ref,
			const Vector1& Vn, Vector1& Un, int num_elem1, int npoint,
			const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
			const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension, 1>& Fb);
  };


  //! class for the computation of radar cross section for Helmholtz equation
  template<class Dimension>
  class VarComputationRCS_Helm : public VarComputationRCS_Base<Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;

    DistributedProblem<Dimension>& var_problem;
    VarSourceProblem_Base& var_source;
    
  public:
    template<class TypeEquation>
    VarComputationRCS_Helm(EllipticProblem<TypeEquation>& var);
    
    void ComputeRCS(const VectReal_wp& U0);
    void ComputeRCS(const VectComplex_wp& U0);
    
    void ComputeIntegralRepresentation(const VectReal_wp& trace_En, const VectReal_wp& trace_Hn,
                                       const MeshInterpolationFEM<Dimension>& mesh, const R_N& pointX,
                                       VectReal_wp& scal_u) const;

    void ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
                                       const MeshInterpolationFEM<Dimension>& mesh, const R_N& pointX,
                                       VectComplex_wp& scal_u) const;
    
  };

  
  //! class for the computation of radar cross section for Helmholtz equation  
  template<class Dimension>
  class VarComputationRCS<HelmholtzEquation<Dimension> >
    : public VarComputationRCS_Helm<Dimension>
  {
  public :
    VarComputationRCS(EllipticProblem<HelmholtzEquation<Dimension> >& var);
    
  };
  

  //! class for transparent condition for Helmholtz equation
  template<class Dimension>
  class TransparencySolver_Helm : public TransparencySolver_Fem<Dimension, 1, 1>
  {    
  protected :
    typedef typename Dimension::R_N R_N;
    DistributedProblem<Dimension>& var_problem;
    VarHelmholtzIndex_Base<Complex_wp, Dimension>& var_helm;
    
  public :
    template<class TypeEquation>    
    TransparencySolver_Helm(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver);
    
    void GetSource(const VectComplex_wp& trace_u, const VectComplex_wp& trace_du_dn, 
		   int n, const Real_wp& omega, const R_N& point, const R_N& normale,
		   Vector<VectComplex_wp>&, int j) const;
    
    void ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				       const MeshInterpolationFEM<Dimension>& surf_mesh,
				       const R_N& pointX, const R_N& normaleX,
				       TinyVector<Complex_wp, 1>& scal_u,
                                       TinyVector<Complex_wp, 1>& scal_du_dn) const;
    
  };

  
  //! class for transparent condition for Helmholtz equation  
  template<class Dimension>
  class TransparencySolver<HelmholtzEquation<Dimension> >
    : public TransparencySolver_Helm<Dimension>
  {    
  public :
    TransparencySolver(EllipticProblem<HelmholtzEquation<Dimension> >& var,
		       All_LinearSolver& solver);
  };


  //! class for the computation of radar cross section for Helmholtz equation
  template<class Dimension>
  class VarComputationRCS<HelmholtzEquationDG<Dimension> >
    : public VarComputationRCS_Helm<Dimension>
  {
  public :
    VarComputationRCS(EllipticProblem<HelmholtzEquationDG<Dimension> >& var);
    
  };


  //! class for transparent condition for Helmholtz equation
  template<class Dimension>
  class TransparencySolver<HelmholtzEquationDG<Dimension> >
    : public TransparencySolver_Helm<Dimension>
  {    
  public :
    TransparencySolver(EllipticProblem<HelmholtzEquationDG<Dimension> >& var,
		       All_LinearSolver& solver);
  };


  //! static condensation for Laplace/Helmholtz equation
  template<class T, class TypeEquation>
  class CondensationBlockSolver_Helm : public CondensationBlockSolver_Fem<T>
  {
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename TypeEquation::Complexe Complexe;

    EllipticProblem<TypeEquation>& var_problem;

    Vector<int> num_u;
    typename Dimension::VectR_N tau;
    Vector<T> mass_pml;

    IVect SharingProcPML; Vector<IVect> SharingRowsPML;
    
  public:    
    CondensationBlockSolver_Helm(EllipticProblem<TypeEquation>& var);

    void Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix);
    void ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    void RecomposeSolution(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;

    void AssemblePML(Vector<T>& x) const;

    size_t GetMemorySize() const;
    
  };


  //! variables needed for condensed solver (a single element)
  template<class T, class Dimension>
  class CondensedBlockSolver_HelmDG
  {
  public:
    bool flag;
    Matrix<T, Symmetric, RowSymPacked> invBtilde;
    Vector<TinyVector<T, Dimension::dim_N> > Sn, SnD;
    Vector<int> NumU;
    Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > invBh;
    Vector<T> Cl;

    size_t GetMemorySize() const;
    
  };

  //! static condensation for Laplace/Helmholtz equation (HDG)
  template<class T, class Complexe, class Dimension>
  class CondensationBlockSolver_HelmDG : public CondensationBlockSolver_Fem<T>
  {
    VarProblem<Dimension>& var_problem;
    VarHelmholtz_Base<Complexe, Dimension>& var_helm;

  protected:
    Vector<CondensedBlockSolver_HelmDG<T, Dimension> > optim_elt;
    
  public:
    template<class TypeEquation>
    CondensationBlockSolver_HelmDG(EllipticProblem<TypeEquation>& var);

    inline Matrix<T, Symmetric, RowSymPacked>& GetBtildeMatrix(int i) { return optim_elt(i).invBtilde; }
    inline const Matrix<T, Symmetric, RowSymPacked>& GetBtildeMatrix(int i) const { return optim_elt(i).invBtilde; }
    
    inline Vector<TinyVector<T, Dimension::dim_N> >& GetSn(int i) { return optim_elt(i).Sn; }
    inline const Vector<TinyVector<T, Dimension::dim_N> >& GetSn(int i) const { return optim_elt(i).Sn; }
    inline Vector<TinyVector<T, Dimension::dim_N> >& GetSnD(int i) { return optim_elt(i).SnD; }
    inline const Vector<TinyVector<T, Dimension::dim_N> >& GetSnD(int i) const { return optim_elt(i).SnD; }

    inline Vector<int>& GetNumU(int i) { return optim_elt(i).NumU; }
    inline const Vector<int>& GetNumU(int i) const { return optim_elt(i).NumU; }
    
    inline Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> >& GetInverseBh(int i) { return optim_elt(i).invBh; }
    inline const Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> >& GetInverseBh(int i) const { return optim_elt(i).invBh; }
    inline Vector<T>& GetCl(int i) { return optim_elt(i).Cl; }
    inline const Vector<T>& GetCl(int i) const { return optim_elt(i).Cl; }
    
    void SetOptimizedCondensation(int i, bool optim);
    void ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    void RecomposeSolution(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;

    size_t GetMemorySize() const;
    
  };


  template<class T, class Dimension>
  class CondensationBlockSolver<T, HelmholtzEquation<Dimension> >
    : public CondensationBlockSolver_Helm<T, HelmholtzEquation<Dimension> >
  {
  public:
    CondensationBlockSolver(EllipticProblem<HelmholtzEquation<Dimension> >& var);
    
  };
  
  
  //! class used to solve Helmholtz equation with nodal elements
  template<class Dimension>
  class EllipticProblem<HelmholtzEquation<Dimension> >
    : public VarHelmholtz_Cplx<HelmholtzEquation<Dimension> >
  {    
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddBoundaryConditionTerms(VirtualMatrix<Real_wp>& mat_sp,
				   const GlobalGenericMatrix<Real_wp>& nat_mat,
				   int offset_row = 0, int offset_col = 0);
    
    void AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>& mat_sp,
				   const GlobalGenericMatrix<Complex_wp>& nat_mat,
				   int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

  };
  
  
  template<class Dimension>
  class IncidentWaveProjector_HelmDG;

  
  template<class T, class Dimension>
  class CondensationBlockSolver<T, HelmholtzEquationDG<Dimension> >
    : public CondensationBlockSolver_HelmDG<T, Complex_wp, Dimension>
  {
  public:
    CondensationBlockSolver(EllipticProblem<HelmholtzEquationDG<Dimension> >& var);
    
  };
  
  
  //! class used to solve Helmholtz equation with DG method
  template<class Dimension>
  class EllipticProblem<HelmholtzEquationDG<Dimension> >
    : public VarHelmholtz_Cplx<HelmholtzEquationDG<Dimension> >
  {
  public :
    EllipticProblem();

    // Inline methods
    IncidentWaveProjector<Complex_wp, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension>& u_inc) const;

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };
  
}

#define MONTJOIE_FILE_VAR_HELMHOLTZ_HXX
#endif

