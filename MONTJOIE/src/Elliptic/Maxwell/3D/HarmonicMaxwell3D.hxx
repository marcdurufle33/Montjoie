#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_3D_HXX

namespace Montjoie
{
  
  //! class to specify the resolution of 3-D Maxwell's equations
  template<class T>
  class MaxwellEquation3D_Base : public GenericEquation<T>
  {
  public :    
    typedef Dimension3 Dimension;
    
    static bool store_dfjm1;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 3, nb_components_hn = 3, nb_unknowns_vec = 3,
          type_element = 2};
    
    static inline bool ComputeDFjm1() { return store_dfjm1; }
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int num_elem, const ElementReference_Dim<Dimension3>& Fb);
    
    template<class TypeEquation, class T0>
    static void GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
				    int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef);
    
  };
  

  //! class to specify the resolution of 3-D time-harmonic Maxwell equations
  /*!
    The time-harmonic maxwell equation solved is
    \f$ -\omega^2 \epsilon \, E + \mbox{curl}( \frac{1}{\mu} \mbox{curl} E ) \, = \, f \f$
  */  
  class HarmonicMaxwellEquation_3D : public MaxwellEquation3D_Base<Complex_wp>
  {
  };

  
  //! base class for solution of Maxwell's equations in 3-D
  class HarmonicMaxwell_3D_Base
  {
  public:
    Vector<Vector<Matrix3_3> > Glob_DFj; // DF matrix for computation of curl
    Vector<Matrix3_3> Glob_DFjm1_Neighbor; // for SIPG formulation
    Vector<VectReal_wp> Glob_invJacobian; //!< inverse of jacobian
    string file_coefficient_q; //!< file where coefficients Q_i are stored

    // for Drude-Lorentz Material
    Vector<DrudeMaterial> ref_drude;
    int nb_dof_drude_vec,nb_dof_drude_vec_all; Vector<int> OffsetDofDrudeV;
    bool use_symm_drude; bool linearize_drude;
    
  private:
    VarProblem<Dimension3>& var_problem;
    VarComputationProblem& var_computation;
    VarBoundaryCondition_Base& var_boundary;
    VarOutputProblem_Dim<Dimension3>& var_output;
#ifdef MONTJOIE_WITH_TRANSMISSION
    VarTransmission_Base<Dimension3>& var_transmission_base;
#endif
    
  public:
    template<class TypeEquation>
    HarmonicMaxwell_3D_Base(EllipticProblem<TypeEquation>& var);
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    inline void ModifyPhysicalCoefPML(Matrix3_3sym& epsilon,
				      Matrix3_3sym& mu, Matrix3_3sym& sigma,
				      int i1, int i) const {}
    
    void ComputeMassMatrix(bool sipg);
    
    bool IsSymmetricProblem(bool eigen = false) const;
    void ComputeDrudeDofs();
    
    template<class T>
    void ComputeEnHnQuadratureHcurl(Vector<Vector<T> >& u_quadrature,
			       Vector<Vector<T> >& curl_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, Vector<T>& En_quad, Vector<T>& Hn_quad) const;
    
    template<class T>
    void ComputeEnHnQuadratureDG(Vector<Vector<T> >& u_quadrature,
                                 Vector<Vector<T> >& grad_quadrature,
                                 int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                                 bool compute_H, Vector<T>& En_quad, Vector<T>& Hn_quad) const;
    
    template<class T>
    void ComputeEnHnNodal(Vector<Vector<T> >& u_nodal, Vector<Vector<T> >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<Vector<T> >& En_nodal, Vector<Vector<T> >& Hn_nodal) const;

    inline int GetDrudeDof(int) const { return -1; }
    inline int GetNbVectorialDofDrude() const { return nb_dof_drude_vec; }
    inline int GetNbVectorialDofDrudeAll() const { return nb_dof_drude_vec_all; }
    
    void ComputeNumberOfDofs();
    void PutOtherGlobalDofs();

    IVect GetScalarDofNumberOnElement(int i) const;
    
    void TreatDirichletCondition();
    
    void WriteCoefQ(Vector<VectComplex_wp> & U0);

    virtual void ComputeCoefficientsQ(const VectComplex_wp& U0, VectReal_wp& coefQ) const = 0;
    
    template<class Vect, class Matrix1, class Complexe>
    void ComputeVariableElementaryMatrix(int iquad, const Vect& B, const Vect& A, bool compute_mass,
					 Matrix1& mat, const ElementReference<Dimension3, 2>& Fb,
					 const GlobalGenericMatrix<Complexe>& nat_mat);

    void FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const;

    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
                                 Vector<int>& intern_node) const;
    
  };


  template<class T>
  class MaxwellGaussianSource;

  template<class Complexe>
  class HarmonicMaxwell3D_PhysGeomInfo
  {
  protected:
    Vector<TinyMatrix<Complexe, Symmetric, 3, 3> > Ah_Bh;
    int offset_Ah, offset_Bh, offset_BhSigma, offset_diagDhSigma, offset_diagDhStiff;
    Vector<Complexe> diag_Dh;    
    Vector<TinyVector<Real_wp, 3> > tauPML;
    Vector<TinyMatrix<Complexe, Symmetric, 2, 2> > BhSurf;
    
  public:
    HarmonicMaxwell3D_PhysGeomInfo();

    size_t GetMemorySize() const;
    
    void ReallocateAhBh(int N1, int N2, int N3);
    void ReallocateBhSurf(int N);

    void SetAh(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M);
    void SetBh(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M);
    void SetBhSigma(int i, const TinyMatrix<Complexe, Symmetric, 3, 3>& M);
    void SetBhSurf(int i, const TinyMatrix<Complexe, Symmetric, 2, 2>& M);
    
    const TinyMatrix<Complexe, Symmetric, 3, 3>& GetBh(int i) const;
    const TinyMatrix<Complexe, Symmetric, 3, 3>& GetBhSigma(int i) const;
    const TinyMatrix<Complexe, Symmetric, 3, 3>& GetAh(int i) const;
    const TinyMatrix<Complexe, Symmetric, 2, 2>& GetBhSurf(int i) const;
    
    void ReallocateDhDiag(int N, int Ns = 0);
    void ReallocateTauPML(int N);

    void SetDh(int i, const Complexe& d);
    void SetDhSigma(int i, const Complexe& d);
    void SetDhStiff(int i, const Complexe& d);
    void SetTauPML(int i, const TinyVector<Real_wp, 3>& t);

    bool IsDiagonalMass() const;
    const Complexe& GetDh(int i) const;
    const Complexe& GetDhSigma(int i) const;
    const Complexe& GetDhStiff(int i) const;
    const TinyVector<Real_wp, 3>& GetTauPML(int i) const;
    
    void ClearBhSigma();
    void ClearDiagonalDh();
    
  };


  //! class to solve 3-D time-harmonic Maxwell equations
  /*!
    The time-harmonic maxwell equation solved is
    \f$ -\omega^2 \epsilon \, E + \mbox{curl}( \frac{1}{\mu} \mbox{curl} E ) \, = \, f \f$
    Edge finite elements are used
  */
  template<class Complexe>
  class HarmonicMaxwell_3D : public HarmonicMaxwell_3D_Base
  {
  public :
    typedef Vector<Complexe> VectComplexe; //!< vector of float or complex floats
    typedef TinyMatrix<Complexe,Symmetric,3,3> Matrix3_3sym_Complexe;
    
    //!< 3x3 "complex" symmetric matrix
    //! dieletric permittivity
    Vector<TensorPhysicalIndice<Dimension3, 3, Complexe> > ref_epsilon;
    Vector<TensorPhysicalIndice<Dimension3, 3, Complexe> > ref_invEpsilon;
    //!< inverse of dieletric permittivity
    Vector<TensorPhysicalIndice<Dimension3, 3, Complexe> > ref_sigma;
    //! magnetic permeability
    Vector<TensorPhysicalIndice<Dimension3, 3, Complexe> > ref_mu;
    Vector<TensorPhysicalIndice<Dimension3, 3, Complexe> > ref_invMu; 
    //!< inverse of magnetic permeability

    //! dieletric permittivity in empty space
    Complexe epsilon0;
    //! magnetic permeability in empty space
    Complexe mu0;
    Complexe invMu0; //!< inverse of magnetic permeability in the vacuum
    Vector<Complexe> coefficient_impedance_absorbing;

    Vector<HarmonicMaxwell3D_PhysGeomInfo<Complexe> > Glob_matMass_elem;
    
  private:
    VarProblem<Dimension3>& var_problem;
    VarComputationProblem& var_computation;
    VarOutputProblem_Dim<Dimension3>& var_output;
    VarBoundaryCondition_Dim<Complexe, Dimension3>& var_boundary;
    VarSourceProblem_Base& var_source_base;
    VarSourceProblem_Cplx<Complexe, Dimension3>& var_source;
    
  public :
    template<class TypeEquation>
    HarmonicMaxwell_3D(EllipticProblem<TypeEquation>& var);
    
    void GetMemoryUsed(map<string, size_t>& var) const;

    void CopyInputData(const HarmonicMaxwell_3D<Complexe>& var_base);
    
    template<class T0, class Prop>
    void GetMassMatrix(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       TinyMatrix<T0, Prop, 3, 3>& Bmass) const;

    template<class T0, class Prop>
    void GetStiffMatrix(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
			TinyMatrix<T0, Prop, 3, 3>& Astiff) const;
    
    void ModifyPhysicalCoefPML(Matrix3_3sym_Complexe& epsilon,
			       Matrix3_3sym_Complexe& mu, Matrix3_3sym_Complexe& sigma,
			       int i1, int i) const;
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complexe>* >& rho_complex,
                           IVect& num_ref, IVect& num_index, IVect& num_component,
                           Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void FinalizeComputationVaryingIndices();
    bool IsVaryingMedia(int i) const;
    
    void AllocateMassMatrices();
    void ComputeCoefficientsQ(const VectComplex_wp& U0, VectReal_wp& coefQ) const;

    template<class T>
    void ComputeElementaryMatrix(int iquad, IVect& num_dof,
				 VirtualMatrix<T>& mat_interac,
				 const GlobalGenericMatrix<T>& nat_mat,
				 const ElementReference<Dimension3, 2>& Fb);

    template<class T>
    void ComputeElementaryMatrixHdg(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				    CondensationBlockSolver_Base<T>& solver,
				    const GlobalGenericMatrix<T>& nat_mat,
				    const ElementReference<Dimension3, 2>& Fb);
    
    template<class T>
    void ModifyCoefficientBC_HDG(const Mesh<Dimension3>& mesh, int ref_boundary, int k,
				 const T& s, const Complexe& coef_v,
				 int ref_domain, int iquad, int num_loc, T& coef_abc) const;

    void ComputeTauCoefficient();
    
    template<class T>
    void AddElementaryFluxesSipg(VirtualMatrix<T>& mat_sp,
				 const GlobalGenericMatrix<T>& nat_mat,
				 int offset_row, int offset_col);
    
    template<class T>
    void AddElementaryFluxSipg(VirtualMatrix<T>& mat_sp,
			       const GlobalGenericMatrix<T>& nat_mat,
			       int iquad, int num_pos1_face, int num_elem2,
			       int num_face, int ref, bool new_face, int& nb_neighbor,
			       int offset_row, int offset_col,
			       const ElementReference<Dimension3, 2>& Fb);
    
  };
  
  
  template<class TypeEquation>
  class HarmonicMaxwell_3D_Eq : public VarHarmonic<TypeEquation>,
				public HarmonicMaxwell_3D<typename TypeEquation::Complexe>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    HarmonicMaxwell_3D_Eq();
    
    void SetFirstOrderFormulation(bool mix_);
    void SetTypeEquation(const string& type_equation);
    
    void ComputeNumberOfDofs();
    void PutOtherGlobalDofs();
    IVect GetScalarDofNumberOnElement(int i) const;
    
    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
                                 Vector<int>& intern_node) const;
    
    void PerformOtherInitializations();
    void ComputeTauCoefficient();
    
    void TreatDirichletCondition();
  
    void WriteDatas(Vector<VectReal_wp>& U0);
    void WriteDatas(Vector<VectComplex_wp>& U0);

    void GetMemoryUsed(map<string, size_t>& var) const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    void CopyInputData(const VarProblem_Base& var_base);

    bool DiagonalScalarMassMatrix(const ElementReference_Dim<Dimension3>&, int i) const;
    bool BlockDiagonalScalarMassMatrix(const ElementReference_Dim<Dimension3>&, int i) const;
    
    bool IsSymmetricProblem(bool eigen = false) const;

    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void FinalizeComputationVaryingIndices();
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;

    void AllocateMassMatrices();    
    void ComputeMassMatrix(bool compute_rho = true, bool delete_points = true);

    VirtualSourceField<Complexe, Dimension3>*
    GetNewVolumeSourceFunction(const IVect&, const VectString& param,
                               int& nb, Vector<Complexe>& polar,
                               VolumetricSource_Base<Complexe, Dimension3>&) const;

    VirtualSourceField<Complexe, Dimension3>*
    GetNewSurfaceSourceFunction(const IVect&, const VectString& param,
				int& nb, Vector<Complexe>& polar,
				VolumetricSource_Base<Complexe, Dimension3>&) const;
    
    ModalSourceBoundary_Dim<Complexe, Dimension3>* GetNewModalSourceEquation() const;

    void FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const;
    
  };
  

  template<class TypeEquation>
  class HarmonicMaxwell_3D_Cplx : public HarmonicMaxwell_3D_Eq<TypeEquation>
  {
  public :
    // Inline methods
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension3, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
  };
    
  
  //! class for Silver-Muller condition and maxwell equations
  template<class T, class Complexe>
  class ImpedanceABC_Maxwell3D
    : public ImpedanceFunction_Base<T, Dimension3>
  {
  private:
    const HarmonicMaxwell_3D<Complexe>& var_maxwell;
    const VarBoundaryCondition_Dim<Complexe, Dimension3>& var_boundary;
    const VarProblem<Dimension3>& var_problem;
    
  public :
    template<class TypeEquation>
    ImpedanceABC_Maxwell3D(const EllipticProblem<TypeEquation>& var);
    
    //! impedance coefficient
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<T>& nat_mat, int ref,
                              const SetPoints<Dimension3>& Pts, const SetMatrices<Dimension3>& Mat);
    
  };


  //! class for Silver-Muller condition and 3-D time-harmonic maxwell equations
  template<>
  class ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_3D>
    : public ImpedanceABC_Maxwell3D<Complex_wp, Complex_wp>
  {
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

  };

  
  //! class to specify the resolution of time-harmonic Maxwell equations with dg method
  /*!
    The maxwell equation solved is
    \f$ - \omega^2 \varepsilon \, E + \mbox{curl}( H) \, = \, f \f$
    \f$ \mu H - \mbox{curl} E \, = \, g \f$
   */
  class HarmonicMaxwellEquation_3D_DG : public GenericEquation<Complex_wp>
  {
  public :
    typedef Dimension3 Dimension;

    static const bool FirstOrderFormulation = true;
    static bool store_dfjm1;
    
    enum {nb_unknowns = 6, nb_unknowns_scal = 3, nb_unknowns_vec = 3,
	  nb_components_en = 3, nb_components_hn = 3, nb_unknowns_hdg=1};

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    template<class TypeEquation>    
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int i, const ElementReference_Dim<Dimension3>&);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& var,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			      int ref, MatMass& mass);

    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref,
				Vector1& U, Vector1& V);
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
				    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc,
				 const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Ephi_grad, MatStiff& Dgrad_phi);
    
    template<class TypeEquation,
	     class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Vn, Vector2& Un);
    
    template<class TypeEquation,
	     class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				     int ref, Vector1& Un, Vector2& Vn);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, R3& normale, int ref, int iquad, int k,
			const GlobalGenericMatrix<T0>& nat_mat, int ref2, const GenericPb& vars,
			const ElementReference<Dimension3, 1>& Fb);

    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, R3& normale, int iquad, int k, int nf,
			   const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dimension3, 1>& Fb);
      
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R3& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
			   int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension3, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    void MltNabc(R3& normale, int ref, const Vector1& Vn, Vector1& Un,
		 int num_elem1, int k, const GlobalGenericMatrix<T0>& nat_mat,
		 int ref2, const EllipticProblem<TypeEquation>& vars,
		 const ElementReference<Dimension3, 1>& Fb);
    
  };
  

  //! base class for computing radar cross section for 3-D Maxwell's equations
  class VarComputationRCS_Maxwell3D : public VarComputationRCS_Base<Dimension3>
  {
  private:
    VarProblem<Dimension3>& var_problem;
    
  public :
    template<class TypeEquation>
    VarComputationRCS_Maxwell3D(EllipticProblem<TypeEquation>& var);
    
    void ComputeRCS(const VectReal_wp& U0);
    void ComputeRCS(const VectComplex_wp& U0);
    
  };
  
  
  //! class for computing radar cross section for 3-D Maxwell's equations
  template<>
  class VarComputationRCS<HarmonicMaxwellEquation_3D>
    : public VarComputationRCS_Maxwell3D
  {
  public :
    template<class TypeEquation>
    VarComputationRCS(EllipticProblem<TypeEquation>& var);

  };

  //! class for computing radar cross section for 3-D Maxwell's equations  
  template<>
  class VarComputationRCS<HarmonicMaxwellEquation_3D_DG>
    : public VarComputationRCS_Maxwell3D
  {
  public :
    template<class TypeEquation>
    VarComputationRCS(EllipticProblem<TypeEquation>& var);
    
  };

  
  //! base class for transparent condition for 3-D Maxwell's equations
  class TransparencySolver_Maxwell3D : public TransparencySolver_Fem<Dimension3, 3, 3>
  { 
  private:
    VarProblem_Base& var_problem;
    
  public :
    //! default constructor
    template<class TypeEquation>
    TransparencySolver_Maxwell3D(EllipticProblem<TypeEquation>& var,
				 All_LinearSolver& solver);
    
    void ComputeIntegralRepresentation(const VectComplex_wp& trace_En,
                                       const VectComplex_wp& trace_Hn,
				       const MeshInterpolationFEM<Dimension3>& mesh_,
				       const R3& pointX, const R3& normaleX,
				       R3_Complex_wp& scal_u, R3_Complex_wp& scal_du_dn) const;
    
    void GetSource(const VectComplex_wp& scal_u,
		   const VectComplex_wp& scal_du_dn,
		   int n, const Real_wp& k_inf, const R3& point, const R3& normale,
		   Vector<VectComplex_wp> &, int j) const;
    
  };
  
  
  //! class for transparent condition for 3-D Maxwell's equations
  template<>
  class TransparencySolver<HarmonicMaxwellEquation_3D>
    : public TransparencySolver_Maxwell3D
  {    
  public :
    TransparencySolver(EllipticProblem<HarmonicMaxwellEquation_3D>& var,
		       All_LinearSolver& solver);

  };


  //! class for transparent condition for 3-D Maxwell's equations (LDG formulation)
  template<>
  class TransparencySolver<HarmonicMaxwellEquation_3D_DG>
    : public TransparencySolver_Maxwell3D
  {    
  public :
    TransparencySolver(EllipticProblem<HarmonicMaxwellEquation_3D_DG>& var,
		       All_LinearSolver& solver);

  };


  //! static condensation for Maxwell's equation
  template<class T, class TypeEquation>
  class CondensationBlockSolver_Maxwell3D : public CondensationBlockSolver_Fem<T>
  {
    typedef typename TypeEquation::Complexe Complexe;

    EllipticProblem<TypeEquation>& var_problem;
  protected:
    Vector<int> num_u;
    Vector<R3> tau;
    Vector<char> coor_u;
    Vector<Complexe> mass_pml, ratio_mass;
    
  public:    
    CondensationBlockSolver_Maxwell3D(EllipticProblem<TypeEquation>& var);

    void Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix);
    void ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    void RecomposeSolution(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    
    size_t GetMemorySize() const;
    
  };


  template<class T>
  class CondensationBlockSolver<T, HarmonicMaxwellEquation_3D>
    : public CondensationBlockSolver_Maxwell3D<T, HarmonicMaxwellEquation_3D>
  {
  public:
    CondensationBlockSolver(EllipticProblem<HarmonicMaxwellEquation_3D>& var);
    
  };

  
  //! class for solving 3-D Maxwell's equations
  template<>
  class EllipticProblem<HarmonicMaxwellEquation_3D>
    : public HarmonicMaxwell_3D_Cplx<HarmonicMaxwellEquation_3D>
  {
  public:
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

    void ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;

    void ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;

    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			       Vector<VectReal_wp>& curl_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;

    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			       Vector<VectComplex_wp>& curl_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;
    
    void ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const;

    void ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const;
        
  };


  template<>
  class FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_3D>
    : public FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquation_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquation_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  //! class for solving 3-D Maxwell's equations (LDG formulation)
  template<>
  class EllipticProblem<HarmonicMaxwellEquation_3D_DG>
    : public HarmonicMaxwell_3D_Cplx<HarmonicMaxwellEquation_3D_DG>
  {
  public:
    IncidentWaveProjector<Complex_wp, Dimension3>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension3>& u_inc) const;

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row, int offset_col);
    
    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row, int offset_col);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;
        
    void ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;

    void ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			     int i, bool u_component, bool grad_component) const;

    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			       Vector<VectReal_wp>& grad_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;

    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			       Vector<VectComplex_wp>& grad_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;
    
    void ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const;

    void ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const;

  };


  //! base class for specifying hdg formulation for Maxwell's equations
  template<class T>
  class MaxwellEquationHdg_3D : public GenericEquation<T>
  {
  public:
    typedef Dimension3 Dimension;

    static const bool FirstOrderFormulation = true;
    static bool store_dfjm1;
    
    enum {nb_unknowns = 2, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 3, nb_components_hn = 3, nb_unknowns_vec = 1,
          type_element = 2};
    
    static inline bool ComputeDFjm1() { return store_dfjm1; }
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var_problem,
				  int iquad, const ElementReference_Dim<Dimension3>& Fb);
    
  };

  
  class HarmonicMaxwellEquationHdg_3D : public MaxwellEquationHdg_3D<Complex_wp>
  {
  };


  //! class for computing radar cross section for 3-D Maxwell's equations  
  template<>
  class VarComputationRCS<HarmonicMaxwellEquationHdg_3D>
    : public VarComputationRCS_Maxwell3D
  {
  public :
    template<class TypeEquation>
    VarComputationRCS(EllipticProblem<TypeEquation>& var);
    
  };


  //! class for transparent condition for 3-D Maxwell's equations (HDG formulation)
  template<>
  class TransparencySolver<HarmonicMaxwellEquationHdg_3D>
    : public TransparencySolver_Maxwell3D
  {    
  public :
    TransparencySolver(EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var,
		       All_LinearSolver& solver);

  };


  //! variables needed for condensed solver (a single element)
  template<class T>
  class CondensedBlockSolver_MaxwellHdg3D
  {
  public:
    bool flag;
    Matrix<T, Symmetric, RowSymPacked> invBtilde;
    Vector<TinyMatrix<T, Symmetric, 3, 3> > invBh;
    Matrix<T, General, ArrayRowSparse> SnD, Mint;
    
    size_t GetMemorySize() const;
    
  };

  
  //! static condensation for Maxwell's equation (HDG formulation)
  template<class T, class Complexe>
  class CondensationBlockSolver_MaxwellHdg3D : public CondensationBlockSolver_Fem<T>
  {
    VarProblem<Dimension3>& var_problem;
    
  protected:
    Vector<CondensedBlockSolver_MaxwellHdg3D<T> > optim_elt;

  public:
    template<class TypeEquation>
    CondensationBlockSolver_MaxwellHdg3D(EllipticProblem<TypeEquation>& var);

    inline Matrix<T, Symmetric, RowSymPacked>& GetBtildeMatrix(int i) { return optim_elt(i).invBtilde; }
    inline const Matrix<T, Symmetric, RowSymPacked>& GetBtildeMatrix(int i) const { return optim_elt(i).invBtilde; }
    
    inline Matrix<T, General, ArrayRowSparse>& GetSnD(int i) { return optim_elt(i).SnD; }
    inline const Matrix<T, General, ArrayRowSparse>& GetSnD(int i) const { return optim_elt(i).SnD; }
    
    inline Vector<TinyMatrix<T, Symmetric, 3, 3> >& GetInverseBh(int i) { return optim_elt(i).invBh; }
    inline const Vector<TinyMatrix<T, Symmetric, 3, 3> >& GetInverseBh(int i) const { return optim_elt(i).invBh; }
    
    inline Matrix<T, General, ArrayRowSparse>& GetMint(int i) { return optim_elt(i).Mint; }
    inline const Matrix<T, General, ArrayRowSparse>& GetMint(int i) const { return optim_elt(i).Mint; }
    
    void SetOptimizedCondensation(int i, bool optim);
    void ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    void RecomposeSolution(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    
    size_t GetMemorySize() const;
    
  };


  template<class T>
  class CondensationBlockSolver<T, HarmonicMaxwellEquationHdg_3D>
    : public CondensationBlockSolver_MaxwellHdg3D<T, Complex_wp>
  {
  public:
    CondensationBlockSolver(EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var);
    
  };
  

  //! class for solving 3-D Maxwell's equations (HDG formulation)
  template<>
  class EllipticProblem<HarmonicMaxwellEquationHdg_3D>
    : public HarmonicMaxwell_3D_Cplx<HarmonicMaxwellEquationHdg_3D>
  {
  public:
    EllipticProblem();
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
  };

  
  template<>
  class FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquationHdg_3D>
    : public FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquationHdg_3D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

  
} // end namespace

#define MONTJOIE_FILE_HARMONIC_MAXWELL_3D_HXX
#endif
