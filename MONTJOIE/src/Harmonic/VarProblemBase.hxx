#ifndef MONTJOIE_FILE_VAR_PROBLEM_BASE_HXX

namespace Montjoie
{
  
  //! base class for all classes of EllipticProblem
  class VarProblem_Base : public VarPhysicalProblem
  {
    friend class VarBoundaryCondition_Base;
    
  public :
    int print_level; //!< variable used to set the level display of the program
    
    // variables that are set directly in the constructor
    const int nb_unknowns_scal, nb_unknowns_vec, nb_unknowns;
    const int nb_components_en, nb_components_hn, nb_unknowns_hdg;
    const int type_element; const Vector<int> other_type_element;

    Vector<int> offset_dof_unknown; //!< offset for numbers of each unknown
    Vector<int> offset_dof_condensed; //!< offset after static condensation
    Vector<int> mesh_num_unknown; //!< which numbering to use for each unknown
    VectString name_other_elements;
    
    //! true if a discontinuous formulation is used
    int dg_formulation; bool sipg_formulation;
    bool compute_dfjm1;
    Vector<int> other_dg_formulation; // if some unknowns are discontinuous and other not
    
    //! penalization parameters for DG methods/ Nedelec's second family
    Real_wp alpha_penalization, delta_penalization;
    
    //! do we use upwind fluxes ?
    bool upwind_fluxes;
    
    //! if true the penalization term is automatically chosen
    bool automatic_choice_penalization;
    VectReal_wp Glob_CoefPenalDG;
    
    Vector<VectString> mesh_data; //!< parameters of the mesh
    string finite_element_name;
    
    //! if a face doesn't have a boundary condition, program is stopped ?
    static bool exit_if_no_boundary_condition;
    
  protected :
    Vector<MeshNumbering_Base<Real_wp>* > all_mesh_num;
    
    //! boundaries of the computational domain
    Real_wp xmin, xmax, ymin, ymax, zmin, zmax;
            
    //! total number of dofs of the problem
    int nodl;
    int default_order; //!< for output on meshes or transparent condition

    int order_over_integration;
    bool dg_exact_integration; //!< exact integration (Gauss points) ?
    bool mixed_formulation; //!< use of mixed formulation ?
    const bool first_order_formulation; //!< only first-order derivatives
    
    //! omega is the pulsation
    Real_wp omega;
    //!  omega2 the square of omega
    Real_wp omega2;
    Real_wp frequency; //!<  frequency = omega/ (2*pi)

    Real_wp wavelength_adim;

  public :    
    //! timer used
    MontjoieTimer var_chrono;

  private:
    
    void InitDefaultValues();

  public :
    template<class TypeEquation>
    VarProblem_Base(EllipticProblem<TypeEquation>&);
    
    virtual ~VarProblem_Base();
    
    // Inline methods    
    int GetNbDof() const;
    void SetNbDof(int);

    int GetNbMeshDof(int n = 0) const;
    int GetNbMeshNumberings() const;
    const MeshNumbering_Base<Real_wp>& GetMeshNumberingBase(int n = 0) const;
    
    int GetOffsetDofUnknown(int) const;
    int GetOffsetDofCondensed(int) const;
    int GetNbDofUnknown(int) const;

    int GetDefaultOrder() const;
    
    int FormulationDG() const;
    bool ComputeDFjm1() const;
    bool FirstOrderFormulation() const;
    bool FirstOrderFormulationDG() const;
    virtual void SetFirstOrderFormulation(bool);
    bool UseExactIntegrationElement() const;
    int GetOverIntegration() const;
        
    const Real_wp& GetXmin() const;
    const Real_wp& GetXmax() const;
    const Real_wp& GetYmin() const;
    const Real_wp& GetYmax() const;
    const Real_wp& GetZmin() const;
    const Real_wp& GetZmax() const;
    void SetComputationalDomain(const Real_wp&, const Real_wp&, const Real_wp&,
				const Real_wp&, const Real_wp&, const Real_wp&);

    const Real_wp& GetSquareOmega() const;
    const Real_wp& GetOmega() const;
    void GetMiomega(Real_wp& ) const;
    void GetMiomega(Complex_wp& ) const;
    void GetMomega2(Real_wp& ) const;
    void GetMomega2(Complex_wp& ) const;
    const Real_wp& GetFrequency() const;
    void SetOmega(const Real_wp&);
    void SetFrequency(const Real_wp&);

    const Real_wp& GetWaveLengthAdim() const;

    // Convenient methods

    void SetInputData(const string& description_field, const VectString& parameters);

    virtual void GetMemoryUsed(map<string, size_t>& var) const;
        
    // Virtual methods
    virtual int GetDimension() const = 0;
    virtual int GetNbComponentsUnknown(int n) const = 0;
    virtual int GetNbComponentsGradient(int n) const = 0;
    
    virtual int GetNbLocalDof(int i, int nm = 0) const = 0;
    virtual int GetNbSurfaceDof(int i, int nm = 0) const = 0;
    virtual int GetNbDofBoundaries(int i, int nm = 0) const = 0;
    virtual int GetNbPointsQuadratureInside(int) const = 0;
    virtual const VectReal_wp& WeightsND(int) const = 0;
    virtual bool ElementInsidePML(int) const = 0;
    virtual const ElementReference_Base& GetReferenceElementBase(int i, int n = 0) const = 0;
    virtual const ElementReference_Base& GetSurfaceElementBase(int i, int n = 0) const = 0;
    
    virtual void WriteMesh(const string&) = 0;
    virtual void SetSameNumberPeriodicDofs() = 0;

    virtual void InitIndices(int n) = 0;
    virtual int GetNbPhysicalIndices() const = 0;
    virtual void SetIndices(int ref, const Vector<string>&) = 0;
    virtual void SetPhysicalIndex(const string&, int, const Vector<string>&) = 0;
    virtual bool IsVaryingMedia(int ref) const = 0;
    virtual bool IsVaryingMedia(int m, int ref) const = 0;
    virtual Real_wp GetCoefficientPenaltyStiffness(int ref) const = 0;
    virtual string GetPhysicalIndexName(int m) const = 0;
    
    template<class Dimension>
    void GetVelocityOnElements(VectReal_wp& velocity, const Mesh<Dimension>&);
    
    virtual Real_wp GetVelocityOfMedia(int ref) const;
    virtual Real_wp GetVelocityOfInfinity() const;

  protected:
    virtual void InitPolarization() = 0;
    
    virtual void UpdateWaveVector() = 0;
    virtual void FillWaveVectorComponents(Real_wp&, Real_wp&, Real_wp&) = 0;

    virtual void ComputeGeometryQuantity(IVect&) = 0;
    virtual void ComputeArraySpecificEquation() = 0;
    virtual void ComputeTauCoefficient();

  public:
    virtual void CopyInputData(const VarProblem_Base& var);
        
    virtual bool IsSymmetricProblem(bool eigen = false) const = 0;
    virtual bool IsSymmetricMassMatrix() const = 0;
    virtual bool IsComplexProblem() const = 0;

    virtual void ComputeMeshAndFiniteElement(const string& name_elt, bool split = true) = 0;
    virtual void PerformOtherInitializations() = 0;

    virtual void SetTypeEquation(const string&);
        
  };
  
  
  //! base class used by AssembleMassMatrix
  class VarComputationProblem_Base
  {
  protected:
    //! threshold used to drop values in finite element matrix
    Real_wp threshold_matrix;
    
  public :
    VarComputationProblem_Base();
    virtual ~VarComputationProblem_Base();
    
    // Inline methods
    Real_wp GetThresholdMatrix() const;
    void SetThresholdMatrix(const Real_wp&);
    
    // Virtual methods
    virtual int GetNbElt() const = 0;
    virtual int GetNbRows() const = 0;
    virtual int GetPrintLevel() const = 0;
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&, CondensationBlockSolver_Base<Real_wp>&,
					 const GlobalGenericMatrix<Real_wp>&) = 0;
    
    virtual void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&, CondensationBlockSolver_Base<Complex_wp>&,
					 const GlobalGenericMatrix<Complex_wp>&) = 0;
    
    virtual void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
					 Vector<int>& intern_node) const;

    virtual CondensationBlockSolver_Base<Real_wp>* GetNewCondensationSolver(const Real_wp&);
    virtual CondensationBlockSolver_Base<Complex_wp>* GetNewCondensationSolver(const Complex_wp&);
    
  };
  
  
  //! class handling computation of finite element 
  class VarComputationProblem : public VarComputationProblem_Base
  {
  private:
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    DistributedProblem_Base& var_comm;
    
    void InitDefaultValues();
    
  protected :
    //! tells if the user wants the matrix to be explicitly computed or not
    /*!
      If equal to MATRIX_STORED, the matrix will be stored as a sparse matrix,
      the "fast" matrix-vector product will not be used.
    */
    int storage_finite_element_matrix;

    //! elementary matrix is symmetric ?
    bool symmetric_elementary_matrix;
    //! global matrix is symmetric ?
    bool symmetric_global_matrix;
    //! sparse elementary matrix
    bool sparse_elementary_matrix;

    //! static condensation prescribed by the leaf class of CondensationBlockSolver
    mutable bool leaf_static_condensation;
    bool light_static_condensation;
    
    bool use_symmetrization_when_possible; //!< symmetrization of matrix when possible
    bool erase_dirichlet_columns; //!< columns associated with Dirichlet dofs are erased
    
  public :
    //! Types of storage for the matrix
    /*!
      MATRIX_AUTO_STORAGE : the program decides if the matrix is stored or not
      MATRIX_STORED : the matrix is effectively stored in a Seldon format
      MATRIX_FREE : the matrix is not stored
    */
    enum {MATRIX_AUTO_STORAGE, MATRIX_STORED, MATRIX_FREE};
        
    template<class TypeEquation>
    VarComputationProblem(EllipticProblem<TypeEquation>&);

    // Inline methods
    bool UseMatrixFreeAlgorithm() const;    
    bool IsSymmetricGlobalMatrix() const;
    int GetStorageFiniteElementMatrix() const;
    void SetStorageFiniteElementMatrix(int type);
    void SetSymmetricElementaryMatrix(bool sym = true);
    void SetLeafStaticCondensation(bool leaf = true) const;
    bool GetLeafStaticCondensation() const;
    bool LightStaticCondensation() const;
    bool GetSymmetrizationUse() const;
    void SetSymmetrizationUse(bool sym = true);
    void SetHomogeneousDirichlet(bool hg_dir);
    bool IsHomogeneousDirichlet() const;
    
    inline VarBoundaryCondition_Base& GetBoundaryConditionProblem() { return var_boundary; }
    inline DistributedProblem_Base& GetDistributedProblem() { return var_comm; }
    inline VarProblem_Base& GetVarProblemBase() { return var_problem; }
    inline const VarProblem_Base& GetVarProblemBase() const { return var_problem; }
    
    int GetNbRows() const;
    int GetPrintLevel() const;
    virtual void SetPrintLevel(int lvl); 
    virtual size_t GetMemorySize() const;
    virtual void ComputeMassMatrix(bool compute_rho = true, bool delete_points = true) = 0;

    // convenient methods
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
				 Vector<int>& intern_node) const;
    
    // Methods related to the computation of finite element matrices        
    template<class T, class MatrixSparse>
    void AddMatrixWithBC(FemMatrixFreeClass_Base<T>&,
			 MatrixSparse& mat_direct, const GlobalGenericMatrix<T>& nat_mat,
			 int offset_row = 0, int offset_col = 0,
			 CondensationBlockSolver_Fem<T>* solver = NULL,
                         bool diag_matrix = false);
    
    template<class T, class Prop, class Storage>
    void AddMatrixWithBC(DistributedMatrix<T, Prop, Storage>& mat_direct,
                         const GlobalGenericMatrix<T>& nat_mat,
                         int offset_row = 0, int offset_col = 0,
			 CondensationBlockSolver_Fem<T>* solver = NULL,
                         bool diag_matrix = false);
    
    template<class T, class Prop>
    void AddMatrixWithBC(Matrix<T, Prop, DiagonalRow>& mat_direct,
                         const GlobalGenericMatrix<T>& nat_mat,
                         int offset_row = 0, int offset_col = 0,
			 CondensationBlockSolver_Fem<T>* solver = NULL,
                         bool diag_matrix = false);
    
    template<class T>
    void AddMatrixWithBC(FemMatrixFreeClass_Base<T>&,
			 const GlobalGenericMatrix<T>& nat_mat,
                         int offset_row = 0, int offset_col = 0,
			 CondensationBlockSolver_Fem<T>* solver = NULL,
                         bool diag_matrix = false);

    template<class T>
    void AddMatrixFEM(FemMatrixFreeClass_Base<T>& A,
		      const GlobalGenericMatrix<T>& nat_mat,
		      int offset_row = 0, int offset_col = 0,
		      CondensationBlockSolver_Fem<T>* solver = NULL,
                      bool diag_matrix = false);
    
    template<class T, class MatrixSparse>
    void AddMatrixFEM(FemMatrixFreeClass_Base<T>& mat_iterative,
		      MatrixSparse& mat_direct, const GlobalGenericMatrix<T>& nat_mat,
		      int offset_row = 0, int offset_col = 0,
		      CondensationBlockSolver_Fem<T>* solver = NULL,
                      bool diag_matrix = false);

    template<class T, class Prop, class Storage>
    void AddMatrixFEM(DistributedMatrix<T, Prop, Storage>& mat_direct,
                      const GlobalGenericMatrix<T>& nat_mat,
                      int offset_row = 0, int offset_col = 0,
		      CondensationBlockSolver_Fem<T>* solver = NULL,
                      bool diag_matrix = false);

    template<class T, class Prop>
    void AddMatrixFEM(Matrix<T, Prop, DiagonalRow>& mat_direct,
                      const GlobalGenericMatrix<T>& nat_mat,
                      int offset_row = 0, int offset_col = 0,
		      CondensationBlockSolver_Fem<T>* solver = NULL,
                      bool diag_matrix = false);
    
    // computation of the diagonal
    template<class T>
    void ComputeDiagonalMatrix(Vector<T>& diagonal,
			       const GlobalGenericMatrix<T>& nat_mat, bool assemble = true);
    
    template<class T>
    void ComputeDiagonalMatrix(Vector<T>& diagonal,
                               const FemMatrixFreeClass_Base<T>& mat,
			       const GlobalGenericMatrix<T>& nat_mat, bool assemble = true);
    
    template<class T, class Prop, class Storage, class Allocator>
    void ComputeDiagonalMatrix(Vector<T>& diagonal,
                               const Matrix<T, Prop, Storage, Allocator>& mat,
			       const GlobalGenericMatrix<T>& nat_mat);
    
    // computation of block-diagonal matrix
    void GetProfilBlockDiagonal(int nblock, Vector<IVect>& num_ddl_blocks);
    
    template<class MatrixDiag, class T>
    void ComputeBlockDiagonalMatrix(MatrixDiag& diagonal, int s,
                                    const FemMatrixFreeClass_Base<T>& mat,
				    const GlobalGenericMatrix<T>& nat_mat);
    
    template<class Matrix1, class MatrixDiag, class T>
    void ComputeBlockDiagonalMatrix(MatrixDiag& diagonal, int s, const Matrix1& mat,
				    const GlobalGenericMatrix<T>& nat_mat);

    void CopyInputData(const VarComputationProblem& var);

    // virtual methods
    virtual bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    virtual bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    virtual bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    virtual bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    virtual bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    virtual bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    
    virtual void GetStaticCondensedRows(IVect& IndexCondensedRows, IVect&, IVect&, IVect&,
					int& nb_scalar_dof, int&,
					IVect&, Vector<IVect>&) const = 0;
    
    virtual void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
				       const GlobalGenericMatrix<Real_wp>& nat_mat,
				       int offset_row = 0, int offset_col = 0);

    virtual void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
				       const GlobalGenericMatrix<Complex_wp>& nat_mat,
				       int offset_row = 0, int offset_col = 0);

    virtual void UpdateShiftAdimensionalization(Real_wp&, Real_wp&);
    virtual void UpdateShiftAdimensionalization(Complex_wp&, Complex_wp&);
    
    virtual void UpdateEigenvaluesAdimensionalization(VectReal_wp&, VectReal_wp&,
						      Matrix<Real_wp, General, ColMajor>&);

    virtual void UpdateEigenvaluesAdimensionalization(VectComplex_wp&, VectComplex_wp&,
						      Matrix<Complex_wp, General, ColMajor>&);

    virtual void FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const;
    
    // other methods
    void ModifySourceSymmetry(VectReal_wp&) const;
    void ModifySourceSymmetry(VectComplex_wp&) const;
    void ModifySourceSymmetry(Matrix<Real_wp, General, ColMajor>&) const;
    void ModifySourceSymmetry(Matrix<Complex_wp, General, ColMajor>&) const;
    
    virtual FemMatrixFreeClass_Base<Real_wp>* GetNewIterativeMatrix(const Real_wp&) const = 0;
    virtual FemMatrixFreeClass_Base<Complex_wp>* GetNewIterativeMatrix(const Complex_wp&) const = 0;
    virtual All_LinearSolver* GetNewLinearSolver() = 0;
    virtual int GetMassMatrixType(Vector<bool>&) const = 0;
    
    virtual All_Preconditioner_Base<Real_wp>* GetNewPreconditioning(const Real_wp&) = 0;
    virtual All_Preconditioner_Base<Complex_wp>* GetNewPreconditioning(const Complex_wp&) = 0;
    
  };
  
}

#define MONTJOIE_FILE_VAR_PROBLEM_BASE_HXX  
#endif
