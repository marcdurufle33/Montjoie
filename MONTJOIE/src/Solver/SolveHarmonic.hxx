#ifndef MONTJOIE_FILE_SOLVE_HARMONIC_HXX

namespace Montjoie
{

  class PreconditioningEnum
  {
  public :
    // type of preconditioners in Montjoie
    enum {IDENTITY, SSOR, DIRECT, JACOBI, BLOCK_JACOBI,
	  LOW_ORDER, MULTIGRID, SUBDOMAIN};
    
  };

  
  template<class T>
  class CondensationBlockSolver_Fem;
  
  //! general class for solution of finite element linear systems
  class All_LinearSolver : public InputDataProblem_Base, public Preconditioner_Base<Complex_wp>
  {
  protected :
    //! iterative or direct resolution ?
    bool iterative_resolution;
    //! choice of solver
    int type_choice_solver_auto;
    //! estimation of the condition number of the matrix ?
    bool estimation_conditioning_number_asked;
    //! number of iterations to know condition number
    /*!
      the used method is a simple iterative power of the matrix
      the initial vector is filled randomly
    */
    int nb_iterations_powers_iterative;
    
    //! direct solver for real numbers
    All_MatrixLU<Real_wp> mat_lu_real;

    //! direct solver for complex numbers
    All_MatrixLU<Complex_wp> mat_lu_cplx;

    // sparse matrices
    DistributedMatrix<Real_wp, General, ArrayRowSparse> mat_unsym_real;
    DistributedMatrix<Complex_wp, General, ArrayRowSparse> mat_unsym_cplx;

    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> mat_sym_real;
    DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> mat_sym_cplx;

    GlobalGenericMatrix<Real_wp> nature_matrix_real;
    GlobalGenericMatrix<Complex_wp> nature_matrix_cplx;
    
    //! iterative matrix for real numbers
    FemMatrixFreeClass_Base<Real_wp>* mat_iterative_real;
    //! iterative matrix for complex numbers
    FemMatrixFreeClass_Base<Complex_wp>* mat_iterative_cplx;

    //! preconditioning for real numbers
    All_Preconditioner_Base<Real_wp>* precond_real;
    //! preconditioning for complex numbers
    All_Preconditioner_Base<Complex_wp>* precond_cplx;

    //! iterative solver with preconditioning for real numbers
    All_IterativeSolver<Real_wp> iterative_solver_real;
    
    //! iterative solver with preconditioning for complex numbers
    All_IterativeSolver<Complex_wp> iterative_solver_cplx;

    //! using csr sparse matrices (to have fast matrix vector product)
    bool use_csr_matrices;
    
    //! attributes for static condensation 
    bool static_condensation;
    bool compress_matrix;
    CondensationBlockSolver_Fem<Real_wp>* solver_real;
    CondensationBlockSolver_Fem<Complex_wp>* solver_cplx;
    
    // scaling of original matrix ?
    bool scaling_matrix;
    
    //! matrix is symmetric ?
    bool sym_matrix;
    
    //! matrix is complex ?
    bool complex_matrix;

    //! use of real preconditioning for a complex matrix ?
    bool use_real_preconditioning;
    //! if true the present preconditioning is not deleted
    bool shallow_copy_real_precond;
    
    // informations to reuse the old pattern for a new factorization
    bool use_old_pattern_Ah, new_pattern_Ah;
    int nnz_old_Ah;
    bool facto_completed;
    bool use_symmetrization;

    //! true if the user wants to store the matrix on the disk
    bool storage_matrix_asked;    
    //! the name of the file where the matrix will be stored
    string file_name_matrix_stored;
    
  private :
    //! object containing datas
    VarProblem_Base& var_problem;
    VarSourceProblem_Base& var_source;
    VarBoundaryCondition_Base& var_boundary;
    VarComputationProblem& var_computation;
    DistributedProblem_Base& var_comm;
    
    void InitDefaultValues();
    
  public :
    bool add_diagonal_part;
    Vector<Real_wp> diagonal_real_part_to_add;
    Vector<Complex_wp> diagonal_complex_part_to_add;
    
    // choice of solver
    enum {AUTOMATIC_SOLVER, DIRECT_SOLVER, ITERATIVE_SOLVER,
	  MULTIGRID_SOLVER, ILUT_SOLVER, USER_CHOICE};
    
    // constructor/destructor
    template<class TypeEquation>
    All_LinearSolver(EllipticProblem<TypeEquation>&);
    
    virtual ~All_LinearSolver();
    
    
    // Inline methods
    int GetTypeSolver() const;
    bool IsIterativeSolver() const;
    
    void SetDirectSolver();
    void SetDirectSolver(const VectString& param);
    void SetIterativeSolver(int type);
    void SetIterativeSolver(int type, int type_precond);
    bool StaticCondensation() const;
    void EnableRealPreconditioning(bool flag = true);
    void SetRealPreconditioning(All_Preconditioner_Base<Real_wp>* prec); 
    
    int GetMaxNumberIteration() const;
    int GetRestart() const;
    Real_wp GetStoppingCriterion() const;
    
    void SetPrintLevel(int p);
    void UseOldPattern(bool use_p = true);
    bool UseCompressedMatrix() const;
    void DisableLocalPreconditioning();
    
    FemMatrixFreeClass_Base<Real_wp>& GetIterativeMatrix(const Real_wp&);
    FemMatrixFreeClass_Base<Complex_wp>& GetIterativeMatrix(const Complex_wp&);

    All_Preconditioner_Base<Real_wp>& GetPreconditioning(const Real_wp&);
    All_Preconditioner_Base<Complex_wp>& GetPreconditioning(const Complex_wp&);

    All_MatrixLU<Real_wp>& GetDirectFactorization(const Real_wp&);
    All_MatrixLU<Complex_wp>& GetDirectFactorization(const Complex_wp&);
    
    All_IterativeSolver<Real_wp>& GetIterativeSolver(const Real_wp&);
    All_IterativeSolver<Complex_wp>& GetIterativeSolver(const Complex_wp&);

    CondensationBlockSolver_Fem<Real_wp>& GetCondensedSolver(const Real_wp&);
    CondensationBlockSolver_Fem<Complex_wp>& GetCondensedSolver(const Complex_wp&);
    
    DistributedMatrix<Real_wp, General, ArrayRowSparse>& GetDirectUnsymmetricMatrix(const Real_wp&);
    DistributedMatrix<Complex_wp, General, ArrayRowSparse>& GetDirectUnsymmetricMatrix(const Complex_wp&);
    
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& GetDirectSymmetricMatrix(const Real_wp&);
    DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>& GetDirectSymmetricMatrix(const Complex_wp&);

    DistributedMatrix<Real_wp, General, RowSparse>& GetDirectUnsymmetricMatrixCSR(const Real_wp&);
    DistributedMatrix<Complex_wp, General, RowSparse>& GetDirectUnsymmetricMatrixCSR(const Complex_wp&);
    
    DistributedMatrix<Real_wp, Symmetric, RowSymSparse>& GetDirectSymmetricMatrixCSR(const Real_wp&);
    DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>& GetDirectSymmetricMatrixCSR(const Complex_wp&);
    
    GlobalGenericMatrix<Real_wp>& GetNatureMatrix(const Real_wp&);
    GlobalGenericMatrix<Complex_wp>& GetNatureMatrix(const Complex_wp&);

    inline bool FactorizationCompleted() const { return facto_completed; }
    
    // selection of the best linear solver
    template<class Dimension, class Prop>
    void SelectOptimalLinearSolver(int order, int N, Dimension dim, Prop property);

    void SetInputData(const string& description_field, const VectString& parameters);
    void SetInputPreconditioning(const string& keyword, const VectString& param);

    void EnableStaticCondensation(bool flag, bool store, bool compress);
    void EnableMatrixStorage(bool flag, const string& file_name);
    
    void CopyParameter(const All_LinearSolver& solver);
    void CopyParameter(const Preconditioner_Base<Complex_wp>&);
    
    // Memory management
    void GetMemoryUsed(map<string, size_t>& var, bool add = false) const;
    size_t GetMemorySize() const;

    void ClearFactorization();
    void Clear();
    
    // General methods    
    template<class T>
    void PerformFactorizationStep(const GlobalGenericMatrix<T>&,
				  bool precond = true);

  protected:
    template<class T>
    void SolveIterativeSystem(Vector<T>& source_rhs, Vector<T>& x_sol);

  public:
    void Solve(const VirtualMatrix<Complex_wp>&, const Vector<Complex_wp>& r, Vector<Complex_wp>& z);
    void TransSolve(const VirtualMatrix<Complex_wp>&, const Vector<Complex_wp>& r, Vector<Complex_wp>& z);

    template<class T>
    void ComputeSolution(Vector<T>& x_sol,
                         bool dirichlet_inhg = true, bool assemble = true);
    
    template<class T>
    void ComputeSolution(Vector<T>& x_sol,
                         const GlobalGenericMatrix<T>& nat_mat,
			 bool dirichlet_inhg = true, bool assemble = true);
    
    template<class T>
    void ComputeSolution(const SeldonTranspose&, Vector<T>& x_sol,
                         const GlobalGenericMatrix<T>& nat_mat,
			 bool dirichlet_inhg = true, bool assemble = true);
    
    template<class T>
    void ComputeSolution(Matrix<T, General, ColMajor>& x_sol,
                         const GlobalGenericMatrix<T>& nat_mat,
			 bool dirichlet_inhg = true, bool assemble = true);

    template<class T>
    void ComputeSolution(const SeldonTranspose& trans, Matrix<T, General, ColMajor>& x_sol,
                         const GlobalGenericMatrix<T>& nat_mat, bool dirichlet_inhg = true,
			 bool assemble = true);
    
    template<class T, class Prop, class Storage>
    void ComputeMatrixAndFactorize(DistributedMatrix<T, Prop, Storage>& Ah,
                                   All_MatrixLU<T>& matLU,
                                   const GlobalGenericMatrix<T>& nature_matrix);
    
    template<class T>
    void ComputeMatrixAndPreconditioning(FemMatrixFreeClass_Base<T>& Ah, 
                                         All_Preconditioner_Base<T>& prec,
					 const GlobalGenericMatrix<T>&,
                                         bool precond = true);
    
  };


  //! for each element, 2x2 decomposition of the finite element matrix
  template<class T, class Prop, class Storage>
  class CondensedBlockSchur
  {
  public:
    IVect num_ddl;
    Matrix<T, Prop, Storage> inv_a22;
    Matrix<T> a12, a21;

    size_t GetMemorySize() const;
  };


  //! object handling static condensation for finite element method
  template<class T>
  class CondensationBlockSolver_Fem : public CondensationBlockSolver_Base<T>
  {
    template<class T0> friend class CondensationBlockSolver_Fem;
    
  protected:
    VarComputationProblem& var_computation;
    DistributedProblem_Base& var_comm;
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    
    IVect* IndexCondensedRows, *GlobalCondensedRows, *OverlappedRows, *OverlappedProcs;
    IVect* SharingProcs; Vector<IVect>* SharingRows;
    int size_compressed_system;
    bool store_block;
    int nb_scalar_dof, nb_global_dof;
    bool optim_condensation, use_symmetrization, compress_system;
    bool symmetric_elem_matrix;
    
    IVect IndexCondensedRows_store, GlobalCondensedRows_store, OverlappedRows_store, OverlappedProcs_store;
    IVect SharingProcs_store; Vector<IVect> SharingRows_store;
    
    mutable Vector<CondensedBlockSchur<T, General, RowMajor> > block_unsym;
    mutable Vector<CondensedBlockSchur<T, Symmetric, RowSymPacked> > block_sym;

    mutable MatrixVectorProductLevel level_volume;

    void InitDefaultValues();
    
  public:    
    template<class TypeEquation>
    CondensationBlockSolver_Fem(EllipticProblem<TypeEquation>& var);

    inline int GetNbRows() const { return size_compressed_system; }
    inline MatrixVectorProductLevel& GetVolumeLevel() { return level_volume; }
    inline const Vector<int>& GetIndexCondensedRow() const { return *IndexCondensedRows; }
    
    inline void SetOptimalCondensation(bool optim = true) { optim_condensation = optim; }
    inline void SetSymmetrization(bool sym = true) { use_symmetrization = sym; }
    inline bool UseSymmetrization() const { return use_symmetrization; }
    inline void SetCompressionSystem(bool comp = true) { compress_system = comp; }
    inline void SetSymmetryElementaryMatrix(bool sym = true) { symmetric_elem_matrix = sym; }
    
    void SetStoreBlock(bool store = true);

    void CopyParameter(const CondensationBlockSolver_Fem<T>& solver);
    
    template<class Prop, class Storage>
    void CompressMatrixGen(DistributedMatrix<T, Prop, Storage>& Ah);

    void CompressMatrixGen(FemMatrixFreeClass_Base<T>& Ah);

    template<class T0>
    void ExtractVector(const Vector<T0>& x_sol, Vector<T0>& xtmp) const;

    template<class T0>
    void ScatterVector(const Vector<T0>& xtmp, Vector<T0>& x_sol) const;

    template<class T0>
    void AddDomains(Vector<T0>& x, int nb_u = -1) const; 
    
#ifdef SELDON_WITH_MPI
    DistributedVector<T>* AllocateDistributedVector(Vector<T>& x);
#endif
    
  protected:
    // internal methods
    template<class Prop0, class Storage0, class Prop, class Storage>
    void GetSchurMatrixGen(int i, int nb_dof_loc, Matrix<T, Prop0, Storage0>& mat_interac,
			   Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a11,
			   Matrix<T, Prop, Storage>& a22, Vector<int>&);

    template<class Prop0, class Storage0, class Prop, class Storage>
    void ModifyRhsGen(const SeldonTranspose& trans, Vector<T>& b_source, int i, const IVect& num_ddl, Matrix<T, Prop0, Storage0>& mat_elt_unsym,
		      Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a22) const;

    template<class Prop0, class Storage0, class Prop, class Storage>
    void RecomposeSolGen(const SeldonTranspose& trans, Vector<T>& b_source, int i, const IVect& num_ddl, Matrix<T, Prop0, Storage0>& mat_elt_unsym,
			 Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a22) const;
    
  public:
    template<class T0>
    void CopyPtr(const CondensationBlockSolver_Fem<T0>& solver);
    
    virtual void Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix);
    virtual void ModifyElementaryMatrix(int i, IVect& num_ddl, VirtualMatrix<T>& mat_interac,
					const GlobalGenericMatrix<T>& nat_mat);
    
    virtual void ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;
    virtual void RecomposeSolution(const SeldonTranspose&, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const;

    virtual size_t GetMemorySize() const;
    
  };


  template<class T, class TypeEquation>
  class CondensationBlockSolver : public CondensationBlockSolver_Fem<T>
  {
  public:
    CondensationBlockSolver(EllipticProblem<TypeEquation>& var);
    
  };
  
}

#define MONTJOIE_FILE_SOLVE_HARMONIC_HXX
#endif
  
