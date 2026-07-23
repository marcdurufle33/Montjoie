#ifndef MONTJOIE_FILE_PRECONDITIONER_HXX

namespace Seldon
{
  
#ifdef SELDON_WITH_HYPRE
  template<>
  class HyprePreconditioner<Complex_wp> : public Preconditioner_Base<Complex_wp>
  {
  public:
    template<class T, class Prop, class Storage>
    void ConstructPreconditioner(DistributedMatrix<T, Prop, Storage>& A,
				 bool keep_matrix = false);
    
    void Solve(const VirtualMatrix<Complex_wp>& A, const Vector<Complex_wp>& b, Vector<Complex_wp>& x);
    void TransSolve(const VirtualMatrix<Complex_wp>& A, const Vector<Complex_wp>& b, Vector<Complex_wp>& x);
    
  };
#endif
  
}

namespace Montjoie
{  
  
  template<class T>
  class IdentityPreconditioner : public Preconditioner_Base<T>
  {
  public :
    IdentityPreconditioner();
    
    void SetInputData(const string& keyword, const VectString& param);
    void SetInputPreconditioning(const string& keyword, const VectString& param);
    
  };
  
#ifdef SELDON_WITH_PRECONDITIONING  
  template<class T>
  class JacobiPreconditioning : public Preconditioner_Base<T>
  {
  protected :
    //! Relaxation parameter
    Real_wp omega;
    //! number of iterations
    int nb_iterations;
    //! 1 for diagonal jacobi, higher for block-diagonal
    int size_block;
    //! diagonal for Jacobi algorithm
    Vector<T> diagonal;
    //! unsymmetric block-diagonal matrix
    DistributedMatrix<T, General, BlockDiagRow> block_diag_unsym;
    //! symmetric block-diagonal matrix
    DistributedMatrix<T, Symmetric, BlockDiagRowSym> block_diag_sym;
    
  public :
    
    JacobiPreconditioning();
    
    T GetRelaxationCoef() const;
    virtual size_t GetMemorySize() const;
    
    virtual void SetInputData(const string& keyword, const VectString& param);
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);
    
    template<class MatrixSparse>
    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
				 const MatrixSparse& mat, const GlobalGenericMatrix<T>& nat_mat,
				 CondensationBlockSolver_Fem<T>& cond_solver);

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x, bool init);

    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x, bool init);
    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x, bool init);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    void CopyParameter(const Preconditioner_Base<T>&);
    void Clear();
    
  };

  
  //! Base class for domain decomposition preconditioning
  template<class T>
  class SubdomainPreconditioning_Base : public Preconditioner_Base<T>
  {
  protected :

    //! parameters for direct solver
    VectString param_direct_solver;    
    T alpha_laplace;
    // operator P
    VectReal_wp OperatorP;
    // solver for the current subdomain
    All_MatrixLU<T> mat_lu;

  private:
    DistributedProblem_Base& var_subdomain;

    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    SubdomainPreconditioning_Base(EllipticProblem<TypeEquation>&);
    
    virtual ~SubdomainPreconditioning_Base();
    
    virtual size_t GetMemorySize() const;
    
    virtual void SetInputData(const string& keyword, const VectString& param);
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);
    void CopyParameter(const Preconditioner_Base<T>&);
    
    virtual void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver) = 0;

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);

    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    
  };
  

  template<class T, class Dimension>
  class SubdomainPreconditioning_Dim : public SubdomainPreconditioning_Base<T>
  {
  private:
    DistributedProblem<Dimension>& var_subdomain;
    VarComputationProblem& var_computation;
    
  public :
    template<class TypeEquation>
    SubdomainPreconditioning_Dim(EllipticProblem<TypeEquation>& var);

    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
				 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature,
				 CondensationBlockSolver_Fem<T>& cond_solver);
    
  };
  
  template<class T, class TypeEquation>
  class SubdomainPreconditioning
    : public SubdomainPreconditioning_Dim<T, typename TypeEquation::Dimension>
  {
  public :
    SubdomainPreconditioning(EllipticProblem<TypeEquation>& var);
    
  };


  //! class for preconditioners depending on equation
  template<class T>
  class LowOrderPreconditioning_Base : public Preconditioner_Base<T>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

  protected :
    
    VectString param_direct_solver;
    //! matching array between dofs of high-order mesh and dofs of low-order mesh
    IVect NumberDof_LowOrder;
    //! scaling between dofs of high-order mesh and dofs of low-order mesh
    VectReal_wp ScalingDof_LowOrder;
    //! solver used for low order problem
    All_MatrixLU<T> mat_lu;

  public :
    
    LowOrderPreconditioning_Base();

    virtual ~LowOrderPreconditioning_Base();
    
    virtual size_t GetMemorySize() const;
    
    virtual void SetInputData(const string& keyword, const VectString& param);
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);
    
    virtual void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver) = 0;

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);

    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    virtual void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x);
    virtual void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x);

  };
  
  template<class T, class Dimension>
  class LowOrderPreconditioning_Dim : public LowOrderPreconditioning_Base<T>
  {
  private:
    DistributedProblem<Dimension>& var_problem;
    
  public :
    template<class TypeEquation>
    LowOrderPreconditioning_Dim(EllipticProblem<TypeEquation>& var);

    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
				 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature,
				 CondensationBlockSolver_Fem<T>& cond_solver);

    virtual void ComputeLocalProlongation(const VectReal_wp&,
					  DistributedProblem<Dimension>& var,
					  DistributedProblem<Dimension>& var_low);
    
  };

  template<class T, class TypeEquation>
  class LowOrderPreconditioning
    : public LowOrderPreconditioning_Dim<T, typename TypeEquation::Dimension>
  {
  public :
    LowOrderPreconditioning(EllipticProblem<TypeEquation>&);
    
  };
  
  
  //! Base class for multigrid preconditioning
  template<class T>
  class MultigridPreconditioning_Base : public Preconditioner_Base<T>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

  protected:
    //! damping parameter for preconditioner
    /*!
      preconditioning based on approximate inverse of
      \f$ -k^2(delta + i beta) - \Delta \f$
      alpha is equal to delta + i beta
    */
    T alpha_laplace;
    //! relaxation in multigrid parameters
    /*!
      for each coarse grid, alpha is multiplied by this coefficient
      so that we can adapt the wavelength to the refinement of the mesh
    */
    T relaxation_laplace;
    
    //! matrix of the restriction operator
    /*!
      the prolongation operator is the transposed of the restriction
    */
    Matrix<Real_wp, General, ArrayRowSparse> Vrestriction;
    //! scaling of restriction operator
    VectReal_wp InvDh;
    //! restriction operator on a reference element
    FiniteElementInterpolator proj_coarse_to_fine;
    //! restriction operator for surface element
    FiniteElementInterpolator surface_coarse_to_fine;
    // direct coarse solver
    All_LinearSolver* coarse_solver;
    //! pointer for the solution on the fine mesh
    All_LinearSolver* fine_solver;
    GlobalGenericMatrix<T> nat_mat_coarse;
    //! parameters for direct solver
    VectString param_direct_solver;
    //! number of steps for the post and pre-smoother
    /*!
      The post and pre-smooter have the same number of steps
      so that the preconditioner is symmetric
    */
    int nb_smoothing_iterations;
    //! \f$ \gamma \f$ of the multigrid
    int nb_subcycles_multigrid;
    bool use_f_cycle;
    int minimal_order_multigrid; //!< order of coarsest grid
    int maximal_order_multigrid;
    bool no_damping_on_first_smoother;
    int nb_iterations_multigrid;
    int type_smoother; //!< smoother used for multigrid iteration
    //! SSOR smoothing
    SorPreconditioner<T> smoother_ssor;
    //! used matrix for SSOR smoothing
    typename GenericStorage<T>::SparseSymMatrix* mat_ssor_sym, mat_stored_sym;
    typename GenericStorage<T>::SparseSymMatrixCSR* matCSR_ssor_sym, matCSR_stored_sym;
    typename GenericStorage<T>::SparseMatrix* mat_ssor_unsym, mat_stored_unsym;
    typename GenericStorage<T>::SparseMatrixCSR* matCSR_ssor_unsym, matCSR_stored_unsym;
    //! Jacobi smoothing
    JacobiPreconditioning<T> smoother_jacobi;
    //! used matrix for Jacobi/SSOR smoothing
    FemMatrixFreeClass_Base<T> *mat_smoothing, *mat_stored_smoothing;
    CondensationBlockSolver_Fem<T>* cond_solver;
    int nb_dof_coarse;
    
  private:
    
    //! fine problem to solve
    VarComputationProblem& var_fine;
    //! coarse problem (order is divided by 2)
    VarComputationProblem* var_coarse;

    void InitDefaultValues();
       
  public :
    // type of smoothers
    enum {SSOR, JACOBI, BLOCK_JACOBI};
    
    template<class TypeEquation>
    MultigridPreconditioning_Base(EllipticProblem<TypeEquation>& var);
    
    virtual ~MultigridPreconditioning_Base();

    virtual size_t GetMemorySize() const;
    
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);
    virtual void SetInputData(const string& keyword, const VectString& param);

    void SetDampingCoefficient(const T& alpha);

    void FinalizePreconditioner(VarComputationProblem* var, int rf, int rc,
				const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nature);
    
    virtual void ConstructPreconditioner(VarComputationProblem&, All_LinearSolver&,
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver) = 0;
    
    virtual void ConstructSmoother(const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nat_mat);
    
    void CopyParameter(const Preconditioner_Base<T>& fine_mg);
    void SetMaximalOrder(int r);

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x, int m = 0, int n = 0);

    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    template<class T0>
    void MultiGrid_Vcycle(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			  Vector<T0>& r_coarse, Vector<T0>& r, bool init_guess_null, int n);

    template<class T0>
    void ApplyPreSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			     Vector<T0>& r, bool init_guess_null);

    template<class T0>
    void ApplyPostSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			      Vector<T0>& r, bool init_guess_null);
    
    virtual void ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
				  Vector<Treal>& r, bool init_guess_null);

    virtual void ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
				  Vector<Tcplx>& r, bool init_guess_null);

    virtual void ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
				   Vector<Treal>& r, bool init_guess_null);
    
    virtual void ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
				   Vector<Tcplx>& r, bool init_guess_null);

    virtual void GetRestriction(const Vector<Treal>& r, Vector<Treal>& r_coarse) const = 0;
    virtual void GetRestriction(const Vector<Tcplx>& r, Vector<Tcplx>& r_coarse) const = 0;
    
    virtual void GetProlongation(const Vector<Treal>& r_coarse, Vector<Treal>& r) const = 0;
    virtual void GetProlongation(const Vector<Tcplx>& r_coarse, Vector<Tcplx>& r) const = 0;
    
    virtual int GetCoarseOrder() const = 0;
    
  };


  template<class T, class Dimension>
  class MultigridPreconditioning_Dim
    : public MultigridPreconditioning_Base<T>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

  protected:
    void ComputeProlongationElement(int i, int rf, int rc, Matrix<Real_wp>& LocalProlongation);

  private:
    DistributedProblem<Dimension>& var_fine;
    DistributedProblem<Dimension>* var_coarse;
    
  public :
    template<class TypeEquation>
    MultigridPreconditioning_Dim(EllipticProblem<TypeEquation>& var);
    
    virtual void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver&, 
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver);

    template<class T0>
    void GetRestrictionGen(const Vector<T0>& r, Vector<T0>& r_coarse) const;

    template<class T0>
    void GetProlongationGen(const Vector<T0>& r_coarse, Vector<T0>& r) const;

    void GetRestriction(const Vector<Treal>& r, Vector<Treal>& r_coarse) const;
    void GetRestriction(const Vector<Tcplx>& r, Vector<Tcplx>& r_coarse) const;
    
    void GetProlongation(const Vector<Treal>& r_coarse, Vector<Treal>& r) const;
    void GetProlongation(const Vector<Tcplx>& r_coarse, Vector<Tcplx>& r) const;

    int GetCoarseOrder() const;
    
  };

  
  template<class T, class TypeEquation>
  class MultigridPreconditioning
    : public MultigridPreconditioning_Dim<T, typename TypeEquation::Dimension>
  {
  public :
    MultigridPreconditioning(EllipticProblem<TypeEquation>& var);
    
  };


  template<class T>
  class LocalPreconditioning_Base;
  
  //! Preconditioners in Montjoie
  template<class T>
  class All_Preconditioner_Base : public Preconditioner_Base<T>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    
  protected :
    //! type of preconditioner
    int type_precond;
    //! LU matrix if direct solver is used as preconditioner
    All_MatrixLU<T> mat_lu;
    Real_wp threshold;
    
    //! SSOR preconditioning
    SorPreconditioner<T> prec_ssor;
    //! used matrix for SSOR iteration
    typename GenericStorage<T>::SparseSymMatrix* mat_ssor_sym, mat_stored_sym;
    typename GenericStorage<T>::SparseMatrix* mat_ssor_unsym, mat_stored_unsym;
    typename GenericStorage<T>::SparseSymMatrixCSR* matCSR_ssor_sym, matCSR_stored_sym;
    typename GenericStorage<T>::SparseMatrixCSR* matCSR_ssor_unsym, matCSR_stored_unsym;

    //! Jacobi preconditioner
    JacobiPreconditioning<T> prec_jacobi;
    //! used matrix for Jacobi iteration
    FemMatrixFreeClass_Base<T>* mat_jacobi, *mat_stored_jacobi;
    //! Subdomain preconditioning
    SubdomainPreconditioning_Base<T>* prec_subdomain;
    //! Multigrid preconditioning
    MultigridPreconditioning_Base<T>* prec_multigrid;
    //! Low-order preconditioning
    LowOrderPreconditioning_Base<T>* prec_low_order;

    //! damping parameter for preconditioner
    /*!
      preconditioning based on approximate inverse of
      \f$ -k^2(delta + i beta) - \Delta \f$
      alpha is equal to delta + i beta
    */
    T alpha_laplace;

    //! if true, local preconditioning is used
    bool local_precond;
    LocalPreconditioning_Base<T>* prec_local;

    //! hypre preconditioning
    Preconditioner_Base<T>* prec_hypre;
    
    void InitDefaultValues();
    
  public :
    // type of preconditioners
    enum {IDENTITY, SSOR, DIRECT, JACOBI,
	  LOW_ORDER, MULTIGRID, SUBDOMAIN, LOCAL, HYPRE};

    // type of mesh splitting (for subdomain preconditioning)
    enum {SUBDIV_DOMAIN_BOXES, SUBDIV_DOMAIN_VERTICES,
	  SUBDIV_DOMAIN_EDGES, SUBDIV_DOMAIN_FACES};
    
    template<class TypeEquation>    
    All_Preconditioner_Base(EllipticProblem<TypeEquation>& var);
    
    virtual ~All_Preconditioner_Base();
    
    virtual void SetInputData(const string& keyword, const VectString& param);
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);

    MultigridPreconditioning_Base<T>& GetMultigridPreconditioning();
    
    void SetDampingCoefficient(const T& alpha);
    void SetPreconditioningType(int type);

    inline void DisableLocalPreconditioning() { local_precond = false; }
    
    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);

    void SolveHypreGen(const SeldonTranspose&, const VirtualMatrix<T>&,
		       const Vector<Complex_wp>& b, Vector<Complex_wp>& x);

    void SolveHypreGen(const SeldonTranspose&, const VirtualMatrix<T>&,
		       const Vector<Real_wp>& b, Vector<Real_wp>& x);

    virtual void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x);
    virtual void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x);

    virtual size_t GetMemorySize() const;
    
    template<class Dimension>
    void GetMaximumSize(T a, Dimension dim, int& nb_max_ilut,
			int& nb_max_iterative, int& nb_max_multigrid,
			int& order_min_multigrid);
    
    void SetOptimalParameters(int order, int N, All_LinearSolver& solver);

    virtual void ConstructPreconditioner(VarComputationProblem&, All_LinearSolver&,
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver);

    void CopyParameter(const Preconditioner_Base<T>&);
    
    template<class Prop, class Storage, class Allocator>
    void SetDirectPreconditioning(DistributedMatrix<T, Prop, Storage, Allocator>& A);
    
    void Clear();
    
  };


  //! base class for local preconditioning
  template<class T>
  class LocalPreconditioning_Base
  {
  protected:
    All_Preconditioner_Base<T>* precond;
    Vector<int> num_elt_precond;
    Vector<int> ListeDof;
    FemMatrixFreeClass_Base<T>* iterative_matrix;
    Vector<T> diagonal;
    CondensationBlockSolver_Fem<T>* local_cond;
    
  public:
    LocalPreconditioning_Base();
    virtual ~LocalPreconditioning_Base();

    void SetInputPreconditioning(const string& keyword, const Vector<string>& param);
    void CopyParameter(const LocalPreconditioning_Base<T>&);

    size_t GetMemorySize() const;
    
    virtual void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
					 const FemMatrixFreeClass_Base<T>& A,
					 const GlobalGenericMatrix<T>& nature,
					 CondensationBlockSolver_Fem<T>& cond_solver) = 0;
    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);

    void FinalizePreconditioner(VarComputationProblem& var,
				VarComputationProblem& var_global, All_LinearSolver& solver,
				const FemMatrixFreeClass_Base<T>& A,
				const GlobalGenericMatrix<T>& nature,
				CondensationBlockSolver_Fem<T>& cond_solver);
    
  };


  //! class for local preconditioning depending on dimension
  template<class T, class Dimension>
  class LocalPreconditioning_Dim : public LocalPreconditioning_Base<T>
  {
  protected:
    DistributedProblem<Dimension>* var_local;

  public:
    LocalPreconditioning_Dim();
    ~LocalPreconditioning_Dim();
    
    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
				 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature,
				 CondensationBlockSolver_Fem<T>& cond_solver);
    
  };


  //! local preconditioning for any equation
  template<class T, class TypeEquation>
  class LocalPreconditioning
    : public LocalPreconditioning_Dim<T, typename TypeEquation::Dimension>
  {
  public :
    LocalPreconditioning(EllipticProblem<TypeEquation>& var);
    
  };

  template<class T, class TypeEquation>
  class All_Preconditioner : public All_Preconditioner_Base<T>
  {
  public :
    All_Preconditioner(EllipticProblem<TypeEquation>& var) : All_Preconditioner_Base<T>(var) {}
  };

#else
  template<class T>
  class All_Preconditioner_Base : public Preconditioner_Base<T>
  {
  protected :
    //! type of preconditioner
    int type_precond;

  public :
    // type of mesh splitting (for subdomain preconditioning)
    enum {SUBDIV_DOMAIN_BOXES, SUBDIV_DOMAIN_VERTICES,
	  SUBDIV_DOMAIN_EDGES, SUBDIV_DOMAIN_FACES};
    
    template<class TypeEquation>
    All_Preconditioner_Base(EllipticProblem<TypeEquation>& var);

    virtual ~All_Preconditioner_Base();
    
    virtual void SetInputData(const string& keyword, const VectString& param);
    virtual void SetInputPreconditioning(const string& keyword, const VectString& param);

    void SetPreconditioningType(int type);

    virtual size_t GetMemorySize() const;

    template<class Dimension>
    void GetMaximumSize(T a, Dimension dim, int& nb_max_ilut,
			int& nb_max_iterative, int& nb_max_multigrid,
			int& order_min_multigrid);
    
    void SetOptimalParameters(int order, int N, All_LinearSolver&);
    
    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver&,
				 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature,
				 CondensationBlockSolver_Fem<T>& cond_solver);
    
    void Clear();
    
  };  

  template<class T, class TypeEquation>
  class All_Preconditioner : public All_Preconditioner_Base<T>
  {
  public :
    All_Preconditioner(EllipticProblem<TypeEquation>& var);
    
  };
#endif
  
}

#define MONTJOIE_FILE_PRECONDITIONER_HXX
#endif
