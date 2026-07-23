#ifndef MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_HXX

namespace Montjoie
{
  
  //! base class for mass matrix and continuous formulation
  /*!
    The mass matrix is equal to \int C varphi_i varphi_j dx
    where C is a physical coefficient (for example density rho)
    This mass matrix may be diagonal, block-diagonal, sparse, etc
    This class implements the different available storages and different
    solvers in order to invert efficiently the mass matrix since this step
    is crucial for time explicit schemes
   */
  class ContinuousMassMatrix_Base : public VirtualMatrix<Real_wp>,
				    public Preconditioner_Base<Real_wp>
  {
    
  protected :
    //! reference to the non-lumped dofs
    IVect& NonLumpedElement;
    IVect& NonLumpedDof;
    IVect& InverseDof;
    IVect& LumpedDof;
    int nodl_mesh;
    
    //! mass matrix stored as a matrix-free object
    FemMatrixFreeClass_Base<Real_wp>* free_mat;    

    //! diagonal if the mass matrix is diagonal
    VectReal_wp diagonal;
    //! inverse of the diagonal part when a part of mass matrix is diagonal
    VectReal_wp invDiagonal;
    
    //! block-diagonal
    DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> block_diagonal;
    //! inverse of the block-diagonal part
    DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> invBlock_diagonal;
    
    //! block-diagonal
    DistributedMatrix<Real_wp, General, BlockDiagRow> block_diagonal_unsym;
    //! inverse of the block-diagonal part
    DistributedMatrix<Real_wp, General, BlockDiagRow> invBlock_diagonal_unsym;
    
    //! diagonal preconditioning
    VectReal_wp precond_diag;
    
    //! case where mass matrix is only a sparse matrix
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> mat_sparse;

    //! case where mass matrix is only a sparse matrix
    DistributedMatrix<Real_wp, General, ArrayRowSparse> mat_sparse_unsym;

    bool iterative_solver;
    All_MatrixLU<Real_wp> mat_lu;
    
    //! current type of this mass matrix
    int type_matrix;
    //! initial guess for iterative solver
    VectReal_wp previous_iterate;
    //! temporary vectors
    VectReal_wp Xn_tmp, Bn_tmp;
    //! stopping criterion used by iterative solver
    Real_wp stopping_criterion;
    int nb_max_iteration, print_level;    

    void InitDefaultValues();
    
  public:
    template<class Dimension>
    ContinuousMassMatrix_Base(VarProblem<Dimension>& var, int type_mat,
			      IVect& non_lumped_elt, IVect& non_lumped_dof, IVect& inverse_dof,
			      IVect& lumped_dof);
    
    ~ContinuousMassMatrix_Base();
    
    int GetNumberOfRows() const;
    int GetMatrixType() const;

    void SetDiagonal(const VectReal_wp& D);
    void SetInverseDiagonal(const VectReal_wp& invD);

    VectReal_wp& GetDiagonal();
    const VectReal_wp& GetDiagonal() const;
    VectReal_wp& GetInverseDiagonal();
    const VectReal_wp& GetInverseDiagonal() const;

    DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& GetBlockDiagonal();
    DistributedMatrix<Real_wp, General, BlockDiagRow>& GetUnsymmetricBlockDiagonal();

    DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& GetInverseBlockDiagonal();
    DistributedMatrix<Real_wp, General, BlockDiagRow>& GetInverseUnsymmetricBlockDiagonal();
    
    DistributedMatrix<Real_wp, General, ArrayRowSparse>& GetSparseMatrix();
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& GetSymSparseMatrix();

    void SetIterativeSolver(bool iterative, int nb_iter, const Real_wp& eps, int lvl);
    void SetPrintLevel(int lvl);
    
    size_t GetMemorySize() const;
    
    void MltMass(const VectReal_wp&, VectReal_wp& Vh);
    void MltMass(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh);

    void SolveSparseMass(const VectReal_wp& B, VectReal_wp& X);
    void SolveSparseMassUnsym(const VectReal_wp& B, VectReal_wp& X);
    
    void SolveMass(VectReal_wp& Y);
    void SolveMass(const Real_wp& alpha, const VectReal_wp&, VectReal_wp& Y);

    template<class Prop, class Storage>
    void ComputePreconditioning(DistributedMatrix<Real_wp, Prop, Storage>& A);
    
    void Invert(bool keep_mass);
    
    // methods for quasi-lumped mass matrices
    void GetSource(const VectReal_wp& Uh, VectReal_wp& B) const;
    
    void AddSolution(const Real_wp& alpha, const VectReal_wp& Uh,
		     const VectReal_wp&, VectReal_wp& X) const;

    template<class Prop, class Storage>
    void AddBlockSolution(const Real_wp& alpha,
			  const DistributedMatrix<Real_wp, Prop, Storage>& M,
			  const VectReal_wp& Uh, const VectReal_wp& X, VectReal_wp& Y) const;
    
    void GetSolution(const VectReal_wp& Uh, VectReal_wp& X) const;

    template<class Prop, class Storage>
    void GetBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
			  const VectReal_wp& Uh, VectReal_wp& X) const;
    
    void MltSolution(const VectReal_wp& Uh, VectReal_wp& X) const;

    template<class Prop, class Storage>
    void ExtractBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
			      const VectReal_wp& U, VectReal_wp& B) const;

    template<class Prop, class Storage>
    void MltBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
			  const VectReal_wp& Uh, const VectReal_wp& U, VectReal_wp& X) const;
    
    void MltAddSolution(const Real_wp& alpha, const VectReal_wp& Uh,
 			const VectReal_wp&, VectReal_wp& X) const;

    template<class Prop, class Storage>
    void MltAddBlockSolution(const Real_wp& alpha, const DistributedMatrix<Real_wp, Prop, Storage>& M,
			     const VectReal_wp& Uh, const VectReal_wp& X, VectReal_wp& Y) const;
    
    void Solve(const VirtualMatrix<Real_wp>& A, const VectReal_wp& r, VectReal_wp& z);
    void TransSolve(const VirtualMatrix<Real_wp>& A, const VectReal_wp& r, VectReal_wp& z);
    
  };


  template<class Dimension>
  class ContinuousMassMatrix : public ContinuousMassMatrix_Base
  {
  private:
    MeshNumbering<Dimension>& mesh_num;
    Mesh<Dimension>& mesh;
    VarProblem<Dimension>& var_problem;
    
  public:
    ContinuousMassMatrix(VarProblem<Dimension>& var, int type_mat,
			 IVect& non_lumped_elt, IVect& non_lumped_dof, IVect& inverse_dof,
			 IVect& lumped_dof);

    void Init();
    
    void MltAddElement(const Real_wp& alpha, int i, int ref, VectReal_wp& Uh, VectReal_wp& Vh,
                       IVect& num_ddl, const VectReal_wp& X, VectReal_wp& Y, const ElementReference_Dim<Dimension>& Fb) const;
    
    void MltAddVector(const Real_wp& alpha, const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y) const;

    void MltAddVector(const Complex_wp& alpha, const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Y) const;
    void MltAddVector(const Real_wp& alpha, const SeldonTranspose&, const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y) const;
    void MltAddVector(const Complex_wp& alpha, const SeldonTranspose&, const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Y) const;
    
    void MltVector(const VectReal_wp& X, VectReal_wp& Y) const;
    void MltVector(const VectComplex_wp& X, VectComplex_wp& Y) const;
    void MltVector(const SeldonTranspose&, const VectReal_wp& X, VectReal_wp& Y) const;
    void MltVector(const SeldonTranspose&, const VectComplex_wp& X, VectComplex_wp& Y) const;
    
  };

} // end namespace

#define MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_HXX
#endif
