#ifndef MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_HXX

namespace Montjoie
{
  
  class VirtualMassMatrix
  {
  private:
    VarInstationary_Base& var_time;
    
    void InitDefaultValues();
    
  public:
    bool compute_Dh_default, compute_Sh_default, compute_invDh_default;
    bool compute_DhMinusdtSh_default, compute_invDhPlusdtSh_default, compute_Bh_default;
    bool compute_invBh_default, compute_BhMinusdtSh_default,
      compute_invBhPlusdtSh_default, compute_ShVec_default;
    
  public:
    template<class TypeEquation>
    VirtualMassMatrix(HyperbolicProblem<TypeEquation>& vars);

    virtual ~VirtualMassMatrix();
    
    void FindMatricesToCompute(bool& compute_Dh, bool& compute_invDh,
			       bool& compute_DhMinusdtSh, bool& compute_invDhPlusdtSh,
                               bool& compute_Sh);
    
    void FindMatricesToComputeVec(bool& compute_Dh, bool& compute_invDh,
				  bool& compute_DhMinusdtSh, bool& compute_invDhPlusdtSh,
                                  bool& compute_Sh);

    virtual size_t GetMemorySize() const = 0;
    virtual const VectReal_wp& GetDiagonalDh() const = 0;
    virtual const VectReal_wp& GetInverseDiagonalDh() const = 0;
    virtual const VectReal_wp& GetDiagonalSh() const = 0;
    virtual const VectReal_wp& GetDiagonalDhMinusdtSh() const = 0;
    virtual const VectReal_wp& GetInverseDiagonalDhPlusdtSh() const = 0;
    
    virtual const VectReal_wp& GetDiagonalBh() const = 0;
    virtual const VectReal_wp& GetDiagonalBhMinusdtSh() const = 0;
    virtual const VectReal_wp& GetInverseDiagonalBhPlusdtSh() const = 0;
    
    virtual void Init(bool compute_time = true) = 0;
    
    virtual void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				 const Real_wp& beta, VectReal_wp& Prod_Uh) = 0;

    virtual void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
					  const Real_wp& beta, VectReal_wp& Prod_Uh) = 0;
    
    virtual void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				 const Real_wp& beta, VectReal_wp& Prod_Uh) = 0;
    
    virtual void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
					  const Real_wp& beta, VectReal_wp& Prod_Uh) = 0;
    
    virtual void SolveOperatorDhPlusdtSh(VectReal_wp& Y) = 0;
    virtual void SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) = 0;
    
    virtual void SolveOperatorDh(VectReal_wp& Y) = 0;
    virtual void SolveMassMatrix(VectReal_wp& Y) = 0;    
    virtual void ApplyMassMatrix(const Real_wp&, const Real_wp& , const VectReal_wp&,
				 const Real_wp&, VectReal_wp& Y) = 0;    
    
    virtual void SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y) = 0;
    
    virtual void SolveOperatorDh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) = 0;

    virtual void ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
				 const Real_wp& beta, VectReal_wp& Vh) = 0;
    
    virtual void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
					  const Real_wp& beta, VectReal_wp& Vh) = 0;
    
    virtual void SolveOperatorBh(VectReal_wp& Y) = 0;
    virtual void SolveOperatorBh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) = 0;
    
    virtual void SolveOperatorBhPlusdtSh(VectReal_wp& Y) = 0;
    virtual void SolveOperatorBhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y) = 0;

    // for HDG
    virtual void SolveOperatorCh(VectReal_wp& Y);
    
  };
  
  
  //! class storing the mass matrix, damping matrix or combinations for time schemes 
  //! and continuous formulation
  /*!
    This class is able to store Dh (mass matrix), Sh (damping matrix)
    and combinations Dh - dt/2 Sh, Dh + dt/2 Sh
    For first-order formulations, it can also store Bh (mass matrix of additional unknowns), ShV
    and combinations Bh - dt/2 ShV, Bh + dt/2 ShV
   */
  class ContinuousUnsteadyMassMatrix_Base : public VirtualMassMatrix
  {
  protected :
    // different mass matrices we need for scalar unknown
    ContinuousMassMatrix_Base* Dh, *DhMinusdtSh, *DhPlusdtSh, *Sh;
    
    // and mass matrices for vectorial unknown
    ContinuousMassMatrix_Base* Bh, *BhMinusdtSh, *BhPlusdtSh, *ShVec;
    
    // dofs/elements which are not mass-lumped
    IVect NonLumpedElement, NonLumpedDof, InverseDof, LumpedDof;

    //! global row numbers for scalar dofs (needed in parallel)
    IVect scalar_rows_;
    DistributedMatrixIntegerArray info_sparse, info_block;
    
  private:
    VarComputationProblem& var_computation;
    VarProblem_Base& var_problem;
    DistributedProblem_Base& var_comm;
    VarBoundaryCondition_Base& var_boundary;
    VarInstationary_Base& var_time;
    
    void InitDefaultValues();

  public :
    template<class TypeEquation>
    ContinuousUnsteadyMassMatrix_Base(HyperbolicProblem<TypeEquation>& vars);

    ~ContinuousUnsteadyMassMatrix_Base();

    void Clear();
    void SolveMassMatrix(VectReal_wp& Y);    
    void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
			 const VectReal_wp& x, const Real_wp& beta, VectReal_wp& Y);    
    
    void Init(bool compute_time = true);

    virtual size_t GetMemorySize() const;

    template<class Prop, class Storage>
    static void ExtractScalarMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
				    int, int nb_global_scalar_dof, Vector<int>& scalar_rows,
				    DistributedMatrix<Real_wp, Prop, Storage>& B);

    static bool IsLocalSymmetricMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
				       int N, Real_wp threshold = 1e-12);

    void ConstructNonLumpedDofs(const Vector<bool>& diag_elt);
    void ConstructLumpedDofs(Vector<int>&);

    void ExtractDiagonal(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
			 const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma,
			 int N, VectReal_wp& Dh_diagonal, VectReal_wp& Sh_diagonal,
			 VectReal_wp& Dh_ass, VectReal_wp& Sh_ass);
    
    void FillDiagonal(const VectReal_wp& Dh_diagonal, const VectReal_wp& Sh_diagonal,
		      const VectReal_wp& Dh_ass, const VectReal_wp& Sh_ass,
		      bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		      bool compute_invDhPlusdtSh, bool compute_Sh);

    void FillBlockDiagonal(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& Dh_block,
			   const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& Sh_block,
			   bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
			   bool compute_invDhPlusdtSh, bool compute_Sh, bool invert = true);

    void FillBlockDiagonal(const DistributedMatrix<Real_wp, General, BlockDiagRow>& Dh_block,
			   const DistributedMatrix<Real_wp, General, BlockDiagRow>& Sh_block,
			   bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
			   bool compute_invDhPlusdtSh, bool compute_Sh, bool invert = true);

    void FillSparseMatrices(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Dh_s,
			    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Sh_s,
			    bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
			    bool compute_invDhPlusdtSh, bool compute_Sh);
    
    void FillSparseMatrices(DistributedMatrix<Real_wp, General, ArrayRowSparse>& Dh_s,
			    DistributedMatrix<Real_wp, General, ArrayRowSparse>& Sh_s,
			    bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
			    bool compute_invDhPlusdtSh, bool compute_Sh);

    void SymmetrizePattern(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass);
    
    virtual void ComputeScalarMassMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
                                         DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma);

    virtual void ComputeVectorialMassMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
                                            DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma);
    
    virtual ContinuousMassMatrix_Base* GetNewMassMatrix(int type) = 0;
    
    const VectReal_wp& GetDiagonalDh() const;
    const VectReal_wp& GetInverseDiagonalDh() const;
    const VectReal_wp& GetDiagonalSh() const;
    const VectReal_wp& GetDiagonalDhMinusdtSh() const;
    const VectReal_wp& GetInverseDiagonalDhPlusdtSh() const;
    const VectReal_wp& GetDiagonalBh() const;
    const VectReal_wp& GetDiagonalBhMinusdtSh() const;
    const VectReal_wp& GetInverseDiagonalBhPlusdtSh() const;

    void InitDiagonalMass();

    void SetDiagonalDh(const VectReal_wp& D);
    void SetDiagonalSh(const VectReal_wp& D);
    void SetInverseDiagonalDhPlusdtSh(const VectReal_wp& D);

    void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);

    void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				    const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void SolveOperatorDhPlusdtSh(VectReal_wp& Y);
    void SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y);
    
    void SolveOperatorDh(VectReal_wp& Y);
    
    void SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y);
    
    void SolveOperatorDh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y);

    void ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                         const Real_wp& beta, VectReal_wp& Vh);
    
    void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                                  const Real_wp& beta, VectReal_wp& Vh);
    
    void SolveOperatorBh(VectReal_wp& Y);
    void SolveOperatorBh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y);
    
    void SolveOperatorBhPlusdtSh(VectReal_wp& Y);
    void SolveOperatorBhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y);
    
  };
  
  
  template<class Dimension>
  class ContinuousUnsteadyMassMatrix_Dim : public ContinuousUnsteadyMassMatrix_Base
  {
  public:
    template<class TypeEquation>
    ContinuousUnsteadyMassMatrix_Dim(HyperbolicProblem<TypeEquation>& vars);
    
  };
  
  
  template<class Dimension>
  class ContinuousUnsteadyMassMatrix : public ContinuousUnsteadyMassMatrix_Dim<Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
    
  public:
    template<class TypeEquation>
    ContinuousUnsteadyMassMatrix(HyperbolicProblem<TypeEquation>& vars);

    ContinuousMassMatrix_Base* GetNewMassMatrix(int type_mat);
    
  };

} // end namespace

#define MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_HXX
#endif
