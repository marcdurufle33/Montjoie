#ifndef MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_HXX

namespace Montjoie
{

  //! base class for a block (of a discontinuous mass matrix)
  class DiscontinuousBlockMassMatrix
  {
  public :
    virtual ~DiscontinuousBlockMassMatrix();

    virtual void MltMass(VectReal_wp& Y) const = 0;
    virtual void SolveMass(VectReal_wp& Y) const = 0;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const = 0;
    virtual size_t GetMemorySize() const = 0;
    
  };


  //! case of a diagonal block (of a discontinuous mass matrix)
  class DiscontinuousDiagonalMassMatrix : public DiscontinuousBlockMassMatrix
  {
  protected:
    // diagonal storage
    Vector<Real_wp> diagonal;
    
  public:
    DiscontinuousDiagonalMassMatrix(const VectReal_wp& d);
    
    virtual void MltMass(VectReal_wp& Y) const;
    virtual void SolveMass(VectReal_wp& Y) const;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const;

    void SetDiagonal(const VectReal_wp& d);

    size_t GetMemorySize() const;
    
  };


  //! case of a block-diagonal block (of a discontinuous mass matrix)
  template<int d>
  class DiscontinuousBlockDiagMassMatrix : public DiscontinuousBlockMassMatrix
  {
  protected:
    // row numbers
    Vector<int> row_num;

    // associated blocks
    Vector<TinyMatrix<Real_wp, Symmetric, d, d> > blocks;

  public:
    DiscontinuousBlockDiagMassMatrix(const Vector<int>& num, const Vector<TinyMatrix<Real_wp, Symmetric, d, d> >& B);
    
    virtual void MltMass(VectReal_wp& Y) const;
    virtual void SolveMass(VectReal_wp& Y) const;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const;

    void SetBlockDiagonal(const Vector<int>& num, const Vector<TinyMatrix<Real_wp, Symmetric, d, d> >& B);
    
    size_t GetMemorySize() const;
    
  };
  
    
  //! case of a sparse block (of a discontinuous mass matrix)
  class DiscontinuousSparseMassMatrix : public DiscontinuousBlockMassMatrix
  {
  protected:
    // sparse storage of Cholesky factor
    Matrix<Real_wp, Symmetric, RowSymSparse> sparse_Mh;
    
  public:
    DiscontinuousSparseMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked>& A, const Real_wp& eps);
    
    virtual void MltMass(VectReal_wp& Y) const;
    virtual void SolveMass(VectReal_wp& Y) const;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const;

    void ConvertDense(Matrix<Real_wp, Symmetric, RowSymPacked>& M, const Real_wp& eps);

    size_t GetMemorySize() const;
    
  };


  //! case of a dense block (of a discontinuous mass matrix)
  class DiscontinuousDenseMassMatrix : public DiscontinuousBlockMassMatrix
  {
  protected:
    // dense storage
    Matrix<Real_wp, Symmetric, RowSymPacked> dense_Mh;
    
  public:
    DiscontinuousDenseMassMatrix(const Matrix<Real_wp, Symmetric, RowSymPacked>& A);
      
    virtual void MltMass(VectReal_wp& Y) const;
    virtual void SolveMass(VectReal_wp& Y) const;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const;

    void SetDense(const Matrix<Real_wp, Symmetric, RowSymPacked>& M);

    size_t GetMemorySize() const;
    
  };


  //! case of a matrix-free block (of a discontinuous mass matrix)
  class DiscontinuousMatrixFreeMassMatrix : public DiscontinuousBlockMassMatrix,
					    public VirtualMatrix<Real_wp>,
					    public Preconditioner_Base<Real_wp>
  {
  protected : 
    Real_wp stopping_criterion;

    // temporary vector
    mutable Vector<Real_wp> xn_tmp;
    
    // diagonal preconditioning
    Vector<Real_wp> diagonal;

    //! reference to finite element class
    ElementReference_Base& Fb;
    
  public:
    DiscontinuousMatrixFreeMassMatrix(ElementReference_Base& elt);
    
    virtual void MltMass(VectReal_wp& Y) const;
    virtual void SolveMass(VectReal_wp& Y) const;

    virtual void SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const;

    void SetStoppingCriterion(const Real_wp& eps);
    Real_wp GetStoppingCriterion() const;
    
    // preconditioning 
    void Solve(const VirtualMatrix<Real_wp>& A,
	       const VectReal_wp& r, VectReal_wp& z);

    size_t GetMemorySize() const;
    
  };


  //! interface class for mass matrix 
  class DiscontinuousMassMatrixProblem
  {
  public :
    // pure virtual methods
    virtual int GetNbElt() const = 0;
    virtual void GetElementInfo(int n, int nu, bool,
				int& ref, bool& affine, bool& varying) const = 0;
    
    virtual void GetJacobianPanel(int n, bool& affine, bool& curved,
				  Real_wp& jacob_affine, VectReal_wp& decomp_jacob,
				  VectReal_wp& jacob_quad) = 0;

    virtual Real_wp GetAffineJacobian(int) const = 0;
    
    virtual const ElementReference_Base& GetReferenceElement(int) const = 0;

    // optional virtual method
    virtual void GetMassDampingCoefficient(int n, int k, int nu, int ref,
					   Real_wp& rho, Real_wp& sigma);

  };


  //! block mass matrix (volume integrals) for unsteady problem
  template<class Dimension>
  class DiscontinuousMassMatrixVol : public DiscontinuousMassMatrixProblem
  {
  private:
    Mesh<Dimension>& mesh;
    VarProblem<Dimension>& var_problem;
    VarInstationary_Base& var_time;
    
  public:
    DiscontinuousMassMatrixVol(VarProblem<Dimension>& var_p,
			       VarInstationary_Base& var_t);
    
    int GetNbElt() const;    
    void GetElementInfo(int i, int num, bool geom_mass,
			int& ref, bool& affine, bool& varying) const;

    void GetJacobianPanel(int n, bool& affine, bool& curved,
			  Real_wp& jacob_affine, VectReal_wp& decomp_jacob,
			  VectReal_wp& jacob_quad);

    void GetMassDampingCoefficient(int n, int k, int nu, int ref,
				   Real_wp& rho, Real_wp& sigma);

    Real_wp GetAffineJacobian(int n) const;

    const ElementReference_Base& GetReferenceElement(int) const;
    
  };


  //! block mass matrix (surface integrals) for unsteady problem
  template<class Dimension>
  class DiscontinuousMassMatrixSurf : public DiscontinuousMassMatrixProblem
  {
  private:
    Mesh<Dimension>& mesh;
    VarProblem<Dimension>& var_problem;
    VarInstationary_Base& var_time;
    
  protected:
    VectReal_wp& ValueDsj;
    IVect& OffsetDsj;
    VectReal_wp coef_surface;
    
  public:
    DiscontinuousMassMatrixSurf(VarProblem<Dimension>& var_p,
				VarInstationary_Base& var_t,
				VectReal_wp& val_ds, IVect& offset_ds, bool geom);

    void ConstructDs();
    void ConstructCoef();
    
    int GetNbElt() const;    
    void GetElementInfo(int i, int num, bool geom_mass,
			int& ref, bool& affine, bool& varying) const;

    void GetJacobianPanel(int n, bool& affine, bool& curved,
			  Real_wp& jacob_affine, VectReal_wp& decomp_jacob,
			  VectReal_wp& jacob_quad);

    void GetMassDampingCoefficient(int n, int k, int nu, int ref,
				   Real_wp& rho, Real_wp& sigma);

    Real_wp GetAffineJacobian(int n) const;
    
    const ElementReference_Base& GetReferenceElement(int) const;

    size_t GetMemorySize() const;
    
  };

  
  //! mass matrix for discontinuous approximation and all the mesh
  class DiscontinuousMassMatrix_Base
  {
  protected :
    // threshold used to convert a dense matrix to a sparse matrix
    Real_wp threshold;
    Real_wp stopping_criterion;

    // diagonal if the mass matrix is diagonal
    VectReal_wp diagonal, invDiagonal;
    
    // block-diagonal
    Matrix<Real_wp, Symmetric, BlockDiagRowSym> block_diagonal;
        
    // mass for each elements for BLOCK_ELT storage
    Vector<DiscontinuousBlockMassMatrix*> block_mass_store;
    
    // pointer for BLOCK_ELT storage
    Vector<Vector<DiscontinuousBlockMassMatrix* > > block_mass;
    IVect block_num_element, block_size;
    Vector<VectReal_wp> coef_block_mass;

    // pointer to geometrical mass matrix
    DiscontinuousMassMatrix_Base* GeomMass;
    
    // elements where the mass matrix is considered
    MatrixVectorProductLevel level_time;
    
    // current type of this mass matrix
    int type_matrix;
    
    // variables used by a block for matrix vector product
    int num_unknown;
    Real_wp coef_mass_, coef_damp_;

    //! object containing geometry and physical information
    DiscontinuousMassMatrixProblem& var_problem;
    
  public:
    DiscontinuousMassMatrix_Base(VarInstationary_Base& var_t,
				 DiscontinuousMassMatrix_Base* var_geom,
				 DiscontinuousMassMatrixProblem&, int type_mat);
    
    virtual ~DiscontinuousMassMatrix_Base();
    
    // inline methods
    void SetThreshold(const Real_wp& eps);
    Real_wp GetThreshold() const;

    void SetStoppingCriterion(const Real_wp& eps);
    Real_wp GetStoppingCriterion() const;

    int GetNbElt() const;
    int GetMatrixType() const;
    void SetUnknownNumber(int m);
    void SetCoefficient(const Real_wp& a, const Real_wp& b);
    
    VectReal_wp& GetDiagonal();
    const VectReal_wp& GetDiagonal() const;
    VectReal_wp& GetInverseDiagonal();
    const VectReal_wp& GetInverseDiagonal() const;

    const IVect& GetOffsetSizeBlocks() const;
    MatrixVectorProductLevel& GetLevelTime();
    
    // other methods    
    void ClearBlocks();

    virtual size_t GetMemorySize() const;
    
    void MltMass(int offset, const VectReal_wp& X, VectReal_wp& Y);
    void MltMass(int offset, const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh);
    
    void SolveMass(int offset, VectReal_wp& Y);
    void SolveMass(int offset, const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y);
    
    void SolveCholesky(int offset, const SeldonTranspose& transA, VectReal_wp& Y);
    
    void ComputeLocalMass(int num_elem, int ref, bool varying, int& nb,
                          bool geom_mass, const ElementReference_Base& Fb);

    void ReallocateStoredMass();
    void SetLocalStoredBlock(int i, DiscontinuousBlockMassMatrix*);
    void CopyBlockCoef(const VectReal_wp& coef, int num_u);

    void ComputeMass(bool, int m1 = 0, int m2 = 1);
    
    void MltLocalMass(int i, int n, int offset, const VectReal_wp& Uh,
		      int&, VectReal_wp& Uloc);
    
    void SolveLocalMass(int i, int n, int offset, const VectReal_wp& Uh,
			int&, VectReal_wp& Uloc);
    
    void SolveLocalCholesky(int i, int offset, const SeldonTranspose& transA,
			    const VectReal_wp& Uh, VectReal_wp& Uloc);

    virtual void SetDiagonalPML(int nb_new, int m1, int m2, int nodl, VectReal_wp& diag_masse);
    
  };


  //! class for unsteady mass matrices
  template<class Dimension>
  class DiscontinuousMassMatrix_Dim : public DiscontinuousMassMatrix_Base
  {
  private:    
    Mesh<Dimension>& mesh;
    MeshNumbering<Dimension>& mesh_num;
    VarBoundaryCondition_Base& var_boundary;
    DistributedProblem<Dimension>& var_problem;
    VarInstationary_Base& var_time;
    
  public:
    DiscontinuousMassMatrix_Dim(DistributedProblem<Dimension>& var,
				VarInstationary_Base& var_t,
				DiscontinuousMassMatrixProblem&, 
				DiscontinuousMassMatrix_Base* var_geom, int type_mat);
    
    void SetDiagonalPML(int nb_new, int m1, int m2, int nodl, VectReal_wp& diag_masse);

  };


  //! leaf class
  template<class Dimension>
  class DiscontinuousMassMatrix : public DiscontinuousMassMatrix_Dim<Dimension>
  {
  protected:
    DiscontinuousMassMatrixVol<Dimension> mass_object;
    
  public:
    DiscontinuousMassMatrix(VarProblem<Dimension>& var,
			    VarInstationary_Fem<Dimension>& var_t,
			    DiscontinuousMassMatrix_Base* var_geom, int type_mat);

    
  };


  //! class for unsteady mass matrices
  template<class Dimension>
  class DiscontinuousSurfaceMassMatrix : public DiscontinuousMassMatrix_Base
  {
  protected:
    VectReal_wp ValueDsj; IVect OffsetDsj;
    DiscontinuousMassMatrixSurf<Dimension> mass_surf;
    
  public:
    DiscontinuousSurfaceMassMatrix(VarProblem<Dimension>& var,
				   VarInstationary_Base& var_t,
				   DiscontinuousMassMatrix_Base* var_geom, int type_mat,
				   VectReal_wp& val_ds, IVect& offset_ds);

    DiscontinuousSurfaceMassMatrix(VarProblem<Dimension>& var,
				   VarInstationary_Base& var_t,
				   DiscontinuousMassMatrix_Base* var_geom, int type_mat);

    VectReal_wp& GetDs();
    IVect& GetOffsetDs();

    size_t GetMemorySize() const;
    
  };
  
} // end namespace

#define MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_HXX
#endif
