#ifndef MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_HXX

namespace Montjoie
{
  
  //! base class storing intermediate variable for the matrix vector product
  template<class Complexe>
  class ExtrapolVariablesProduct_Base
  {
  public :

#ifdef SELDON_WITH_MPI
    Vector<Vector<MPI_Request> > request;
    Vector<Vector<Vector<Complexe> > > Xneighbor, Xsend;
    Vector<Vector<Vector<int64_t> > > Xneighbor_tmp, Xsend_tmp;
#endif
    
  public:
    // Inline methods
    void ReallocateExchangeVector(int n);

    // Other methods
    template<class Dimension, class T0>
    static void
    GetFaceOrientationAndNumber(const DistributedProblem<Dimension>& var,
				int num_face, int num_elem1, int num_pos1_face,
				int offset_face1, int & ref_line, int& cond,
				bool& new_face, int& num_elem2, int& num_pos2_face,
				int& offset_face2, int& rot,
				T0& phase, T0& phase_conj);
    
  };
  
  
  //! variables used for H1 matrix-vector product
  template<class Complexe, class TypeEquation>
  class ExtrapolVariablesProductFEM
    : public ExtrapolVariablesProduct_Base<Complexe>
  {
  public :
    enum {nb_unknowns = TypeEquation::nb_unknowns};
    
    TinyVector<Vector<Complexe>, nb_unknowns>
    extrapolU, extrapolDun, extrapolCgradU;
    
    Vector<TinyVector<Vector<Complexe>, nb_unknowns> > Un_quad;
    TinyVector<Vector<Complexe>, nb_unknowns> Uloc, Uquad;
    TinyVector<Vector<Complexe>, nb_unknowns>
    Uneighbor, DunNeighbor, CgradUNeighbor;
    
    TinyVector<bool, nb_unknowns> unknown_to_derive, fct_test_to_derive;
    
  public:
    // Inline methods
    template<class Dimension, int type>
    TinyVector<Vector<Complexe>,TypeEquation::nb_unknowns>&
    GetUnQuad(int i, const ElementReference<Dimension, type>& Fb);

    size_t GetMemorySize() const;
    
  };


  // file to define "black-box" matrices, where the matrix
  // may not be stored, and matrix-vector product is done "on the fly"
  // particularly useful for hexahedral elements
  
  //! Base structure for finite element matrices
  /*!
    It is a matrix class, that is used for matrix-free
    algorithms, or in the case where the matrix is stored.
    Functions MltAdd are overloaded depending of finite element
  */
  template<class T>
  class FemMatrixFreeClass_Base : public VirtualMatrix<T>
  {
  public :
    //! additional interactions coming from boundary conditions
    typename GenericStorage<T>::SparseSymMatrix mat_boundary_sym;
    typename GenericStorage<T>::SparseMatrix mat_boundary_unsym;
    
    //! sparse matrix if storage is required (with triangles for example)
    typename GenericStorage<T>::SparseSymMatrix mat_iterative_sym;
    typename GenericStorage<T>::SparseMatrix mat_iterative_unsym;
    
    //! additional interactions coming from boundary conditions
    typename GenericStorage<T>::SparseSymMatrixCSR matCSR_boundary_sym;
    typename GenericStorage<T>::SparseMatrixCSR matCSR_boundary_unsym;
    
    //! sparse matrix if storage is required (with triangles for example)    
    typename GenericStorage<T>::SparseSymMatrixCSR matCSR_iterative_sym;
    typename GenericStorage<T>::SparseMatrixCSR matCSR_iterative_unsym;

  private :    

    //! reference to the object containing information about Dirichlet dofs
    const VarBoundaryCondition_Base& var_dir;
    
    //! reference to the object VarProblem
    const VarProblem_Base& var_problem;
    
    //! reference to the object handling communications
    const DistributedProblem_Base& var_comm;

    void InitDefaultValues();
    
  protected :
    //! if true rows and columns are scaled
    bool row_scaling, column_scaling;
    
    //! scaling for rows and colums
    VectReal_wp* row_scale, *col_scale;
    //! if true, the matrix vector product is performed ignoring Dirichlet dofs
    /*!
      If true, it should be similar to a matrix-vector product
      with Neumann condition (=> no boundary terms)
    */
    bool ignore_dirichlet_dof;
    
    //! symmetry of the matrix ?
    bool sym_matrix, sym_dirichlet_condition;    
    
    //! coefficient on the diagonal for each dirichlet dof
    Real_wp coef_dirichlet;
    
    //! object containing the elements to consider for the matrix vector product
    MatrixVectorProductLevel* prod_level;
    
    //! inhomogeneous Dirichlet condition
    Vector<Vector<Real_wp> > source_real_inhg_dirichlet;
    Vector<Vector<Complex_wp> > source_cplx_inhg_dirichlet;
    
    //! columns of matrix corresponding to Dirichlet dofs
    DistributedMatrix<T, General, RowSparse> column_dirichlet;

    //! coefficients used in the finite element matrix
    GlobalGenericMatrix<T> nature_matrix;

    //! solver used to remove internal nodes
    CondensationBlockSolver_Fem<T>* condensed_solver;
    bool use_symmetrization;
    
  public :
    
    template<class TypeEquation>
    FemMatrixFreeClass_Base(const EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    void SetCoefficientDirichlet(const Real_wp& coef);    
    void SetCoefficientMatrix(const GlobalGenericMatrix<T>& nat_mat);
    
    T GetCoefMass() const;

    bool IsSymmetric() const;
    int FormulationDG() const;
    void SetCondensedSolver(CondensationBlockSolver_Fem<T>* solver);
    
    const T operator()(int i, int j) const;

    bool DirichletDofIgnored() const;
    void IgnoreDirichletDof();

    void SetScaling(VectReal_wp& diagonal_scale_left,
                    VectReal_wp& diagonal_scale_right);

    bool SucceedInAffectingPointer(typename GenericStorage<T>::SparseSymMatrix*&,
				   typename GenericStorage<T>::SparseSymMatrixCSR*&);
    
    bool SucceedInAffectingPointer(typename GenericStorage<T>::SparseMatrix*&,
				   typename GenericStorage<T>::SparseMatrixCSR*&);

    inline void SetSymmetrizationUse(bool flag) { use_symmetrization = flag; }
    inline bool GetSymmetrizationUse() const { return use_symmetrization; }
    inline const DistributedProblem_Base& GetDistributedProblem() const { return var_comm; }
    
    // other methods
    void Reallocate(int m, int n);
    void Clear();
    
    size_t GetMemorySize() const;
    
    void InitSymmetricMatrix();
    void InitUnsymmetricMatrix();

    template<class Vector1>
    void ApplyRightScaling(const Vector1& B2, Vector1& C2, Vector1& B,
			   Vector1& C) const;
    
    template<class Vector1>
    void ApplyLeftScaling(const Vector1& B2, Vector1& C2, Vector1& B,
			  Vector1& C) const;
        
    void CompressMatrix();
    void WriteText(const string& file_out) const;    
    
    void AddExtraBoundaryTerms(const Real_wp& alpha,
                               const Vector<Real_wp>&, Vector<Real_wp>&) const;

    void AddExtraBoundaryTerms(const Complex_wp& alpha,
                               const Vector<Complex_wp>&, Vector<Complex_wp>&) const;
    
    // Dirichlet stuff
    void SetNbDirichletCondition(int k);
    
    void ApplyDirichletCondition(const SeldonTranspose&, Vector<Real_wp>& b_rhs, int k = 0, bool assemble = false);    
    void ApplyDirichletCondition(const SeldonTranspose&, Vector<Complex_wp>& b_rhs, int k = 0, bool assemble = false);
    
    template<class T0>
    void MltAddHetereogeneousDirichlet(const T0& alpha, const SeldonTranspose& trans,
                                       const Vector<T0>& B, Vector<T0>& C) const;
    
    template<class MatrixSparse>
    void SetDirichletCondition(MatrixSparse& mat_sp, int offset_row = 0, int offset_col = 0,
                               bool erase_col = false);
    
    void SetDirichletCondition(int offset_row = 0, int offset_col = 0, bool erase_col = false);
    
    template<class Prop>
    void SetDirichletCondition(Matrix<T, Prop, DiagonalRow>& A, int offset_row = 0,
                               int offset_col = 0, bool erase_col = false);

    template<class Prop>
    void SetDirichletCondition(DistributedMatrix<T, Prop, BlockDiagRow>& A, int offset_row = 0,
                               int offset_col = 0, bool erase_col = false);

    template<class Prop>
    void SetDirichletCondition(DistributedMatrix<T, Prop, BlockDiagRowSym>& A, int offset_row = 0,
                               int offset_col = 0, bool erase_col = false);
    
    void InitDirichletCondition(VectReal_wp&, int k = 0);
    void InitDirichletCondition(VectComplex_wp&, int k = 0);

    void ImposeDirichletCondition(const SeldonTranspose&, VectReal_wp&, int k = 0);
    void ImposeDirichletCondition(const SeldonTranspose&, VectComplex_wp&, int k = 0);    

    // other functions
    void AddRowSum(VectReal_wp&);
    void AddRowColSum(VectReal_wp&, VectReal_wp&);

    // matrix vector products
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    void MltAddVector(const Treal&, const Vector<Treal>&, const Treal&,
		      Vector<Treal>&, bool assemble) const;

    void MltAddVector(const Tcplx&, const Vector<Tcplx>&, const Tcplx&,
		      Vector<Tcplx>&, bool assemble) const;

    void MltAddVector(const Treal&, const Vector<Treal>&, const Treal&, Vector<Treal>&) const;
    void MltAddVector(const Tcplx&, const Vector<Tcplx>&, const Tcplx&, Vector<Tcplx>&) const;

    void MltAddVector(const Treal&, const SeldonTranspose&, const Vector<Treal>&,
		      const Treal&, Vector<Treal>&) const;
    
    void MltAddVector(const Tcplx&, const SeldonTranspose&, const Vector<Tcplx>&,
		      const Tcplx&, Vector<Tcplx>&) const;

    void MltAddVector(const Treal&, const SeldonTranspose&, const Vector<Treal>&, const Treal&,
		      Vector<Treal>&, bool assemble) const;

    void MltAddVector(const Tcplx&, const SeldonTranspose&, const Vector<Tcplx>&, const Tcplx&,
		      Vector<Tcplx>&, bool assemble) const;

    void MltVector(const Vector<Treal>&, Vector<Treal>&) const;
    void MltVector(const Vector<Tcplx>&, Vector<Tcplx>&) const;

    void MltVector(const SeldonTranspose&, const Vector<Treal>&, Vector<Treal>&) const;
    void MltVector(const SeldonTranspose&, const Vector<Tcplx>&, Vector<Tcplx>&) const;

    void MltVector(const SeldonTranspose&, const Vector<Treal>&, Vector<Treal>&, bool assemble) const;
    void MltVector(const SeldonTranspose&, const Vector<Tcplx>&, Vector<Tcplx>&, bool assemble) const;

    template<class T0>
    void MltVectorGen(const SeldonTranspose& trans, const Vector<T0>& X, Vector<T0>& Y, bool assemble) const;
    
    template<class T0>
    void MltAddVectorGen(const T0& alpha, const SeldonTranspose& trans, const Vector<T0>& X,
			 const T0& beta, Vector<T0>& Y, bool assemble) const;
    
    virtual void MltAddFree(const GlobalGenericMatrix<T>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Treal>& X, Vector<Treal>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<T>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Tcplx>& X, Vector<Tcplx>& Y) const;

    virtual void ApplySor(const SeldonTranspose&, Vector<Treal>& x, const Vector<Treal>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    virtual void ApplySor(const SeldonTranspose&, Vector<Tcplx>& x, const Vector<Tcplx>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    template<class T0>
    friend void CopyReal(const FemMatrixFreeClass_Base<complex<T0> >& A, FemMatrixFreeClass_Base<T0>& B);
    
  };

  template<class T>
  void CopyReal(const FemMatrixFreeClass_Base<complex<T> >& A, FemMatrixFreeClass_Base<T>& B);

  template<class T>
  void CopyReal(const FemMatrixFreeClass_Base<T>& A, FemMatrixFreeClass_Base<T>& B);

  // Inline functions
  template<class T>
  bool IsSymmetricMatrix(const FemMatrixFreeClass_Base<T>& A);
  
  template<class T>
  bool IsComplexMatrix(const FemMatrixFreeClass_Base<T>& A);  

  template<class T>
  bool IsComplexMatrix(const FemMatrixFreeClass_Base<complex<T> >& A);  
  
  //! Finite element matrix
  template<class T, class TypeEquation>
  class FemMatrixFreeClass_Eq : public FemMatrixFreeClass_Base<T>
  {
  public :
    //! object containing datas on solved problem
    const EllipticProblem<TypeEquation>& var;
    
  protected :
    //! object containing additional arrays needed to perform the matrix-vector product
    ExtrapolVariablesProductFEM<T, TypeEquation> var_extra;

  public :    
    FemMatrixFreeClass_Eq(const EllipticProblem<TypeEquation>& var_);

    // inline methods
    ExtrapolVariablesProductFEM<T, TypeEquation>& GetExtrapolVariables();

    size_t GetMemorySize() const;
    
  };


  template<class T, class TypeEquation>
  class FemMatrixFreeClass : public FemMatrixFreeClass_Eq<T, TypeEquation>
  {
  public:
    FemMatrixFreeClass(const EllipticProblem<TypeEquation>& var_);
    
  };

  
  //! base class for mass matrix
  class FemMassMatrix
  {
  public :
    // differents type of matrix (diagonal, block-diagonal, sparse, matrix-free, etc)
    enum { DIAGONAL, BLOCK_DIAGONAL, BLOCK_ELT, MATRIX_SPARSE, MATRIX_FREE,
	   BLOCK_DIAGONAL_UNSYM, MATRIX_SPARSE_UNSYM, 
	   DIAG_SPARSE, BLOCK_DIAG_SPARSE, DIAG_MATRIX_FREE, BLOCK_DIAG_MATRIX_FREE,
           DIAG_SPARSE_UNSYM, BLOCK_DIAG_SPARSE_UNSYM};

    static int GetSymmetryType(int type);
    
  };
    
}

#define MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_HXX
#endif
