#ifndef SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_HXX

namespace Seldon
{
  //! Parameters for Slepc package
  class SlepcParamPep
  {
  protected:
    //! which solver ?
    int type_solver;

  public :
    enum {TOAR, STOAR, QARNOLDI, LINEAR, JD};

    SlepcParamPep();
    
    int GetEigensolverType() const;
    void SetEigensolverType(int type);

  };

  
  //! Base class for polynomial eigenvalue solver
  template<class T>
  class PolynomialEigenProblem_Base : public GeneralEigenProblem<T>
  {
  public:
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

  protected:
    //! if true teta is searched (lambda = shift + 1/teta)
    bool use_spectral_transfo;
    
    //! polynomial degree
    int pol_degree;
    
    //! Slepc parameters
    SlepcParamPep slepc_param;
    FeastParam feast_param;
    
    //! mass diagonal ?
    bool diagonal_mass;
    Vector<T> invDiag;
    
  public:
    PolynomialEigenProblem_Base();

    bool UseSpectralTransformation() const;
    void SetSpectralTransformation(bool t = true);
    
    SlepcParamPep& GetSlepcParameters();
    FeastParam& GetFeastParameters();
    
    int GetPolynomialDegree() const;
    
    void SetDiagonalMass(bool diag = true);
    bool DiagonalMass();

    // virtual methods overloaded in leaf classes
    virtual void ComputeOperator(int num, const Vector<T>& coef);
    virtual void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    virtual void FactorizeMass();
    virtual void SolveMass(const SeldonTranspose&, const Vector<T>& x, Vector<T>& y);
    
    virtual void FactorizeOperator(const Vector<T>& coef);
    virtual void SolveOperator(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };

  //! matrix-free implementation of polynomial eigenvalue solver
  template<class T>
  class PolynomialEigenProblem : public PolynomialEigenProblem_Base<T>
  {
  protected:
    Vector<VirtualMatrix<T>* > list_op;
    Vector<Vector<T> > list_coef;
    
  public:
    void InitMatrix(const Vector<VirtualMatrix<T>* >& op, int n = -1);
    
    void ComputeOperator(int num, const Vector<T>& coef);
    void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;

  };

  //! implementation of polynomial eigenvalue solver for dense problem
  template<class T, class Prop, class Storage>
  class PolynomialDenseEigenProblem : public PolynomialEigenProblem<T>
  {
  protected:
    Vector<Matrix<T, Prop, Storage>* > list_mat;
    Matrix<T, Prop, Storage> mat_lu;
    Vector<int> pivot;
    
  public:
    void InitMatrix(const Vector<Matrix<T, Prop, Storage>* >& op);
    
    void FactorizeMass();
    void SolveMass(const SeldonTranspose&, const Vector<T>& x, Vector<T>& y);

    void FactorizeOperator(const Vector<T>& coef);
    void SolveOperator(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };

  //! implementation of sparse polynomial eigenvalue solver
#ifdef SELDON_WITH_MPI
  template<class T, class MatStiff,
           class MatMass = DistributedMatrix<T, Symmetric, ArrayRowSymSparse> >
#else
  template<class T, class MatStiff,
           class MatMass = Matrix<T, Symmetric, ArrayRowSymSparse> >
#endif
  class PolynomialSparseEigenProblem : public PolynomialEigenProblem<T>
  {
  protected:

#ifdef SELDON_WITH_MPI
    //! LU factorization of sparse matrix
    SparseDistributedSolver<T> mat_lu;
#else
    //! LU factorization of sparse matrix
    SparseDirectSolver<T> mat_lu;
#endif
    
    // stiffness matrices and mass matrix
    Vector<MatStiff*> Kh;
    MatMass* Mh;
    
    // parallel stuff
    bool distributed; int nloc;
    Vector<int> local_col_numbers;
    Vector<int>* ProcSharingRows;
    Vector<Vector<int> >* SharingRowNumbers;
    int nodl_scalar_, nb_unknowns_scal_;

    int RetrieveLocalNumbers(MatStiff& K);
    
  public:
    PolynomialSparseEigenProblem();
    
    void InitMatrix(Vector<MatStiff*>&, MatMass&);
    
    void FactorizeMass();
    void SolveMass(const SeldonTranspose&, const Vector<T>& x, Vector<T>& y);
    
    void FactorizeOperator(const Vector<T>& coef);
    void SolveOperator(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);

    void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);

    void DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec);
    
  };

  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(PolynomialEigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver = TypeEigenvalueSolver::DEFAULT);
}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_HXX
#endif
