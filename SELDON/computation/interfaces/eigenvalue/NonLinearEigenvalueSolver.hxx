#ifndef SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_HXX

#ifdef SELDON_WITH_SLEPC
#include <slepcnep.h>
#endif

namespace Seldon
{
  
#ifdef SELDON_WITH_SLEPC
  
  //! Parameters for Slepc package
  class SlepcParamNep
  {
  protected:
    //! which solver ?
    int type_solver;

    //! do we use parameters of the command line
    bool use_command_line_options;

    //! region of the spectrum
    double a_sp, b_sp, c_sp, d_sp;

    // parameters for nleigs
    bool full_basis, locking_variant;
    int degree_interpol;
    double tol_interpol, restart;

    // do we use default solver proposed by Petsc ?
    bool use_default_petsc_solver;

    // shifts for nleigs
    Vector<Petsc_Scalar> shifts;
    
  public :
    enum {RII, SLP, NARNOLDI, CISS, INTERPOL, NLEIGS};
    
    SlepcParamNep();
    
    int GetEigensolverType() const;
    void SetEigensolverType(int type);

    double GetLrMin() const;
    double GetLrMax() const;
    double GetLiMin() const;
    double GetLiMax() const;
    
    void SetIntervalRegion(double, double, double, double);
    bool InsideRegion(const complex<double>& z) const;
    
    void EnableCommandLineOptions(bool yes = true);
    bool UseCommandLineOptions() const;

    void SetDefaultPetscSolver(bool yes = true);
    bool UseDefaultPetscSolver() const;
    
    void SetFullBasis(bool yes = true);
    bool FullBasis() const;

    void SetLockingVariant(bool yes = true);
    bool LockingVariant() const;

    void SetInterpolationDegree(int);
    int GetInterpolationDegree() const;

    void SetInterpolationTolerance(double);
    double GetInterpolationTolerance() const;
    
    void SetRestartNleigs(double);
    double GetRestartNleigs() const;

    void SetRKShifts(const Vector<Petsc_Scalar>& s);
    const Vector<Petsc_Scalar>& GetRKShifts() const;
    
  };

  //! Base class for non-linear eigenvalue solver
  template<class T>
  class NonLinearEigenProblem_Base : public GeneralEigenProblem<T>
  {
  public:
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

  protected:        
    //! Slepc parameters
    SlepcParamNep slepc_param;

    //! singularities
    Vector<T> singular_points;

    bool exact_preconditioning;
    bool explicit_matrix;
    int nb_split_matrix;
    Vector<Vector<T> > numer_pol_split, denom_pol_split;
    bool rational_problem;
    Vector<FN> non_rational_fct;
    
  public:
    NonLinearEigenProblem_Base();

    SlepcParamNep& GetSlepcParameters();
    
    void SetExactPreconditioning(bool yes = true);
    bool ExactPreconditioning() const;
    
    const Vector<T>& GetSingularities() const;
    void SetSingularities(const Vector<T>& s);

    void SetExplicitMatrix(bool yes = true);
    bool ExplicitMatrix() const;
    
    void SetSplitMatrices(int);
    bool UseSplitMatrices() const;
    int GetNbSplitMatrices() const;

    void SetNumeratorSplitFct(int i, const Vector<T>&);
    void SetDenominatorSplitFct(int i, const Vector<T>&);
    const Vector<T>& GetNumeratorSplitFct(int i);
    const Vector<T>& GetDenominatorSplitFct(int i);

    bool RationalSplitFunctions();
    const Vector<FN>& GetSlepcFunction_NonRational() const;
    void SetSlepcFunction_NonRational(Vector<FN>& fct);
    
    void CheckValueL(const T& L);
    
    virtual void ComputeOperator(const T& L);
    virtual void MltOperator(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    virtual void ComputeOperatorExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A);
    
    virtual void ComputeJacobian(const T& L);
    virtual void MltJacobian(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    virtual void ComputeJacobianExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A);
    
    virtual void ComputePreconditioning(const T& L);
    virtual void ComputePreconditioning(const Vector<T>& L, const Vector<T>& coef);
    virtual void ComputeSplitPreconditioning(const Vector<int>& numL, const Vector<T>& coef);
    virtual void ComputeExplicitPreconditioning(DistributedMatrix<T, General, ArrayRowSparse>& A);
    virtual void ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);

    virtual void ComputeOperatorSplitExplicit(int i, DistributedMatrix<T, General, ArrayRowSparse>& A);
    virtual void MltOperatorSplit(int i, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
  };

  //! class solving non-linear eigenvalue with a split form
  /*!
    We consider a non-linear eigenvalue problem T(L) X = 0
    where T(L) can be written in the following split form :
    T(L) = \sum f_i(L) A_i
    where A_i are sparse matrices and f_i rational functions
   */
  template<class T, class Prop, class Storage>
  class SplitSparseNonLinearEigenProblem : public NonLinearEigenProblem_Base<T>
  {
    bool distributed; int nloc;
    Vector<int> local_col_numbers;
    Vector<int>* ProcSharingRows;
    Vector<Vector<int> >* SharingRowNumbers;
    int nodl_scalar_, nb_unknowns_scal_;
    
    
  protected:
    //! matrices A_i
    Vector<DistributedMatrix<T, Prop, Storage> >* vec_Ai;

#ifdef SELDON_WITH_MPI
    //! LU factorization of sparse matrix
    SparseDistributedSolver<T> mat_lu;
#else
    //! LU factorization of sparse matrix
    SparseDirectSolver<T> mat_lu;
#endif
    
    void CheckPresenceMatrix() const;
    
  public:
    SplitSparseNonLinearEigenProblem();

    void InitMatrix(Vector<DistributedMatrix<T, Prop, Storage> >& A,
                    const Vector<Vector<T> >& numer, const Vector<Vector<T> >& denom);

#ifdef SELDON_WITH_MPI
    SparseDistributedSolver<T>& GetDirectSolver();
#else
    SparseDirectSolver<T>& GetDirectSolver();
#endif
    
    bool IsSymmetricProblem() const;
    
    void ComputeSplitPreconditioning(const Vector<int>& numL, const Vector<T>& coef);
    void ComputeExplicitPreconditioning(DistributedMatrix<T, General, ArrayRowSparse>& A);
    
    void ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);

    void ComputeOperatorSplitExplicit(int i, DistributedMatrix<T, General, ArrayRowSparse>& A);
    
    void MltOperatorSplit(int i, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);

    void DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec);
    
  };

  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(NonLinearEigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver = TypeEigenvalueSolver::DEFAULT);

#endif
  
}
 
#define SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_HXX
#endif
