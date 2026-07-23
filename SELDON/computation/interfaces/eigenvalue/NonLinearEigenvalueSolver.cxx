#ifndef SELDON_FILE_NON_LINEAR_EIGENVALUE_SOLVER_CXX

#include "NonLinearEigenvalueSolver.hxx"

namespace Seldon
{

#ifdef SELDON_WITH_SLEPC

  SlepcParamNep::SlepcParamNep()
  {
    type_solver = NLEIGS;
    use_command_line_options = false;
    a_sp = 0.1; b_sp = 10.0; c_sp = 0.1; d_sp = 10.0;
    full_basis = false;
    locking_variant = true;
    restart = 0.5;
    tol_interpol = -1.0;
    degree_interpol = -1;
    use_default_petsc_solver = false;
  }

  int SlepcParamNep::GetEigensolverType() const
  {
    return type_solver;
  }

  void SlepcParamNep::SetEigensolverType(int type)
  {
    type_solver = type;
  }

  double SlepcParamNep::GetLrMin() const
  {
    return a_sp;
  }
  
  double SlepcParamNep::GetLrMax() const
  {
    return b_sp;
  }
  
  double SlepcParamNep::GetLiMin() const
  {
    return c_sp;
  }
  
  double SlepcParamNep::GetLiMax() const
  {
    return d_sp;
  }
  
  void SlepcParamNep::SetIntervalRegion(double a, double b, double c, double d)
  {
    a_sp = a; b_sp = b; c_sp = c; d_sp = d;
  }

  bool SlepcParamNep::InsideRegion(const complex<double>& z) const
  {
    if ((real(z) >= a_sp) && (real(z) <= b_sp)
        && (imag(z) >= c_sp) && (imag(z) <= d_sp))
      return true;

    return false;
  }
  
  void SlepcParamNep::EnableCommandLineOptions(bool yes)
  {
    use_command_line_options = yes;
  }
  
  bool SlepcParamNep::UseCommandLineOptions() const
  {
    return use_command_line_options;
  }
  
  void SlepcParamNep::SetDefaultPetscSolver(bool yes)
  {
    use_default_petsc_solver = yes;
  }

  bool SlepcParamNep::UseDefaultPetscSolver() const
  {
    return use_default_petsc_solver;
  }

  void SlepcParamNep::SetFullBasis(bool yes)
  {
    full_basis = yes;
  }
  
  bool SlepcParamNep::FullBasis() const
  {
    return full_basis;
  }

  void SlepcParamNep::SetLockingVariant(bool yes)
  {
    locking_variant = yes;
  }
  
  bool SlepcParamNep::LockingVariant() const
  {
    return locking_variant;
  }

  void SlepcParamNep::SetInterpolationDegree(int p)
  {
    degree_interpol = p;
  }
  
  int SlepcParamNep::GetInterpolationDegree() const
  {
    return degree_interpol;
  }

  void SlepcParamNep::SetInterpolationTolerance(double eps)
  {
    tol_interpol = eps;
  }
  
  double SlepcParamNep::GetInterpolationTolerance() const
  {
    return tol_interpol;
  }
    
  void SlepcParamNep::SetRestartNleigs(double r)
  {
    restart = r;
  }
  
  double SlepcParamNep::GetRestartNleigs() const
  {
    return restart;
  }

  void SlepcParamNep::SetRKShifts(const Vector<Petsc_Scalar>& s)
  {
    shifts = s;
  }
  
  const Vector<Petsc_Scalar>& SlepcParamNep::GetRKShifts() const
  {
    return shifts;
  }

  
  /******************************
   * NonLinearEigenProblem_Base *
   ******************************/

  
  //! default constructor
  template<class T>
  NonLinearEigenProblem_Base<T>::NonLinearEigenProblem_Base()
  {
    exact_preconditioning = false;
    explicit_matrix = false;
    nb_split_matrix = 0;
    this->automatic_selection_arnoldi_vectors = false; // we let Slepc decide
    this->nb_arnoldi_vectors = 0;
    rational_problem = true;
  }

  
  //! returns object storing slepc parameters
  template<class T>
  SlepcParamNep& NonLinearEigenProblem_Base<T>::GetSlepcParameters()
  {
    return slepc_param;
  }
  

  //! specifies that the preconditioning is exact (solves exactly T(L) x = b)
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetExactPreconditioning(bool yes)
  {
    exact_preconditioning = yes;
  }


  //! returns true if the preconditioning is exact
  template<class T>
  bool NonLinearEigenProblem_Base<T>::ExactPreconditioning() const
  {
    return exact_preconditioning;
  }
  

  //! returns the poles of T(L)
  template<class T>
  const Vector<T>& NonLinearEigenProblem_Base<T>::GetSingularities() const
  {
    return singular_points;
  }

  
  //! sets the poles of T(L)
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetSingularities(const Vector<T>& s)
  {
    singular_points = s;
  }

  
  //! specifies that the matrix T(L) must be stored
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetExplicitMatrix(bool yes)
  {
    explicit_matrix = yes;
  }

  
  //! returns true if the matrix T(L) is explicitely computed
  template<class T>
  bool NonLinearEigenProblem_Base<T>::ExplicitMatrix() const
  {
    return explicit_matrix;
  }

  
  //! sets the number of split matrices (T(L) = \sum_k f_k(L) A_k)
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetSplitMatrices(int n)
  {
    nb_split_matrix = n;
    numer_pol_split.Reallocate(n);
    denom_pol_split.Reallocate(n);
  }

  
  //! returns true if the split form T(L) = \sum_k f_k(L) A_k is used
  template<class T>
  bool NonLinearEigenProblem_Base<T>::UseSplitMatrices() const
  {
    if (nb_split_matrix > 1)
      return true;

    return false;
  }

  
  //! returns the number of split matrices A_k (T(L) = \sum_k f_k(L) A_k)
  template<class T>
  int NonLinearEigenProblem_Base<T>::GetNbSplitMatrices() const
  {
    return nb_split_matrix;
  }

  
  //! sets the numerator of rational function f_i(L)
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetNumeratorSplitFct(int i, const Vector<T>& coef)
  {
    rational_problem = true;
    numer_pol_split(i) = coef;
  }

  
  //! sets the denominator of rational function f_i(L)
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetDenominatorSplitFct(int i, const Vector<T>& coef)
  {
    rational_problem = true;
    denom_pol_split(i) = coef;
  }

  
  //! returns the numerator of rational function f_i(L)
  template<class T>
  const Vector<T>& NonLinearEigenProblem_Base<T>::GetNumeratorSplitFct(int i)
  {
    return numer_pol_split(i);
  }

  
  //! returns the denominator of rational function f_i(L)
  template<class T>
  const Vector<T>& NonLinearEigenProblem_Base<T>::GetDenominatorSplitFct(int i)
  {
    return denom_pol_split(i);
  }


  //! returns true if the problem is rational
  template<class T>
  bool NonLinearEigenProblem_Base<T>::RationalSplitFunctions()
  {
    return rational_problem;
  }


  //! returns the vector of functions to use for non-rational problems
  template<class T>
  const Vector<FN>& NonLinearEigenProblem_Base<T>::GetSlepcFunction_NonRational() const
  {
    return non_rational_fct;
  }


  //! sets the vector of functions to use for non-rational problems
  template<class T>
  void NonLinearEigenProblem_Base<T>::SetSlepcFunction_NonRational(Vector<FN>& fct)
  {
    rational_problem = false;
    non_rational_fct = fct;
  }
  
  
  //! checks if the value of L is valid
  template<class T>
  void NonLinearEigenProblem_Base<T>::CheckValueL(const T& L)
  {
    if (isinf(abs(L)) || isnan(abs(L)))
      {
        cout << "Invalid value for L : " << L << endl;
        abort();
      }
  }  
  

  //! computes the operator T(L) (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeOperator(const T& L)
  {
    cout << "ComputeOperator not overloaded" << endl;
    abort();
  }

  
  //! computes Y = T(L) X or T(L)^T X (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::MltOperator(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltOperator not overloaded" << endl;
    abort();
  }

  
  //! computes explicitely T(L) stored in A (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeOperatorExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    cout << "ComputeOperatorExplicit not overloaded" << endl;
    abort();

  }
  
  //! computes the operator T'(L) (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeJacobian(const T& L)
  {
    cout << "ComputeJacobian not overloaded" << endl;
    abort();
  }
  
  //! computes Y = T'(L) X or T'(L)^T X (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::MltJacobian(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltJacobian not overloaded" << endl;
    abort();
  }

  //! computes explicitely T'(L) stored in A (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeJacobianExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    cout << "ComputeOperatorExplicit not overloaded" << endl;
    abort();

  }

  //! computes the preconditioning used to solve T(L) x = b (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputePreconditioning(const T& L)
  {
    cout << "ComputePreconditioning not overloaded" << endl;
    abort();
  }

  //! computes the preconditioning used to solve \sum_k coef(k) T(L_k) x = b (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputePreconditioning(const Vector<T>& L, const Vector<T>& coef)
  {
    cout << "ComputePreconditioning not overloaded" << endl;
    abort();
  }


  //! computes the preconditioning used to solve \sum_k coef(k) A_{num(k)} x = b (to overload)
  //! where A_k are matrices of the split formulation
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeSplitPreconditioning(const Vector<int>& num, const Vector<T>& coef)
  {
    cout << "ComputeSplitPreconditioning not overloaded" << endl;
    abort();
  }

  
  //! computes the preconditioning used to solve A x = b (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ComputeExplicitPreconditioning(DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    cout << "ComputeExplicitPreconditioning not overloaded" << endl;
    abort();
  }
  
  //! applies preconditioning for the solution of A x = b (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>::ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "ApplyPreconditioning not overloaded" << endl;
    abort();
  }

  //! Computes matrix A_i of the split formulation T(L) = \sum f_i(L) A_i (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>
  ::ComputeOperatorSplitExplicit(int i, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    cout << "ComputeOperatorSplitExplicit not overloaded" << endl;
    abort();
  }

  //! Computes Y = A_i X where A_i is a matrix of the split formulation T(L) = \sum f_i(L) A_i (to overload)
  template<class T>
  void NonLinearEigenProblem_Base<T>
  ::MltOperatorSplit(int i, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltOperatorSplit not overloaded" << endl;
    abort();
  }

  
  /************************************
   * SplitSparseNonLinearEigenProblem *
   ************************************/

  //! checks if matrices Ai have been provided
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>::CheckPresenceMatrix() const
  {
    if (vec_Ai == NULL)
      {
        cout << "Matrices A_i not initialized" << endl;
        cout << "Did you call InitMatrix ? " << endl;
        abort();
      }
  }


  //! default constructor
  template<class T, class Prop, class Storage>
  SplitSparseNonLinearEigenProblem<T, Prop, Storage>::SplitSparseNonLinearEigenProblem()
  {
    distributed = false;
    vec_Ai = NULL; nloc = 0;
    this->exact_preconditioning = true; // direct solver => preconditioning is exact
  }
  
  
  //! inits the problem by providing matrices A_i and functions f_i
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::InitMatrix(Vector<DistributedMatrix<T, Prop, Storage> >& A,
               const Vector<Vector<T> >& numer, const Vector<Vector<T> >& denom)
  {
    this->SetSplitMatrices(A.GetM());
    this->vec_Ai = &A;
    this->numer_pol_split = numer;
    this->denom_pol_split = denom;
    this->mat_lu.Clear();
    
    int nb_proc(1);
#ifdef SELDON_WITH_MPI
    this->SetCommunicator(A(0).GetCommunicator());
    MPI_Comm_size(this->comm, &nb_proc);
#endif
    
    if (nb_proc <= 1)
      {
        distributed = false;
        local_col_numbers.Clear();
        this->Init(A(0).GetM());
      }
    else
      {
#ifdef SELDON_WITH_MPI
        distributed = true;
        // checking that all matrices have the same rows
        for (int i = 1; i < A.GetM(); i++)
          if (!A(0).SameDistributedRows(A(i)))
            {
              cout << "Matrices A_i must have the same distributed rows" << endl;
              cout << "Matrix " << i << " is different" << endl;
              abort();
            }

        // counting the number of local columns
        int m = A(0).GetLocalM(); nloc = m;
        const IVect& OverlapRow = A(0).GetOverlapRowNumber();
	int noverlap = OverlapRow.GetM();
	int n = m - noverlap;

        // filling local_col_numbers
        local_col_numbers.Reallocate(n);
	Vector<bool> OverlappedRow(m); OverlappedRow.Fill(false);
	for (int i = 0; i < noverlap; i++)
	  OverlappedRow(OverlapRow(i)) = true;
	
	int ncol = 0;
	for (int i = 0; i < m; i++)
	  if (!OverlappedRow(i))
	    local_col_numbers(ncol++) = i;

        // other data needed to assemble vectors
        ProcSharingRows = &A(0).GetProcessorSharingRows();
	SharingRowNumbers = &A(0).GetSharingRowNumbers();
	nodl_scalar_ = A(0).GetNodlScalar();
	nb_unknowns_scal_ = A(0).GetNbScalarUnknowns();
        
        this->Init(n);
#endif
      }
  }


#ifdef SELDON_WITH_MPI
  //! returns the direct solver used as preconditioning
  template<class T, class Prop, class Storage>
  SparseDistributedSolver<T>& SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::GetDirectSolver()
  {
    return mat_lu;
  }
#else
  //! returns the direct solver used as preconditioning
  template<class T, class Prop, class Storage>
  SparseDirectSolver<T>& SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::GetDirectSolver()
  {
    return mat_lu;
  }
#endif

  //! returns true if all the matrices A_i are symmetric
  template<class T, class Prop, class Storage>
  bool SplitSparseNonLinearEigenProblem<T, Prop, Storage>::IsSymmetricProblem() const
  {
    CheckPresenceMatrix();
    for (int i = 0; i < vec_Ai->GetM(); i++)
      if (!(*vec_Ai)(i).IsSymmetric())
        return false;

    return true;
  }

  
  //! computes preconditioning to solve \sum coef_i A_(numL(i)) x = b
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::ComputeSplitPreconditioning(const Vector<int>& numL, const Vector<T>& coef)
  {
    DistributedMatrix<T, Prop, Storage> A;
    A = (*vec_Ai)(numL(0)); Mlt(coef(0), A);
    for (int i = 1; i < numL.GetM(); i++)
      Add(coef(i), (*vec_Ai)(numL(i)), A);

    if (this->print_level >= 2)
      mat_lu.ShowMessages();
    
    mat_lu.Factorize(A);
    mat_lu.HideMessages();
  }


  //! computes preconditioning to solve A x = b
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::ComputeExplicitPreconditioning(DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    if (distributed)
      {
        cout << "not implemented" << endl;
        abort();
      }

    if (this->print_level >= 2)
      mat_lu.ShowMessages();

    mat_lu.Factorize(A);
    mat_lu.HideMessages();
  }


  //! applies preconditioning to vector X
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::ApplyPreconditioning(const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (distributed)
      {
        Vector<T> Xd(nloc); Xd.Zero();
        for (int i = 0; i < this->n_; i++)
	  Xd(local_col_numbers(i)) = X(i);

        mat_lu.Solve(trans, Xd, false);
        
        for (int i = 0; i < this->n_; i++)
	  Y(i) = Xd(local_col_numbers(i));
      }
    else
      {
        Y = X;
        mat_lu.Solve(trans, Y);
      }
  }


  //! Computes T(L) explicitely in matrix A
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::ComputeOperatorSplitExplicit(int i, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    if (distributed)
      {
        cout << "not implemented" << endl;
        abort();
      }

    Copy((*vec_Ai)(i), A);
  }


  //! Computes Y = A_i X or Y = A_i^T X
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>
  ::MltOperatorSplit(int i, const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
        Vector<T> Xd(nloc), Yd(nloc); Xd.Zero();
        for (int i = 0; i < this->n_; i++)
	  Xd(local_col_numbers(i)) = X(i);

        AssembleVector(Xd, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
                       this->comm, nodl_scalar_, nb_unknowns_scal_, 18);	    
        
        Mlt(trans, (*vec_Ai)(i), Xd, Yd);
        
        for (int i = 0; i < this->n_; i++)
	  Y(i) = Yd(local_col_numbers(i));

#endif
      }
    else
      Mlt(trans, (*vec_Ai)(i), X, Y);
  }

  
  //! changes final eigenvectors if needed
  template<class T, class Prop, class Storage>
  void SplitSparseNonLinearEigenProblem<T, Prop, Storage>::DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec)
  {
#ifdef SELDON_WITH_MPI
    if (distributed)
      {
        this->AssembleEigenvectors(eigen_vec, local_col_numbers,  ProcSharingRows,
                                   SharingRowNumbers, nloc, nodl_scalar_, nb_unknowns_scal_);
      }
#endif
  }

  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(NonLinearEigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver)
  {
    if (type_solver == TypeEigenvalueSolver::DEFAULT)
      type_solver = TypeEigenvalueSolver::GetDefaultNonLinearSolver();

    if (type_solver == TypeEigenvalueSolver::SLEPC)
      {
#ifdef SELDON_WITH_SLEPC
	T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
      }
    else
      {
        cout << "Recompile with eigenvalue solver" << endl;
        abort();
      }
    
  }

#endif
  
}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX
#endif
