#ifndef SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX

#include "PolynomialEigenvalueSolver.hxx"

namespace Seldon
{
  SlepcParamPep::SlepcParamPep()
  {
    type_solver = TOAR;
  }

  int SlepcParamPep::GetEigensolverType() const
  {
    return type_solver;
  }

  void SlepcParamPep::SetEigensolverType(int type)
  {
    type_solver = type;
  }


  /*******************************
   * PolynomialEigenProblem_Base *
   *******************************/

  //! default constructor
  template<class T>
  PolynomialEigenProblem_Base<T>::PolynomialEigenProblem_Base()
  {
    use_spectral_transfo = false;
    pol_degree = 0;    
    diagonal_mass = false;
  }


  //! returns true if a spectral transformation has to be used
  template<class T>
  bool PolynomialEigenProblem_Base<T>::UseSpectralTransformation() const
  {
    return use_spectral_transfo;
  }

  
  //! enables a spectral transformation
  template<class T>
  void PolynomialEigenProblem_Base<T>::SetSpectralTransformation(bool t)
  {
    use_spectral_transfo = t;
  }

  
  //! returns object storing slepc parameters
  template<class T>
  SlepcParamPep& PolynomialEigenProblem_Base<T>::GetSlepcParameters()
  {
    return slepc_param;
  }

  
  //! returns object storing Feast parameters
  template<class T>
  FeastParam& PolynomialEigenProblem_Base<T>::GetFeastParameters()
  {
    return feast_param;

  }

  //! returns the polynomial degree
  template<class T>
  int PolynomialEigenProblem_Base<T>::GetPolynomialDegree() const
  {
    return pol_degree;
  }
  
  
  //! sets a diagonal mass
  template<class T>    
  void PolynomialEigenProblem_Base<T>::SetDiagonalMass(bool diag)
  {
    diagonal_mass = diag;
  }

  
  //! returns true if the mass is diagonal
  template<class T>
  bool PolynomialEigenProblem_Base<T>::DiagonalMass()
  {
    return diagonal_mass;
  }
  

  //! to overload
  template<class T>
  void PolynomialEigenProblem_Base<T>::ComputeOperator(int num, const Vector<T>& coef)
  {
    cout << "ComputeOperator not overloaded" << endl;
    abort();
  }

  
  //! to overload
  template<class T>
  void PolynomialEigenProblem_Base<T>::MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "MltOperator not overloaded" << endl;
    abort();
  }

  
  //! to overload
  template<class T>    
  void PolynomialEigenProblem_Base<T>::FactorizeMass()
  {
    cout << "FactorizeMass not overloaded" << endl;
    abort();
  }

  
  //! to overload for non-diagonal mass
  template<class T>    
  void PolynomialEigenProblem_Base<T>
  ::SolveMass(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      {
        cout << "SolveMass not overloaded" << endl;
        abort();
      }
  }

  
  template<class T>
  void PolynomialEigenProblem_Base<T>::FactorizeOperator(const Vector<T>& coef)
  {
    cout << "FactorizeOperator not overloaded" << endl;
    abort();
  }

  
  template<class T>
  void PolynomialEigenProblem_Base<T>
  ::SolveOperator(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    cout << "SolveOperator not overloaded" << endl;
    abort(); 
  }

  
  /**************************
   * PolynomialEigenProblem *
   **************************/

  
  //! inits the operators of the polynomial
  template<class T>
  void PolynomialEigenProblem<T>::InitMatrix(const Vector<VirtualMatrix<T>* >& op, int n)
  {
    if (n == -1)
      this->Init(op(0)->GetM());
    else
      this->Init(n);
    
    list_op = op;
    this->pol_degree = list_op.GetM()-1;
  }
  
  //! computes the operator with coefficients stored in coef
  template<class T>
  void PolynomialEigenProblem<T>::ComputeOperator(int num, const Vector<T>& coef)
  {
    // default implementation : we store the coefficients, no matrix stored
    if (num >= list_coef.GetM())
      list_coef.Resize(num+1);
    
    list_coef(num) = coef;
  }
  
  //! Computes Y = A X where A is the operator num
  template<class T>
  void PolynomialEigenProblem<T>::MltOperator(int num, const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (list_op.GetM() <= 0)
      return;
    
    if (list_coef.GetM() == 0)
      {
        // no coeffients, the operator num is stored in list_op(num)
        list_op(num)->MltVector(trans, X, Y);
        return;
      }

    // coefficients, the operator num is a linear combination of stored operators
    T zero; SetComplexZero(zero);
    T one; SetComplexOne(one);
    if (list_coef(num)(0) != zero)
      {
        list_op(0)->MltVector(trans, X, Y);
        Mlt(list_coef(num)(0), Y);
      }
    else
      Y.Zero();
    
    for (int i = 1; i < list_op.GetM(); i++)
      if (list_coef(num)(i) != zero)
        list_op(i)->MltAddVector(list_coef(num)(i), trans, X, one, Y);
  }


  //! returns true if the problem is symmetric
  template<class T>
  bool PolynomialEigenProblem<T>::IsSymmetricProblem() const
  {
    for (int i = 0; i < this->list_op.GetM(); i++)
      if (!this->list_op(i)->IsSymmetric())
        return false;

    return true;
  }
  

  //! returns true if the problem is hermitian
  template<class T>
  bool PolynomialEigenProblem<T>::IsHermitianProblem() const
  {
    // polynomial eigenvalue problem => hermitian not possible
    return false;
  }
  
  
  /*******************************
   * PolynomialDenseEigenProblem *
   *******************************/

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::InitMatrix(const Vector<Matrix<T, Prop, Storage>* >& op)
  {
    Vector<VirtualMatrix<T>* > op0(op.GetM());
    for (int i = 0; i < op0.GetM(); i++)
      op0(i) = op(i);
    
    PolynomialEigenProblem<T>::InitMatrix(op0);
    list_mat = op;
  }
  
  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::FactorizeMass()
  {    
    if (this->DiagonalMass())
      {
        T one; SetComplexOne(one);
        this->invDiag.Reallocate(this->n_);
        Matrix<T, Prop, Storage>& M = *list_mat(this->pol_degree);
        for (int i = 0; i < this->n_; i++)
          this->invDiag(i) = one / M(i, i);
      }
    else
      {
        mat_lu = *list_mat(this->pol_degree);
        GetLU(mat_lu, pivot);
      }
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>
  ::SolveMass(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (this->DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      {
        if (!trans.NoTrans())
          {
            cout << "Not implemented" << endl;
            abort();
          }
        
        y = x;
        SolveLU(mat_lu, pivot, y);
      }
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>::FactorizeOperator(const Vector<T>& coef)
  {
    T zero; SetComplexZero(zero);
    if (coef(0) == zero)
      {
        mat_lu.Reallocate(this->n_, this->n_);
        mat_lu.Zero();
      }
    else
      {
        mat_lu = *list_mat(0);
        Mlt(coef(0), mat_lu);
      }

    for (int k = 1; k < list_mat.GetM(); k++)
      if (coef(k) != zero)
        Add(coef(k), *list_mat(k), mat_lu);
    
    GetLU(mat_lu, pivot);
  }

  template<class T, class Prop, class Storage>
  void PolynomialDenseEigenProblem<T, Prop, Storage>
  ::SolveOperator(const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (!trans.NoTrans())
      {
        cout << "Not implemented" << endl;
        abort();
      }
    
    Y = X;
    SolveLU(mat_lu, pivot, Y);
  }
  
  
  /*********************************
   *  PolynomialSparseEigenProblem *
   *********************************/

  template<class T, class MatStiff, class MatMass>
  PolynomialSparseEigenProblem<T, MatStiff, MatMass>::PolynomialSparseEigenProblem()
  {
    Mh = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = nb_unknowns_scal_ = 0;
    nloc = 0;
  }
    
  template<class T, class MatStiff, class MatMass>
  int PolynomialSparseEigenProblem<T, MatStiff, MatMass>::RetrieveLocalNumbers(MatStiff& K)
  {
    nloc = K.GetM();
    
#ifdef SELDON_WITH_MPI
    try
      {
	DistributedMatrix_Base<typename MatStiff::entry_type>& A
	  = dynamic_cast<DistributedMatrix_Base<typename MatStiff::entry_type>& >(K);

	MPI_Comm comm = A.GetCommunicator();
	this->SetCommunicator(comm);
	int nb_proc;
	MPI_Comm_size(comm, &nb_proc);

	// only one processor => sequential case
	if (nb_proc <= 1)
	  return -1;
        
	// parallel case
	int m = A.GetLocalM();
	const IVect& OverlapRow = A.GetOverlapRowNumber();
	int noverlap = OverlapRow.GetM();
	int n = m - noverlap;
	local_col_numbers.Reallocate(n);
	Vector<bool> OverlappedRow(m); OverlappedRow.Fill(false);
	for (int i = 0; i < noverlap; i++)
	  OverlappedRow(OverlapRow(i)) = true;
	
	int ncol = 0;
	for (int i = 0; i < m; i++)
	  if (!OverlappedRow(i))
	    local_col_numbers(ncol++) = i;

	ProcSharingRows = &A.GetProcessorSharingRows();
	SharingRowNumbers = &A.GetSharingRowNumbers();
	nodl_scalar_ = A.GetNodlScalar();
	nb_unknowns_scal_ = A.GetNbScalarUnknowns();
	
	return n;
      }
    catch (const std::bad_cast&)
      {
	// a sequential matrix has been provided
	this->SetCommunicator(MPI_COMM_SELF);

	local_col_numbers.Clear();
      }
#endif

    return -1;
  }

  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::InitMatrix(Vector<MatStiff*>& K, MatMass& M)
  {
    int n = RetrieveLocalNumbers(*K(0));
    distributed = false;
    if (n >= 0)
      distributed = true;
    
    Vector<VirtualMatrix<T>* > list_op(K.GetM() + 1);
    for (int i = 0; i < K.GetM(); i++)
      list_op(i) = K(i);
    
    list_op(K.GetM()) = &M;
    PolynomialEigenProblem<T>::InitMatrix(list_op, n);
    
    Kh.Reallocate(K.GetM());
    for (int i = 0; i < K.GetM(); i++)
      Kh(i) = K(i);
    
    Mh = &M;
  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::FactorizeMass()
  {
    if (this->DiagonalMass())
      {
        T one; SetComplexOne(one);
        this->invDiag.Reallocate(nloc);
        for (int i = 0; i < nloc; i++)
          this->invDiag(i) = (*Mh)(i, i);
        
	// D is assembled for distributed matrices
	if (distributed)
	  {
#ifdef SELDON_WITH_MPI
	    Vector<T> M(this->invDiag);
	    AssembleVector(M, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
			   this->comm, nodl_scalar_, nb_unknowns_scal_, 15);
            
	    this->invDiag.Reallocate(local_col_numbers.GetM());
	    for (int i = 0; i < local_col_numbers.GetM(); i++)
	      this->invDiag(i) = one / M(local_col_numbers(i));
#endif
	  }
        else
          for (int i = 0; i < nloc; i++)
            this->invDiag(i) = one / this->invDiag(i);
      }
    else
      {
        mat_lu.Factorize(*Mh);
      }
  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>
  ::SolveMass(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (this->DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      this->SolveOperator(trans, x, y);

  }
  
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::FactorizeOperator(const Vector<T>& coef)
  {
    MatStiff A;
    T zero; SetComplexZero(zero);
    A = *Kh(0);
    Mlt(coef(0), A);
    for (int k = 1; k < Kh.GetM(); k++)
      if (coef(k) != zero)
        Add(coef(k), *Kh(k), A);
    
    Add(coef(Kh.GetM()), *Mh, A);
    
    mat_lu.Factorize(A);
  }

  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>
  ::SolveOperator(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> X(nloc);
    if (distributed)
      {
        X.Zero();
        for (int i = 0; i < this->n_; i++)
          X(local_col_numbers(i)) = x(i);                
      }
    else
      X = x;
    
    mat_lu.Solve(trans, X, false);  
    if (distributed)
      {
        for (int i = 0; i < this->n_; i++)
          y(i) = X(local_col_numbers(i));
      }
    else
      y = X;
  }

  //! Computes Y = A X where A is the operator num
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::MltOperator(int num, const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (Kh.GetM() <= 0)
      return;
    
    Vector<T> X, Y;
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
        X.Reallocate(nloc); Y.Reallocate(nloc);
        X.Zero();
        for (int i = 0; i < this->n_; i++)
          X(local_col_numbers(i)) = x(i);
        
        AssembleVector(X, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
          this->comm, nodl_scalar_, nb_unknowns_scal_, 17);	    
        
        if (this->list_coef.GetM() == 0)
          this->list_op(num)->MltVector(trans, X, Y);
        else
          {
            // coefficients, the operator num is a linear combination of stored operators
            T zero; SetComplexZero(zero);
            T one; SetComplexOne(one);
            if (this->list_coef(num)(0) != zero)
              {
                this->list_op(0)->MltVector(trans, X, Y);
                Mlt(this->list_coef(num)(0), Y);
              }
            else
              Y.Zero();
            
            for (int i = 1; i < this->list_op.GetM(); i++)
              if (this->list_coef(num)(i) != zero)
                this->list_op(i)->MltAddVector(this->list_coef(num)(i), trans, X, one, Y);
          }
        
        for (int i = 0; i < this->n_; i++)
          y(i) = Y(local_col_numbers(i));
#endif
      }
    else
      PolynomialEigenProblem<T>::MltOperator(num, trans, x, y);
  }


  //! changes final eigenvectors if needed
  template<class T, class MatStiff, class MatMass>
  void PolynomialSparseEigenProblem<T, MatStiff, MatMass>::DistributeEigenvectors(Matrix<T, General, ColMajor>& eigen_vec)
  {
    if (distributed)
      {
#ifdef SELDON_WITH_MPI
        this->AssembleEigenvectors(eigen_vec, local_col_numbers,  ProcSharingRows,
                                   SharingRowNumbers, nloc, nodl_scalar_, nb_unknowns_scal_);
#endif
      }
  }


  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(PolynomialEigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec,
				  int type_solver)
  {
    if (type_solver == TypeEigenvalueSolver::DEFAULT)
      type_solver = TypeEigenvalueSolver::GetDefaultPolynomialSolver();

    if (type_solver == TypeEigenvalueSolver::FEAST)
      {
#ifdef SELDON_WITH_FEAST
        T zero; SetComplexZero(zero);
        Matrix<T, General, ColMajor> eigen_old;
        FindEigenvaluesFeast(var_eig, lambda, lambda_imag, eigen_old);
        
        // eigenvalues are sorted by ascending order
        SortEigenvalues(lambda, lambda_imag, eigen_old,
                        eigen_vec, var_eig.LARGE_EIGENVALUES,
                        var_eig.GetTypeSorting(), zero, zero);
#else
        cout << "Recompile with MKL or Feast" << endl;
        abort();
#endif
      }
    else if (type_solver == TypeEigenvalueSolver::SLEPC)
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

}
 
#define SELDON_FILE_POLYNOMIAL_EIGENVALUE_SOLVER_CXX
#endif
