#ifndef MONTJOIE_FILE_EIGENVALUES_HARMONIC_CXX

#include "EigenvaluesHarmonic.hxx"

namespace Montjoie
{
 
  /************************
   * EigenProblemMontjoie *
   ************************/
  

  //! Sets values of attributes to default values
  template<class T>
  void EigenProblemMontjoie<T>::InitDefaultValues()
  {
    eigenvalue_solver_used = TypeEigenvalueSolver::GetDefaultSolver();
    file_name_eigenval = "Omega.dat";
    type_solver_chol = -1;
  }
  
  
  //! Basic initializations before eigenvalue computation
  template<class T>
  void EigenProblemMontjoie<T>::InitComputation()
  {
    Vector<bool> diag_elt;
    int type_mass_matrix = var_computation.GetMassMatrixType(diag_elt);

    bool diag_m = true;
    if (!IsComplexNumber(T(0)))
      {
	if (this->shift_imag != T(0))
	  diag_m = false;
      }
    
    if (diag_m && (type_mass_matrix == FemMassMatrix::DIAGONAL))
      this->SetDiagonalMass();
    
    if (!var_problem.IsSymmetricMassMatrix() && !this->DiagonalMass())
      {
        // for non-symmetric mass matrix, only invert mode should be used
        this->eigenvalue_computation_mode = this->INVERT_MODE;
      }
    
    this->print_level = var_problem.print_level/2;     
    
    // constructing local_col_numbers
    int m = var_problem.GetNbDof();
    int ndir = var_boundary.GetNbDirichletDof();
    
#ifdef SELDON_WITH_MPI
    this->SetCommunicator(var_comm.comm_group_mode);    
    
    int noverlap = var_comm.GetNbOverlappedDof();
    Vector<bool> OverlappedRow(m); OverlappedRow.Fill(false);
    for (int i = 0; i < noverlap; i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;
    
    int noverlap2 = noverlap;
    if (ndir > 0)
      for (int i = 0; i < m; i++)
        if (var_boundary.IsDofDirichlet(i) && OverlappedRow(i))
          noverlap2--;
#else
    int noverlap2 = 0;
#endif
    
    int n = m - noverlap2 - ndir;
    local_col_numbers.Reallocate(n);

    int ncol = 0;
#ifdef SELDON_WITH_MPI    
    if (ndir > 0)
      {
        for (int i = 0; i < m; i++)
          if (!var_boundary.IsDofDirichlet(i) && !OverlappedRow(i))
            local_col_numbers(ncol++) = i;
      }
    else
      {
        for (int i = 0; i < m; i++)
          if (!OverlappedRow(i))
            local_col_numbers(ncol++) = i;
      }
#else
    if (ndir > 0)
      {
        for (int i = 0; i < m; i++)
          if (!var_boundary.IsDofDirichlet(i))
            local_col_numbers(ncol++) = i;
      }
#endif
    
    this->complex_system = false;    
    this->Init(n);
    
    /*T one, zero;
    SetComplexZero(zero); SetComplexOne(one);
    
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(one);
    nat_mat.SetCoefStiffness(zero);
    nat_mat.SetCoefDamping(zero);

    DistributedMatrix<T, General, ArrayRowSparse> K, M;
    var_computation.AddMatrixWithBC(M, nat_mat);
    
    nat_mat.SetCoefMass(zero);
    nat_mat.SetCoefStiffness(one);
    nat_mat.SetCoefDamping(one);
    var_computation.AddMatrixWithBC(K, nat_mat);
    
    K.WriteText("Kh.dat");
    M.WriteText("Mh.dat");*/
  }   
    

  template<class T>
  size_t EigenProblemMontjoie<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille = local_col_numbers.GetMemorySize();
    taille += chol_facto_mass_matrix.GetMemorySize();
    taille = Xchol_real.GetMemorySize();
    taille = Xchol_imag.GetMemorySize();
    return taille;
  }


  //! Sets eigenproblem to solve
  template<class T>
  void EigenProblemMontjoie<T>
  ::InitMatrix(FemMatrixFreeClass_Base<T>& K, FemMatrixFreeClass_Base<T>& M)
  {
    mat_mass = &M;
    mat_stiff = &K;
    this->SetMatrix(K, M);
    InitComputation();
  }
  
  
  //! modifying parameters with data file
  template<class T>
  void EigenProblemMontjoie<T>
  ::SetInputData(const string& description_field, const Vector<string>& parameters)
  {
    if (!description_field.compare("EigenvalueTolerance"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "EigenvalueTolerance needs more parameters, for instance :" << endl;
	    cout << "EigenvalueTolerance = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        this->stopping_criterion = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("EigenvalueMaxNumberIterations"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "EigenvalueMaxNumberIterations needs more parameters, for instance :" << endl;
	    cout << "EigenvalueMaxNumberIterations = n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        this->nb_maximum_iterations = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("UseCholeskyForEigenvalue"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "UseCholeskyForEigenvalue needs more parameters, for instance :" << endl;
	    cout << "UseCholeskyForEigenvalue = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          this->SetCholeskyFactoForMass(true);
        else
          this->SetCholeskyFactoForMass(false);
      }
    else if (!description_field.compare("FileEigenvalue"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "FileEigenvalue needs more parameters, for instance :" << endl;
	    cout << "FileEigenvalue = file_name_eigenval" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        file_name_eigenval = parameters(0);
      }
    else if (!description_field.compare("EigenvalueSolver"))
      {
	if (parameters(0) == "Arpack")
	  eigenvalue_solver_used = TypeEigenvalueSolver::ARPACK;
	else if (parameters(0) == "Anasazi")
	  eigenvalue_solver_used = TypeEigenvalueSolver::ANASAZI;
	else if (parameters(0) == "Feast")
	  eigenvalue_solver_used = TypeEigenvalueSolver::FEAST;
	else if (parameters(0) == "Slepc")
	  eigenvalue_solver_used = TypeEigenvalueSolver::SLEPC;
	else
	  {
	    cout << "Unknown eigenvalue solver " << parameters(0) << endl;
	    abort();
	  }
      }
    else if (!description_field.compare("Eigenvalue"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "Eigenvalue needs more parameters, for instance :" << endl;
	    cout << "Eigenvalue = YES SHIFTED SMALL n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("NO"))
          {
            this->nb_eigenvalues_wanted = 0;
          }
        else
          {
            int nb = 2;
            if (!parameters(1).compare("SHIFTED"))
              this->eigenvalue_computation_mode = this->SHIFTED_MODE;
	    else if (!parameters(1).compare("COMPLEX_SHIFT"))
	      this->eigenvalue_computation_mode = this->IMAG_SHIFTED_MODE;           
            else if (!parameters(1).compare("BUCKLING"))
              this->eigenvalue_computation_mode = this->BUCKLING_MODE;
            else if (!parameters(1).compare("CAYLEY"))
              this->eigenvalue_computation_mode = this->CAYLEY_MODE;
            else if (!parameters(1).compare("INVERT"))
              this->eigenvalue_computation_mode = this->INVERT_MODE;
            else
              this->eigenvalue_computation_mode = this->REGULAR_MODE;
            
            if (!parameters(2).compare("SMALL"))
              {
                this->type_spectrum_wanted = this->SMALL_EIGENVALUES;
                this->shift = 0; nb++;
              }
            else if (!parameters(2).compare("LARGE"))
              {
                this->type_spectrum_wanted = this->LARGE_EIGENVALUES;
                nb++;
              }
            else if (!parameters(2).compare("COMPLEX"))
              {
                this->type_spectrum_wanted = this->CENTERED_EIGENVALUES;
                nb++;
                this->shift = 2.0*pi_wp*to_num<T>(parameters(nb++));
                this->shift_imag = 2.0*pi_wp*to_num<T>(parameters(nb++));
		
		if (parameters.GetM() <= 4)
		  {
		    cout << "In SetInputData of EigenProblemMontjoie" << endl;
		    cout << "Eigenvalue needs more parameters, for instance :" << endl;
		    cout << "Eigenvalue = YES COMPLEX a b n" << endl;
		    cout << "Current parameters are : " << endl << parameters << endl;
		    abort();
		  }
              }
            else if (!parameters(2).compare("CENTERED"))
              {
                this->type_spectrum_wanted = this->CENTERED_EIGENVALUES;
                nb++;
                
                if (IsComplexNumber(T(0)))
                  {
                    this->shift = to_num<T>(parameters(nb++));
                    SetComplexZero(this->shift_imag);
                  }
                else
                  {
                    // frequency is given
                    // we deduce the pulsation of square of pulsation
                    // (for second-order formulation)
                    T omega = 2.0*pi_wp*to_num<T>(parameters(nb++));
                    
                    if (var_problem.FirstOrderFormulation())
                      {
                        this->shift = 0;
                        this->shift_imag = omega;
                      }
                    else
                      {
                        this->shift = omega*omega;
                        this->shift_imag = 0;
                      }
                  }
              }
	    else
	      {
		cout << "Unknown spectrum" << endl;
                abort();
	      }

            // number of eigenvalues
            to_num(parameters(nb++), this->nb_eigenvalues_wanted);
          }
      }
  }
  

  //! computation of mass matrix (when diagonal)
  template<class T>
  void EigenProblemMontjoie<T>::ComputeDiagonalMass()
  {
    T one, zero; SetComplexOne(one); SetComplexZero(zero);
    
    // assembling the matrix with only the diagonal
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(one);
    nat_mat.SetCoefStiffness(zero);
    nat_mat.SetCoefDamping(zero);
    
    Vector<T> D;
    FemMatrixFreeClass_Base<T>* Mh;
    Mh = var_computation.GetNewIterativeMatrix(zero);
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    var_computation.ComputeDiagonalMatrix(D, *Mh, nat_mat);

    this->sqrt_diagonal_mass.Reallocate(local_col_numbers.GetM());
    for (int i = 0; i < local_col_numbers.GetM(); i++)
      this->sqrt_diagonal_mass(i) = D(local_col_numbers(i));
    
    delete Mh;
  }
  
    
  //! Prepares computation before Cholesky factorization
  template<class T>
  void EigenProblemMontjoie<T>::ComputeMassForCholesky()
  {
    // everything is done in FactorizeCholeskyMass
  }
  
  
  //! Cholesky factorisation of mass matrix
  template<class T>
  void EigenProblemMontjoie<T>::FactorizeCholeskyMass()
  {
    GetFactoCholeskyMass(T(1.0));
  }


  template<class T>
  void EigenProblemMontjoie<T>::SetCholeskySolver(int solver)
  {
    type_solver_chol = solver;
  }
  
 
  //! Cholesky factorisation of mass matrix
  template<>
  void EigenProblemMontjoie<Complex_wp>::GetFactoCholeskyMass(const Real_wp&)
  {
    cout << "Not implemented for complex numbers" << endl;
    abort();
  }
  
  
  //! Cholesky factorisation of mass matrix
  template<>
  void EigenProblemMontjoie<Real_wp>::GetFactoCholeskyMass(const Real_wp&)
  {
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Mh;
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0);
    nat_mat.SetCoefStiffness(0.0);
    nat_mat.SetCoefDamping(0.0);
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    
    var_computation.AddMatrixWithBC(Mh, nat_mat);

    if (this->type_solver_chol >= 0)
      this->chol_facto_mass_matrix.SelectDirectSolver(this->type_solver_chol);
    else
      {
        if (this->chol_facto_mass_matrix.GetDirectSolver() <= 0)
          {
            cout << "Efficient sparse Cholesky solver not available" << endl;
            cout << "Recompile with Cholmod or Pastix" << endl;
            cout << "You can also use invert mode to compute eigenvalues in order to solve a non-symmetric eigenproblem (no Cholesky solver needed)" << endl;
            abort();
          }
      }
    
    if (this->print_level > 0)
      this->chol_facto_mass_matrix.ShowMessages();

    this->chol_facto_mass_matrix.Factorize(Mh);

    if (this->print_level < 2)
      this->chol_facto_mass_matrix.HideMessages();

    Xchol_real.Reallocate(this->n_);
    Xchol_imag.Reallocate(this->n_);
  }
  
  
  //! unused function
  template<class T>
  void EigenProblemMontjoie<T>::GetFactoCholeskyMass(const Complex_wp&)
  {
    // nothing
    cout << "Cholesky facto not available for complex numbers" << endl;
    abort();
  }


  template<class T> template<class T0>
  void EigenProblemMontjoie<T>::ExpandVector(const Vector<T0>& X, Vector<T0>& Y,
                                             bool assemble)
  {
    Y.Reallocate(var_problem.GetNbDof());
    Y.Zero();
    for (int i = 0; i < local_col_numbers.GetM(); i++)
      Y(local_col_numbers(i)) = X(i);

    if (assemble)
      var_comm.AddDomains(Y);
  }
  

  template<class T> template<class T0>
  void EigenProblemMontjoie<T>::ContractVector(const Vector<T0>& X, Vector<T0>& Y)
  {
    for (int i = 0; i < local_col_numbers.GetM(); i++)
      Y(i) = X(local_col_numbers(i));
  }

  
  //! Solves L y = x or L^T y = x where M = L L^T
  template<class T>
  void EigenProblemMontjoie<T>
  ::SolveCholeskyMass(const SeldonTranspose& transA, Vector<Real_wp>& X)
  {
    VectReal_wp Y;
    ExpandVector(X, Y);
    
    chol_facto_mass_matrix.Solve(transA, Y);

    ContractVector(Y, X);
  }
  
  
  //! Solves L y = x or L^T y = x where M = L L^T
  template<class T>
  void EigenProblemMontjoie<T>
  ::SolveCholeskyMass(const SeldonTranspose& transA, Vector<Complex_wp>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }

    VectReal_wp Yreal, Yimag;
    ExpandVector(Xchol_real, Yreal);
    ExpandVector(Xchol_imag, Yimag);
    
    chol_facto_mass_matrix.Solve(transA, Yreal);
    chol_facto_mass_matrix.Solve(transA, Yimag);
    
    ContractVector(Yreal, Xchol_real);
    ContractVector(Yimag, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));
  }
  

  //! computes y = L x or y = L^T x where M = L L^T
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltCholeskyMass(const SeldonTranspose& transA, Vector<Real_wp>& X)
  {
    VectReal_wp Y;
    ExpandVector(X, Y);

    chol_facto_mass_matrix.Mlt(transA, Y);

    ContractVector(Y, X);
  }
  

  //! computes y = L x or y = L^T x where M = L L^T
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltCholeskyMass(const SeldonTranspose& transA, Vector<Complex_wp>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }

    VectReal_wp Yreal, Yimag;
    ExpandVector(Xchol_real, Yreal);
    ExpandVector(Xchol_imag, Yimag);
    
    chol_facto_mass_matrix.Mlt(transA, Yreal);
    chol_facto_mass_matrix.Mlt(transA, Yimag);

    ContractVector(Yreal, Xchol_real);
    ContractVector(Yimag, Xchol_imag);
    
    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<Treal>(Xchol_real(i), Xchol_imag(i));    
  }

  
  //! computation of mass matrix
  template<class T>
  void EigenProblemMontjoie<T>::ComputeMassMatrix()
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    ComputeStiffnessMatrix(one, zero, *mat_mass, true);    
  }
  

  //! matrix vector product with mass matrix, Y = M X
  template<class T> template<class T0>
  void EigenProblemMontjoie<T>::MltMassGen(const SeldonTranspose& trans, const Vector<T0>& X, Vector<T0>& Y)
  {
    Vector<T0> X2, Y2(var_problem.GetNbDof());
    ExpandVector(X, X2, true);
    
    mat_mass->MltVector(trans, X2, Y2);
    
    ContractVector(Y2, Y);
  }


  //! matrix vector product with mass matrix, Y = M X
  template<class T>
  void EigenProblemMontjoie<T>::MltMass(const Vector<Treal>& X, Vector<Treal>& Y)  
  {
    MltMassGen(SeldonNoTrans, X, Y);
  }


  //! matrix vector product with mass matrix, Y = M X
  template<class T>
  void EigenProblemMontjoie<T>::MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y)  
  {
    MltMassGen(SeldonNoTrans, X, Y);
  }


  //! matrix vector product with mass matrix, Y = M X
  template<class T>
  void EigenProblemMontjoie<T>::MltMass(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y)  
  {
    MltMassGen(trans, X, Y);
  }


  //! matrix vector product with mass matrix, Y = M X
  template<class T>
  void EigenProblemMontjoie<T>::MltMass(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y)  
  {
    MltMassGen(trans, X, Y);
  }


  //! computation of stiffness matrix
  template<class T>
  void EigenProblemMontjoie<T>::ComputeStiffnessMatrix()
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    ComputeStiffnessMatrix(zero, one, *mat_stiff, true);
  }
  
    
  //! computation of  a M + b K, where M and K are mass and stiffness matrices
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeStiffnessMatrix(const T& a, const T& b)
  {
    ComputeStiffnessMatrix(a, b, *mat_stiff, true);
  }
  
  
  //! computation of a M + b K
  /*!
    \param[in] iterative if false, the matrix is stored
   */
  template<class T> template<class T0>
  void EigenProblemMontjoie<T>
  ::ComputeStiffnessMatrix(const T0& a, const T0& b, FemMatrixFreeClass_Base<T0>& A, bool iterative)
  {
    GlobalGenericMatrix<T0> nat_mat;
    nat_mat.SetCoefMass(a);
    nat_mat.SetCoefStiffness(b);
    nat_mat.SetCoefDamping(b);
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    
    int type_storage = var_computation.GetStorageFiniteElementMatrix();
    if (!iterative)
      var_computation.SetStorageFiniteElementMatrix(var_computation.MATRIX_STORED);

    A.Clear();
    var_computation.AddMatrixWithBC(A, nat_mat);
    
    if (!iterative)
      var_computation.SetStorageFiniteElementMatrix(type_storage);
  }
  
  
  //! multiplication by the stiffness matrix Y = K_h*X
  template<class T> template<class T0>
  void EigenProblemMontjoie<T>
  ::MltStiffnessGen(const SeldonTranspose& trans, const Vector<T0>& X, Vector<T0>& Y)
  {
    Vector<T0> X2, Y2(var_problem.GetNbDof());
    ExpandVector(X, X2, true);
    
    mat_stiff->MltVector(trans, X2, Y2);
    
    ContractVector(Y2, Y);
  }
  

  //! computation of Y = K X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y)
  {
    MltStiffnessGen(SeldonNoTrans, X, Y);
  }


  //! computation of Y = K X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    MltStiffnessGen(SeldonNoTrans, X, Y);
  }


  //! computation of Y = (a M + b K) X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const T& a, const T& b,
		 const Vector<Treal>& X, Vector<Treal>& Y)
  {
    MltStiffnessGen(SeldonNoTrans, X, Y);
  }


  //! computation of Y = (a M + b K) X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const T& a, const T& b,
		 const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    MltStiffnessGen(SeldonNoTrans, X, Y);
  }


  //! computation of Y = K^T X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const SeldonTranspose& trans,
                 const Vector<Treal>& X, Vector<Treal>& Y)
  {
    MltStiffnessGen(trans, X, Y);
  }


  //! computation of Y = K^T X
  template<class T>
  void EigenProblemMontjoie<T>
  ::MltStiffness(const SeldonTranspose& trans,
                 const Vector<Tcplx>& X, Vector<Tcplx>& Y)
  {
    MltStiffnessGen(trans, X, Y);
  }

  
  //! computation and factorization of the matrix a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeAndFactorizeStiffnessMatrix(const Real_wp& a, const Real_wp& b, int which)
  {
    ComputeAndFactoRealMatrix(T(0), a, b, which);
  }
  
  
  //! computation and factorization of the matrix a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>::
  ComputeAndFactorizeStiffnessMatrix(const Complex_wp& a,
                                     const Complex_wp& b, int which)
  {
    this->selected_part = which;
    
    nature_matrix_complex.SetCoefMass(a);
    nature_matrix_complex.SetCoefDamping(b);
    nature_matrix_complex.SetCoefStiffness(b);
    if (var_problem.print_level > 2)
      glob_solver.SetPrintLevel(var_problem.print_level);
    else
      glob_solver.SetPrintLevel(0);
    
    this->complex_system = true;

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    glob_solver.PerformFactorizationStep(nature_matrix_complex);
    if ((var_problem.print_level < 6) || (!glob_solver.IsIterativeSolver()))
      glob_solver.SetPrintLevel(0);
    
  }


  //! computation and factorization of the matrix a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeAndFactoRealMatrix(const Real_wp&, const Real_wp& a, const Real_wp& b, int which)
  {
    this->selected_part = which;

    nature_matrix_real.SetCoefMass(a);
    nature_matrix_real.SetCoefDamping(b);
    nature_matrix_real.SetCoefStiffness(b);
    if (var_problem.print_level > 2)
      glob_solver.SetPrintLevel(var_problem.print_level);
    else 
      glob_solver.SetPrintLevel(0);
    
    this->complex_system = false;

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    glob_solver.PerformFactorizationStep(nature_matrix_real);
    if ((var_problem.print_level < 6) || (!glob_solver.IsIterativeSolver()))
      glob_solver.SetPrintLevel(0);
  }


  //! computation and factorization of the matrix a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeAndFactoRealMatrix(const Complex_wp&, const Real_wp& a, const Real_wp& b, int which)
  {
    cout << "Incompatibles types" << endl;
    cout << "Provide coefficients a and b of the same type as T" << endl;
    abort();
  }


  //! solving by a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeSolution(const Vector<Real_wp>& X, Vector<Real_wp>& Y)
  {
    ComputeSolution(SeldonNoTrans, X, Y);
  }


  //! solving by a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeSolution(const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
  {
    ComputeSolution(SeldonNoTrans, X, Y);
  }
  
  
  //! solving by a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeSolution(const SeldonTranspose& transA,
		    const Vector<Real_wp>& X, Vector<Real_wp>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    cout << "The result can not be a real vector" << endl;
	    abort();
	  }
	
	Vector<Complex_wp> Xcplx(this->n_);
	for (int i = 0; i < this->n_; i++)
	  Xcplx(i) = Complex_wp(X(i), 0);
	
	Solve(transA, nature_matrix_complex, Xcplx);
	
	if (this->selected_part == EigenProblem_Base<T>::IMAG_PART)
	  for (int i = 0; i < this->n_; i++)
	    Y(i) = imag(Xcplx(i));
	else
	  for (int i = 0; i < this->n_; i++)
	    Y(i) = real(Xcplx(i));
      }
    else
      {
        Copy(X, Y);
        Solve(transA, nature_matrix_real, Y);
      }
  }

  //! solving by a*M + b*K
  template<class T>
  void EigenProblemMontjoie<T>
  ::ComputeSolution(const SeldonTranspose& transA,
		    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
  {
    if (this->complex_system)
      {
	if (this->selected_part == EigenProblem_Base<T>::COMPLEX_PART)
	  {
	    Copy(X, Y);
	    Solve(transA, nature_matrix_complex, Y);
	  }
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
      }
    else
      {
        cout << "not implemented" << endl;
        abort();
      }
  }


  //! Internal function 
  template<class T> template<class Complexe>
  void EigenProblemMontjoie<T>::Solve(const SeldonTranspose& transA,
				      GlobalGenericMatrix<Complexe>& nat_mat, Vector<Complexe>& Y)
  {
    Vector<Complexe> Y2;
    ExpandVector(Y, Y2);
    
    glob_solver.ComputeSolution(transA, Y2, nat_mat);
    
    ContractVector(Y2, Y);
  }
  
  
  //! memory is released
  template<class T>
  void EigenProblemMontjoie<T>::Clear()
  {
    glob_solver.ClearFactorization();
  }
  
  
  //! computation of eigenmodes (eigenvalues and eigenvectors)
  template<class T>
  void EigenProblemMontjoie<T>::ComputeEigenModes()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    var_computation.UpdateShiftAdimensionalization(this->shift, this->shift_imag);
    if (rank_proc == 0)
      {DISP(this->shift); DISP(this->shift_imag); DISP(var_problem.GetOmega()); }

    FemMatrixFreeClass_Base<T>* Kh;
    FemMatrixFreeClass_Base<T> *Mh;
    
    Mh = var_computation.GetNewIterativeMatrix(T(0));
    Kh = var_computation.GetNewIterativeMatrix(T(0));
    
    InitMatrix(*Kh, *Mh);
    
    Vector<T> eigen_values, lambda_imag;
    Matrix<T, General, ColMajor> eigen_vectors;
    GetEigenvaluesEigenvectors(*this, eigen_values, lambda_imag, eigen_vectors, eigenvalue_solver_used);
    
    var_computation.UpdateEigenvaluesAdimensionalization(eigen_values, lambda_imag, eigen_vectors);
    int n = var_problem.GetNbDof();

    if (rank_proc == 0)
      cout << "Eigenvalues have been obtained with " << this->GetNbMatrixVectorProducts() << " linear solves" << endl;

    if (rank_proc == 0)
      {
        DISP(eigen_values); DISP(lambda_imag);
      }
      
    if (IsComplexNumber(T(0)))
      {
        // for harmonic problems, we multiply with -omega to obtain eigenpulsations
        // associated with M_h^{-1} K_h
        Mlt(-var_problem.GetOmega(), eigen_values);
        eigen_values.Write(file_name_eigenval);
      }
    else
      {
        // for stationary problems, lambda = 1j*omega for first order formulation
        // and lambda = omega**2 for second order formulation
        Vector<Complex_wp> omega(eigen_values.GetM());
        if (var_problem.FirstOrderFormulation())
          {
            for (int i = 0; i < eigen_values.GetM(); i++)
              omega(i) = -Iwp*(eigen_values(i) + Iwp*lambda_imag(i));
          }
        else
          {
            for (int i = 0; i < eigen_values.GetM(); i++)
              omega(i) = sqrt(eigen_values(i) + Iwp*lambda_imag(i));
          }
        
        omega.Write(file_name_eigenval);
      }
    
    this->Clear();
    delete Kh; delete Mh;
    
    // the modes are written in files
    Vector<T> eigen_mode(n);
    Vector<Vector<T> > eigen_mode_v;
    eigen_mode_v.SetData(1, &eigen_mode);
    for (int i = 0; i < this->nb_eigenvalues_wanted; i++)
      {
        eigen_mode.Zero();
        for (int j = 0; j < eigen_vectors.GetM(); j++)
          eigen_mode(local_col_numbers(j)) = eigen_vectors(j, i);
        
        var_comm.AddDomains(eigen_mode);
        
        string numero = NumberToString(i);
        for (int j = 0; j < var_output.output_mesh_param.GetM(); j++)
          var_output.output_mesh_param(j).
            SetFileName(2, GetBaseString(var_output.output_mesh_param(j).GetTotalFieldFile())
                        + "." + numero.substr(numero.size()-3,3) + string(".bb"));
					 
        for (int j = 0; j < var_output.output_grid_param.GetM(); j++)
          var_output.output_grid_param(j).
            SetFileName(2, GetBaseString(var_output.output_grid_param(j).GetTotalFieldFile())
                        + numero + string(".dat"));
					 
        var_output.WriteOutputFile(eigen_mode_v, 2);
      }
				 
    eigen_mode_v.Nullify();
  }

  
  template<class T>
  bool EigenProblemMontjoie<T>::IsSymmetricProblem() const
  {
    return var_problem.IsSymmetricProblem(true);
  }
  

  template<class T>
  bool EigenProblemMontjoie<T>::IsHermitianProblem() const
  {
    if (IsComplexNumber(T(0)))
	return false;
    
    return var_problem.IsSymmetricProblem(true);
  }
  
  
  /**********************************
   * PolynomialEigenProblemMontjoie *
   **********************************/
  
  
  //! modifies the object with a line of the data file
  template<class T>
  void PolynomialEigenProblemMontjoie<T>
  ::SetInputData(const string& description_field, const Vector<string>& parameters)
  {
    if (!description_field.compare("EigenvalueTolerance"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "EigenvalueTolerance needs more parameters, for instance :" << endl;
	    cout << "EigenvalueTolerance = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        this->stopping_criterion = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("EigenvalueMaxNumberIterations"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "EigenvalueMaxNumberIterations needs more parameters, for instance :" << endl;
	    cout << "EigenvalueMaxNumberIterations = n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        this->nb_maximum_iterations = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("FileEigenvalue"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "FileEigenvalue needs more parameters, for instance :" << endl;
	    cout << "FileEigenvalue = file_name_eigenval" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        file_name_eigenval = parameters(0);
      }
    else if (!description_field.compare("PolynomialEigenvalue"))
      {
	if (parameters.GetM() < 3)
	  {
	    cout << "In SetInputData of EigenProblemMontjoie" << endl;
	    cout << "Eigenvalue needs more parameters, for instance :" << endl;
	    cout << "Eigenvalue = YES n SHIFTED shift" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (parameters(0) == "NO")
          {
            this->nb_eigenvalues_wanted = 0;
          }
        else
          {
            this->nb_eigenvalues_wanted = to_num<int>(parameters(1));
            if (parameters(2) == "SHIFTED")
              {
                this->type_spectrum_wanted = this->CENTERED_EIGENVALUES;
                this->shift = to_num<T>(parameters(3));
                this->use_spectral_transfo = true;
              }
            else
              {
                this->type_spectrum_wanted = this->LARGE_EIGENVALUES;
                this->use_spectral_transfo = false;
              }
          }
      }
  }
  
  
  //! returns row numbers involved in the computation of eigenvalues
  template<class T>
  const IVect& PolynomialEigenProblemMontjoie<T>::GetLocalColumnNumbers() const
  {
    return local_col_numbers;
  }
  
  
  //! inits the eigenvalue computation
  template<class T>
  void PolynomialEigenProblemMontjoie<T>::InitComputation()
  {
    this->pol_degree = 2;
    Vector<bool> diag_elt;
    int type_mass_matrix = var_computation.GetMassMatrixType(diag_elt);
    if (type_mass_matrix == FemMassMatrix::DIAGONAL)
      this->SetDiagonalMass(true);
    else
      this->SetDiagonalMass(false);
    
    // constructing local_col_numbers
    int m = var_problem.GetNbDof();
    int ndir = var_boundary.GetNbDirichletDof();
    this->print_level = var_problem.print_level/2;
    
#ifdef SELDON_WITH_MPI
    this->SetCommunicator(var_comm.comm_group_mode);    
    
    int noverlap = var_comm.GetNbOverlappedDof();
    Vector<bool> OverlappedRow(m); OverlappedRow.Fill(false);
    for (int i = 0; i < noverlap; i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;
    
    int noverlap2 = noverlap;
    if (ndir > 0)
      for (int i = 0; i < m; i++)
        if (var_boundary.IsDofDirichlet(i) && OverlappedRow(i))
          noverlap2--;
#else
    int noverlap2 = 0;
#endif
    
    int n = m - noverlap2 - ndir;
    local_col_numbers.Reallocate(n);

    int ncol = 0;
#ifdef SELDON_WITH_MPI    
    if (ndir > 0)
      {
        for (int i = 0; i < m; i++)
          if (!var_boundary.IsDofDirichlet(i) && !OverlappedRow(i))
            local_col_numbers(ncol++) = i;
      }
    else
      {
        for (int i = 0; i < m; i++)
          if (!OverlappedRow(i))
            local_col_numbers(ncol++) = i;
      }
#else
    if (ndir > 0)
      {
        for (int i = 0; i < m; i++)
          if (!var_boundary.IsDofDirichlet(i))
            local_col_numbers(ncol++) = i;
      }
#endif
    
    this->Init(n);
  }
    

  template<class T>
  void PolynomialEigenProblemMontjoie<T>::ExpandVector(const Vector<T>& X, Vector<T>& Y,
                                                       bool assemble)
  {
    Y.Reallocate(var_problem.GetNbDof());
    Y.Zero();
    for (int i = 0; i < local_col_numbers.GetM(); i++)
      Y(local_col_numbers(i)) = X(i);

    if (assemble)
      var_comm.AddDomains(Y);
  }
  

  template<class T>
  void PolynomialEigenProblemMontjoie<T>::ContractVector(const Vector<T>& X, Vector<T>& Y)
  {
    for (int i = 0; i < local_col_numbers.GetM(); i++)
      Y(i) = X(local_col_numbers(i));
  }

  
  //! computes the finite element matrix with coef
  template<class T>
  void PolynomialEigenProblemMontjoie<T>::ComputeOperator(int num, const Vector<T>& coef)
  {
    mat_op(num) = var_computation.GetNewIterativeMatrix(T(0));

    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(coef(2));
    nat_mat.SetCoefDamping(coef(1));
    nat_mat.SetCoefStiffness(coef(0));
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    
    var_computation.AddMatrixWithBC(*mat_op(num), nat_mat);
  }
  
  
  template<class T>
  void PolynomialEigenProblemMontjoie<T>
  ::MltOperator(int num, const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    Vector<T> X2, Y2(var_problem.GetNbDof());
    ExpandVector(X, X2, true);

    mat_op(num)->MltVector(trans, X2, Y2);

    ContractVector(Y2, Y);
  }
    
  
  template<class T>
  void PolynomialEigenProblemMontjoie<T>::FactorizeMass()
  {
    Vector<T> coef(3);
    coef(0) = Real_wp(1); coef(1) = Real_wp(0); coef(2) = Real_wp(0);
    ComputeOperator(0, coef);

    coef(0) = Real_wp(0); coef(1) = Real_wp(1); coef(2) = Real_wp(0);
    ComputeOperator(1, coef);
    
    coef(0) = Real_wp(0); coef(1) = Real_wp(0); coef(2) = Real_wp(1);
    ComputeOperator(2, coef);
    
    if (this->DiagonalMass())
      {
        this->invDiag.Reallocate(this->n_);
        
        // assembling the matrix with only the diagonal
        T one, zero; SetComplexOne(one); SetComplexZero(zero);
        GlobalGenericMatrix<T> nat_mat;
        nat_mat.SetCoefMass(one);
        nat_mat.SetCoefStiffness(zero);
        nat_mat.SetCoefDamping(zero);
        var_boundary.SetCoefficientDirichlet(Real_wp(1));

        Vector<T> D;
        FemMatrixFreeClass_Base<T>* Mh;
        Mh = var_computation.GetNewIterativeMatrix(zero);
        var_computation.ComputeDiagonalMatrix(D, *Mh, nat_mat);
        
        for (int i = 0; i < this->n_; i++)
          this->invDiag(i) = one / D(local_col_numbers(i));

        delete Mh;
      }
    else
      {        
        coef(0) = Real_wp(0); coef(1) = Real_wp(0); coef(2) = Real_wp(1);
        FactorizeOperator(coef);    
      }
  }
  
  
  template<class T>
  void PolynomialEigenProblemMontjoie<T>
  ::SolveMass(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y)
  {
    if (this->DiagonalMass())
      for (int i = 0; i < x.GetM(); i++)
        y(i) = x(i)*this->invDiag(i);
    else
      SolveOperator(trans, x, y);
  }

  
  template<class T>
  void PolynomialEigenProblemMontjoie<T>::FactorizeOperator(const Vector<T>& coef)
  {
    /*
     {
      T one, zero;
      SetComplexZero(zero); SetComplexOne(one);
      
      GlobalGenericMatrix<T> nat_mat;
      nat_mat.SetCoefMass(one);
      nat_mat.SetCoefStiffness(zero);
      nat_mat.SetCoefDamping(zero);

      DistributedMatrix<T, General, ArrayRowSparse> K, S, M;
      var_computation.AddMatrixWithBC(M, nat_mat);

      nat_mat.SetCoefMass(zero);
      nat_mat.SetCoefStiffness(zero);
      nat_mat.SetCoefDamping(one);
      var_computation.AddMatrixWithBC(S, nat_mat);

      nat_mat.SetCoefMass(zero);
      nat_mat.SetCoefStiffness(one);
      nat_mat.SetCoefDamping(zero);
      var_computation.AddMatrixWithBC(K, nat_mat);

      K.WriteText("Kh.dat");
      S.WriteText("Sh.dat");
      M.WriteText("Mh.dat");
     }
    */
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(coef(2));
    nat_mat.SetCoefDamping(coef(1));
    nat_mat.SetCoefStiffness(coef(0));
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    
    nat_mat_solver = nat_mat;
    glob_solver.SetPrintLevel(2);
    glob_solver.PerformFactorizationStep(nat_mat);
    glob_solver.SetPrintLevel(0);
  }
  
  
  template<class T>
  void PolynomialEigenProblemMontjoie<T>
  ::SolveOperator(const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    Vector<T> Y2;
    ExpandVector(X, Y2);

    glob_solver.ComputeSolution(trans, Y2, nat_mat_solver);
    
    ContractVector(Y2, Y);
  }
    
  
  // main method called to retrieve all the eigenvalues and eigenvectors
  template<class T>
  void PolynomialEigenProblemMontjoie<T>::ComputeEigenModes()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // the computation is initialized
    this->InitComputation();
    
    // computation of eigenvalues with Slepc
    Vector<T> eigen_values, lambda_imag;
    Matrix<T, General, ColMajor> eigen_vectors;
#ifdef SELDON_WITH_SLEPC
    FindEigenvaluesSlepc(*this, eigen_values, lambda_imag, eigen_vectors);
#else
    cout << "Recompile with Slepc" << endl;
#endif
    
    if (rank_proc == 0)
      cout << "Eigenvalues have been obtained with " << this->nb_linear_solves << " linear solves" << endl;
    
    // eigenvalues are multiplied by omega to obtain eigenpulsations
    for (int i = 0; i < eigen_values.GetM(); i++)
      eigen_values(i) *= var_problem.GetOmega();

    // pulsations are written
    eigen_values.Write(file_name_eigenval);

    // the modes are written in files
    int n = var_problem.GetNbDof();
    Vector<T> eigen_mode(n);
    Vector<Vector<T> > eigen_mode_v;
    eigen_mode_v.SetData(1, &eigen_mode);
    for (int i = 0; i < this->nb_eigenvalues_wanted; i++)
      {
        eigen_mode.Zero();
        for (int j = 0; j < eigen_vectors.GetM(); j++)
          eigen_mode(local_col_numbers(j)) = eigen_vectors(j, i);
        
        var_comm.AddDomains(eigen_mode);
        
        string numero = NumberToString(i);
        for (int j = 0; j < var_output.output_mesh_param.GetM(); j++)
          var_output.output_mesh_param(j).
            SetFileName(2, GetBaseString(var_output.output_mesh_param(j).GetTotalFieldFile())
                        + "." + numero.substr(numero.size()-3,3) + string(".bb"));
					 
        for (int j = 0; j < var_output.output_grid_param.GetM(); j++)
          var_output.output_grid_param(j).
            SetFileName(2, GetBaseString(var_output.output_grid_param(j).GetTotalFieldFile())
                        + numero + string(".dat"));
					 
        var_output.WriteOutputFile(eigen_mode_v, 2);
      }
				 
    eigen_mode_v.Nullify();
  }
  

  /*******************
   * Other functions *
   *******************/
  

  //! returns spectral radius of M^-1 K
  template<class T>
  Real_wp GetSpectralRadius(VirtualMatrix<T>& mat_stiff, VirtualMatrix<Real_wp>& mat_mass)
  {
    VirtualEigenProblem<T, T, Real_wp> var_eig;
    
    var_eig.SetStoppingCriterion(1e-11);
    var_eig.SetNbAskedEigenvalues(10);
    
    // finding large eigenvalues of M^-1 K
    var_eig.InitMatrix(mat_stiff, mat_mass);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
    
    Vector<T> lambda, lambda_imag;
    Matrix<T, General, ColMajor> eigen_vec;

    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    Real_wp rho = 0;
    for (int i = 0; i < lambda.GetM(); i++)
      rho = max(rho, abs(lambda(i) + Iwp*lambda_imag(i)));
    
    return rho;
  }


  //! returns spectral radius of K
  template<class T>
  Real_wp GetSpectralRadius(VirtualMatrix<T>& mat_stiff)
  {
    VirtualEigenProblem<T> var_eig;
    
    var_eig.SetStoppingCriterion(1e-11);
    var_eig.SetNbAskedEigenvalues(10);
    
    // finding large eigenvalues of K
    var_eig.InitMatrix(mat_stiff);
    var_eig.SetComputationalMode(var_eig.REGULAR_MODE);
    var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0, var_eig.SORTED_MODULUS);
    
    Vector<T> lambda, lambda_imag;
    Matrix<T, General, ColMajor> eigen_vec;

    GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
    
    Real_wp rho = 0;
    for (int i = 0; i < lambda.GetM(); i++)
      rho = max(rho, abs(lambda(i) + Iwp*lambda_imag(i)));
    
    return rho;
  }
    
}
			 
#define MONTJOIE_FILE_EIGENVALUES_HARMONIC_CXX
#endif
