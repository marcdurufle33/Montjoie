#ifndef MONTJOIE_FILE_SOLVE_HARMONIC_INLINE_CXX

namespace Montjoie
{
  
  /********************
   * All_LinearSolver *
   ********************/
  
  
  template<class TypeEquation>
  inline All_LinearSolver::All_LinearSolver(EllipticProblem<TypeEquation>& var)
    : mat_iterative_real(new FemMatrixFreeClass<Real_wp, TypeEquation>(var)),
      mat_iterative_cplx(new FemMatrixFreeClass<Complex_wp, TypeEquation>(var)),
      precond_real(new All_Preconditioner<Real_wp, TypeEquation>(var)),
      precond_cplx(new All_Preconditioner<Complex_wp, TypeEquation>(var)),
      iterative_solver_real(*precond_real), iterative_solver_cplx(*this),
      solver_real(new CondensationBlockSolver<Real_wp, TypeEquation>(var)),
      solver_cplx(new CondensationBlockSolver<Complex_wp, TypeEquation>(var)),
      var_problem(var), var_source(var), var_boundary(var),
      var_computation(var), var_comm(var)
  {    
    InitDefaultValues();
  }
  
  
  //! returns linear solver used
  inline int All_LinearSolver::GetTypeSolver() const
  {
    return type_choice_solver_auto;
  }
  
  
  //! returns true if an iterative solver is used
  inline bool All_LinearSolver::IsIterativeSolver() const
  {
    return iterative_resolution;
  }


  //! forces the use of a direct solver
  inline void All_LinearSolver::SetDirectSolver()
  {
    iterative_resolution = false;
  }


  //! forces the use of a direct solver
  inline void All_LinearSolver::SetDirectSolver(const VectString& param)
  {
    iterative_resolution = false;
    mat_lu_real.SelectDirectSolver(param);
    mat_lu_cplx.SelectDirectSolver(param);
  }

  
  //! sets the iterative solver to use
  inline void All_LinearSolver::SetIterativeSolver(int type)
  {
    iterative_resolution = true;
    iterative_solver_real.SetSolverType(type);
    iterative_solver_cplx.SetSolverType(type);
  }
  

  //! returns true if a static condensation is performed
  inline bool All_LinearSolver::StaticCondensation() const
  {
    return static_condensation;
  }
  

  //! sets a real preconditioning
  inline void All_LinearSolver::EnableRealPreconditioning(bool flag)
  {
    use_real_preconditioning = flag;
  }


  //! changes the real preconditioning
  inline void All_LinearSolver::SetRealPreconditioning(All_Preconditioner_Base<Real_wp>* prec)
  {
    precond_real = prec;
    shallow_copy_real_precond = true;
  }

  
  //! sets the iterative solver to use
  inline void All_LinearSolver::SetIterativeSolver(int type, int type_precond)
  {
    iterative_resolution = true;
    iterative_solver_real.SetSolverType(type);
    iterative_solver_cplx.SetSolverType(type);
    precond_real->SetPreconditioningType(type_precond);
    precond_cplx->SetPreconditioningType(type_precond);
  }

  
  //! returns the maximal number of iterations for iterative solver
  inline int All_LinearSolver::GetMaxNumberIteration() const
  {
    if (complex_matrix)
      return iterative_solver_cplx.GetMaxNumberIteration(); 
    
    return iterative_solver_real.GetMaxNumberIteration(); 
  }
  
  
  //! returns the restart number for iterative solver
  inline int All_LinearSolver::GetRestart() const
  {
    if (complex_matrix)
      return iterative_solver_cplx.GetRestart();
    
    return iterative_solver_real.GetRestart(); 
  }
  
  
  //! returns the stopping criterion used by the iterative solver
  inline Real_wp All_LinearSolver::GetStoppingCriterion() const
  {
    if (complex_matrix)
      return iterative_solver_cplx.GetStoppingCriterion(); 
    
    return iterative_solver_real.GetStoppingCriterion(); 
  }
    

  //! sets the verbose level
  inline void All_LinearSolver::SetPrintLevel(int p)
  {
    mat_lu_real.SetPrintLevel(p);
    mat_lu_cplx.SetPrintLevel(p);
    iterative_solver_real.SetPrintLevel(p);
    iterative_solver_cplx.SetPrintLevel(p);
  }


  inline void All_LinearSolver::UseOldPattern(bool use_p)
  {
    use_old_pattern_Ah = use_p;
    new_pattern_Ah = true;
  }


  inline bool All_LinearSolver::UseCompressedMatrix() const
  {
    return compress_matrix;
  }


  inline void All_LinearSolver::DisableLocalPreconditioning()
  {
    precond_real->DisableLocalPreconditioning();
    precond_cplx->DisableLocalPreconditioning();
  }

  
  //! returns a reference to the iterative matrix
  inline FemMatrixFreeClass_Base<Real_wp>& All_LinearSolver
  ::GetIterativeMatrix(const Real_wp&)
  {
    return *mat_iterative_real;
  }


  //! returns a reference to the iterative matrix
  inline FemMatrixFreeClass_Base<Complex_wp>& All_LinearSolver
  ::GetIterativeMatrix(const Complex_wp&)
  {
    return *mat_iterative_cplx;
  }


  //! returns a reference to the iterative matrix
  inline All_Preconditioner_Base<Real_wp>& All_LinearSolver
  ::GetPreconditioning(const Real_wp&)
  {
    return *precond_real;
  }


  //! returns a reference to the iterative matrix
  inline All_Preconditioner_Base<Complex_wp>& All_LinearSolver
  ::GetPreconditioning(const Complex_wp&)
  {
    return *precond_cplx;
  }

  
  //! returns a reference to the direct matrix
  inline DistributedMatrix<Real_wp, General, ArrayRowSparse>& 
  All_LinearSolver::GetDirectUnsymmetricMatrix(const Real_wp&)
  {
    return mat_unsym_real;
  }
  

  //! returns a reference to the direct matrix
  inline DistributedMatrix<Complex_wp, General, ArrayRowSparse>&
  All_LinearSolver::GetDirectUnsymmetricMatrix(const Complex_wp&)
  {
    return mat_unsym_cplx;
  }


  //! returns a reference to the direct matrix
  inline DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& 
  All_LinearSolver::GetDirectSymmetricMatrix(const Real_wp&)
  {
    return mat_sym_real;
  }
  
  
  //! returns a reference to the direct matrix
  inline DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&
  All_LinearSolver::GetDirectSymmetricMatrix(const Complex_wp&)
  {
    return mat_sym_cplx;
  }
  
  
  //! returns a reference to the iterative solver
  inline All_IterativeSolver<Real_wp>& All_LinearSolver::GetIterativeSolver(const Real_wp&)
  {
    return iterative_solver_real;
  }
  
  
  //! returns a reference to the iterative solver
  inline All_IterativeSolver<Complex_wp>& All_LinearSolver::GetIterativeSolver(const Complex_wp&)
  {
    return iterative_solver_cplx;
  }
  

  //! returns a reference to the iterative matrix
  inline All_MatrixLU<Real_wp>& All_LinearSolver::GetDirectFactorization(const Real_wp&)
  {
    return mat_lu_real;
  }


  //! returns a reference to the iterative matrix
  inline All_MatrixLU<Complex_wp>& All_LinearSolver::GetDirectFactorization(const Complex_wp&)
  {
    return mat_lu_cplx;
  }


  //! returns a reference to the solver handling static condensation
  inline CondensationBlockSolver_Fem<Real_wp>& All_LinearSolver::GetCondensedSolver(const Real_wp&)
  {
    return *solver_real;
  }


  //! returns a reference to the solver handling static condensation
  inline CondensationBlockSolver_Fem<Complex_wp>& All_LinearSolver::GetCondensedSolver(const Complex_wp&)
  {
    return *solver_cplx;
  }
  
  
  inline GlobalGenericMatrix<Real_wp>& All_LinearSolver
  ::GetNatureMatrix(const Real_wp&)
  {
    return nature_matrix_real;
  }


  inline GlobalGenericMatrix<Complex_wp>& All_LinearSolver
  ::GetNatureMatrix(const Complex_wp&)
  {
    return nature_matrix_cplx;
  }


  template<class T> template<class TypeEquation>
  inline CondensationBlockSolver_Fem<T>::CondensationBlockSolver_Fem(EllipticProblem<TypeEquation>& var)
    : var_computation(var), var_comm(var), var_problem(var), var_boundary(var)
  {
    InitDefaultValues();
  }
  

  template<class T, class TypeEquation>
  inline CondensationBlockSolver<T, TypeEquation>::CondensationBlockSolver(EllipticProblem<TypeEquation>& var)
    : CondensationBlockSolver_Fem<T>(var)
  {
  }

} // namespace Montjoie

#define MONTJOIE_FILE_SOLVE_HARMONIC_INLINE_CXX
#endif
