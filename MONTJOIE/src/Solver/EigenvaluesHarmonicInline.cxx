#ifndef MONTJOIE_FILE_EIGENVALUES_HARMONIC_INLINE_CXX

#include "EigenvaluesHarmonic.hxx"

namespace Montjoie
{
  
  //! constructor with a generic problem and a given linear solver
  template<class T> template<class TypeEquation>
  inline EigenProblemMontjoie<T>::
  EigenProblemMontjoie(EllipticProblem<TypeEquation>& var,
                       All_LinearSolver& solver) 
    : glob_solver(solver), var_problem(var), var_computation(var), var_boundary(var),
      var_comm(var), var_output(var)
  { 
    InitDefaultValues();
  }

  
  //! constructor with a generic problem and a given linear solver
  template<class T> template<class TypeEquation>
  inline PolynomialEigenProblemMontjoie<T>::
  PolynomialEigenProblemMontjoie(EllipticProblem<TypeEquation>& var,
                                 All_LinearSolver& solver) 
    : glob_solver(solver), mat_op(3), var_problem(var), var_computation(var), 
      var_boundary(var), var_comm(var), var_output(var)
  { 
    file_name_eigenval = "Omega.dat";
  }

}

#define MONTJOIE_FILE_EIGENVALUES_HARMONIC_INLINE_CXX
#endif

  
