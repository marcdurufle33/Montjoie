#include "Solver/MontjoieSolverHeader.hxx"
#include "Solver/MontjoieSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/NonLinearEquations.cxx"
#include "Solver/NonLinearLeastSquares.cxx"
#include "Solver/NewtonSolver.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template int SolveMinpack(NonLinearEquations_Base<Real_wp, R2, Matrix2_2>&,
					  R2& x_sol, R2& fvec, Matrix2_2& fjac, R2& scale_eqn,
					  const IVect& Control, const VectReal_wp& RControl);

  SELDON_EXTERN template int SolveMinpack(NonLinearEquations_Base<Real_wp, R3, Matrix3_3>&,
					  R3& x_sol, R3& fvec, Matrix3_3& fjac, R3& scale_eqn,
					  const IVect& Control, const VectReal_wp& RControl);

  SELDON_EXTERN template int SolveMinpack(NonLinearEquations_Base<Real_wp, VectReal_wp, Matrix<Real_wp> >&,
					  VectReal_wp&, VectReal_wp&, Matrix<Real_wp>&, VectReal_wp&,
					  const IVect& Control, const VectReal_wp& RControl);
  
  SELDON_EXTERN template class VirtualMinimizedFunction<double>;
  SELDON_EXTERN template class VirtualLeastSquaresFunction<Real_wp>;

#ifdef MONTJOIE_WITH_MULTIPLE
  SELDON_EXTERN template class VirtualMinimizedFunction<Real_wp>;
#endif

#ifdef MONTJOIE_WITH_GSL
  SELDON_EXTERN template double MinimizeParametersGsl(VirtualMinimizedFunction<double>&, Vector<double>&, double, unsigned);
#ifdef MONTJOIE_WITH_MULTIPLE
  SELDON_EXTERN template Real_wp MinimizeParametersGsl(VirtualMinimizedFunction<Real_wp>&, Vector<Real_wp>&, double, unsigned);
#endif
#endif

#ifdef MONTJOIE_WITH_ALGLIB
  SELDON_EXTERN template double MinimizeParametersAlglib(VirtualMinimizedFunction<double>&, Vector<double>&, double, unsigned);
#ifdef MONTJOIE_WITH_MULTIPLE
  SELDON_EXTERN template Real_wp MinimizeParametersAlglib(VirtualMinimizedFunction<Real_wp>&, Vector<Real_wp>&, double, unsigned);
#endif
#endif

#ifdef SELDON_WITH_MKL
  SELDON_EXTERN template Real_wp SolveLeastSquaresMkl(VirtualLeastSquaresFunction<Real_wp>&, VectReal_wp&, Real_wp, unsigned);
#endif

  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp>&, VectReal_wp&, Real_wp, Real_wp, unsigned, int, Real_wp);

  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp, R2, R2, Matrix2_2, Matrix2_2sym>&, R2&, Real_wp, Real_wp, unsigned, int, Real_wp);

  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp, R3, R3, Matrix3_3, Matrix3_3sym>&, R3&, Real_wp, Real_wp, unsigned, int, Real_wp);
  
  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp, VectReal_wp, VectReal_wp, Matrix<Real_wp> >&, VectReal_wp&, Real_wp, Real_wp, unsigned, int, Real_wp);

  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp, TinyVector<Real_wp, 1>, R2, TinyMatrix<Real_wp, General, 2, 1>, TinyMatrix<Real_wp, Symmetric, 1, 1> >&, TinyVector<Real_wp, 1>&, Real_wp, Real_wp, unsigned, int, Real_wp);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template Real_wp SolveLeastSquaresLvm(VirtualLeastSquaresFunction<Real_wp, VectReal_wp, VectReal_wp, DistributedMatrix<Real_wp, General, RowMajor> >&, VectReal_wp&, Real_wp, Real_wp, unsigned, int, Real_wp);
#endif

  SELDON_EXTERN template class NewtonSolver<Real_wp>;
  SELDON_EXTERN template class NewtonSolver<Complex_wp>;
  
  SELDON_EXTERN template class NewtonSolver<Real_wp, R2, R2>;
  SELDON_EXTERN template class NewtonSolver<Real_wp, R3, R3>;

}
