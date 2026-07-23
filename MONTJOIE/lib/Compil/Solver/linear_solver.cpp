#include "MontjoieFlag.hxx"

#include "Solver/MontjoieSolverHeader.hxx"
#include "Solver/MontjoieSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/SolveSystem.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class All_MatrixLU<Real_wp>;
  SELDON_EXTERN template class All_MatrixLU<Complex_wp>;

  SELDON_EXTERN template void All_MatrixLU<Real_wp>::SelectOptimalDirectSolver(int, Real_wp, Dimension2&);
  SELDON_EXTERN template void All_MatrixLU<Real_wp>::SelectOptimalDirectSolver(int, Real_wp, Dimension3&);

  SELDON_EXTERN template void All_MatrixLU<Complex_wp>::SelectOptimalDirectSolver(int, Complex_wp, Dimension2&);
  SELDON_EXTERN template void All_MatrixLU<Complex_wp>::SelectOptimalDirectSolver(int, Complex_wp, Dimension3&);

  SELDON_EXTERN template class VirtualMatrixTranspose<Real_wp>;
  SELDON_EXTERN template class VirtualMatrixTranspose<Complex_wp>;
  
  SELDON_EXTERN template class All_IterativeSolver<Real_wp>;
  SELDON_EXTERN template class All_IterativeSolver<Complex_wp>;
  
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const VirtualMatrix<Real_wp>&, VectReal_wp&, const VectReal_wp&);
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, VectReal_wp&, const VectReal_wp&);
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const VirtualMatrix<Real_wp>&, VectReal_wp&, const VectReal_wp&, Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);

  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const VirtualMatrix<Complex_wp>&, VectComplex_wp&, const VectComplex_wp&);
  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&, VectComplex_wp&, const VectComplex_wp&);
  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const VirtualMatrix<Complex_wp>&, VectComplex_wp&, const VectComplex_wp&, Preconditioner_Base<Complex_wp>&, Iteration<Real_wp>&);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&);
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&);
  SELDON_EXTERN template int All_IterativeSolver<Real_wp>::Solve(const VirtualMatrix<Real_wp>&, DistributedVector<Real_wp>&, const DistributedVector<Real_wp>&, Preconditioner_Base<Real_wp>&, Iteration<Real_wp>&);

  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const VirtualMatrix<Complex_wp>&, DistributedVector<Complex_wp>&, const DistributedVector<Complex_wp>&);
  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&, DistributedVector<Complex_wp>&, const DistributedVector<Complex_wp>&);
  SELDON_EXTERN template int All_IterativeSolver<Complex_wp>::Solve(const VirtualMatrix<Complex_wp>&, DistributedVector<Complex_wp>&, const DistributedVector<Complex_wp>&, Preconditioner_Base<Complex_wp>&, Iteration<Real_wp>&);
#endif
  
  SELDON_EXTERN template void All_IterativeSolver<Real_wp>::SelectOptimalIterativeSolver(int, Real_wp, Dimension2&, General);
  SELDON_EXTERN template void All_IterativeSolver<Real_wp>::SelectOptimalIterativeSolver(int, Real_wp, Dimension2&, Symmetric);
  SELDON_EXTERN template void All_IterativeSolver<Real_wp>::SelectOptimalIterativeSolver(int, Real_wp, Dimension3&, General);
  SELDON_EXTERN template void All_IterativeSolver<Real_wp>::SelectOptimalIterativeSolver(int, Real_wp, Dimension3&, Symmetric);

  SELDON_EXTERN template void All_IterativeSolver<Complex_wp>::SelectOptimalIterativeSolver(int, Complex_wp, Dimension2&, General);
  SELDON_EXTERN template void All_IterativeSolver<Complex_wp>::SelectOptimalIterativeSolver(int, Complex_wp, Dimension2&, Symmetric);
  SELDON_EXTERN template void All_IterativeSolver<Complex_wp>::SelectOptimalIterativeSolver(int, Complex_wp, Dimension3&, General);
  SELDON_EXTERN template void All_IterativeSolver<Complex_wp>::SelectOptimalIterativeSolver(int, Complex_wp, Dimension3&, Symmetric);

}
