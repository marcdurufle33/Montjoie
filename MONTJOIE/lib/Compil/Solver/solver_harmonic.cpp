#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/SolveHarmonic.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template void All_LinearSolver::SelectOptimalLinearSolver(int, int, Dimension2, General);
  SELDON_EXTERN template void All_LinearSolver::SelectOptimalLinearSolver(int, int, Dimension3, General);
  SELDON_EXTERN template void All_LinearSolver::SelectOptimalLinearSolver(int, int, Dimension2, Symmetric);
  SELDON_EXTERN template void All_LinearSolver::SelectOptimalLinearSolver(int, int, Dimension3, Symmetric);

  SELDON_EXTERN template void All_LinearSolver::PerformFactorizationStep(const GlobalGenericMatrix<Real_wp>&, bool);
  SELDON_EXTERN template void All_LinearSolver::PerformFactorizationStep(const GlobalGenericMatrix<Complex_wp>&, bool);

  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(VectReal_wp&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(VectComplex_wp&, bool, bool);

  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(VectReal_wp&, const GlobalGenericMatrix<Real_wp>&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(VectComplex_wp&, const GlobalGenericMatrix<Complex_wp>&, bool, bool);

  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(const SeldonTranspose&, VectReal_wp&, const GlobalGenericMatrix<Real_wp>&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(const SeldonTranspose&, VectComplex_wp&, const GlobalGenericMatrix<Complex_wp>&, bool, bool);

  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(Matrix<Real_wp, General, ColMajor>&, const GlobalGenericMatrix<Real_wp>&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(const SeldonTranspose&, Matrix<Real_wp, General, ColMajor>&, const GlobalGenericMatrix<Real_wp>&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(Matrix<Complex_wp, General, ColMajor>&, const GlobalGenericMatrix<Complex_wp>&, bool, bool);
  SELDON_EXTERN template void All_LinearSolver::ComputeSolution(const SeldonTranspose&, Matrix<Complex_wp, General, ColMajor>&, const GlobalGenericMatrix<Complex_wp>&, bool, bool);

  SELDON_EXTERN template class CondensationBlockSolver_Fem<Real_wp>;
  SELDON_EXTERN template class CondensationBlockSolver_Fem<Complex_wp>;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::AddDomains(VectReal_wp&, int) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::AddDomains(VectComplex_wp&, int) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::AddDomains(VectReal_wp&, int) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::AddDomains(VectComplex_wp&, int) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::ExtractVector(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::ExtractVector(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::ExtractVector(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::ExtractVector(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::ScatterVector(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::ScatterVector(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::ScatterVector(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::ScatterVector(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::CompressMatrixGen(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>::CompressMatrixGen(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::CompressMatrixGen(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::CompressMatrixGen(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::CompressMatrixGen(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>::CompressMatrixGen(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                 Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                 Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                 Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                    Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                    Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Real_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectReal_wp&, int, const IVect&,
                    Matrix<Real_wp, Symmetric, RowSymPacked>&, Matrix<Real_wp>&, Matrix<Real_wp>&, Matrix<Real_wp, Symmetric, RowSymPacked>&) const;
  
  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                 Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                 Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::ModifyRhsGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                 Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                    Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                    Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&) const;

  SELDON_EXTERN template void CondensationBlockSolver_Fem<Complex_wp>
  ::RecomposeSolGen(const SeldonTranspose&, VectComplex_wp&, int, const IVect&,
                    Matrix<Complex_wp, Symmetric, RowSymPacked>&, Matrix<Complex_wp>&, Matrix<Complex_wp>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&) const;
  
}
