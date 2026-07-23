#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Computation/AssembleMatrix.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class GlobalGenericMatrix<Real_wp>;
  SELDON_EXTERN template class GlobalGenericMatrix<Complex_wp>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const GlobalGenericMatrix<Complex_wp>&);

  SELDON_EXTERN template class CondensationBlockSolver_Base<Real_wp>;
  SELDON_EXTERN template class CondensationBlockSolver_Base<Complex_wp>;

  SELDON_EXTERN template void AssembleMatrix(VirtualMatrix<Real_wp>&, VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, VarComputationProblem_Base&, CondensationBlockSolver_Base<Real_wp>&, int, int);
  SELDON_EXTERN template void AssembleMatrix(VirtualMatrix<Complex_wp>&, VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, VarComputationProblem_Base&, CondensationBlockSolver_Base<Complex_wp>&, int, int);
  
  SELDON_EXTERN template void GetRowSum(VectReal_wp&, FemMatrixFreeClass_Base<Real_wp>&);
  SELDON_EXTERN template void GetRowSum(VectReal_wp&, FemMatrixFreeClass_Base<Complex_wp>&);

  SELDON_EXTERN template void GetRowColSum(VectReal_wp&, VectReal_wp&, FemMatrixFreeClass_Base<Real_wp>&);
  SELDON_EXTERN template void GetRowColSum(VectReal_wp&, VectReal_wp&, FemMatrixFreeClass_Base<Complex_wp>&);

  SELDON_EXTERN template void ScaleMatrix(FemMatrixFreeClass_Base<Real_wp>&, VectReal_wp&, VectReal_wp&);
  SELDON_EXTERN template void ScaleMatrix(FemMatrixFreeClass_Base<Complex_wp>&, VectReal_wp&, VectReal_wp&);
  
}
