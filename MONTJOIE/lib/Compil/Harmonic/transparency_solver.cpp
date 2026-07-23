#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/TransparencyCondition.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarComputationRCS_Base<Dimension2>;
  SELDON_EXTERN template class VarComputationRCS_Base<Dimension3>;

  SELDON_EXTERN template class TransparencySolver_Dim<Dimension2>;
  SELDON_EXTERN template class TransparencySolver_Dim<Dimension3>;

  SELDON_EXTERN template class TransparencySolver_Fem<Dimension2, 1, 1>;
  SELDON_EXTERN template class TransparencySolver_Fem<Dimension3, 1, 1>;

  SELDON_EXTERN template class TransparencySolver_Fem<Dimension2, 1, 2>;
  SELDON_EXTERN template class TransparencySolver_Fem<Dimension3, 1, 2>;

  SELDON_EXTERN template class TransparencySolver_Fem<Dimension2, 2, 2>;
  SELDON_EXTERN template class TransparencySolver_Fem<Dimension3, 2, 2>;

  SELDON_EXTERN template class TransparencySolver_Fem<Dimension2, 3, 3>;
  SELDON_EXTERN template class TransparencySolver_Fem<Dimension3, 3, 3>;
  
  //SELDON_EXTERN template void VarComputationProblem::AddMatrixWithBC(FemMatrixFreeClass_Base<Real_wp>&);
  
}
