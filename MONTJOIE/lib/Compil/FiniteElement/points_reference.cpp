#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/PointsReference.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class SetPoints<Dimension2>;
  SELDON_EXTERN template class SetPoints<Dimension3>;
  SELDON_EXTERN template class SetMatrices<Dimension2>;
  SELDON_EXTERN template class SetMatrices<Dimension3>;

  SELDON_EXTERN template class FjInverseProblem<Dimension2>;
  SELDON_EXTERN template class FjInverseProblem<Dimension3>;

}

