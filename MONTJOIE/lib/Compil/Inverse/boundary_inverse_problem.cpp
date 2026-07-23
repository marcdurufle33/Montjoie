#include "MontjoieFlag.hxx"

#include "Inverse/MontjoieInverseHeader.hxx"
#include "Inverse/MontjoieInverseInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Inverse/BoundaryInverseProblem.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class BoundaryInverseProblem_Base<Dimension2>;

  SELDON_EXTERN template class BoundaryInverseProblem<HelmholtzEquation<Dimension2> >;

  SELDON_EXTERN template class BoundaryInverseProblem<HarmonicElasticEquation<Dimension2> >;

  SELDON_EXTERN template class BoundaryInverseProblem_Base<Dimension3>;

  SELDON_EXTERN template class BoundaryInverseProblem<HelmholtzEquation<Dimension3> >;
  
  SELDON_EXTERN template class BoundaryInverseProblem<HarmonicElasticEquation<Dimension3> >;
}
