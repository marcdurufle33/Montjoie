#include "MontjoieFlag.hxx"

#include "Hyperbolic/Elastic/MontjoieTimeElasticHeader.hxx"
#include "Hyperbolic/Elastic/MontjoieTimeElasticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/TimeElastic.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class TimeElastic_Base<TimeElasticEquation<Dimension2> >;
  SELDON_EXTERN template class TimeElastic_Base<TimeElasticEquation<Dimension3> >;

  SELDON_EXTERN template class ContinuousMassMatrixElastic<Dimension2>;
  SELDON_EXTERN template class ContinuousMassMatrixElastic<Dimension3>;

  SELDON_EXTERN template class HyperbolicProblem<TimeElasticEquation<Dimension2> >;
  SELDON_EXTERN template class HyperbolicProblem<TimeElasticEquation<Dimension3> >;

}
