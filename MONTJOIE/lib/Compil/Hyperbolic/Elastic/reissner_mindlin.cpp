#include "MontjoieFlag.hxx"

#include "Hyperbolic/Elastic/MontjoieTimeElasticHeader.hxx"
#include "Hyperbolic/Elastic/MontjoieTimeElasticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/ReissnerMindlin.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VarHarmonic<ReissnerMindlinEquation>;
  SELDON_EXTERN template class VarHarmonic<HarmonicReissnerMindlinEquation>;

  SELDON_EXTERN template class ReissnerMindlin_Eq<ReissnerMindlinEquation>;
  SELDON_EXTERN template class ReissnerMindlin_Eq<HarmonicReissnerMindlinEquation>;

  SELDON_EXTERN template class EllipticProblem<ReissnerMindlinEquation>;
  SELDON_EXTERN template class EllipticProblem<HarmonicReissnerMindlinEquation>;

}
