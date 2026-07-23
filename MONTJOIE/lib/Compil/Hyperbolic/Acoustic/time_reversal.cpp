#include "MontjoieFlag.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Acoustic/TimeReversal.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class TimeReversalSource<Dimension2>;
  SELDON_EXTERN template class TimeReversalSource<Dimension3>;

  SELDON_EXTERN template class TimeAcousticReversal<Dimension2>;
  SELDON_EXTERN template class TimeAcousticReversal<Dimension3>;

}
