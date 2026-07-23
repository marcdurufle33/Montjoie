#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/LocalTimeSchemes.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class LocalTimeSteppingRK<Real_wp>;
  SELDON_EXTERN template class LocalTimeSteppingRK<Complex_wp>;

  SELDON_EXTERN template class LocalImexRK<Real_wp>;
  SELDON_EXTERN template class LocalImexRK<Complex_wp>;

}
