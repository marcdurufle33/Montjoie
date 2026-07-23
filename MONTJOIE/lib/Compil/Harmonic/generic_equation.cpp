#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/GenericEquation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class GenericEquation_Base<Real_wp>;
  SELDON_EXTERN template class GenericEquation_Base<Complex_wp>;
    
}
