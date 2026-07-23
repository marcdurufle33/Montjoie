#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/VarHarmonicBase.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarHarmonic_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarHarmonic_Base<Real_wp, Dimension3>;

  SELDON_EXTERN template class VarHarmonic_Base<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarHarmonic_Base<Complex_wp, Dimension3>;
   
}
