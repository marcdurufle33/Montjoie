#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/TransmissionModel.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarTransmission_Base<Dimension2>;
  SELDON_EXTERN template class VarTransmission_Base<Dimension3>;

}
