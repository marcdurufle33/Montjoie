#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/NonLinearOpticsProblem.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void RamanEffect::StoreNewField(const Real_wp&, const VectReal_wp&);
  SELDON_EXTERN template void RamanEffect::StoreNewField(const Real_wp&, const VectComplex_wp&);
}
