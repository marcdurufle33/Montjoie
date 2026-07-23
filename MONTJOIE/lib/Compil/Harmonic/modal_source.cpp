#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Source/ModalSource.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class ModalSourceBoundary_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class ModalSourceBoundary_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class ModalSourceBoundary_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class ModalSourceBoundary_Dim<Complex_wp, Dimension3>;

}
