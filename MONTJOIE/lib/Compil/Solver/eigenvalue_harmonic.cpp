#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/EigenvaluesHarmonic.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class EigenProblemMontjoie<Real_wp>;
  SELDON_EXTERN template class EigenProblemMontjoie<Complex_wp>;

  SELDON_EXTERN template void EigenProblemMontjoie<Real_wp>::ExpandVector(Vector<Real_wp> const&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void EigenProblemMontjoie<Real_wp>::ExpandVector(Vector<Complex_wp> const&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template class PolynomialEigenProblemMontjoie<Real_wp>;
  SELDON_EXTERN template class PolynomialEigenProblemMontjoie<Complex_wp>;
  
}
