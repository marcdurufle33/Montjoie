#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"
#include "Instationary/MontjoieTimeHeader.hxx"
#include "Instationary/MontjoieTimeInline.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.hxx"
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void KerrProblem::ApplyRamanEffect(const VectReal_wp&, VectReal_wp&, Real_wp);
  SELDON_EXTERN template void KerrProblem::ApplyRamanEffect(const VectReal_wp&, VectReal_wp&, VectReal_wp&, Real_wp);
  SELDON_EXTERN template void KerrProblem::ApplyRamanEffect(const VectReal_wp&, VectComplex_wp&, Real_wp);
  SELDON_EXTERN template void KerrProblem::ApplyRamanEffect(const VectReal_wp&, VectComplex_wp&, VectComplex_wp&, Real_wp);

  SELDON_EXTERN template void KerrProblem::ExtractNonZeroComponents(Real_wp&, Real_wp&, VectReal_wp&, const Real_wp&);
  SELDON_EXTERN template void KerrProblem::ExtractNonZeroComponents(Real_wp&, Real_wp&, VectComplex_wp&, const Real_wp&);

}
