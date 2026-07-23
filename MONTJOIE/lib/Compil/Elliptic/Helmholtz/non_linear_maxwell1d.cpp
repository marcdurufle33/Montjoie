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
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<2>;
  SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<4>;
  SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<6>;
  SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<8>;
  SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<10>;
  //SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<12>;
  //SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<14>;
  //SELDON_EXTERN template class TinyStiffnessOperatorMaxwell<16>;  
}
