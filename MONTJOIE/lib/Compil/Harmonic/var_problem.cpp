#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/VarProblem.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarProblem<Dimension2>;
  SELDON_EXTERN template class VarProblem<Dimension3>;
  
}
