#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/VarGeometryProblem.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarGeometryProblem<Dimension2>;
  SELDON_EXTERN template class VarGeometryProblem<Dimension3>;
  
}
