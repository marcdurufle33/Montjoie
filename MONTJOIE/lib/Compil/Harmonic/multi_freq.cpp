#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/MultiFrequencyProblem.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void MultiFrequencyProblem::GatherVectorFrequency(Vector<int>&);
  SELDON_EXTERN template void MultiFrequencyProblem::GatherVectorFrequency(Vector<Real_wp>&);
  SELDON_EXTERN template void MultiFrequencyProblem::GatherVectorFrequency(Vector<Complex_wp>&);
  
}

