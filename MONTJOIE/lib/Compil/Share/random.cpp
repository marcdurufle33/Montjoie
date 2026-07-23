#include "MontjoieFlag.hxx"

#include "Share/MontjoieCommonHeader.hxx"
#include "Share/MontjoieCommonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Share/RandomGenerator.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template void EvaluateContinuedFraction(const Vector<int>&, Real_wp&);
  SELDON_EXTERN template void DecomposeContinuedFraction(const Real_wp&, const Real_wp&, Vector<int>&);
  SELDON_EXTERN template void GetNumeratorDenominator(const Vector<int>&, int&, int&);

  SELDON_EXTERN template void EvaluateContinuedFraction(const Vector<int64_t>&, Real_wp&);
  SELDON_EXTERN template void DecomposeContinuedFraction(const Real_wp&, const Real_wp&, Vector<int64_t>&);
  SELDON_EXTERN template void GetNumeratorDenominator(const Vector<int64_t>&, int64_t&, int64_t&);

}
