#include "MontjoieFlag.hxx"

#include "Output/MontjoieOutputHeader.hxx"
#include "Output/MontjoieOutputInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/MeshInterpolation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class MeshInterpolation_Base<Dimension2>;
  SELDON_EXTERN template class MeshInterpolation_Base<Dimension3>;

  SELDON_EXTERN template class MeshInterpolation<Dimension2>;
  SELDON_EXTERN template class MeshInterpolation<Dimension3>;
  
}
