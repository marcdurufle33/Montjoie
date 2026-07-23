#include "MontjoieFlag.hxx"

#include "Inverse/MontjoieInverseHeader.hxx"
#include "Inverse/MontjoieInverseInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/MultiDimPhysicalIndex.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VariableParameterND_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class VariableParameterND_Base<Complex_wp, Dimension2>;
  
  SELDON_EXTERN template class SameMeshIndex<Real_wp, Dimension2>;
  SELDON_EXTERN template class SameMeshIndex<Complex_wp, Dimension2>;

  SELDON_EXTERN template class UniformParameterIndex<Real_wp, Dimension2>;
  SELDON_EXTERN template class UniformParameterIndex<Complex_wp, Dimension2>;
  
  SELDON_EXTERN template class RectangularBSplineIndex<Real_wp, Dimension2>;
  SELDON_EXTERN template class RectangularBSplineIndex<Complex_wp, Dimension2>;

  SELDON_EXTERN template class VariableParameterND_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class VariableParameterND_Base<Complex_wp, Dimension3>;
  
  SELDON_EXTERN template class SameMeshIndex<Real_wp, Dimension3>;
  SELDON_EXTERN template class SameMeshIndex<Complex_wp, Dimension3>;

  SELDON_EXTERN template class UniformParameterIndex<Real_wp, Dimension3>;
  SELDON_EXTERN template class UniformParameterIndex<Complex_wp, Dimension3>;
  
  SELDON_EXTERN template class RectangularBSplineIndex<Real_wp, Dimension3>;
  SELDON_EXTERN template class RectangularBSplineIndex<Complex_wp, Dimension3>;
  
}
