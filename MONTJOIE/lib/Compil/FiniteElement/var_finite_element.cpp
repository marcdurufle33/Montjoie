#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/VarFiniteElement.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension1, 1>;

  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension2, 1>;
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension2, 2>;
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension2, 3>;
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension2, 4>;
  
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension3, 1>;
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension3, 2>;
  SELDON_EXTERN template class VarFiniteElementEnum_Base<Dimension3, 3>;

  SELDON_EXTERN template class VarFiniteElementEnum<Dimension1, 1>;
  
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension2, 1>;
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension2, 2>;
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension2, 3>;
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension2, 4>;
  
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension3, 1>;
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension3, 2>;
  SELDON_EXTERN template class VarFiniteElementEnum<Dimension3, 3>;

  SELDON_EXTERN template class VarFiniteElement<Dimension1>;
  SELDON_EXTERN template class VarFiniteElement<Dimension2>;
  SELDON_EXTERN template class VarFiniteElement<Dimension3>;

}
