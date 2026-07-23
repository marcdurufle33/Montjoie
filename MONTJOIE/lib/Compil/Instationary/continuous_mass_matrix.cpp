#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/ContinuousMassMatrix.cxx"
#include "Instationary/ContinuousUnsteadyMassMatrix.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class ContinuousMassMatrix<Dimension2>;
  SELDON_EXTERN template class ContinuousMassMatrix<Dimension3>;

  SELDON_EXTERN template class ContinuousUnsteadyMassMatrix_Dim<Dimension2>;
  SELDON_EXTERN template class ContinuousUnsteadyMassMatrix_Dim<Dimension3>;

  SELDON_EXTERN template class ContinuousUnsteadyMassMatrix<Dimension2>;
  SELDON_EXTERN template class ContinuousUnsteadyMassMatrix<Dimension3>;
  
}
