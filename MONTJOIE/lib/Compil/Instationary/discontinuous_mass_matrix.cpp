#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/DiscontinuousMassMatrix.cxx"
#include "Instationary/DiscontinuousUnsteadyMassMatrix.cxx"
#include "vector/Vector.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void Vector<Montjoie::DiscontinuousBlockMassMatrix*>::Resize(size_t);  
}

namespace Montjoie
{
  SELDON_EXTERN template class DiscontinuousBlockDiagMassMatrix<2>;
  SELDON_EXTERN template class DiscontinuousBlockDiagMassMatrix<3>;
  
  SELDON_EXTERN template class DiscontinuousMassMatrix_Dim<Dimension2>;
  SELDON_EXTERN template class DiscontinuousMassMatrix_Dim<Dimension3>;

  SELDON_EXTERN template class DiscontinuousMassMatrix<Dimension2>;

  SELDON_EXTERN template class DiscontinuousMassMatrix<Dimension3>;

  SELDON_EXTERN template class DiscontinuousUnsteadyMassMatrix<Dimension2>;

  SELDON_EXTERN template class DiscontinuousUnsteadyMassMatrix<Dimension3>;

  SELDON_EXTERN template class DiscontinuousMassMatrixVol<Dimension2>;
  SELDON_EXTERN template class DiscontinuousMassMatrixSurf<Dimension2>;

  SELDON_EXTERN template class DiscontinuousSurfaceMassMatrix<Dimension2>;
  SELDON_EXTERN template class DiscontinuousSurfaceMassMatrix<Dimension3>;

}
