#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/VarInstationary.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void VarInstationary_Base::ExtractStiffnessMatrix(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void VarInstationary_Base::ExtractStiffnessMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template class VarInstationary_Dim<Dimension2>;
  SELDON_EXTERN template class VarInstationary_Dim<Dimension3>;

  SELDON_EXTERN template class VarInstationary_Fem<Dimension2>;
  SELDON_EXTERN template class VarInstationary_Fem<Dimension3>;

}
