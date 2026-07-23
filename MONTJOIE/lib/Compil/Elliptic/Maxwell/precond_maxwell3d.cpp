#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/3D/PrecondMaxwell.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class PreconditionerHelmholtzMaxwell<Real_wp>;
  SELDON_EXTERN template class PreconditionerHelmholtzMaxwell<Complex_wp>;

  SELDON_EXTERN template class LowOrderPreconditioning_Maxwell<Real_wp>;
  SELDON_EXTERN template class LowOrderPreconditioning_Maxwell<Complex_wp>;

  SELDON_EXTERN template class MultigridPreconditioning_Maxwell<Real_wp>;
  SELDON_EXTERN template class MultigridPreconditioning_Maxwell<Complex_wp>;

  SELDON_EXTERN template class All_Preconditioner_Maxwell<Real_wp>;
  SELDON_EXTERN template class All_Preconditioner_Maxwell<Complex_wp>;
  
}
