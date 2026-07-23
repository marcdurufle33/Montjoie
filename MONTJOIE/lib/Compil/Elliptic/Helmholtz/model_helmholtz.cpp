#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/ThinSlotHelmholtzModel.cxx"
#include "Elliptic/Helmholtz/TransmissionModelHelmholtz.cxx"
#include "Elliptic/Helmholtz/ImpedanceHelmholtz.cxx"
#include "vector/Vector.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void Vector<Montjoie::SlotModelParameters<Montjoie::Dimension2> >::Resize(size_t);

}

namespace Montjoie
{
  
  SELDON_EXTERN template class SlotModelParameters<Dimension2>;

  SELDON_EXTERN template class VarTransmission_Helm<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarTransmission_Helm<Complex_wp, Dimension3>;
  SELDON_EXTERN template class VarTransmission_Helm<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarTransmission_Helm<Real_wp, Dimension3>;

  SELDON_EXTERN template class ImpedanceABC_Helm<Dimension2>;
  SELDON_EXTERN template class ImpedanceABC_Helm<Dimension3>;

  SELDON_EXTERN template class ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension3> >;

  SELDON_EXTERN template class VarGeneralizedImpedance_Helm<Dimension2>;
  SELDON_EXTERN template class VarGeneralizedImpedance_Helm<Dimension3>;
    
}
