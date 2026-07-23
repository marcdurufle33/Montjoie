#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class MaxwellGaussianSource<Real_wp>;
  SELDON_EXTERN template class MaxwellGaussianSource<Complex_wp>;

  SELDON_EXTERN template class CircularCoaxialModeField<Real_wp>;
  SELDON_EXTERN template class CircularCoaxialModeField<Complex_wp>;

  SELDON_EXTERN template class VolumetricSource_MaxwellHdg3D<Real_wp>;  
  SELDON_EXTERN template class VolumetricSource_MaxwellHdg3D<Complex_wp>;

}
