#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/DefineSourceHelmholtz.cxx"
#include "Elliptic/Helmholtz/LaplacianModalSource.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class ModalSourceLaplace_Dim<Dimension2>;
  SELDON_EXTERN template class ModalSourceLaplace_Dim<Dimension3>;

  SELDON_EXTERN template class ModalSourceLaplace<Real_wp, Dimension2>;
  SELDON_EXTERN template class ModalSourceLaplace<Complex_wp, Dimension2>;
  SELDON_EXTERN template class ModalSourceLaplace<Real_wp, Dimension3>;
  SELDON_EXTERN template class ModalSourceLaplace<Complex_wp, Dimension3>;

  SELDON_EXTERN template class DiffractedWaveSource_Helm<Dimension2>;
  SELDON_EXTERN template class DiffractedWaveSource_Helm<Dimension3>;

  SELDON_EXTERN template class TotalWaveSource_Helm<Dimension2>;
  SELDON_EXTERN template class TotalWaveSource_Helm<Dimension3>;

  SELDON_EXTERN template class DiffractedWaveSource_HelmDG<Dimension2>;
  SELDON_EXTERN template class DiffractedWaveSource_HelmDG<Dimension3>;

  SELDON_EXTERN template class TotalWaveSource_HelmDG<Dimension2>;
  SELDON_EXTERN template class TotalWaveSource_HelmDG<Dimension3>;

  SELDON_EXTERN template class IncidentWaveProjector_HelmDG<Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector_HelmDG<Dimension3>;

  SELDON_EXTERN template class VolumetricSource_HelmDG<Real_wp, Dimension2>;
  SELDON_EXTERN template class VolumetricSource_HelmDG<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VolumetricSource_HelmDG<Real_wp, Dimension3>;
  SELDON_EXTERN template class VolumetricSource_HelmDG<Complex_wp, Dimension3>;

}
