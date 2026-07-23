#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/TimeSchemes.cxx"
#include "Instationary/ImplicitTimeSchemes.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VirtualOdeSystem<Real_wp>;
  SELDON_EXTERN template class VirtualOdeSystem<Complex_wp>;

  SELDON_EXTERN template class VirtualTimeScheme<Real_wp>;
  SELDON_EXTERN template class VirtualTimeScheme<Complex_wp>;

  SELDON_EXTERN template class All_TimeScheme<Real_wp>;
  SELDON_EXTERN template class All_TimeScheme<Complex_wp>;
  
  SELDON_EXTERN template class TetaScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class TetaScheme_Iterator<Complex_wp>;

  SELDON_EXTERN template class SdirkScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class SdirkScheme_Iterator<Complex_wp>;
  
  SELDON_EXTERN template class DirkScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class DirkScheme_Iterator<Complex_wp>;

  SELDON_EXTERN template class GaussRungeKutta_Iterator<Real_wp>;
  SELDON_EXTERN template class GaussRungeKutta_Iterator<Complex_wp>;

  SELDON_EXTERN template class AdamsImplicit_Iterator<Real_wp>;
  SELDON_EXTERN template class AdamsImplicit_Iterator<Complex_wp>;

  SELDON_EXTERN template class MilneSimpson_Iterator<Real_wp>;
  SELDON_EXTERN template class MilneSimpson_Iterator<Complex_wp>;

  SELDON_EXTERN template class BackwardDifferentiation_Iterator<Real_wp>;
  SELDON_EXTERN template class BackwardDifferentiation_Iterator<Complex_wp>;

  SELDON_EXTERN template class ImplicitSymmetricMultistep_Iterator<Real_wp>;
  SELDON_EXTERN template class ImplicitSymmetricMultistep_Iterator<Complex_wp>;

  SELDON_EXTERN template class PadeScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class PadeScheme_Iterator<Complex_wp>;

  SELDON_EXTERN template class LinearSdirkScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class LinearSdirkScheme_Iterator<Complex_wp>;
  
}
