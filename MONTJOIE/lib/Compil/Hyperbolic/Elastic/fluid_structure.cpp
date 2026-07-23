#include "MontjoieFlag.hxx"

#include "Hyperbolic/Elastic/MontjoieTimeElasticHeader.hxx"
#include "Hyperbolic/Elastic/MontjoieTimeElasticInline.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#include "Hyperbolic/Elastic/FluidStructureInteraction.hxx"
#include "Hyperbolic/Elastic/FluidStructureInteractionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/FluidStructureInteraction.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class FluidStructureInteraction_Base<Dimension2>;
  SELDON_EXTERN template class FluidStructureInteraction_Base<Dimension3>;

  SELDON_EXTERN template class FluidStructureInteraction_Elas<Real_wp, Dimension2>;
  SELDON_EXTERN template class FluidStructureInteraction_Elas<Complex_wp, Dimension2>;
  SELDON_EXTERN template class FluidStructureInteraction_Elas<Real_wp, Dimension3>;
  SELDON_EXTERN template class FluidStructureInteraction_Elas<Complex_wp, Dimension3>;

  SELDON_EXTERN template class FluidStructureInteractionHarmonic<HelmholtzEquation<Dimension2>, HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class FluidStructureInteractionHarmonic<HelmholtzEquation<Dimension3>, HarmonicElasticEquation<Dimension3> >;
  
  SELDON_EXTERN template class FluidStructureInteraction<AcousticEquation<Dimension2>, TimeElasticEquation<Dimension2> >;
  SELDON_EXTERN template class FluidStructureInteraction<AcousticEquation<Dimension3>, TimeElasticEquation<Dimension3> >;

  SELDON_EXTERN template void FluidStructureLeapFrogScheme::SetInitialCondition(const Real_wp&, const Real_wp&, FluidStructureInteraction<AcousticEquation<Dimension3>, TimeElasticEquation<Dimension3> >&);

  SELDON_EXTERN template void FluidStructureLeapFrogScheme::Advance(const Real_wp&, int, FluidStructureInteraction<AcousticEquation<Dimension3>, TimeElasticEquation<Dimension3> >&);


  SELDON_EXTERN template void FluidStructureLeapFrogScheme::SetInitialCondition(const Real_wp&, const Real_wp&, FluidStructureInteraction<AcousticEquation<Dimension2>, TimeElasticEquation<Dimension2> >&);

  SELDON_EXTERN template void FluidStructureLeapFrogScheme::Advance(const Real_wp&, int, FluidStructureInteraction<AcousticEquation<Dimension2>, TimeElasticEquation<Dimension2> >&);


  SELDON_EXTERN template class AcousticAcousticInteraction<Dimension2>;
  SELDON_EXTERN template class AcousticAcousticInteraction<Dimension3>;

}
