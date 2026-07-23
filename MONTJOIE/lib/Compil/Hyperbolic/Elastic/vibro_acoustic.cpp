#include "MontjoieFlag.hxx"

#include "Hyperbolic/Elastic/MontjoieTimeElasticHeader.hxx"
#include "Hyperbolic/Elastic/MontjoieTimeElasticInline.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#include "Hyperbolic/Elastic/VibroAcousticSession.hxx"
#include "Hyperbolic/Elastic/VibroAcousticSessionInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Elastic/VibroAcousticSession.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> >;

  SELDON_EXTERN template void VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> >::ComputeCbarre(Matrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void VibroAcousticLeapFrogScheme::SetInitialCondition(const Real_wp&, const Real_wp&, VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> >&);

  SELDON_EXTERN template void VibroAcousticLeapFrogScheme::Advance(const Real_wp&, int, VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> >&);
  
}
