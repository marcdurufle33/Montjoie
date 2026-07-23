#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Header.hxx"
#include "Instationary/MontjoieUnsteadyHeader.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Inline.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#include "Hyperbolic/Maxwell/3D/TimeMaxwell3D.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Maxwell/3D/TimeMaxwell3D.cxx"
#include "Hyperbolic/Maxwell/3D/ProdMatVectMaxwell3D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class TimeMaxwell_3D<TimeMaxwellEquation_3D>;
  SELDON_EXTERN template class TimeMaxwell_3D<TimeMaxwellEquation_3D_DG>;

  SELDON_EXTERN template class HyperbolicProblem<TimeMaxwellEquation_3D>;
  SELDON_EXTERN template class HyperbolicProblem<TimeMaxwellEquation_3D_DG>;

}
