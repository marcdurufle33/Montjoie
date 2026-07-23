#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwellAxiHeader.hxx"
#include "Instationary/MontjoieUnsteadyHeader.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwellAxiInline.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#include "Hyperbolic/Maxwell/Axi/TimeMaxwellAxi.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Maxwell/Axi/TimeMaxwellAxi.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class TimeMaxwell_Axi<TimeMaxwellEquation_Axi_DG>;
  SELDON_EXTERN template class TimeMaxwell_Axi<TimeMaxwellEquation_HcurlAxi>;

  SELDON_EXTERN template class HyperbolicProblem<TimeMaxwellEquation_Axi_DG>;
  SELDON_EXTERN template class HyperbolicProblem<TimeMaxwellEquation_HcurlAxi>;

}
