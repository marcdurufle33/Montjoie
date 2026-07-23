#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell2D_Header.hxx"
#include "Instationary/MontjoieUnsteadyHeader.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell2D_Inline.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#include "Hyperbolic/Maxwell/2D/TimeMaxwell2D.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Maxwell/2D/TimeMaxwell2D.cxx"
#include "Hyperbolic/Maxwell/2D/ProdMatVectMaxwell2D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class TimeMaxwell_2D<TimeMaxwellEquation_2D>;
  SELDON_EXTERN template class HyperbolicProblem<TimeMaxwellEquation_2D>;

}
