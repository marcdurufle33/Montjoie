#include "MontjoieFlag.hxx"

#include "Hyperbolic/Acoustic/MontjoieAcousticHeader.hxx"
#include "Hyperbolic/Acoustic/MontjoieAcousticInline.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"
#include "Hyperbolic/Acoustic/AxiSymAcoustic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Hyperbolic/Acoustic/TimeAcoustic.cxx"
#include "Hyperbolic/Acoustic/AxiSymAcoustic.cxx"
#include "Hyperbolic/Acoustic/ProdMatVectAcoustic.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class TimeAcoustic_Base<Dimension2>;
  SELDON_EXTERN template class TimeAcoustic_Base<Dimension3>;

  SELDON_EXTERN template class DiscontinuousMassMatrixAcoustic<Dimension2>;
  SELDON_EXTERN template class DiscontinuousMassMatrixAcoustic<Dimension3>;

  SELDON_EXTERN template class ContinuousMassMatrixAcoustic<Dimension2>;
  SELDON_EXTERN template class ContinuousMassMatrixAcoustic<Dimension3>;

  SELDON_EXTERN template class TimeAcoustic_Eq<AcousticEquation<Dimension2> >;
  SELDON_EXTERN template class TimeAcoustic_Eq<AcousticEquation<Dimension3> >;

  SELDON_EXTERN template class TimeAcoustic_Eq<AcousticEquationDG<Dimension2> >;
  SELDON_EXTERN template class TimeAcoustic_Eq<AcousticEquationDG<Dimension3> >;

  SELDON_EXTERN template class HyperbolicProblem<AcousticEquation<Dimension2> >;
  SELDON_EXTERN template class HyperbolicProblem<AcousticEquation<Dimension3> >;

  SELDON_EXTERN template class HyperbolicProblem<AcousticEquationDG<Dimension2> >;
  SELDON_EXTERN template class HyperbolicProblem<AcousticEquationDG<Dimension3> >;
  
}
