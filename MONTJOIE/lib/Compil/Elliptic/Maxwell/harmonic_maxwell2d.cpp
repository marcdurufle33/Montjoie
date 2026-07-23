#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell2D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell2D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/2D/HarmonicMaxwell2D.cxx"
#include "Elliptic/Maxwell/2D/StaticMaxwell2D.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"

#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ElementaryMatrixHcurl2D.cxx"
#include "Computation/ProdMatVectHcurl2D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class HarmonicMaxwell_2D<Real_wp>;
  SELDON_EXTERN template class HarmonicMaxwell_2D<Complex_wp>;

  SELDON_EXTERN template class VarHarmonic<HarmonicMaxwellEquation_2D>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_2D>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_2D_DG>;

  SELDON_EXTERN template class HarmonicMaxwell_2D_Eq<HarmonicMaxwellEquation_2D>;
  SELDON_EXTERN template class HarmonicMaxwell_2D_Eq<StaticMaxwellEquation_2D>;
  SELDON_EXTERN template class HarmonicMaxwell_2D_Eq<StaticMaxwellEquation_2D_DG>;

  SELDON_EXTERN template class HarmonicMaxwell_2D_Cplx<HarmonicMaxwellEquation_2D>;
  SELDON_EXTERN template class StaticMaxwell_2D<StaticMaxwellEquation_2D>;
  SELDON_EXTERN template class StaticMaxwell_2D<StaticMaxwellEquation_2D_DG>;

  SELDON_EXTERN template class ImpedanceABC<Real_wp, StaticMaxwellEquation_2D>;

  SELDON_EXTERN template class EllipticProblem<HarmonicMaxwellEquation_2D>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_2D>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_2D_DG>;

  SELDON_EXTERN template void StaticMaxwellEquation_2D_DG::GetTensorMass(EllipticProblem<StaticMaxwellEquation_2D_DG> const&, int, int, GlobalGenericMatrix<Real_wp> const&, int, TinyMatrix<Real_wp, General, 3, 3>&);
}
