#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/3D/HarmonicMaxwell3D.cxx"
#include "Elliptic/Maxwell/3D/StaticMaxwell3D.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"

#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectHcurl3D.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class HarmonicMaxwell_3D<Real_wp>;
  SELDON_EXTERN template class HarmonicMaxwell_3D<Complex_wp>;

  SELDON_EXTERN template class VarHarmonic<HarmonicMaxwellEquation_3D>;
  SELDON_EXTERN template class VarHarmonic<HarmonicMaxwellEquation_3D_DG>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_3D>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_3D_DG>;

  SELDON_EXTERN template class HarmonicMaxwell_3D_Eq<HarmonicMaxwellEquation_3D>;
  SELDON_EXTERN template class HarmonicMaxwell_3D_Eq<HarmonicMaxwellEquation_3D_DG>;
  SELDON_EXTERN template class HarmonicMaxwell_3D_Eq<StaticMaxwellEquation_3D>;
  SELDON_EXTERN template class HarmonicMaxwell_3D_Eq<StaticMaxwellEquation_3D_DG>;

  SELDON_EXTERN template class HarmonicMaxwell_3D_Cplx<HarmonicMaxwellEquation_3D>;
  SELDON_EXTERN template class HarmonicMaxwell_3D_Cplx<HarmonicMaxwellEquation_3D_DG>;
  SELDON_EXTERN template class StaticMaxwell_3D<StaticMaxwellEquation_3D>;
  SELDON_EXTERN template class StaticMaxwell_3D<StaticMaxwellEquation_3D_DG>;

  SELDON_EXTERN template class EllipticProblem<HarmonicMaxwellEquation_3D>;
  SELDON_EXTERN template class EllipticProblem<HarmonicMaxwellEquation_3D_DG>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_3D>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_3D_DG>;

  SELDON_EXTERN template void GetExtrapolationU(EllipticProblem<StaticMaxwellEquation_3D_DG> const&, int, TinyVector<Vector<Real_wp>, 3> const&, TinyVector<Vector<Real_wp>, 3>&, TinyVector<Vector<Real_wp>, 3>&, ElementReference<Dimension3, 1> const&, TinyVector<Vector<Real_wp>, 3>&, bool);  

  SELDON_EXTERN template class CondensationBlockSolver_Maxwell3D<Real_wp, HarmonicMaxwellEquation_3D>;
  SELDON_EXTERN template class CondensationBlockSolver_Maxwell3D<Complex_wp, HarmonicMaxwellEquation_3D>;

  SELDON_EXTERN template class CondensationBlockSolver_Maxwell3D<Real_wp, StaticMaxwellEquation_3D>;
  SELDON_EXTERN template class CondensationBlockSolver_Maxwell3D<Complex_wp, StaticMaxwellEquation_3D>;
  
}
