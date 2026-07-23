#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#include "Elliptic/Helmholtz/HelmholtzHdiv.hxx"
#include "Elliptic/Helmholtz/HelmholtzHdivInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/HelmholtzHdiv.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class HelmholtzEquationHdiv_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class HelmholtzEquationHdiv_Base<Complex_wp, Dimension2>;
  SELDON_EXTERN template class HelmholtzEquationHdiv_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class HelmholtzEquationHdiv_Base<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarHelmholtzHdiv_Base<Dimension2>;
  SELDON_EXTERN template class VarHelmholtzHdiv_Base<Dimension3>;

  SELDON_EXTERN template class VarHarmonic<HelmholtzEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HelmholtzEquationHdiv<Dimension3> >;
  SELDON_EXTERN template class VarHarmonic<LaplaceEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<LaplaceEquationHdiv<Dimension3> >;

  SELDON_EXTERN template class VarHelmholtzHdiv_Eq<HelmholtzEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtzHdiv_Eq<HelmholtzEquationHdiv<Dimension3> >;
  SELDON_EXTERN template class VarHelmholtzHdiv_Eq<LaplaceEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtzHdiv_Eq<LaplaceEquationHdiv<Dimension3> >;

  SELDON_EXTERN template class DiffractedWaveSource_HelmHdiv<Dimension2>;
  SELDON_EXTERN template class DiffractedWaveSource_HelmHdiv<Dimension3>;

  SELDON_EXTERN template class IncidentWaveProjector_HelmHdiv<Dimension2>;
  SELDON_EXTERN template class IncidentWaveProjector_HelmHdiv<Dimension3>;

  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HelmholtzEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HelmholtzEquationHdiv<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HelmholtzEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HelmholtzEquationHdiv<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<LaplaceEquationHdiv<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<LaplaceEquationHdiv<Dimension3> >;

}
