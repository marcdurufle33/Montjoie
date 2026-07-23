#include "MontjoieFlag.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticHeader.hxx"
#include "Elliptic/Aeroacoustic/MontjoieAeroacousticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Aeroacoustic/HarmonicGalbrun.cxx"
#include "Elliptic/Aeroacoustic/HarmonicGalbrunEquation.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VarGalbrunIndex_Base<Dimension2>;
  SELDON_EXTERN template class VarGalbrunIndex_Base<Dimension3>;

  SELDON_EXTERN template class VarGalbrun_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarGalbrun_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarGalbrun_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarGalbrun_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquation<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquation<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquation<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquation<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquation<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquationDG<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquationDG<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquationDG<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationDG<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquationDG<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquationDG<Dimension3> >;

}
