#include "MontjoieFlag.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticHeader.hxx"
#include "Elliptic/Aeroacoustic/MontjoieAeroacousticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Aeroacoustic/HarmonicGalbrunEquation.cxx"
#include "Elliptic/Aeroacoustic/SymmetricHarmonicGalbrun.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VarGalbrunIndex_Base<Dimension2>;
  SELDON_EXTERN template class VarGalbrunIndex_Base<Dimension3>;

  SELDON_EXTERN template class ImpedanceABC<Real_wp, HarmonicGalbrunEquationSipg<Dimension2> >;
  SELDON_EXTERN template class ImpedanceABC<Real_wp, HarmonicGalbrunEquationSipg<Dimension3> >;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HarmonicGalbrunEquationSipg<Dimension2> >;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HarmonicGalbrunEquationSipg<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationSipg<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationSipg<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquationSipg<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquationSipg<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationH1>;
  SELDON_EXTERN template class VarHarmonic<GalbrunStationaryEquationH1>;
  
  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquationSipg<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquationSipg<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquationSipg<Dimension2> >;
  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquationSipg<Dimension3> >;

  SELDON_EXTERN template class VarGalbrun_Eq<HarmonicGalbrunEquationH1>;
  SELDON_EXTERN template class VarGalbrun_Eq<GalbrunStationaryEquationH1>;

  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationSipg<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationSipg<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquationSipg<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquationSipg<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationH1>;
  SELDON_EXTERN template class EllipticProblem<GalbrunStationaryEquationH1>;
  
  SELDON_EXTERN template class VolumetricSource_GalbrunH1<Real_wp>;
  SELDON_EXTERN template class VolumetricSource_GalbrunH1<Complex_wp>;

  SELDON_EXTERN template class ImpedanceABC<Real_wp, HarmonicGalbrunEquationH1>;
  SELDON_EXTERN template class ImpedanceABC<Complex_wp, HarmonicGalbrunEquationH1>;

}
