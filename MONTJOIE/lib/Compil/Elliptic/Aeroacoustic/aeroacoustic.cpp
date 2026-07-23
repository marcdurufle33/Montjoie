#include "MontjoieFlag.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticHeader.hxx"
#include "Elliptic/Aeroacoustic/MontjoieAeroacousticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Aeroacoustic/AeroAcoustic.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class AeroAcoustic_Base<Dimension2>;
  SELDON_EXTERN template class AeroAcoustic_Base<Dimension3>;

  SELDON_EXTERN template class VarHarmonic<HarmonicLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicLinearizedEulerEquation<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<StationaryLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<StationaryLinearizedEulerEquation<Dimension3> >;
  
  SELDON_EXTERN template class VarHarmonic<HarmonicAeroEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<HarmonicAeroEquation<Dimension3> >;

  SELDON_EXTERN template class VarHarmonic<AeroStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class VarHarmonic<AeroStationaryEquation<Dimension3> >;
  
  SELDON_EXTERN template class AeroAcoustic_Eq<HarmonicLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class AeroAcoustic_Eq<HarmonicLinearizedEulerEquation<Dimension3> >;

  SELDON_EXTERN template class AeroAcoustic_Eq<StationaryLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class AeroAcoustic_Eq<StationaryLinearizedEulerEquation<Dimension3> >;
  
  SELDON_EXTERN template class AeroAcoustic_Eq<HarmonicAeroEquation<Dimension2> >;
  SELDON_EXTERN template class AeroAcoustic_Eq<HarmonicAeroEquation<Dimension3> >;

  SELDON_EXTERN template class AeroAcoustic_Eq<AeroStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class AeroAcoustic_Eq<AeroStationaryEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HarmonicLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicLinearizedEulerEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<StationaryLinearizedEulerEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<StationaryLinearizedEulerEquation<Dimension3> >;
  
  SELDON_EXTERN template class EllipticProblem<HarmonicAeroEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HarmonicAeroEquation<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<AeroStationaryEquation<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<AeroStationaryEquation<Dimension3> >;

}
