#include "MontjoieFlag.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticAxisymHeader.hxx"
#include "Elliptic/Aeroacoustic/MontjoieAeroacousticAxisymInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Aeroacoustic/LinearizedEulerEquationAxisym.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  //SELDON_EXTERN template class VarAeroacoustic_Axi<Real_wp>;
  SELDON_EXTERN template class VarAeroacoustic_Axi<Complex_wp>;

  SELDON_EXTERN template class VarHarmonic<HarmonicLinearizedEulerEquationAxi>;
  SELDON_EXTERN template class VarHarmonic<HarmonicAeroacousticAxiEquation>;
  //SELDON_EXTERN template class VarHarmonic<StationaryAeroacousticAxiEquation>;
  //SELDON_EXTERN template class VarHarmonic<StationaryLinearizedEulerEquationAxi>;

  SELDON_EXTERN template class VarAeroacousticAxi_Eq<HarmonicLinearizedEulerEquationAxi>;
  SELDON_EXTERN template class VarAeroacousticAxi_Eq<HarmonicAeroacousticAxiEquation>;
  //SELDON_EXTERN template class VarAeroacousticAxi_Eq<StationaryAeroacousticAxiEquation>;
  //SELDON_EXTERN template class VarAeroacousticAxi_Eq<StationaryLinearizedEulerEquationAxi>;

  SELDON_EXTERN template class EllipticProblem<HarmonicLinearizedEulerEquationAxi>;
  SELDON_EXTERN template class EllipticProblem<HarmonicAeroacousticAxiEquation>;
  //SELDON_EXTERN template class EllipticProblem<StationaryAeroacousticAxiEquation>;
  //SELDON_EXTERN template class EllipticProblem<StationaryLinearizedEulerEquationAxi>;
  
}
