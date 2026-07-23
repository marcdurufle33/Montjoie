#include "MontjoieFlag.hxx"

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticAxisymHeader.hxx"
#include "Elliptic/Aeroacoustic/MontjoieAeroacousticAxisymInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Aeroacoustic/AxiSymGalbrun.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"
#include "Computation/ProdMatVectH1.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VarGalbrun_Axi<Real_wp>;
  SELDON_EXTERN template class VarGalbrun_Axi<Complex_wp>;

  SELDON_EXTERN template class VarHarmonic<HarmonicGalbrunEquationAxi>;
//SELDON_EXTERN template class VarHarmonic<StationaryGalbrunEquationAxi>;

  SELDON_EXTERN template class VarGalbrunAxi_Eq<HarmonicGalbrunEquationAxi>;
//SELDON_EXTERN template class VarGalbrunAxi_Eq<StationaryGalbrunEquationAxi>;

  SELDON_EXTERN template class EllipticProblem<HarmonicGalbrunEquationAxi>;
//SELDON_EXTERN template class EllipticProblem<StationaryGalbrunEquationAxi>;

}
