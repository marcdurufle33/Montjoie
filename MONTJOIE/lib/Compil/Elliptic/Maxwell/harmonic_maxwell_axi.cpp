#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwellAxiHeader.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwellAxiInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/Axi/MaxwellAxiSymHarmonic.cxx"
#include "Elliptic/Maxwell/Axi/AxiSymHcurlMaxwell.cxx"

#include "Elliptic/Maxwell/Axi/StaticMaxwellAxi.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class HarmonicMaxwell_Axi<Real_wp>;
  SELDON_EXTERN template class HarmonicMaxwell_Axi<Complex_wp>;

  SELDON_EXTERN template class VarHarmonic<HarmonicMaxwellEquation_HcurlAxi>;
  //SELDON_EXTERN template class VarHarmonic<HarmonicMaxwellEquationAxiDG>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_HcurlAxi>;
  SELDON_EXTERN template class VarHarmonic<StaticMaxwellEquation_Axi_DG>;

  SELDON_EXTERN template class HarmonicMaxwellAxi_Eq<HarmonicMaxwellEquation_HcurlAxi>;
  //SELDON_EXTERN template class HarmonicMaxwellAxi_Eq<HarmonicMaxwellEquationAxiDG>;
  SELDON_EXTERN template class HarmonicMaxwellAxi_Eq<StaticMaxwellEquation_HcurlAxi>;
  SELDON_EXTERN template class HarmonicMaxwellAxi_Eq<StaticMaxwellEquation_Axi_DG>;

  SELDON_EXTERN template class HarmonicMaxwellAxi_Hcurl<HarmonicMaxwellEquation_HcurlAxi>;
  SELDON_EXTERN template class HarmonicMaxwellAxi_Hcurl<StaticMaxwellEquation_HcurlAxi>;

  SELDON_EXTERN template class EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>;
  //SELDON_EXTERN template class EllipticProblem<HarmonicMaxwellEquationAxiDG>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_HcurlAxi>;
  SELDON_EXTERN template class EllipticProblem<StaticMaxwellEquation_Axi_DG>;
  
}
