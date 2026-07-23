#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/GeneralizedImpedanceModel.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template void VarGeneralizedImpedance_Base::FindDofsGibc(VarProblem<Dimension2>&, const IVect&);
  SELDON_EXTERN template void VarGeneralizedImpedance_Base::FindDofsGibc(VarProblem<Dimension3>&, const IVect&);
  
}
