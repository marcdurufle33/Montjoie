#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/EvaluationCflWave.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class MeshGeometryCFL<Dimension2>;
  SELDON_EXTERN template class MeshGeometryCFL<Dimension3>;
  
  SELDON_EXTERN template void ComputeLocalTimeStep(VarInstationary_Dim<Dimension2>& var, VectReal_wp&);
  SELDON_EXTERN template void ComputeLocalTimeStep(VarInstationary_Dim<Dimension3>& var, VectReal_wp&);
  
}
