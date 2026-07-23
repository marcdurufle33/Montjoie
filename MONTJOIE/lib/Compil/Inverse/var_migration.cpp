#include "MontjoieFlag.hxx"

#include "Inverse/MontjoieInverseHeader.hxx"
#include "Inverse/MontjoieInverseInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Inverse/VarMigration.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarMigration_Base<Dimension2>;
  SELDON_EXTERN template class VarMigration_Base<Dimension3>;
  
  SELDON_EXTERN template void VarMigration_Base<Dimension2>::ProjectExperimentData(VarHarmonic_Base<Complex_wp, Dimension2>&, const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);

  SELDON_EXTERN template void VarMigration_Base<Dimension2>::ProjectSimulationData(VarHarmonic_Base<Complex_wp, Dimension2>&, const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);

  SELDON_EXTERN template void VarMigration_Base<Dimension3>::ProjectExperimentData(VarHarmonic_Base<Complex_wp, Dimension3>&, const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);

  SELDON_EXTERN template void VarMigration_Base<Dimension3>::ProjectSimulationData(VarHarmonic_Base<Complex_wp, Dimension3>&, const Matrix<Complex_wp, General, ColMajor>&, Matrix<Complex_wp, General, ColMajor>&);
  
}
