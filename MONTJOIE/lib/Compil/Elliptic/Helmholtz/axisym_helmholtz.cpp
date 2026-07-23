#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ProdMatVectScalarH1.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarHelmholtz_Axi<Real_wp>;
  SELDON_EXTERN template class VarHelmholtz_Axi<Complex_wp>;

  SELDON_EXTERN template void VarHelmholtz_Axi<Complex_wp>
  ::ComputeElementaryMatrixHelmAxi(int, IVect&, VirtualMatrix<Complex_wp>&,
				   const GlobalGenericMatrix<Complex_wp>& nat_mat,
				   const ElementReference<Dimension2, 1>& Fb);
  
  SELDON_EXTERN template class VarHelmholtzAxi_Eq<HelmholtzEquationAxi>;
  SELDON_EXTERN template class VarHelmholtzAxi_Eq<LaplaceEquationAxi>;

  SELDON_EXTERN template class VarHelmholtzAxi_Cplx<HelmholtzEquationAxi>;
  SELDON_EXTERN template class VarHelmholtzAxi_Real<LaplaceEquationAxi>;

  SELDON_EXTERN template class ImpedanceABC<Real_wp, LaplaceEquationAxi>;
  
  SELDON_EXTERN template class ImpedanceGenericAxiHelm<Real_wp>;
  SELDON_EXTERN template class ImpedanceGenericAxiHelm<Complex_wp>;

  SELDON_EXTERN template class VolumetricSource_AxiHelm<Real_wp>;
  SELDON_EXTERN template class VolumetricSource_AxiHelm<Complex_wp>;

  SELDON_EXTERN template class EllipticProblem<HelmholtzEquationAxi>;
  SELDON_EXTERN template class EllipticProblem<LaplaceEquationAxi>;

}
