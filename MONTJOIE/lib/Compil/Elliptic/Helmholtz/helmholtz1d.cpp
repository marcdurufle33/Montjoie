#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class VarHelmholtz_1D<Real_wp>;
  SELDON_EXTERN template class VarHelmholtz_1D<Complex_wp>;

  SELDON_EXTERN template void VarHelmholtz_1D<Real_wp>::ComputeElementaryMatrixGen(int, IVect&, VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarHelmholtz_1D<Real_wp>::ComputeElementaryMatrixGen(int, IVect&, VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarHelmholtz_1D<Complex_wp>::ComputeElementaryMatrixGen(int, IVect&, VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&);

  SELDON_EXTERN template void VarHelmholtz_1D<Real_wp>::AddBoundaryTerms(VirtualMatrix<Real_wp>&, GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarHelmholtz_1D<Real_wp>::AddBoundaryTerms(VirtualMatrix<Complex_wp>&, GlobalGenericMatrix<Complex_wp>&);
  SELDON_EXTERN template void VarHelmholtz_1D<Complex_wp>::AddBoundaryTerms(VirtualMatrix<Complex_wp>&, GlobalGenericMatrix<Complex_wp>&);

}
