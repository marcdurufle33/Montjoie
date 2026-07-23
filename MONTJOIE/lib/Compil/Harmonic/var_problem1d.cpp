#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/VarProblem1D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void VarProblem_1D::SetPeriodicCondition(Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void VarProblem_1D::SetPeriodicCondition(Matrix<Complex_wp, General, BandedCol>&);

  SELDON_EXTERN template void VarProblem_1D::SetDirichletCondition(Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void VarProblem_1D::SetDirichletCondition(Matrix<Complex_wp, General, BandedCol>&);

  SELDON_EXTERN template void VarProblem_1D::AddMatrixFEM(VirtualMatrix<Real_wp>&, GlobalGenericMatrix<Real_wp>&);
  SELDON_EXTERN template void VarProblem_1D::AddMatrixFEM(VirtualMatrix<Complex_wp>&, GlobalGenericMatrix<Complex_wp>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 3, 3>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 3, 3>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 5, 5>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 5, 5>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 7, 7>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 7, 7>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 9, 9>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 9, 9>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 11, 11>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 11, 11>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 13, 13>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 13, 13>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 15, 15>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 15, 15>&);

  SELDON_EXTERN template void VarProblem_1D::GetStiffnessMatrix(TinyMatrix<Real_wp, General, 17, 17>&);
  SELDON_EXTERN template void VarProblem_1D::GetGradientMatrix(TinyMatrix<Real_wp, General, 17, 17>&);

  SELDON_EXTERN template void VarProblem_1D::AddDomains(VectReal_wp&, int) const;
  SELDON_EXTERN template void VarProblem_1D::AddDomains(VectComplex_wp&, int) const;

  SELDON_EXTERN template Real_wp VarProblem_1D::GetInterpolate(const VectReal_wp&, int, const Real_wp&, VectReal_wp&, bool, bool) const;
  SELDON_EXTERN template Complex_wp VarProblem_1D::GetInterpolate(const VectComplex_wp&, int, const Real_wp&, VectReal_wp&, bool, bool) const;

  SELDON_EXTERN template void VarProblem_1D::ComputeInterpolationU(const VectReal_wp&, const GridInterpolation<Dimension1>&, VectReal_wp&, bool, bool) const;
  SELDON_EXTERN template void VarProblem_1D::ComputeInterpolationU(const VectComplex_wp&, const GridInterpolation<Dimension1>&, VectComplex_wp&, bool, bool) const;

  SELDON_EXTERN template void VarProblem_1D::WriteDatas(VectReal_wp&);
  SELDON_EXTERN template void VarProblem_1D::WriteDatas(VectComplex_wp&);
  
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension1, Real_wp>;
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension1, Complex_wp>;

}
