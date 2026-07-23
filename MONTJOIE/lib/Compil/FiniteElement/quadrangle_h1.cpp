#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/Quadrangle/QuadrangleGauss.cxx"
#include "FiniteElement/Quadrangle/QuadrangleLobatto.cxx"
#include "FiniteElement/Quadrangle/QuadrangleRadau.cxx"
#include "FiniteElement/Quadrangle/QuadrangleHierarchic.cxx"
#include "FiniteElement/Quadrangle/QuadrangleDgGauss.cxx"
#include "FiniteElement/Quadrangle/QuadrangleDgOrtho.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template void QuadrangleHierarchic::ComputeProjectionDofGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleHierarchic::ComputeProjectionDofGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

  SELDON_EXTERN template void QuadrangleGauss::ApplyChGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleGauss::ApplyChGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleGauss::ApplyChTransposeGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleGauss::ApplyChTransposeGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhTransposeGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhTransposeGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhBoundaryGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhBoundaryGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhBoundaryTransposeGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void QuadrangleLobatto::ApplyRhBoundaryTransposeGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
}
