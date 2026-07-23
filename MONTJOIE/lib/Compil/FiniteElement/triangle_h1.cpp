#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/Triangle/TriangleClassical.cxx"
#include "FiniteElement/Triangle/TriangleQuasiLumped.cxx"
#include "FiniteElement/Triangle/TriangleHierarchic.cxx"
#include "FiniteElement/Triangle/TriangleDgMassLumped.cxx"
#include "FiniteElement/Triangle/TriangleDgOrtho.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template void TriangleClassical::ApplyNablaShTransposeGen(int, const VectReal_wp& feval, VectReal_wp& contrib, int) const;
  SELDON_EXTERN template void TriangleClassical::ApplyNablaShTransposeGen(int, const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

  SELDON_EXTERN template void TriangleHierarchic::ComputeProjectionDofGen(const VectReal_wp& feval, VectReal_wp& contrib) const;
  SELDON_EXTERN template void TriangleHierarchic::ComputeProjectionDofGen(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
  
}
