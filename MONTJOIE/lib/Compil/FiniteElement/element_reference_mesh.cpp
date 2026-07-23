#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#include "FiniteElement/FiniteElementH1.hxx"
#include "FiniteElement/FiniteElementHcurl2D.hxx"
#include "FiniteElement/FiniteElementHcurl3D.hxx"
#include "FiniteElement/FiniteElementHdiv.hxx"

#include "FiniteElement/ElementReference.hxx"
#include "FiniteElement/FaceReference.hxx"
#include "FiniteElement/VolumeReference.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/ElementReferenceBase.cxx"

namespace Montjoie
{
  // we fill with empty methods (not used for mesh routines)
  size_t ElementReferenceType<Dimension2, 2>::GetMemorySize() const
  {
    return 0;
  }
  
  void ElementReferenceType<Dimension2, 2>::GetValuePhiOnQuadraturePoint(int k, VectR2& phi) const {}
  void ElementReferenceType<Dimension2, 2>::GetCurlPhiOnQuadraturePoint(int k, VectReal_wp& phi) const {}

  size_t ElementReferenceType<Dimension3, 2>::GetMemorySize() const
  {
    return 0;
  }

  void ElementReferenceType<Dimension3, 2>::GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const {}
  void ElementReferenceType<Dimension3, 2>::GetCurlPhiOnQuadraturePoint(int k, VectR3& phi) const
  {
  }

}

#endif

namespace Montjoie
{

  /*SELDON_EXTERN template class ElementReferenceType<Dimension2, 1>;
  SELDON_EXTERN template class ElementReferenceType<Dimension3, 1>;

  SELDON_EXTERN template void ElementReferenceType<Dimension2, 1>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 1>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension2, 1>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 1>&, const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 1>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 1>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 1>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 1>&, const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template class ElementReferenceType<Dimension2, 2>;
  SELDON_EXTERN template class ElementReferenceType<Dimension3, 2>;

  SELDON_EXTERN template void ElementReferenceType<Dimension2, 2>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 2>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension2, 2>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 2>&, const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 2>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 2>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 2>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 2>&, const VectComplex_wp&, VectComplex_wp&) const;
  
  SELDON_EXTERN template class ElementReferenceType<Dimension2, 3>;
  SELDON_EXTERN template class ElementReferenceType<Dimension3, 3>;

  SELDON_EXTERN template void ElementReferenceType<Dimension2, 3>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 3>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension2, 3>::ProjectQuadratureToDofGen(const ElementReference<Dimension2, 3>&, const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 3>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 3>&, const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReferenceType<Dimension3, 3>::ProjectQuadratureToDofGen(const ElementReference<Dimension3, 3>&, const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template class ElementReference_Dim<Dimension2>;
  SELDON_EXTERN template class ElementReference_Dim<Dimension3>;

  SELDON_EXTERN template void ElementReference_Dim<Dimension2>::ComputeValuesPhiQuadratureRef(const R2&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference_Dim<Dimension3>::ComputeValuesPhiQuadratureRef(const R3&, VectReal_wp&) const;
  */
  SELDON_EXTERN template void ElementReference_Base::ComputeIntegralRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference_Base::ComputeIntegralRef(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void ElementReference_Base::ComputeIntegralSurfaceRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ComputeIntegralSurfaceRef(const VectComplex_wp&, VectComplex_wp&, int) const;
  
  SELDON_EXTERN template void ElementReference_Base::ComputeGaussIntegralSurfaceGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ComputeGaussIntegralSurfaceGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyRhQuadratureGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyRhQuadratureGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyRhQuadratureTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyRhQuadratureTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyShTransposeGen(int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyShTransposeGen(int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyShGen(const Real_wp&, int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyShGen(const Complex_wp&, int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyShQuadratureTransposeGen(int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyShQuadratureTransposeGen(int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyShQuadratureGen(const Real_wp&, int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyShQuadratureGen(const Complex_wp&, int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShTransposeGen(int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShTransposeGen(int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShGen(const Real_wp&, int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShGen(const Complex_wp&, int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShQuadratureTransposeGen(int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShQuadratureTransposeGen(int, const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShQuadratureGen(const Real_wp&, int, const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference_Base::ApplyNablaShQuadratureGen(const Complex_wp&, int, const VectComplex_wp&, VectComplex_wp&, int) const;
  
}
