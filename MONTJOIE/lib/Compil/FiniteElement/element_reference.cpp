#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/ElementReference.cxx"
#include "FiniteElement/ElementReferenceBase.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class ElementReferenceType<Dimension2, 1>;
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
