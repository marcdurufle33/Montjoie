#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/VolumeReference.cxx"
#include "FiniteElement/Tetrahedron/TetrahedronReference.cxx"
#include "FiniteElement/Pyramid/PyramidReference.cxx"
#include "FiniteElement/Wedge/WedgeReference.cxx"
#include "FiniteElement/Hexahedron/HexahedronReference.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VolumeReference<1>;
  SELDON_EXTERN template class VolumeReference<2>;
  SELDON_EXTERN template class VolumeReference<3>;

  SELDON_EXTERN template class ElementReference<Dimension3, 1>;
  SELDON_EXTERN template class ElementReference<Dimension3, 2>;
  SELDON_EXTERN template class ElementReference<Dimension3, 3>;

  SELDON_EXTERN template class TetrahedronReference<1>;
  SELDON_EXTERN template class TetrahedronReference<2>;
  SELDON_EXTERN template class TetrahedronReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const TetrahedronReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const TetrahedronReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const TetrahedronReference<3>&);

  SELDON_EXTERN template class PyramidReference<1>;
  SELDON_EXTERN template class PyramidReference<2>;
  SELDON_EXTERN template class PyramidReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const PyramidReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const PyramidReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const PyramidReference<3>&);

  SELDON_EXTERN template class WedgeReference<1>;
  SELDON_EXTERN template class WedgeReference<2>;
  SELDON_EXTERN template class WedgeReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const WedgeReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const WedgeReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const WedgeReference<3>&);

  SELDON_EXTERN template class HexahedronReference<1>;
  SELDON_EXTERN template class HexahedronReference<2>;
  SELDON_EXTERN template class HexahedronReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const HexahedronReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const HexahedronReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const HexahedronReference<3>&);

  SELDON_EXTERN template void VolumeReference<1>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void VolumeReference<1>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void VolumeReference<2>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void VolumeReference<2>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void VolumeReference<3>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void VolumeReference<3>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;

  // H^1 elements
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeGradientBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeGradientBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhQuadratureSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhQuadratureSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhGen(const VectReal_wp&, VectReal_wp& ) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhTransposeGen(const VectReal_wp&, VectReal_wp& ) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ApplyConstantRhTransposeGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddConstantMassMatrixGen(int, int, const Real_wp&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddConstantMassMatrixGen(int, int, const Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddConstantElemMatrixGen(int, int, const Real_wp&, const TinyMatrix<Real_wp, General, 3, 3>&,
											const R3&, const R3&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddConstantElemMatrixGen(int, int, const Complex_wp&, const TinyMatrix<Complex_wp, General, 3, 3>&,
											const R3_Complex_wp&, const R3_Complex_wp&, const TinyVector<bool, 4>&,
											VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrixGen(int, int, const VectReal_wp&, const Vector<TinyMatrix<Real_wp, General, 3, 3> >&,
			     const Vector<R3>&, const Vector<R3>&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrixGen(int, int, const VectComplex_wp&, const Vector<TinyMatrix<Complex_wp, General, 3, 3> >&,
			     const Vector<R3_Complex_wp>&, const Vector<R3_Complex_wp>&, const TinyVector<bool, 4>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddVariableMassMatrixGen(int, int, const Vector<Real_wp>&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddVariableMassMatrixGen(int, int, const Vector<Complex_wp>&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrixOpt(int, int, const VectReal_wp&, const Vector<TinyMatrix<Real_wp, General, 3, 3> >&,
			     const Vector<R3>&, const Vector<R3>&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>
  ::AddVariableElemMatrixOpt(int, int, const VectComplex_wp&, const Vector<TinyMatrix<Complex_wp, General, 3, 3> >&,
			     const Vector<R3_Complex_wp>&, const Vector<R3_Complex_wp>&, const TinyVector<bool, 4>&, VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddVariableMassMatrixOpt(int, int, const Vector<Real_wp>& A, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::AddVariableMassMatrixOpt(int, int, const Vector<Complex_wp>&, VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeIntegralGradientRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeIntegralGradientRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeIntegralSurfaceGradientRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 1>::ComputeIntegralSurfaceGradientRef(const VectComplex_wp&, VectComplex_wp&, int) const;

  // H(curl) elements
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeCurlBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeCurlBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddConstantMassMatrixGen(int, int, const Matrix3_3&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddConstantMassMatrixGen(int, int, const Matrix3_3_Complex_wp&, VirtualMatrix<Complex_wp>&) const;
        
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddConstantStiffnessMatrixGen(int, int, const Matrix3_3&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddConstantStiffnessMatrixGen(int, int, const Matrix3_3_Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddVariableStiffnessMatrixGen(int, int, const Vector<Matrix3_3>&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddVariableStiffnessMatrixGen(int, int, const Vector<Matrix3_3_Complex_wp>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddVariableMassMatrixGen(int, int, const Vector<Matrix3_3>& A, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::AddVariableMassMatrixGen(int, int, const Vector<Matrix3_3_Complex_wp>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeIntegralCurlRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeIntegralCurlRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeIntegralSurfaceCurlRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 2>::ComputeIntegralSurfaceCurlRef(const VectComplex_wp&, VectComplex_wp&, int) const;

  // VolumeHcurlReference
  SELDON_EXTERN template void VolumeHcurlReference::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void VolumeHcurlReference::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  // H(div) elements
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeDivBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeDivBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyNablaShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ApplyNablaShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;
  
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::AddConstantMassMatrixGen(int, int, const Matrix3_3&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::AddConstantMassMatrixGen(int, int, const Matrix3_3_Complex_wp&, VirtualMatrix<Complex_wp>&) const;
        
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::AddConstantStiffnessMatrixGen(int, int, const Real_wp&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::AddConstantStiffnessMatrixGen(int, int, const Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeIntegralDivRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeIntegralDivRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeIntegralSurfaceDivRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension3, 3>::ComputeIntegralSurfaceDivRef(const VectComplex_wp&, VectComplex_wp&, int) const;
  
}
