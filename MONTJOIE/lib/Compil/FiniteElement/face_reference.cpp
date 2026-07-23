#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/FaceReference.cxx"
#include "FiniteElement/Triangle/TriangleReference.cxx"
#include "FiniteElement/Quadrangle/QuadrangleReference.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class FaceReference<1>;
  SELDON_EXTERN template class FaceReference<2>;
  SELDON_EXTERN template class FaceReference<3>;

  SELDON_EXTERN template class ElementReference<Dimension2, 1>;
  SELDON_EXTERN template class ElementReference<Dimension2, 2>;
  SELDON_EXTERN template class ElementReference<Dimension2, 3>;

  SELDON_EXTERN template class TriangleReference<1>;
  SELDON_EXTERN template class TriangleReference<2>;
  SELDON_EXTERN template class TriangleReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const TriangleReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const TriangleReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const TriangleReference<3>&);

  SELDON_EXTERN template class QuadrangleReference<1>;
  SELDON_EXTERN template class QuadrangleReference<2>;
  SELDON_EXTERN template class QuadrangleReference<3>;

  SELDON_EXTERN template ostream& operator<<(ostream&, const QuadrangleReference<1>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const QuadrangleReference<2>&);
  SELDON_EXTERN template ostream& operator<<(ostream&, const QuadrangleReference<3>&);

  SELDON_EXTERN template void FaceReference<1>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void FaceReference<1>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void FaceReference<2>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void FaceReference<2>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;
  SELDON_EXTERN template void FaceReference<3>::ComputeNodalValuesGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void FaceReference<3>::ComputeNodalValuesGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeGradientBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeGradientBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhQuadratureSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhQuadratureSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhGen(const VectReal_wp&, VectReal_wp& ) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhSplitGen(const VectReal_wp&, VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhSplitGen(const VectComplex_wp&, VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhTransposeGen(const VectReal_wp&, VectReal_wp& ) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ApplyConstantRhTransposeGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddConstantMassMatrixGen(int, int, const Real_wp&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddConstantMassMatrixGen(int, int, const Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddConstantElemMatrixGen(int, int, const Real_wp&, const TinyMatrix<Real_wp, General, 2, 2>&,
											const R2&, const R2&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddConstantElemMatrixGen(int, int, const Complex_wp&, const TinyMatrix<Complex_wp, General, 2, 2>&,
											const R2_Complex_wp&, const R2_Complex_wp&, const TinyVector<bool, 4>&,
											VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrixGen(int, int, const VectReal_wp&, const Vector<TinyMatrix<Real_wp, General, 2, 2> >&,
			     const Vector<R2>&, const Vector<R2>&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrixGen(int, int, const VectComplex_wp&, const Vector<TinyMatrix<Complex_wp, General, 2, 2> >&,
			     const Vector<R2_Complex_wp>&, const Vector<R2_Complex_wp>&, const TinyVector<bool, 4>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddVariableMassMatrixGen(int, int, const Vector<Real_wp>& A, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddVariableMassMatrixGen(int, int, const Vector<Complex_wp>&, VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrixOpt(int, int, const VectReal_wp&, const Vector<TinyMatrix<Real_wp, General, 2, 2> >&,
			     const Vector<R2>&, const Vector<R2>&, const TinyVector<bool, 4>&, VirtualMatrix<Real_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrixOpt(int, int, const VectComplex_wp&, const Vector<TinyMatrix<Complex_wp, General, 2, 2> >&,
			     const Vector<R2_Complex_wp>&, const Vector<R2_Complex_wp>&, const TinyVector<bool, 4>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddVariableMassMatrixOpt(int, int, const Vector<Real_wp>& A, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::AddVariableMassMatrixOpt(int, int, const Vector<Complex_wp>&, VirtualMatrix<Complex_wp>&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeIntegralGradientRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeIntegralGradientRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeIntegralSurfaceGradientRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 1>::ComputeIntegralSurfaceGradientRef(const VectComplex_wp&, VectComplex_wp&, int) const;

  // H(curl) elements
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeCurlBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeCurlBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyShTransposeGen(int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyShTransposeGen(int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyNablaShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyNablaShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyNablaShTransposeGen(int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ApplyNablaShTransposeGen(int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddConstantMassMatrixGen(int, int, const Matrix2_2&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddConstantMassMatrixGen(int, int, const Matrix2_2_Complex_wp&, VirtualMatrix<Complex_wp>&) const;
        
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddConstantStiffnessMatrixGen(int, int, const Real_wp&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddConstantStiffnessMatrixGen(int, int, const Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddVariableStiffnessMatrixGen(int, int, const Vector<Real_wp>&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddVariableStiffnessMatrixGen(int, int, const Vector<Complex_wp>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddVariableMassMatrixGen(int, int, const Vector<Matrix2_2>& A, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::AddVariableMassMatrixGen(int, int, const Vector<Matrix2_2_Complex_wp>&, VirtualMatrix<Complex_wp>&) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeIntegralCurlRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeIntegralCurlRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeIntegralSurfaceCurlRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 2>::ComputeIntegralSurfaceCurlRef(const VectComplex_wp&, VectComplex_wp&, int) const;

  // H(div) elements
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeValueBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeValueBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeDivBoundaryGen(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeDivBoundaryGen(const VectComplex_wp&, VectComplex_wp&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyChGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyChGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyChTransposeGen(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyChTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyNablaShGen(const Real_wp&, int, const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ApplyNablaShGen(const Complex_wp&, int, const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;
  
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::AddConstantMassMatrixGen(int, int, const Matrix2_2&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::AddConstantMassMatrixGen(int, int, const Matrix2_2_Complex_wp&, VirtualMatrix<Complex_wp>&) const;
        
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::AddConstantStiffnessMatrixGen(int, int, const Real_wp&, VirtualMatrix<Real_wp>&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::AddConstantStiffnessMatrixGen(int, int, const Complex_wp&, VirtualMatrix<Complex_wp>&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeIntegralDivRef(const VectReal_wp&, VectReal_wp&) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeIntegralDivRef(const VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeIntegralSurfaceDivRef(const VectReal_wp&, VectReal_wp&, int) const;
  SELDON_EXTERN template void ElementReference<Dimension2, 3>::ComputeIntegralSurfaceDivRef(const VectComplex_wp&, VectComplex_wp&, int) const;

}
