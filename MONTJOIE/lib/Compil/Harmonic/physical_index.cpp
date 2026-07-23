#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/OneDimPhysicalIndex.cxx"
#include "Elliptic/PhysicalProperty.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VariableParameter1D_Base<Real_wp>;
  SELDON_EXTERN template class VariableParameter1D_Base<Complex_wp>;

  SELDON_EXTERN template class SplineInterpolationIndex<Real_wp>;
  SELDON_EXTERN template class SplineInterpolationIndex<Complex_wp>;

  SELDON_EXTERN template class BSplineInterpolationIndex<Real_wp>;
  SELDON_EXTERN template class BSplineInterpolationIndex<Complex_wp>;

  SELDON_EXTERN template class PolynomialInterpolationIndex<Real_wp>;
  SELDON_EXTERN template class PolynomialInterpolationIndex<Complex_wp>;

  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension1, Real_wp>;
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension1, Complex_wp>;

  SELDON_EXTERN template class RectangleCutOff<Dimension2, Real_wp>;
  SELDON_EXTERN template class RectangleCutOff<Dimension3, Real_wp>;
  
  SELDON_EXTERN template void RectangleCutOff<Dimension2, Real_wp>::ApplyCutOff(const R2&, Real_wp&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension2, Real_wp>::ApplyCutOff(const R2&, Complex_wp&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension3, Real_wp>::ApplyCutOff(const R3&, Real_wp&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension3, Real_wp>::ApplyCutOff(const R3&, Complex_wp&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension2, Real_wp>::ApplyCutOff(const R2&, Real_wp&, TinyVector<Real_wp, 2>&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension2, Real_wp>::ApplyCutOff(const R2&, Complex_wp&, TinyVector<Complex_wp, 2>&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension3, Real_wp>::ApplyCutOff(const R3&, Real_wp&, TinyVector<Real_wp, 3>&) const;
  SELDON_EXTERN template void RectangleCutOff<Dimension3, Real_wp>::ApplyCutOff(const R3&, Complex_wp&, TinyVector<Complex_wp, 3>&) const;
  
  SELDON_EXTERN template class PhysicalSinusoidalMedia<Dimension2, Real_wp>;
  SELDON_EXTERN template class PhysicalSinusoidalMedia<Dimension2, Complex_wp>;
  SELDON_EXTERN template class PhysicalSinusoidalMedia<Dimension3, Real_wp>;
  SELDON_EXTERN template class PhysicalSinusoidalMedia<Dimension3, Complex_wp>;

  SELDON_EXTERN template class PhysicalRegularMedia<Dimension2, Real_wp>;
  SELDON_EXTERN template class PhysicalRegularMedia<Dimension2, Complex_wp>;
  SELDON_EXTERN template class PhysicalRegularMedia<Dimension3, Real_wp>;
  SELDON_EXTERN template class PhysicalRegularMedia<Dimension3, Complex_wp>;

  SELDON_EXTERN template class PhysicalMeshMedia<Dimension2, Real_wp>;
  SELDON_EXTERN template class PhysicalMeshMedia<Dimension2, Complex_wp>;
  SELDON_EXTERN template class PhysicalMeshMedia<Dimension3, Real_wp>;
  SELDON_EXTERN template class PhysicalMeshMedia<Dimension3, Complex_wp>;

  //SELDON_EXTERN template class RadialVaryingMedia<Dimension1, Real_wp>;
  //SELDON_EXTERN template class RadialVaryingMedia<Dimension1, Complex_wp>;
  SELDON_EXTERN template class RadialVaryingMedia<Dimension2, Real_wp>;
  SELDON_EXTERN template class RadialVaryingMedia<Dimension2, Complex_wp>;
  SELDON_EXTERN template class RadialVaryingMedia<Dimension3, Real_wp>;
  SELDON_EXTERN template class RadialVaryingMedia<Dimension3, Complex_wp>;

  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension2, Real_wp>;
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension2, Complex_wp>;
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension3, Real_wp>;
  SELDON_EXTERN template class PhysicalVaryingMedia<Dimension3, Complex_wp>;
  
  SELDON_EXTERN template class GenericPhysicalIndice<Dimension2, Real_wp>;
  SELDON_EXTERN template class GenericPhysicalIndice<Dimension2, Complex_wp>;
  SELDON_EXTERN template class GenericPhysicalIndice<Dimension3, Real_wp>;
  SELDON_EXTERN template class GenericPhysicalIndice<Dimension3, Complex_wp>;

  SELDON_EXTERN template class ScalarPhysicalIndice<Dimension2, Real_wp>;
  SELDON_EXTERN template class ScalarPhysicalIndice<Dimension2, Complex_wp>;
  SELDON_EXTERN template class ScalarPhysicalIndice<Dimension3, Real_wp>;
  SELDON_EXTERN template class ScalarPhysicalIndice<Dimension3, Complex_wp>;

  SELDON_EXTERN template class VectorPhysicalIndice<Dimension2, 2, Real_wp>;
  SELDON_EXTERN template class VectorPhysicalIndice<Dimension2, 2, Complex_wp>;
  SELDON_EXTERN template class VectorPhysicalIndice<Dimension2, 3, Real_wp>;
  SELDON_EXTERN template class VectorPhysicalIndice<Dimension2, 3, Complex_wp>;
  SELDON_EXTERN template class VectorPhysicalIndice<Dimension3, 3, Real_wp>;
  SELDON_EXTERN template class VectorPhysicalIndice<Dimension3, 3, Complex_wp>;

  SELDON_EXTERN template class TensorPhysicalIndice<Dimension2, 2, Real_wp>;
  SELDON_EXTERN template class TensorPhysicalIndice<Dimension2, 2, Complex_wp>;
  SELDON_EXTERN template class TensorPhysicalIndice<Dimension2, 3, Real_wp>;
  SELDON_EXTERN template class TensorPhysicalIndice<Dimension2, 3, Complex_wp>;
  SELDON_EXTERN template class TensorPhysicalIndice<Dimension3, 3, Real_wp>;
  SELDON_EXTERN template class TensorPhysicalIndice<Dimension3, 3, Complex_wp>;

  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 2, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 2>&, TinyVector<Real_wp, 2>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 2, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 2>&, TinyVector<Complex_wp, 2>&) const;
    SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 2, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Complex_wp, 2>&, TinyVector<Complex_wp, 2>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 2, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Complex_wp, 2>&, TinyVector<Complex_wp, 2>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 2, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 2>&, TinyVector<Complex_wp, 2>&) const;

  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 3, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Real_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 3, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 3, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Complex_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension2, 3, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Complex_wp, 3>&) const;

  SELDON_EXTERN template void TensorPhysicalIndice<Dimension3, 3, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Real_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension3, 3, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
    SELDON_EXTERN template void TensorPhysicalIndice<Dimension3, 3, Real_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Complex_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension3, 3, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Complex_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
  SELDON_EXTERN template void TensorPhysicalIndice<Dimension3, 3, Complex_wp>::MltMatrix(const VarPhysicalProblem&, int, int, const TinyVector<Real_wp, 3>&, TinyVector<Complex_wp, 3>&) const;
  
  SELDON_EXTERN template class ElasticPhysicalIndice<Dimension2, 2, Real_wp>;
  SELDON_EXTERN template class ElasticPhysicalIndice<Dimension2, 2, Complex_wp>;
  SELDON_EXTERN template class ElasticPhysicalIndice<Dimension3, 3, Real_wp>;
  SELDON_EXTERN template class ElasticPhysicalIndice<Dimension3, 3, Complex_wp>;

  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Real_wp>::MltMatrixVar(int, int, const TinyVector<int, 2>&, const TinyVector<R2, 2>&, TinyVector<R2, 2>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Real_wp>::MltMatrixVar(int, int, const TinyVector<int, 2>&, const TinyVector<R2_Complex_wp, 2>&, TinyVector<R2_Complex_wp, 2>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Real_wp>::MltMatrixVar(int, int, const TinyVector<int, 3>&, const TinyVector<R3, 3>&, TinyVector<R3, 3>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Real_wp>::MltMatrixVar(int, int, const TinyVector<int, 3>&, const TinyVector<R3_Complex_wp, 3>&, TinyVector<R3_Complex_wp, 3>&) const;

  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Complex_wp>::MltMatrixVar(int, int, const TinyVector<int, 2>&, const TinyVector<R2_Complex_wp, 2>&, TinyVector<R2_Complex_wp, 2>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Complex_wp>::MltMatrixVar(int, int, const TinyVector<int, 3>&, const TinyVector<R3_Complex_wp, 3>&, TinyVector<R3_Complex_wp, 3>&) const;

  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Real_wp>::MltVectorVar(int, int, const TinyVector<int, 2>&, const TinyVector<Real_wp, 4>&, TinyVector<Real_wp, 4>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Real_wp>::MltVectorVar(int, int, const TinyVector<int, 3>&, const TinyVector<Real_wp, 9>&, TinyVector<Real_wp, 9>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Complex_wp>::MltVectorVar(int, int, const TinyVector<int, 2>&, const TinyVector<Complex_wp, 4>&, TinyVector<Complex_wp, 4>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Complex_wp>::MltVectorVar(int, int, const TinyVector<int, 3>&, const TinyVector<Complex_wp, 9>&, TinyVector<Complex_wp, 9>&) const;

  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Real_wp>::MltVectorPoint(int, const Vector<Real_wp>&, const TinyVector<Real_wp, 2>&, const TinyVector<Real_wp, 4>&, TinyVector<Real_wp, 4>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Real_wp>::MltVectorPoint(int, const Vector<Real_wp>&, const TinyVector<Real_wp, 3>&, const TinyVector<Real_wp, 9>&, TinyVector<Real_wp, 9>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Complex_wp>::MltVectorPoint(int, const Vector<Real_wp>&, const TinyVector<Real_wp, 2>&, const TinyVector<Complex_wp, 4>&, TinyVector<Complex_wp, 4>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Complex_wp>::MltVectorPoint(int, const Vector<Real_wp>&, const TinyVector<Real_wp, 3>&, const TinyVector<Complex_wp, 9>&, TinyVector<Complex_wp, 9>&) const;

  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Real_wp>::GetCoefficient(const VarPhysicalProblem&, int, int, TinyMatrix<TinyMatrix<Real_wp, General, 2, 2>, General, 2, 2>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension2, 2, Complex_wp>::GetCoefficient(const VarPhysicalProblem&, int, int, TinyMatrix<TinyMatrix<Complex_wp, General, 2, 2>, General, 2, 2>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Real_wp>::GetCoefficient(const VarPhysicalProblem&, int, int, TinyMatrix<TinyMatrix<Real_wp, General, 3, 3>, General, 3, 3>&) const;
  SELDON_EXTERN template void ElasticPhysicalIndice<Dimension3, 3, Complex_wp>::GetCoefficient(const VarPhysicalProblem&, int, int, TinyMatrix<TinyMatrix<Complex_wp, General, 3, 3>, General, 3, 3>&) const;
}
