#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/ProjectionOperator.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<TinyVector<Real_wp, 1> >&, Vector<TinyVector<Real_wp, 1> >&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<R2>&, Vector<R2>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<R3>&, Vector<R3>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<TinyVector<Complex_wp, 1> >&, Vector<TinyVector<Complex_wp, 1> >&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<R2_Complex_wp>&, Vector<R2_Complex_wp>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<R3_Complex_wp>&, Vector<R3_Complex_wp>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<Matrix2_2>&, Vector<Matrix2_2>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<Matrix3_3>&, Vector<Matrix3_3>&) const;

  SELDON_EXTERN template void FiniteElementProjector::TransposeProject(const Vector<R2>&, Vector<R2>&) const;
  SELDON_EXTERN template void FiniteElementProjector::TransposeProject(const Vector<R3>&, Vector<R3>&) const;

  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<Real_wp>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void FiniteElementProjector::Project(const Vector<Complex_wp>&, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void FiniteElementProjector::TransposeProject(const Vector<Real_wp>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void FiniteElementProjector::TransposeProject(const Vector<Complex_wp>&, Vector<Complex_wp>&) const;

  SELDON_EXTERN template class TensorizedProjector_Base<Dimension2>;
  SELDON_EXTERN template class TensorizedProjector_Base<Dimension3>;

  SELDON_EXTERN template class DenseProjector<Dimension1>;
  SELDON_EXTERN template class DenseProjector<Dimension2>;
  SELDON_EXTERN template class DenseProjector<Dimension3>;
  
  SELDON_EXTERN template void TensorizedProjector_Base<Dimension2>::SetChOperator(int m, const Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void TensorizedProjector_Base<Dimension3>::SetChOperator(int m, const Matrix<Real_wp, General, ArrayRowSparse>&);


  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(const Vector<ElementGeomReference<Dimension2>* >&, const Vector<VectR2>&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(const Vector<ElementGeomReference<Dimension2>* >&, const VectReal_wp&, const Vector<VectR2>&);

  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(const Vector<ElementGeomReference<Dimension3>* >&, const Vector<VectR3>&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(const Vector<ElementGeomReference<Dimension3>* >&, const VectReal_wp&, const Vector<VectR3>&);

  SELDON_EXTERN template void FiniteElementInterpolator::Project(const Vector<R2>&, Vector<R2>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::Project(const Vector<R3>&, Vector<R3>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::Project(const Vector<Matrix2_2>&, Vector<Matrix2_2>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::Project(const Vector<Matrix3_3>&, Vector<Matrix3_3>&, int) const;

  SELDON_EXTERN template void FiniteElementInterpolator::TransposeProject(const Vector<R2>&, Vector<R2>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::TransposeProject(const Vector<R3>&, Vector<R3>&, int) const;

  SELDON_EXTERN template void FiniteElementInterpolator::ProjectScalar(const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::ProjectScalar(const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::TransposeProjectScalar(const Vector<Real_wp>&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void FiniteElementInterpolator::TransposeProjectScalar(const Vector<Complex_wp>&, Vector<Complex_wp>&, int) const;

  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension1> const* > const&, VectReal_wp const&, Vector<VectReal_wp> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension1> const* > const&, TinyVector<VectReal_wp, 2> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension2> const* > const&, VectReal_wp const&, Vector<VectR2> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension2>* > const&, TinyVector<VectR2, 4> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension2> const* > const&, TinyVector<VectR2, 2> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension3> const* > const&, VectReal_wp const&, Vector<VectR3> const&);
  SELDON_EXTERN template void FiniteElementInterpolator::InitProjection(Vector<ElementGeomReference<Dimension3>* > const&, TinyVector<VectR3, 4> const&);

}

