#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/ElementGeomReference.cxx"
#include "FiniteElement/FaceGeomReference.cxx"
#include "FiniteElement/Triangle/TriangleGeomReference.cxx"
#include "FiniteElement/Quadrangle/QuadrangleGeomReference.cxx"
#include "vector/Vector.cxx"
#include "vector/Functions_Arrays.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template void Vector<TinyVector<int, 2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<int, 3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<int, 4> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<int, 8> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<Montjoie::Real_wp, 2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<Montjoie::Real_wp, 3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<Montjoie::Real_wp, 4> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyVector<Montjoie::Real_wp, 5> >::Resize(size_t);

  SELDON_EXTERN template void Vector<Montjoie::Matrix2_2>::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Matrix3_3>::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyMatrix<Montjoie::Real_wp, Symmetric, 3, 3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<TinyMatrix<Montjoie::Complex_wp, Symmetric, 3, 3> >::Resize(size_t);
  
  SELDON_EXTERN template void Vector<Vector<Montjoie::R2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Montjoie::R3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Montjoie::Real_wp, VectSparse> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Montjoie::Complex_wp, VectSparse> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<string> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<int> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Real_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Complex_wp> > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Vector<Real_wp> > > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Vector<Vector<Complex_wp> > > >::Resize(size_t);
  SELDON_EXTERN template void Vector<Vector<Matrix<Montjoie::Real_wp, General, RowSparse> > >::Resize(size_t);

  SELDON_EXTERN template void Vector<Matrix<bool> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Matrix<int> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Matrix<Montjoie::Real_wp, General, RowSparse> >::Resize(size_t);

  SELDON_EXTERN template void Vector<Array3D<int> >::Resize(size_t);  

  SELDON_EXTERN template void Vector<Montjoie::Edge<Montjoie::Dimension1> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Edge<Montjoie::Dimension2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Edge<Montjoie::Dimension3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Face<Montjoie::Dimension2> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Face<Montjoie::Dimension3> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::Volume>::Resize(size_t);  

  SELDON_EXTERN template void Vector<Montjoie::R2>::ReadText(std::string);
  SELDON_EXTERN template void Vector<Montjoie::R3>::ReadText(std::string);

  SELDON_EXTERN template void Vector<Montjoie::R2>::Read(std::string, bool);
  SELDON_EXTERN template void Vector<Montjoie::R3>::Read(std::string, bool);
  SELDON_EXTERN template void Vector<Montjoie::R2>::Read(istream&, bool);
  SELDON_EXTERN template void Vector<Montjoie::R3>::Read(istream&, bool);

  SELDON_EXTERN template void Vector<Montjoie::R2>::WriteText(std::string) const;
  SELDON_EXTERN template void Vector<Montjoie::R3>::WriteText(std::string) const;

  SELDON_EXTERN template void Vector<Montjoie::R2>::Write(std::string, bool) const;
  SELDON_EXTERN template void Vector<Montjoie::R3>::Write(std::string, bool) const;
  SELDON_EXTERN template void Vector<Montjoie::R2>::Write(ostream&, bool) const;
  SELDON_EXTERN template void Vector<Montjoie::R3>::Write(ostream&, bool) const;

  SELDON_EXTERN template void Vector<TinyMatrix<Montjoie::Real_wp, Symmetric, 2, 2> >::Write(std::string, bool) const;
  
  SELDON_EXTERN template void Sort(long, Vector<int>&, Vector<Complex_wp>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::R2>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::R3>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<Montjoie::R2>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<Montjoie::R3>&, Vector<int>&);
  SELDON_EXTERN template void Sort(Vector<Montjoie::R2>&, Vector<Montjoie::R2>&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::Edge<Montjoie::Dimension2> >&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::Edge<Montjoie::Dimension3> >&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::Face<Montjoie::Dimension2> >&, Vector<int>&);
  SELDON_EXTERN template void Sort(long, Vector<Montjoie::Face<Montjoie::Dimension3> >&, Vector<int>&);

  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<Montjoie::R2>&, Vector<int>&);
  SELDON_EXTERN template void RemoveDuplicate(long&, Vector<TinyVector<int, 2> >&, Vector<int>&);
  
}

namespace Montjoie
{
  SELDON_EXTERN template class ElementGeomReference_Base<Dimension2>;
  SELDON_EXTERN template class ElementGeomReference_Base<Dimension3>;

  SELDON_EXTERN template class ElementGeomReferenceContainer_Base<Dimension2>;
  SELDON_EXTERN template class ElementGeomReferenceContainer_Base<Dimension3>;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::SetMesh(Mesh<Dimension2>&);
  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::SetMesh(Mesh<Dimension3>&);
  SELDON_EXTERN template void ElementGeomReference_Base<Dimension3>::SetMesh(Mesh<Dimension3>&);

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::ComputeNodalGradient(const SetMatrices<Dimension2>&, const VectReal_wp&, VectR2&) const;
  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::ComputeNodalGradient(const SetMatrices<Dimension2>&, const VectComplex_wp&, VectR2_Complex_wp&) const;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension3>::ComputeNodalGradient(const SetMatrices<Dimension3>&, const VectReal_wp&, VectR3&) const;
  SELDON_EXTERN template void ElementGeomReference_Base<Dimension3>::ComputeNodalGradient(const SetMatrices<Dimension3>&, const VectComplex_wp&, VectR3_Complex_wp&) const;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::
  ComputeNodalGradient(SetMatrices<Dimension2> const&, Vector<Vector<Complex_wp> > const&, Vector<Vector<Complex_wp> >&) const;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension2>::
  ComputeNodalGradient(SetMatrices<Dimension2> const&, Vector<Vector<Real_wp> > const&, Vector<Vector<Real_wp> >&) const;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension3>::
  ComputeNodalGradient(SetMatrices<Dimension3> const&, Vector<Vector<Complex_wp> > const&, Vector<Vector<Complex_wp> >&) const;

  SELDON_EXTERN template void ElementGeomReference_Base<Dimension3>::
  ComputeNodalGradient(SetMatrices<Dimension3> const&, Vector<Vector<Real_wp> > const&, Vector<Vector<Real_wp> >&) const;

  SELDON_EXTERN template void TriangleGeomReference::FjElemNodalCurve(const VectR3&, SetPoints<Dimension3>&, const Mesh<Dimension3>&, int, const Face<Dimension3>&) const;

  SELDON_EXTERN template void QuadrangleGeomReference::FjElemNodalCurve(const VectR3&, SetPoints<Dimension3>&, const Mesh<Dimension3>&, int, const Face<Dimension3>&) const;
  
}

