#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/NumberMap.cxx"
#include "Mesh/NumberMesh.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class MeshNumbering_Base<Real_wp>;
  
  SELDON_EXTERN template class MeshNumbering_Dim<Dimension2, Real_wp>;
  SELDON_EXTERN template class MeshNumbering_Dim<Dimension3, Real_wp>;
  
  SELDON_EXTERN template class MeshNumbering<Dimension2, Real_wp>;
  SELDON_EXTERN template class MeshNumbering<Dimension3, Real_wp>;

  SELDON_EXTERN template void InitStaticMeshData(const Dimension2& dim, const Real_wp&);
  SELDON_EXTERN template void InitStaticMeshData(const Dimension3& dim, const Real_wp&);

  SELDON_EXTERN template IVect NumberMap::GetDofNumberOnElement(const MeshNumbering<Dimension2>&, int) const;
  SELDON_EXTERN template IVect NumberMap::GetDofNumberOnElement(const MeshNumbering<Dimension3>&, int) const;
  SELDON_EXTERN template int NumberMap::GetNbDofElement(int, const Face<Dimension2>&) const;
  SELDON_EXTERN template int NumberMap::GetNbDofElement(int, const Face<Dimension3>&) const;

  SELDON_EXTERN template int NumberMap::GetGmshEntityNumber(const Edge<Dimension2>&, int);
  SELDON_EXTERN template int NumberMap::GetGmshEntityNumber(const Edge<Dimension3>&, int);
  SELDON_EXTERN template int NumberMap::GetGmshEntityNumber(const Face<Dimension2>&, int);
  SELDON_EXTERN template int NumberMap::GetGmshEntityNumber(const Face<Dimension3>&, int);

  SELDON_EXTERN template void NumberMap::GetLocalUnknownVector(MeshNumbering<Dimension2> const&, const VectReal_wp&, int, Vector<Vector<Real_wp> >&) const;
  SELDON_EXTERN template void NumberMap::GetLocalUnknownVector(MeshNumbering<Dimension2> const&, const VectComplex_wp&, int, Vector<Vector<Complex_wp> >&) const;

  SELDON_EXTERN template void NumberMap::GetLocalUnknownVector(MeshNumbering<Dimension3> const&, const VectReal_wp&, int, Vector<Vector<Real_wp> >&) const;
  SELDON_EXTERN template void NumberMap::GetLocalUnknownVector(MeshNumbering<Dimension3> const&, const VectComplex_wp&, int, Vector<Vector<Complex_wp> >&) const;

  // 1 unknown
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 1>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 1>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 1>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 1>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 1>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 1>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 1>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 1>&, int, VectComplex_wp&) const;


  // 2 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 2>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 2>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 2>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 2>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 2>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 2>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 2>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 2>&, int, VectComplex_wp&) const;


  // 3 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 3>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 3>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 3>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 3>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 3>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 3>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 3>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 3>&, int, VectComplex_wp&) const;


  // 4 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 4>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 4>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 4>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 4>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 4>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 4>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 4>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 4>&, int, VectComplex_wp&) const;


  // 5 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 5>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 5>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 5>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 5>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 5>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 5>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 5>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 5>&, int, VectComplex_wp&) const;


  // 6 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 6>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 6>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 6>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 6>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 6>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 6>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 6>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 6>&, int, VectComplex_wp&) const;


  // 7 unknowns
  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 7>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectReal_wp&, int, TinyVector<VectReal_wp, 7>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension2>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 7>&) const;

  SELDON_EXTERN template void NumberMap::
  GetLocalUnknownVector(const MeshNumbering<Dimension3>&, const VectComplex_wp&, int, TinyVector<VectComplex_wp, 7>&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Real_wp&, TinyVector<VectReal_wp, 7>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Real_wp&, TinyVector<VectReal_wp, 7>&, int, VectReal_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension2>&, const Complex_wp&, TinyVector<VectComplex_wp, 7>&, int, VectComplex_wp&) const;

  SELDON_EXTERN template void NumberMap::
  AddLocalUnknownVector(const MeshNumbering<Dimension3>&, const Complex_wp&, TinyVector<VectComplex_wp, 7>&, int, VectComplex_wp&) const;

  
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Complex_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<TinyVector<Real_wp, 1> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<TinyVector<Real_wp, 2> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<TinyVector<Real_wp, 3> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Complex_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<TinyVector<Real_wp, 1> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<TinyVector<Real_wp, 3> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension2> const&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension2> const&, Vector<Complex_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension2> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension2> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension3> const&, Vector<Real_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension3> const&, Vector<Complex_wp>&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension3> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalComponentVector(MeshNumbering<Dimension3> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::GetGlobalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::GetGlobalUnknownVector(MeshNumbering<Dimension2> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::GetGlobalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Vector<Real_wp> >&, int) const;
  SELDON_EXTERN template void NumberMap::GetGlobalUnknownVector(MeshNumbering<Dimension3> const&, Vector<Vector<Complex_wp> >&, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<Real_wp>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<TinyVector<Real_wp, 1> >&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<TinyVector<Real_wp, 3> >&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<R2>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<Complex_wp>&, int, int, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<Real_wp>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<TinyVector<Real_wp, 1> >&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<R3>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalRowMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<Complex_wp>&, int, int, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalColumnMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<Real_wp>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalColumnMatrix(MeshNumbering<Dimension2> const&, VirtualMatrix<Complex_wp>&, int, int, int) const;

  SELDON_EXTERN template void NumberMap::ModifyLocalColumnMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<Real_wp>&, int, int, int) const;
  SELDON_EXTERN template void NumberMap::ModifyLocalColumnMatrix(MeshNumbering<Dimension3> const&, VirtualMatrix<Complex_wp>&, int, int, int) const;

  SELDON_EXTERN template void NumberMap::ApplyOperatorFace(const SeldonTranspose&, int, int, const Edge<Dimension2>&, const Vector<Real_wp>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void NumberMap::ApplyOperatorFace(const SeldonTranspose&, int, int, const Face<Dimension3>&, const Vector<Real_wp>&, Vector<Real_wp>&) const;
  SELDON_EXTERN template void NumberMap::ApplyOperatorFace(const SeldonTranspose&, int, int, const Edge<Dimension2>&, const Vector<Complex_wp>&, Vector<Complex_wp>&) const;
  SELDON_EXTERN template void NumberMap::ApplyOperatorFace(const SeldonTranspose&, int, int, const Face<Dimension3>&, const Vector<Complex_wp>&, Vector<Complex_wp>&) const;
}
