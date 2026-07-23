#ifndef MONTJOIE_FILE_MESH3D_BOUNDARIES_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor
  template<class T>
  inline MeshBoundaries<Dimension3, T>::MeshBoundaries()
    : MeshBoundaries_Base<Dimension3, T>()
  {
  }
  
  
  //! returns a nodal point of a curved edge (real position of the point)
  /*!
    \param[in] ne referenced edge number
    \param[in] k local number of the point in the edge
   */
  template<class T>
  inline const TinyVector<T, 3>& MeshBoundaries<Dimension3, T>
  ::GetPointInsideEdge(int ne, int k) const
  {
    return PointsEdgeRef(ne, k);
  }

  
  //! returns the number of nodal points strictly inside the face ne
  template<class T>
  inline int MeshBoundaries<Dimension3, T>::GetNbPointsInsideFace(int ne) const
  {
    return PointsFaceRef(ne).GetM();
  }


  //! returns a nodal point of a curved face (real position of the point)
  /*!
    \param[in] ne referenced face number
    \param[in] k local number of the point in the face
   */  
  template<class T>
  inline const TinyVector<T, 3>& MeshBoundaries<Dimension3, T>
  ::GetPointInsideFace(int ne, int k) const
  {
    return PointsFaceRef(ne)(k);
  }


  //! returns nodal points of a curved face (real position of the point)
  template<class T>
  inline const Vector<TinyVector<T, 3> >& MeshBoundaries<Dimension3, T>
  ::GetPointInsideFace(int ne) const
  {
    return PointsFaceRef(ne);
  }

  
  //! sets a nodal point of a curved edge (real position of the point)
  /*!
    \param[in] ne referenced edge number
    \param[in] k local number of the point in the edge
    \param[in] pt position of the point
   */
  template<class T>
  inline void MeshBoundaries<Dimension3, T>
  ::SetPointInsideEdge(int ne, int k, const TinyVector<T, 3>& pt)
  {
    PointsEdgeRef(ne, k) = pt;
  }

  
  //! sets a nodal point of a curved face (real position of the point)
  /*!
    \param[in] ne referenced face number
    \param[in] k local number of the point in the face
    \param[in] pt position of the point
   */  
  template<class T>
  inline void MeshBoundaries<Dimension3, T>
  ::SetPointInsideFace(int ne, int k, const TinyVector<T, 3>& pt)
  {
    PointsFaceRef(ne)(k) = pt;
  }

  
  //! sets nodal points of a curved face (real position of the points)
  /*!
    \param[in] ne referenced face number
    \param[in] k local number of the point in the face
    \param[in] pt position of the points
   */  
  template<class T>
  inline void MeshBoundaries<Dimension3, T>
  ::SetPointInsideFace(int ne, const Vector<TinyVector<T, 3> >& pt)
  {
    PointsFaceRef(ne) = pt;
  }
  
  
  //! returns triangular finite element
  template<class T>
  inline const TriangleGeomReference& MeshBoundaries<Dimension3, T>::GetSurfaceFiniteElement() const
  {
    return triangle_reference;
  }
  
  
  //! returns quadrilateral finite element
  template<class T>
  inline const QuadrangleGeomReference& MeshBoundaries<Dimension3, T>::GetSurfaceFiniteElement2() const
  {
    return quadrangle_reference;
  }
  
  
  //! Adds boundary faces
  template<class T>  
  inline void MeshBoundaries<Dimension3, T>::AddBoundaryEdges()
  {
    AddBoundaryFaces();
  }
  
  
  //! finding faces around a vertex,
  //! but only if the vertex has the lowest number among other vertices of each face
  template<class T>
  inline void MeshBoundaries<Dimension3, T>
  ::GetNeighboringBoundariesAroundVertex(int i, IVect& num) const
  {
    num.Clear();
    for (int e = this->GetLeafClass().head_faces_minv(i);
         e != -1; e = this->GetLeafClass().next_face(e))
      num.PushBack(e);
  }

  
  //! returns the nodal point i on the unit triangle
  template<class T>
  inline const R2& MeshBoundaries<Dimension3, T>::GetNodalPointInsideTriangle(int i) const
  {
    return triangle_reference.PointsNodalND(3*triangle_reference.GetGeometryOrder()+i);
  }
  
}

#define MONTJOIE_FILE_MESH3D_BOUNDARIES_INLINE_CXX
#endif
