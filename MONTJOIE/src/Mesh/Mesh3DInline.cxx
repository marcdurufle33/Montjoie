#ifndef MONTJOIE_FILE_MESH3D_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  template<class T>
  inline Mesh<Dimension3, T>::Mesh() : Mesh_Base<Dimension3, T>(),
				       MeshBoundaries<Dimension3, T>()
  {
    regular_mesh_with_slight_modification = 0;
    this->SetGeometryOrder(1);
  }
  
  
  //! returns the referenced edge i
  template<class T>
  inline Edge<Dimension3>& Mesh<Dimension3, T>::EdgeRef(int i)
  {
    return this->edges_ref(i);
  }
  
  
  //! returns the referenced edge i
  template<class T>
  inline const Edge<Dimension3>& Mesh<Dimension3, T>::EdgeRef(int i) const
  {
    return this->edges_ref(i);
  }
  
  
  //! returns the referenced face i
  template<class T>
  inline Face<Dimension3>& Mesh<Dimension3, T>::BoundaryRef(int i)
  {
    return this->faces_ref(i);
  }
  
  
  //! returns the referenced face i
  template<class T>
  inline const Face<Dimension3>& Mesh<Dimension3, T>::BoundaryRef(int i) const
  {
    return this->faces_ref(i);
  }
  

  //! returns the face i
  template<class T>
  inline Face<Dimension3>& Mesh<Dimension3, T>::Boundary(int i)
  {
    return this->faces(i);
  }
  
  
  //! returns the face i
  template<class T>
  inline const Face<Dimension3>& Mesh<Dimension3, T>::Boundary(int i) const
  {
    return this->faces(i);
  }
  
  
  //! returns the number of faces in the mesh
  template<class T>
  inline int Mesh<Dimension3, T>::GetNbBoundary() const
  {
    return this->faces.GetM();
  }


  //! returns the global face number of local face i
  template<class T>
  inline int Mesh<Dimension3, T>::GetGlobalBoundaryNumber(int i) const
  {
    return this->GlobFaceNumber_Subdomain(i);
  }
  
  
  //! returns the number of referenced faces in the mesh
  template<class T>
  inline int Mesh<Dimension3, T>::GetNbBoundaryRef() const
  {
    return this->faces_ref.GetM();
  }
  
  
  //! returns the number of referenced edges in the mesh
  template<class T>
  inline int Mesh<Dimension3, T>::GetNbEdgesRef() const
  {
    return this->edges_ref.GetM();
  }
  
  
  //! returns the number of faces in the mesh
  template<class T>
  inline int Mesh<Dimension3, T>::GetNbFaces() const
  {
    return this->faces.GetM();
  }
  

  template<class T>
  inline const TriangleGeomReference& Mesh<Dimension3, T>
  ::GetTriangleReferenceElement() const
  {
    return this->triangle_reference; 
  }
  
  
  template<class T>
  inline const QuadrangleGeomReference& Mesh<Dimension3, T>
  ::GetQuadrilateralReferenceElement() const
  {
    return this->quadrangle_reference; 
  }  
  
  
  //! fills an array containing reference elements for edges
  template<class T>
  inline void Mesh<Dimension3, T>
  ::GetReferenceElementSurface(Vector<const ElementGeomReference<Dimension2>* >& ref) const
  {
    ref.Reallocate(2);
    ref(0) = &this->triangle_reference; ref(1) = &this->quadrangle_reference;
  }
  
  
  //! default constructor
  template<class T>
  inline SurfacicMesh<Dimension3, T>::SurfacicMesh() : SurfacicMesh_Base<Dimension3, T>()
  {
  }
  
  
}

#define MONTJOIE_FILE_MESH3D_INLINE_CXX
#endif
