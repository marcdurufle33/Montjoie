#ifndef MONTJOIE_FILE_MESH2D_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class T>
  inline Mesh<Dimension2, T>::Mesh() : Mesh_Base<Dimension2, T>(),
				       MeshBoundaries<Dimension2, T>()
  {
    SetGeometryOrder(1);
  }
  
  
  //! returns an edge of the mesh
  template<class T>
  inline Edge<Dimension2>& Mesh<Dimension2, T>::Boundary(int i)
  {
    return this->edges(i);
  }
  
  
  //! returns an edge of the mesh
  template<class T>
  inline const Edge<Dimension2>& Mesh<Dimension2, T>::Boundary(int i) const
  {
    return this->edges(i);
  }

  
  //! returns the number of edges in the mesh
  template<class T>
  inline int Mesh<Dimension2, T>::GetNbBoundary() const
  {
    return this->edges.GetM();
  }


  //! returns the global edge number of local edge i
  template<class T>
  inline int Mesh<Dimension2, T>::GetGlobalBoundaryNumber(int i) const
  {
    return this->GlobEdgeNumber_Subdomain(i);
  }
  

  //! returns the referenced edge i
  template<class T>
  inline Edge<Dimension2>& Mesh<Dimension2, T>::BoundaryRef(int i)
  {
    return this->edges_ref(i);
  }
  
  
  //! returns the referenced edge i
  template<class T>
  inline const Edge<Dimension2>& Mesh<Dimension2, T>::BoundaryRef(int i) const
  {
    return this->edges_ref(i);
  }
  
  
  //! returns the number of referenced edges in the mesh
  template<class T>
  inline int Mesh<Dimension2, T>::GetNbBoundaryRef() const
  {
    return this->edges_ref.GetM();
  }


  //! returns the number of referenced edges in the mesh
  template<class T>
  inline int Mesh<Dimension2, T>::GetNbEdgesRef() const
  {
    return this->edges_ref.GetM();
  }

  
  template<class T>
  inline const TriangleGeomReference& Mesh<Dimension2, T>
  ::GetTriangleReferenceElement() const
  {
    return triangle_reference; 
  }
  

  template<class T>
  inline const QuadrangleGeomReference& Mesh<Dimension2, T>
  ::GetQuadrilateralReferenceElement() const
  {
    return quadrangle_reference;
  }
  
  
  //! fills an array containing reference elements for edges
  template<class T>
  inline void Mesh<Dimension2, T>
  ::GetReferenceElementSurface(Vector<const ElementGeomReference<Dimension1>* >& ref) const
  {
    ref.Reallocate(1);
    ref(0) = &this->interval_reference;
  }
  

  //! returns true if the mesh contains both triangles and quadrilaterals
  template<class T>
  inline bool Mesh<Dimension2, T>::HybridMesh() const
  {
    if ((GetNbTriangles() > 0)&&(GetNbQuadrangles() > 0))
      return true;
    
    return false;
  }
  
  
  //! returns true if mesh is only made of triangles
  template<class T>
  inline bool Mesh<Dimension2, T>::IsOnlyTriangular() const
  {
    if (this->GetNbQuadrangles() == 0)
      return true;
    
    return false;
  }
  
  
  //! returns true if mesh is only made of quadrilaterals
  template<class T>
  inline bool Mesh<Dimension2, T>::IsOnlyQuadrilateral() const
  {
    if (this->GetNbTriangles() == 0)
      return true;
    
    return false;
  }
  
  
  //! no hexahedra in 2-D mesh
  template<class T>
  inline int Mesh<Dimension2, T>::GetNbHexahedra() const
  {
    return 0;
  }
  
  
  //! adds an element in the mesh (internal function)
  template<class T>
  inline void Mesh<Dimension2, T>::AddElement(int& nb, int& nb_edges, const Face<Dimension2>& e)
  {
    this->Element(nb) = e;
    for (int j = 0; j < e.GetNbEdges(); j++)
      this->AddEdgeOfElement(nb_edges, nb, j);
    
    nb++;
  }
  
  
  //! default constructor
  template<class T>
  inline SurfacicMesh<Dimension2, T>::SurfacicMesh()
    : SurfacicMesh_Base<Dimension2, T>()
  {
  }


}

#define MONTJOIE_FILE_MESH2D_INLINE_CXX
#endif
