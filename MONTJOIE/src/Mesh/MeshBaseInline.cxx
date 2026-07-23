#ifndef MONTJOIE_FILE_MESH_BASE_INLINE_CXX

namespace Montjoie
{
    
  /****************
   * MESH_BASE    *
   ****************/
  
  
  //! returns reference to leaf class
  template<class Dimension, class T>
  inline Mesh<Dimension, T>& Mesh_Base<Dimension, T>::GetLeafClass()
  {
    return static_cast<Mesh<Dimension, T>& >(*this);
  }
  
  
  //! returns reference to leaf class
  template<class Dimension, class T>
  inline const Mesh<Dimension, T>& Mesh_Base<Dimension, T>::GetLeafClass() const
  {
    return static_cast<const Mesh<Dimension, T>& >(*this);
  }
  

  //! returns the number of elements
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNbElt() const
  {
    return this->elements.GetM();
  }
  
    
  //! returns the number of vertices
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNbVertices() const
  {
    return this->Vertices.GetM();
  }
  
  
  //! return the number of edges
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNbEdges() const
  {
    return this->edges.GetM();
  }
  
  
  //! returns the original references of faces at interface between processors
  template<class Dimension, class T>
  inline const IVect& Mesh_Base<Dimension, T>::GetOriginalNeighborReference() const
  {
    return reference_original_neighbor;
  }
  
  
  //! sets the original references of faces at interface between processors
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::SetOriginalNeighborReference(const IVect& ref)
  {
    reference_original_neighbor = ref;
  }
  
  
  //! sets directory name where meshes are searched when calling ConstructMesh
  template <class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::SetPathName(const string& path)
  {
    mesh_path = path;
  }
  
  
  //! returns directory name where meshes are searched when calling ConstructMesh
  template <class Dimension, class T>
  const string& Mesh_Base<Dimension, T>::GetPathName() const
  {
    return mesh_path;
  }

  
  //! returns true if the element i is affine
  template<class Dimension, class T>
  inline bool Mesh_Base<Dimension, T>::IsElementAffine(int i) const
  {
    return affine_element(i);
  }
  
  
  //! changes the size of array containing elements, previous elements are lost
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::ReallocateElements(int i)
  {
    elements.Reallocate(i);
  }

  
  //! changes the size of array containing elements, previous elements are kept
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::ResizeElements(int i)
  {
    elements.Resize(i);
  }


  //! removes all the elements contained in the mesh
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::ClearElements()
  {
    elements.Clear();
  }
  
  
  //! gives access to an element of the mesh
  template<class Dimension, class T>
  inline const typename Dimension::Element_Rn& Mesh_Base<Dimension, T>::Element(int i) const
  {
    return this->elements(i);
  }

  
  //! gives access to an element of the mesh
  template<class Dimension, class T>
  inline typename Dimension::Element_Rn& Mesh_Base<Dimension, T>::Element(int i)
  {
    return this->elements(i);
  }
  
  
  //! return type of element num_elt
  /*!
    \param[in] num_elt element number
   */
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetTypeElement(int num_elt) const
  {
    return this->elements(num_elt).GetHybridType();
  }

  
  //! gives access to an edge of the mesh
  template<class Dimension, class T>
  inline Edge<Dimension>& Mesh_Base<Dimension, T>::GetEdge(int i)
  {
    return this->edges(i);
  }
  
  
  //! gives access to an edge of the mesh
  template<class Dimension, class T>
  inline const Edge<Dimension>& Mesh_Base<Dimension, T>::GetEdge(int i) const
  {
    return this->edges(i);
  
  }
  
  
  //! removes all the vertices
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::ClearVertices()
  {
    Vertices.Clear();
  }
  
  
  //! gives access to a vertex of the mesh
  template<class Dimension, class T>
  inline typename Dimension::R_N& Mesh_Base<Dimension, T>::Vertex(int i)
  {
    return this->Vertices(i);
  }
  
  
  //! gives access to a vertex of the mesh
  template<class Dimension, class T>
  inline const typename Dimension::R_N& Mesh_Base<Dimension, T>::Vertex(int i) const
  {
    return this->Vertices(i);
  }

  
  //! gives access to vertices
  template<class Dimension, class T>
  inline Vector<typename Dimension::R_N>& Mesh_Base<Dimension, T>::Vertex()
  {
    return this->Vertices;
  }


  //! gives access to vertices
  template<class Dimension, class T>
  inline const Vector<typename Dimension::R_N>& Mesh_Base<Dimension, T>::Vertex() const
  {
    return this->Vertices;
  }


  //! returns the minimum of x-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetXmin() const
  {
    return xmin_;
  }
  
  
  //! returns the maximum of x-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetXmax() const
  {
    return xmax_;
  }
  

  //! returns the minimum of y-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetYmin() const
  {
    return ymin_;
  }
  
  
  //! returns the maximum of y-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetYmax() const
  {
    return ymax_;
  }
  
  
  //! returns the minimum of z-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetZmin() const
  {
    return zmin_;
  }
  
  
  //! returns the maximum of z-coordinate among all vertices of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetZmax() const
  {
    return zmax_;
  }

  
  template<class Dimension, class T>
  inline const ElementGeomReference<Dimension>& Mesh_Base<Dimension, T>
  ::GetReferenceElement(int i) const
  {
    return *reference_element(i);
  }
  
  
  template<class Dimension, class T>
  inline const Vector<ElementGeomReference<Dimension>* >&
  Mesh_Base<Dimension, T>::GetReferenceElement() const
  {
    return reference_element;
  }
  

  /**********************************************************************
   * Methods to use mesh for localization of points (GridInterpolation) *
   **********************************************************************/
  
  
  //! computes nodal points in the real space (after transformation Fj)
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>
  ::FjElemNodal(const VectR_N& s, SetPoints<Dimension>& pts,
		const Mesh<Dimension>& mesh, int i) const
  {
    reference_element(i)->FjElemNodal(s, pts, mesh, i);
  }
    
  
  //! computes jacobian matrices DFj for nodal points of the elements
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>
  ::DFjElemNodal(const VectR_N& s, const SetPoints<Dimension>& pts,
		 SetMatrices<Dimension>& mat,
		 const Mesh<Dimension>& mesh, int i) const
  {
    reference_element(i)->DFjElemNodal(s, pts, mat, mesh, i);
  }
  
  
  //! applies transform Fj to a single point of an element
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>
  ::Fj(const VectR_N& s, const SetPoints<Dimension>& pts,
       const R_N& pt_loc, R_N& pt_glob, const Mesh<Dimension>& mesh, int i) const
  {
    reference_element(i)->Fj(s, pts, pt_loc, pt_glob, mesh, i);
  }
  
  
  //! computes jacobian matrix DFj for a single point of an element
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>
  ::DFj(const VectR_N& s, const SetPoints<Dimension>& pts,
        const R_N& pt_loc, MatrixN_N& dfj, const Mesh<Dimension>& mesh, int i) const
  {
    reference_element(i)->DFj(s, pts, pt_loc, dfj, mesh, i);
  }
  
  
  //! returns the number of nodal points on element i
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNbPointsNodalElt(int i) const
  {
    return reference_element(i)->GetNbPointsNodalElt();
  }
  

  //! returns the number of nodal points on an edge
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNbPointsNodalBoundary(int i, int num_loc) const
  {
    return reference_element(i)->GetNbNodalBoundary(num_loc);
  }
  
  
  //! returns the number of a nodal point of an element
  /*!
    \param[in] i element number
    \param[in] num_loc local position of the edge
    \param[in] k local number of the nodal point on the edge
   */
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::GetNodalNumber(int i, int num_loc, int k) const
  {
    return reference_element(i)->GetNodalNumber(num_loc, k);
  }
 
  
  //! computes the normale and lineic element ds from DFj^{-1}
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::GetNormale(const MatrixN_N& dfjm1, R_N& normale,
						  Real_wp& dsj, int i, int n) const
  {
    reference_element(i)->GetNormale(dfjm1, normale, dsj, n);
  }
  
  
  //! returns the distance of a point of the reference element to the border
  template<class Dimension, class T>
  inline Real_wp Mesh_Base<Dimension, T>::GetDistanceToBoundary(const R_N& pt, int i) const
  {
    return reference_element(i)->GetDistanceToBoundary(pt);
  }
    
  
  //! projects a point outside the reference element to the border of this element
  template<class Dimension, class T>
  inline int Mesh_Base<Dimension, T>::ProjectPointOnBoundary(R_N& pt, int i) const
  {
    return reference_element(i)->ProjectPointOnBoundary(pt);
  }
  
  
  //! wrapping polygon
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::GetBoundingBox(const VectR_N& s,
						      const Real_wp& coef, VectR_N& box,
						      TinyVector<R_N, 2>& enveloppe) const
  {
    return ElementGeomReference<Dimension>::GetBoundingBox(s, coef, box, enveloppe);
  }
  
  
  //! evaluates nodal shape functions on a point of element i
  template<class Dimension, class T>
  inline void Mesh_Base<Dimension, T>::ComputeValuesPhiNodalRef(int i, const R_N& pt_loc, VectReal_wp& phi) const
  {
    return reference_element(i)->ComputeValuesPhiNodalRef(pt_loc, phi);
  }
  

  //! computation of grad(u) on nodal points
  template<class Dimension, class T> template<class Vector1, class Vector2>
  inline void Mesh_Base<Dimension, T>
  ::ComputeNodalGradient(const SetMatrices<Dimension>& mat,
                         const Vector1& Uloc, Vector2& Unode, int i) const
  {
    reference_element(i)->ComputeNodalGradient(mat, Uloc, Unode);
  }
  

  /*********************
   * SurfacicMesh_Base *
   *********************/
  
  
  //! default constructor
  template<class Dimension, class T>
  inline SurfacicMesh_Base<Dimension, T>::SurfacicMesh_Base() : Mesh<Dimension>() 
  { 
  }
        
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BASE_INLINE_CXX
#endif


