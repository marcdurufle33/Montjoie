#ifndef MONTJOIE_FILE_NUMBER_MESH_INLINE_CXX

namespace Montjoie
{

  /**********************
   * MeshNumbering_Base *
   **********************/
  
  
  //! returns a reference to the numbering of element i
  template<class T>
  inline ElementNumbering& MeshNumbering_Base<T>::Element(int i)
  {
    return dof_element(i);
  }
  
  
  //! returns a reference to the numbering of element i
  template<class T>
  inline const ElementNumbering& MeshNumbering_Base<T>::Element(int i) const
  {
    return dof_element(i);
  }
  
  
  //! returns the number of neighboring elements
  template<class T>  
  inline int MeshNumbering_Base<T>::GetNbNeighborElt() const 
  {
    return order_neighbor_elt.GetM();
  }
  
  
  //! returns the order of approximation for neighboring element i
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrderNeighborElement(int i) const
  {
    return order_neighbor_elt(i);
  }
  
  
  //! returns the element type of neighboring element i
  template<class T>
  inline int MeshNumbering_Base<T>::GetTypeNeighborElement(int i) const
  {
    return type_neighbor_elt(i);
  }
  
  
  //! returns the edge number associated with the neighboring element i
  template<class T>
  inline int MeshNumbering_Base<T>::GetEdgeNeighborElement(int i) const
  {
    return num_edge_neighbor_elt(i);
  }
  

  template<class T>
  inline int MeshNumbering_Base<T>::GetLocalEdgeNumberNeighborElement(int i) const
  {
    return inv_num_edge_neighbor_elt(i);
  }
  

  //! sets the order of approximation for neighboring element i  
  template<class T>
  inline void MeshNumbering_Base<T>::SetOrderNeighborElement(int i, int r)
  {
    order_neighbor_elt(i) = r;
  }
  

  //! sets the element type of neighboring element i  
  template<class T>
  inline void MeshNumbering_Base<T>::SetTypeNeighborElement(int i, int type)
  {
    type_neighbor_elt(i) = type;
  }
  

  //! sets the edge number associated with the neighboring element i  
  template<class T>
  inline void MeshNumbering_Base<T>::SetEdgeNeighborElement(int i, int ne)
  {
    num_edge_neighbor_elt(i) = ne;
  }
  
  
  //! returns the type of periodicity for dof i
  template<class T>
  inline int MeshNumbering_Base<T>::GetPeriodicityTypeForDof(int i) const
  {
    return type_periodicity_dof(i);
  }


  //! sets the type of periodicity for dof i
  template<class T>
  inline void MeshNumbering_Base<T>::SetPeriodicityTypeForDof(int i, int n)
  {
    type_periodicity_dof(i) = n;
  }
  
  
  //! returns the type of periodicity for the edge/face i
  template<class T>
  inline int MeshNumbering_Base<T>::GetPeriodicityTypeForBoundary(int i) const
  {
    return type_periodicity_boundary(i);
  }


  //! sets the type of periodicity for the edge/face i
  template<class T>
  inline void MeshNumbering_Base<T>::SetPeriodicityTypeForBoundary(int i, int n)
  {
    type_periodicity_boundary(i) = n;
  }

  
  //! returns the number of periodic dofs
  template <class T>
  inline int MeshNumbering_Base<T>::GetNbPeriodicDof() const
  {
    return PeriodicDofs.GetM();
  }
  
  
  //! returns the global dof number of i-th periodic dof
  template <class T>
  inline int MeshNumbering_Base<T>::GetPeriodicDof(int i) const
  {
    return PeriodicDofs(i);
  }

  
  //! sets the global dof number of i-th periodic dof
  template <class T>
  inline void MeshNumbering_Base<T>::SetPeriodicDof(int i, int n)
  {
    PeriodicDofs(i) = n;
  }
  
  
  //! returns the global dof number of the original dof of i-th periodic dof
  template <class T>
  inline int MeshNumbering_Base<T>::GetOriginalPeriodicDof(int i) const
  {
    return PeriodicDofs_Original(i);
  }

  
  //! sets the global dof number of the original dof of i-th periodic dof
  template <class T>
  inline void MeshNumbering_Base<T>::SetOriginalPeriodicDof(int i, int n)
  {
    PeriodicDofs_Original(i) = n;
  }


  //! returns the number of the opposite boundary
  template <class T>
  inline int MeshNumbering_Base<T>::GetPeriodicBoundary(int i) const
  {
    return PeriodicityBoundary(i);
  }

  
  template <class T>
  inline void MeshNumbering_Base<T>::SetPeriodicBoundary(int i, int nf)
  {
    PeriodicityBoundary(i) = nf;
  }


  template<class T>
  inline void MeshNumbering_Base<T>::SetSameNumberPeriodicDofs()
  {
    type_formulation_periodic = SAME_PERIODIC_DOFS;
  }


  template<class T>
  inline void MeshNumbering_Base<T>::SetFormulationForPeriodicCondition(int type)
  {
    type_formulation_periodic = type;
  }


  template<class T>
  inline int MeshNumbering_Base<T>::GetFormulationForPeriodicCondition() const
  {
    return type_formulation_periodic;
  }

  
  template<class T>
  inline bool MeshNumbering_Base<T>::UseSameNumberForPeriodicDofs() const
  {
    return (type_formulation_periodic == SAME_PERIODIC_DOFS);
  }

  
  //! returns order of approximation when uniform in the mesh
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrder() const
  {
    return order;
  }
  
  
  //! returns the number of degrees of freedom in the mesh
  template<class T>
  inline int MeshNumbering_Base<T>::GetNbDof() const
  {
    return nodl;
  }

  
  //! sets the number of degrees of freedom in the mesh
  template<class T>
  inline void MeshNumbering_Base<T>::SetNbDof(int nb_dof)
  {
    nodl = nb_dof;
  }
  
  
  //! returns the number of dofs on the i-th element
  template<class T>
  inline int MeshNumbering_Base<T>::GetNbLocalDof(int i) const
  {
    return nb_dof_element(i);
  }
  
  
  //! returns the number of dofs contained in PML layers
  template<class T>
  inline int MeshNumbering_Base<T>::GetNbDofPML() const
  {
    return nb_dof_pml;
  }
  
  
  //! returns the global dof number of a dof in PML layers
  template<class T>
  inline int MeshNumbering_Base<T>::GetDofPML(int i) const
  {
    return IndexDofPML(i);
  }
  
  
  //! allocates arrays to store n dofs in PML layers
  template<class T>
  inline void MeshNumbering_Base<T>::ReallocateDofPML(int n)
  {
    IndexDofPML.Reallocate(nodl);
    IndexDofPML.Fill(-1);
    nb_dof_pml = n;
  }
  
  
  //! sets the global dof number of a dof in PML layers
  template<class T>
  inline void MeshNumbering_Base<T>::SetDofPML(int i, int j)
  {
    IndexDofPML(i) = j;
  }

  
  //! sets an uniform order on all the elements
  template<class T>
  inline void MeshNumbering_Base<T>::SetOrder(int r)
  {
    order = r;
    type_variable_order = CONSTANT_ORDER;
  }
  
  
  //! returns true if the order is variable
  template<class T>
  inline bool MeshNumbering_Base<T>::IsOrderVariable() const
  {
    return (type_variable_order != CONSTANT_ORDER);
  }
  
  
  //! returns the algorithm used to determine the order for each edge, face, element
  template<class T>
  inline int MeshNumbering_Base<T>::GetVariableOrder() const
  {
    return type_variable_order;
  }
  
  
  //! sets the order to use in the mesh (uniform, variable)
  template<class T>
  inline void MeshNumbering_Base<T>::SetVariableOrder(int type)
  {
    type_variable_order = type;
  }
  
  
  //! sets a variable order defined by the edge length
  /*!
    For each edge, if the length is comprised between step(i) and step(i+1),
    the order is equal to i
    The order of faces and elements is taken as the maximal order of edges belonging to the entity
   */
  template<class T>
  inline void MeshNumbering_Base<T>::SetMeshSizeVariableOrder(const Vector<T>& step)
  {
    mesh_size_variable_order = step;
  }
  
  
  //! you can specify a coefficient to compute the edge length in order
  //! to take into account physical hetereogeneities
  template<class T>
  inline void MeshNumbering_Base<T>::SetCoefficientVariableOrder(const Vector<T>& coef)
  {
    coef_elt_variable_order = coef;
  }

  
  //! returns the order of element num_elt
  /*!
    \param[in] num_elt element number
   */
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrderElement(int num_elt) const
  {
    if (type_variable_order != CONSTANT_ORDER)
      return order_element(num_elt);
   
    return order;
  }
  
  
  //! returns the order of edge num_edge
  /*!
    \param[in] num_edge edge number
   */
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrderEdge(int num_edge) const
  {
    if (type_variable_order != CONSTANT_ORDER)
      return order_edge(num_edge);
   
    return order;
  }

  
  //! returns the order of face num_face
  /*!
    \param[in] num_face face number
   */
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrderFace(int num_face) const
  {
    if (type_variable_order != CONSTANT_ORDER)
      return order_face(num_face);
   
    return order;
  }

  
  //! returns the order associated with the interior of the element num_elt
  template<class T>
  inline int MeshNumbering_Base<T>::GetOrderInside(int num_elt) const
  {
    if (type_variable_order != CONSTANT_ORDER)
      return order_inside(num_elt);
   
    return order;
  }
  
  
  //! sets the order associated with the interior of the edge num_edge
  template<class T>
  inline void MeshNumbering_Base<T>::SetOrderEdge(int num_edge, int r)
  {
    order_edge(num_edge) = r;
  }
  
  
  //! sets the order associated with the interior of the face num_face
  template<class T>
  inline void MeshNumbering_Base<T>::SetOrderFace(int num_face, int r)
  {
    order_face(num_face) = r;
  }
  
  
  //! sets the order associated with the interior of the element num_elt
  template<class T>
  inline void MeshNumbering_Base<T>::SetOrderInside(int num_elt, int r)
  {
    order_inside(num_elt) = r;
  }
  
  
  //! returns the number of quadrature points on the boundary nf
  template<class T>
  inline int MeshNumbering_Base<T>::GetNbPointsQuadratureBoundary(int nf) const
  {
    return nb_points_quadrature_on_boundary(nf);
  }

  
  /*********************
   * MeshNumbering_Dim *
   *********************/


  template <class Dimension, class T>
  inline int MeshNumbering_Dim<Dimension, T>::GetOrder() const
  {
    return this->order;
  }


  template <class Dimension, class T>
  inline const Mesh<Dimension, T>& MeshNumbering_Dim<Dimension, T>::GetMesh() const
  {
    return this->mesh;
  }

  
  template <class Dimension, class T>
  inline void MeshNumbering_Dim<Dimension, T>::InitPeriodicBoundary()
  {
    this->PeriodicityBoundary.Reallocate(mesh.GetNbBoundaryRef());
    this->type_periodicity_boundary.Reallocate(mesh.GetNbBoundaryRef());
    this->PeriodicityBoundary.Fill(-1);
    this->type_periodicity_boundary.Fill(-1);
  }
    

  //! returns the translation vector between the original and periodic dof
  template <class Dimension, class T>
  inline const TinyVector<T, Dimension::dim_N>& MeshNumbering_Dim<Dimension, T>::
  GetTranslationPeriodicDof(int i) const
  {
    return TranslationPeriodicDof(i);
  }


  //! sets the translation vector between the original and periodic dof
  template <class Dimension, class T>
  inline void MeshNumbering_Dim<Dimension, T>::
  SetTranslationPeriodicDof(int i, const TinyVector<T, Dimension::dim_N>& v)
  {
    TranslationPeriodicDof(i) = v;
  }
  
  
  //! returns the translation vector between the original and periodic boundary
  template <class Dimension, class T>
  inline const TinyVector<T, Dimension::dim_N>& MeshNumbering_Dim<Dimension, T>::
  GetTranslationPeriodicBoundary(int i) const
  {
    return TranslationPeriodicBoundary(i);
  }


  //! returns the translation vector between the original and periodic boundary
  template <class Dimension, class T>
  inline Vector<TinyVector<T, Dimension::dim_N> >& MeshNumbering_Dim<Dimension, T>::
  GetTranslationPeriodicBoundary()
  {
    return TranslationPeriodicBoundary;
  }

  
  //! returns leaf class
  template<class Dimension, class T>
  inline MeshNumbering<Dimension, T>& MeshNumbering_Dim<Dimension, T>::GetLeafClass()
  {
    return static_cast<MeshNumbering<Dimension, T>& >(*this);
  }


  //! returns leaf class
  template<class Dimension, class T>
  inline const MeshNumbering<Dimension, T>& MeshNumbering_Dim<Dimension, T>::GetLeafClass() const
  {
    return static_cast<const MeshNumbering<Dimension, T>& >(*this);
  }
  
    
  /**************
   * Dimension2 *
   **************/
  

  template<class T>
  inline int MeshNumbering<Dimension2, T>::GetLocalVertexBoundary(int type, int num_loc, int i)
  {
    if (i == 0)
      return FirstExtremityEdge(type, num_loc);
    else
      return SecondExtremityEdge(type, num_loc);
  }

  
  //! returns the order of quadrature rule to use on edge ne
  template<class T>
  inline int MeshNumbering<Dimension2, T>::GetOrderQuadrature(int ne) const
  {
    if (this->type_variable_order != this->CONSTANT_ORDER)
      return this->order_boundary(ne);
    
    return this->order;
  }
  
  
  //! sets the order of quadrature rule for edge/face ne
  template<class T>
  inline void MeshNumbering<Dimension2, T>::SetOrderQuadrature(int ne, int r)
  {
    this->order_boundary(ne) = r;
  }

  
  //! returns quadrature points for unit interval
  template<class T>
  inline const VectReal_wp& MeshNumbering<Dimension2, T>::GetReferenceQuadrature(int r, int type_elt) const
  {
    return this->number_map.GetEdgeQuadrature(r);
  }
  
  
  //! sets the order associated with the element num_elem
  template<class T>
  inline void MeshNumbering<Dimension2, T>::SetOrderElement(int num_elem, int r)
  {
    if ((this->order != r) && (this->order_element.GetM() == 0))
      this->SetVariableOrder(this->USER_ORDER);
    
    if (this->type_variable_order != this->CONSTANT_ORDER)
      this->order_element(num_elem) = r;    
  }
  
  
  //! return the difference of orientation between two faces
  template<class T>
  inline int MeshNumbering<Dimension2, T>::GetRotationFace(int rot1, int rot2, int nv)
  {
    if (rot1 == rot2)
      return 0;
    
    return 1;
  }
  
  
  //! returns the opposite rotation
  template<class T>
  inline int MeshNumbering<Dimension2, T>::
  GetOppositeOrientationFace(int rot, const Edge<Dimension2>&)
  {
    return rot;
  }


#ifdef MONTJOIE_WITH_THREE_DIM  
  //! returns the order of quadrature rule associated with the edge/face ne
  template<class T>
  inline int MeshNumbering<Dimension3, T>::GetOrderQuadrature(int ne) const
  {
    if (this->type_variable_order != this->CONSTANT_ORDER)
      return this->order_boundary(ne);
    
    return this->order;
  }
  
  
  //! sets the order of quadrature rule associated with the edge/face ne
  template<class T>
  inline void MeshNumbering<Dimension3, T>::SetOrderQuadrature(int ne, int r)
  {
    this->order_boundary(ne) = r;
  }
  
  
  //! returns quadrature points for unit interval
  template<class T>
  inline const VectR2& MeshNumbering<Dimension3, T>::GetReferenceQuadrature(int r, int type_elt) const
  {
    if (type_elt == 0)
      return this->number_map.GetTriangleQuadrature(r);
    
    return this->number_map.GetQuadrangleQuadrature(r);
  }


  //! sets the order associated with the element num_elem
  template<class T>
  inline void MeshNumbering<Dimension3, T>::SetOrderElement(int num_elem, int r)
  {
    if ((this->order != r) && (this->order_element.GetM() == 0))
      this->SetVariableOrder(this->USER_ORDER);
    
    if (this->type_variable_order != this->CONSTANT_ORDER)
      this->order_element(num_elem) = r;
  }


  //! returns rotation of face 2 compared to face 1
  /*!
    \param[in] way_face1 way of the first face
    \param[in] way_face2 way of the second face
    \return rotation between the two faces
   */
  template<class T>
  inline int MeshNumbering<Dimension3, T>::GetRotationFace(int way_face1,
							  int way_face2, int nv)
  {
    if (nv == 3)
      return RotationFaceTri(way_face1, way_face2);
    
    return RotationFaceQuad(way_face1, way_face2);
  }  

  
  //! returns opposite orientation of a face
  template<class T>
  inline int MeshNumbering<Dimension3, T>::
  GetOppositeOrientationFace(int rot, bool triangular_face)
  {  
    if (triangular_face)
      {
        if ( rot == 2)
          return 1;
        else if ( rot == 1)
          return 2;
      }
    else
      {
        if ( rot == 3)
          return 1;
        else if ( rot == 1)
          return 3;
      }
    
    // other orientations are not modified
    return rot;
  }

  
  //! returns opposite orientation of a face
  template<class T>
  inline int MeshNumbering<Dimension3, T>::
  GetOppositeOrientationFace(int rot, const Face<Dimension3>& f)
  {
    return GetOppositeOrientationFace(rot, f.IsTriangular());
  }
  
  
  //! returns the vertex number of local vertex num of local face nf
  /*!
    type is the hybrid type of an element (0: tet, 1: pyramid, 2: prism, 3: hex)
   */
  template<class T>
  inline int MeshNumbering<Dimension3, T>::
  GetVertexNumberOfFace(int type, int nf, int num)
  {
    switch (num)
      {
      case 0 : return FirstExtremityFace(type, nf);
      case 1 : return SecondExtremityFace(type, nf);
      case 2 : return ThirdExtremityFace(type, nf);
      case 3 : return FourthExtremityFace(type, nf);
      }
    
    return -1;
  }
  
  
  template<class T>
  inline int MeshNumbering<Dimension3, T>::GetLocalVertexBoundary(int type, int nf, int num)
  {
    return GetVertexNumberOfFace(type, nf, num);
  }

  
  //! returns the edge number of local edge num in local face nf
  /*!
    type is the hybrid type of an element (0: tet, 1: pyramid, 2: prism, 3: hex)
   */
  template<class T>
  inline int MeshNumbering<Dimension3, T>::
  GetEdgeNumberOfFace(int type, int nf, int num)
  {
    switch (num)
      {
      case 0 : return FirstEdgeFace(type, nf);
      case 1 : return SecondEdgeFace(type, nf);
      case 2 : return ThirdEdgeFace(type, nf);
      case 3 : return FourthEdgeFace(type, nf);
      }
    
    return -1;
  }
  
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_NUMBER_MESH_INLINE_CXX
#endif
