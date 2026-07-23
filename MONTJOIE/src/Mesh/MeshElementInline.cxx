#ifndef MONTJOIE_FILE_MESH_ELEMENT_INLINE_CXX

namespace Montjoie
{
  
  /*************
   * Edge_Base *
   *************/
  
  
  //! default constructor
  inline Edge_Base::Edge_Base()
  {
    ref = 0;
    num_vertex(0) = -1;
    num_vertex(1) = -1;
  }
  
  
  //! initialization of edge 
  inline void Edge_Base::Init(int n1, int n2)
  {
    num_vertex(0) = n1;
    num_vertex(1) = n2;
  }
  
  
  //! initialization of edge
  inline void Edge_Base::Init(int n1, int n2, int ref_)
  {
    num_vertex(0) = n1;
    num_vertex(1) = n2;
    ref = ref_;
  }
    
  
  //! returns reference of the edge
  inline int Edge_Base::GetReference() const
  {
    return ref;
  }
  
  
  //! modifies the reference
  inline void Edge_Base::SetReference(int ref_)
  {
    ref = ref_;
  }
    
  
  //! returns vertex ne of the edge (ne is 0 or 1) 
  inline int Edge_Base::numVertex(int ne) const
  {
    return num_vertex(ne);
  }
  
  
  //! 2 vertices
  inline int Edge_Base::GetNbVertices()
  {
    return 2;
  }
  
  
  //! nothing to be done
  inline void Edge_Base::ClearConnectivity()
  {
  }
  
  
  /********************
   * Edge<Dimension1> *
   ********************/
  
  
  //! default constructor
  inline Edge<Dimension1>::Edge() : Edge_Base()
  {
    pml_elt = false;
  }
  
  
  //! returns memory used in bytes
  inline int Edge<Dimension1>::GetMemorySize() const
  {
    return sizeof(Edge<Dimension1>);
  }
  
  
  //! the element is an element in PML layers
  inline void Edge<Dimension1>::SetPML()
  {
    pml_elt = true;
  }
  

  //! the element is not an element in PML layers
  inline void Edge<Dimension1>::UnsetPML()
  {
    pml_elt = false;
  }
  
  
  //! returns true if the element is PML
  inline bool Edge<Dimension1>::IsPML() const
  {
    return pml_elt;
  }
  
  
  //! Displays informations about 1-D edge
  inline ostream& operator<<(ostream& out, const Edge<Dimension1>& e)
  {
    out << static_cast<const Edge_Base&>(e);
    return out;
  }
  
  
  /********************
   * Edge<Dimension2> *
   ********************/

  
  //! default constructor
  inline Edge<Dimension2>::Edge() : Edge_Base()
  {
    num_element[0] = -1;
    num_element[1] = -1;
    geom_ref = 0;
  }
  

  //! returns memory used in bytes
  inline int Edge<Dimension2>::GetMemorySize() const
  {
    return sizeof(Edge<Dimension2>);
  }

  
  //! returns number of elements adjacent to the edge
  inline int Edge<Dimension2>::GetNbElements() const
  {
    if (num_element[0] == -1)
      return 0;
    
    if (num_element[1] == -1)
      return 1;
    
    return 2;
  }
  
  
  //! returns global element number of the element ne adjacent to the edge (ne = 0 or 1)
  inline int Edge<Dimension2>::numElement(int ne) const
  {
    return num_element[ne];
  }
  

  //! adds the face nf to the list of elements adjacent to the edge
  inline void Edge<Dimension2>::AddElement(int nf)
  {
    if (num_element[0] == -1)
      num_element[0] = nf;
    else
      num_element[1] = nf;
  }
  

  //! sets the list of elements adjacents to (nf0, nf1)
  inline void Edge<Dimension2>::SetElement(int nf0, int nf1)
  {
    num_element[0] = nf0;
    num_element[1] = nf1;
  }
  
  
  //! returns 0
  inline int Edge<Dimension2>::GetHybridType()
  {
    return 0;
  }
  
  
  //! returns reference of the edge
  inline int Edge<Dimension2>::GetGeometryReference() const
  {
    return geom_ref;
  }
  
  
  //! modifies the reference
  inline void Edge<Dimension2>::SetGeometryReference(int ref_)
  {
    geom_ref = ref_;
  }
  
  
  //! clears connectivity
  inline void Edge<Dimension2>::ClearConnectivity()
  {
    num_element[0] = -1;
    num_element[1] = -1;
  }
    
  
  /********************
   * Edge<Dimension3> *
   ********************/
  

  //! default constructor
  inline Edge<Dimension3>::Edge()
  {
    nb_faces = 0;
  }
  

  //! returns memory used in bytes
  inline int Edge<Dimension3>::GetMemorySize() const
  {
    return sizeof(Edge<Dimension3>) + 4*num_elt_face.GetM();
  }

  
  //! returns number of faces adjacent to the edge
  inline int Edge<Dimension3>::GetNbFaces() const
  {
    return nb_faces;
  }
  
  
  //! returns number of volumes adjacent to the edge
  inline int Edge<Dimension3>::GetNbElements() const
  {
    return num_elt_face.GetM()-nb_faces;
  }
  

  //! returns the face number of the face ne adjacent to the edge
  inline int Edge<Dimension3>::numFace(int ne) const
  {
    return num_elt_face(ne);
  }
  

  //! returns the volume number of the volume ne adjacent to the edge
  inline int Edge<Dimension3>::numElement(int ne) const
  {
    return num_elt_face(ne+nb_faces);
  }
  

  //! clears connectivity
  inline void Edge<Dimension3>::ClearConnectivity()
  {
    nb_faces = 0;;
    num_elt_face.Clear();
  }
  

  //! returns true if the edges u and v are the same
  template<class Dimension>
  inline bool operator==(const Edge<Dimension>& u, const Edge<Dimension> & v)
  {
    return (u.num_vertex == v.num_vertex);
  }
  
  
  //! returns true if the edges u and v are different
  template<class Dimension>
  inline bool operator!=(const Edge<Dimension> & u, const Edge<Dimension>& v)
  {
    return (u.num_vertex != v.num_vertex);
  }
  
  
  //! returns true if the vertices of edge u has lower numbers than those of edge v
  template<class Dimension>
  inline bool operator<(const Edge<Dimension> & u, const Edge<Dimension>& v)
  {
    return (u.num_vertex < v.num_vertex);
  }
  
  
  //! returns true if the vertices of edge u has greater numbers than those of edge v
  template<class Dimension>
  inline bool operator>(const Edge<Dimension> & u, const Edge<Dimension>& v)
  {
    return (u.num_vertex > v.num_vertex);
  }
  
  
  /*************
   * Face_Base *
   *************/
  
  
  //! default constructor
  inline Face_Base::Face_Base()
  {
    ref = 0;
    hybrid_type = 0;
    nb_vertices = 0;
  }
  
  
  //! returns true if the face is a triangle 
  inline bool Face_Base::IsTriangular() const
  {
    if (nb_vertices == 3)
      return true;
    
    return false;
  }
  
  
  //! returns true if the face is a quadrilateral
  inline bool Face_Base::IsQuadrangular() const
  {
    if (nb_vertices == 4)
      return true;
    
    return false;
  }
  

  //! returns the type of face
  inline int Face_Base::GetHybridType() const
  {
    return hybrid_type;
  }
  
    
  //! returns the orientation between local edge and global edge
  inline bool Face_Base::GetOrientationEdge(int ne) const
  {
    return way_edge(ne);
  }
    

  //! returns the orientation between local edge and global edge
  inline int Face_Base::GetOrientationBoundary(int ne) const
  {
    return way_edge(ne);
  }
  
  
  //! returns the global vertex number of the local vertex ne of the face 
  inline int Face_Base::numVertex(int ne) const
  {
    return num_vertex(ne);
  }
  
  
  //! returns the global edge number of the local edge ne of the face
  inline int Face_Base::numEdge(int ne) const
  {
    return num_edge(ne);
  }
    

  //! returns reference of the face
  inline int Face_Base::GetReference() const
  {
    return ref;
  }
  

  //! returns the number of vertices of the face (3 or 4)
  inline int Face_Base::GetNbVertices() const
  {
    return nb_vertices;
  }
  

  //! returns the number of found edges of the face
  inline int Face_Base::GetNbEdges() const
  {
    return nb_vertices;
  }
  

  //! modifies reference of the face
  inline void Face_Base::SetReference(int ref_)
  { 
    ref = ref_;
  } 
  
  
  //! sets an edge of the element
  inline void Face_Base::SetEdge(int num_loc, int ne)
  {
    num_edge(num_loc) = ne;
    
  }
  
  
  //! returns local edge number of the global edge ne (-1 is the edge is not found)
  inline int Face_Base::GetPositionBoundary(int ne) const
  {
    for (int i = 0; i < nb_vertices; i++)
      if (num_edge(i) == ne)
	return i;
    
    return -1;
  }
  
  
  //! not used
  inline int Face_Base::GetOrientationFace(int ne) const
  {
    return -1;
  }
  
  
  //! clears the connectivity
  inline void Face_Base::ClearConnectivity()
  {
    num_edge.Fill(-1);
  }
  
  
  //! clears edges
  inline void Face_Base::ClearEdges()
  {
    num_edge.Fill(-1);
  }
  
  
  //! returns true if the faces u and v are different
  template<class DimensionB>
  inline bool operator !=(const Face<DimensionB> & u, const Face<DimensionB>& v)
  {
    return !(u == v);
  }
    
  
  /********************
   * Face<Dimension2> *
   ********************/  
  
  
  //! default constructor
  inline Face<Dimension2>::Face()
  {
    curved_elt = false;
    pml_num = 0;
    pml_type = 0;
  }
  

  //! returns memory used in bytes
  inline int Face<Dimension2>::GetMemorySize() const
  {
    return sizeof(Face<Dimension2>);
  }

  
  //! the element is an element in PML layers
  inline void Face<Dimension2>::SetPML(int num, int type)
  {
    pml_type = type;
    pml_num = num+1;
  }
  

  //! the element is not an element in PML layers
  inline void Face<Dimension2>::UnsetPML()
  {
    pml_num = 0;
    pml_type = 0;
  }
  
  
  //! returns true if the element is PML
  inline bool Face<Dimension2>::IsPML() const
  {
    return (pml_num > 0);
  }


  //! returns number of PML
  inline int Face<Dimension2>::GetTypePML() const
  {
    return pml_type;
  }

  
  //! returns type of PML
  inline int Face<Dimension2>::GetNumberPML() const
  {
    return pml_num-1;
  }

  
  //! the element is set to be curved
  inline void Face<Dimension2>::SetCurved()
  {
    curved_elt = true;
  }

  
  //! the element is set to be curved
  inline void Face<Dimension2>::UnsetCurved()
  {
    curved_elt = false;
  }
  
  
  //! returns true if the element is curved
  inline bool Face<Dimension2>::IsCurved() const
  {
    return curved_elt;
  } 

  
  //! return the number of edges of the face
  inline int Face<Dimension2>::GetNbBoundary() const
  {
    return nb_vertices;
  }

  
  //! returns the global edge number of the local edge ne
  inline int Face<Dimension2>::numBoundary(int ne) const
  {
    return num_edge(ne);
  }
  
  
  /********************
   * Face<Dimension3> *
   ********************/  

  
  //! default constructor
  inline Face<Dimension3>::Face()
  {
    num_element[0] = -1;
    num_element[1] = -1;
    geom_ref = 0;
  }
  
  
  //! returns memory used in bytes
  inline int Face<Dimension3>::GetMemorySize() const
  {
    return sizeof(Face<Dimension3>);
  }

  
  //! returns number of adjacent volumes to the face (1 or 2)
  inline int Face<Dimension3>::GetNbElements() const
  {
    if (num_element[0] == -1)
      return 0;
    
    if (num_element[1] == -1)
      return 1;
    
    return 2;
  }
  
  
  //! returns global volume number of the adjacent volume (ne is 0 or 1)
  inline int Face<Dimension3>::numElement(int ne) const
  {
    return num_element[ne];
  }
  

  //! adds volume ne to the list of volumes
  inline void Face<Dimension3>::AddElement(int ne)
  {
    if (num_element[0] == -1)
      num_element[0] = ne;
    else
      num_element[1] = ne;
  }
  
  
  //! sets list of volumes to (n1, n2)
  inline void Face<Dimension3>::SetElement(int n1,int n2)
  {
    num_element[0] = n1;
    num_element[1] = n2;
  }
    
  
  //! returns geometric reference of the face
  inline int Face<Dimension3>::GetGeometryReference() const
  {
    return geom_ref;
  }
  
  
  //! modifies the geometric reference
  inline void Face<Dimension3>::SetGeometryReference(int ref_)
  {
    geom_ref = ref_;
  }
  
  
  //! clears connectivity
  inline void Face<Dimension3>::ClearConnectivity()
  {
    Face_Base::ClearConnectivity();
    num_element[0] = -1;
    num_element[1] = -1;
  }
  
  
  /**********
   * Volume *
   **********/
  
  
  //! returns memory used in bytes
  inline int Volume::GetMemorySize() const
  {
    return sizeof(Volume);
  }

  
  //! returns reference of the volume
  inline int Volume::GetReference() const
  {
    return ref;
  }
  
  
  //! returns type of the volume
  inline int Volume::GetHybridType() const
  {
    return hybrid_type;
  }
  
  
  //! sets the volume to a PML volume
  inline void Volume::SetPML(int num, int type)
  {
    pml_num = num+1;
    pml_type = type;
  }
  
  
  //! sets the volume to a non-PML volume
  inline void Volume::UnsetPML()
  {
    pml_num = 0;
    pml_type = 0;
  }
  
  
  //! returns true if the volume is PML
  inline bool Volume::IsPML() const
  {
    return (pml_num > 0);
  }


  //! returns type of PML
  inline int Volume::GetTypePML() const
  {
    return pml_type;
  }


  //! returns number of PML
  inline int Volume::GetNumberPML() const
  {
    return pml_num-1;
  }

  
  //! sets the volume to a curved volume
  inline void Volume::SetCurved()
  {
    curved_elt = true;
  }
  
  
  //! sets the volume to a curved volume
  inline void Volume::UnsetCurved()
  {
    curved_elt = false;
  }
  
  
  //! returns true if the volume is curved
  inline bool Volume::IsCurved() const
  {
    return curved_elt;
  }
  
  
  //! returns the number of vertices
  inline int Volume::GetNbVertices() const
  {
    return nb_vertices;
  }
  

  //! returns the number of edges 
  inline int Volume::GetNbEdges() const
  {
    return nb_edges;
  }
  
  
  //! returns the number of faces
  inline int Volume::GetNbFaces() const
  {
    return nb_faces;
  }
  
  
  //! returns the number of faces
  inline int Volume::GetNbBoundary() const
  {
    return nb_faces;
  }
  
  
  //! returns global vertex number of the local vertex ne
  inline int Volume::numVertex(int ne) const
  {
    return num_vertex(ne);
  }
  
  
  //! returns global edge number of the local edge ne
  inline int Volume::numEdge(int ne) const
  {
    return num_edge(ne);
  }
  
  
  //! returns global face number of the local face nf
  inline int Volume::numFace(int nf) const
  {
    return num_face(nf);
  }
   

  //! returns global face number of the local face nf
  inline int Volume::numBoundary(int nf) const
  {
    return num_face(nf);
  }
  
  
  //! returns the orientation between local edge and global edge
  inline bool Volume::GetOrientationEdge(int ne) const
  {
    return way_edge(ne);
  }

  
  //! sets the orientation between local face and global face
  inline void Volume::SetOrientationFace(int nf, int way)
  {
    way_face(nf) = way;
  }
  

  //! returns the orientation between local face and global face
  inline int Volume::GetOrientationFace(int nf) const
  {
    return way_face(nf);
  }
    
  
  //! returns the orientation between local face and global face
  inline int Volume::GetOrientationBoundary(int nf) const
  {
    return way_face(nf);
  }
  
    
  //! modifies the reference
  inline void Volume::SetReference(int ref_)
  {
    ref = ref_;
  }
    
  
  //! tetrahedral volume
  inline void Volume::InitTetrahedral(const IVect& num, int ref_)
  {
    InitTetrahedral(num(0), num(1), num(2), num(3), ref_);
  }
  
  
  //! pyramidal volume
  inline void Volume::InitPyramidal(const IVect& num, int ref_)
  {
    InitPyramidal(num(0), num(1), num(2), num(3), num(4), ref_);
  }
  
    
  //! wedge
  inline void Volume::InitWedge(const IVect& num, int ref_)
  {
    InitWedge(num(0), num(1), num(2), num(3), num(4), num(5), ref_);
  }
  
  
  //! hexahedral volume
  inline void Volume::InitHexahedral(const IVect& num, int ref_)
  {
    InitHexahedral(num(0), num(1), num(2), num(3), num(4),
                   num(5), num(6), num(7), ref_);
  }
  
  
  //! sets an edge
  inline void Volume::SetEdge(int num_loc, int ne)
  {
    num_edge(num_loc) = ne;
  }
  
  
  //! sets a face
  inline void Volume::SetFace(int num_loc, int nf)
  {
    num_face(num_loc) = nf;
  }
    
  
  //! sets a face
  inline void Volume::SetFace(int num_loc, int nf, int way)
  {
    num_face(num_loc) = nf;
    way_face(num_loc) = way;
  }
  
  
  //! returns local number of the global edge ne
  inline int Volume::GetPositionEdge(int ne) const
  {
    for (int i = 0; i < nb_edges; i++)
      if (num_edge(i) == ne)
	return i;
    
    return -1;
  }
  
  
  //! returns local number of the global face nf
  inline int Volume::GetPositionBoundary(int nf) const
  {
    for (int i = 0; i < nb_faces; i++)
      if (num_face(i) == nf)
	return i;
    
    return -1;
  }
  
  
  //! clears the connectivity
  inline void Volume::ClearConnectivity()
  {
    num_edge.Fill(-1);
    num_face.Fill(-1);
    way_face.Fill(true);
  }
  
  
  //! clears edges
  inline void Volume::ClearEdges()
  {
    num_edge.Fill(-1);
  }
  
  
  //! clears faces
  inline void Volume::ClearFaces()
  {
    num_face.Fill(-1);
    way_face.Fill(true);
  }
  
  
  /********************
   * ElementNumbering *
   ********************/
  
  
  //! returns memory used by the object in bytes
  inline int ElementNumbering::GetMemorySize() const
  {
    int taille = sizeof(ElementNumbering) + sizeof(int)*(Nodle.GetM() + negative_dof_num.GetM());
    return taille;
  }
  
    
  //! returns number of dofs of the element
  inline int ElementNumbering::GetNbDof() const
  {
    return Nodle.GetM();
  }
  
  
  //! returns global dof number
  inline int ElementNumbering::GetNumberDof(int i) const
  {
    return Nodle(i);
  }
  
  
  //! returns the global dof numbers
  inline const IVect& ElementNumbering::GetNodle() const
  {
    return Nodle;
  }
  
  
  //! returns local dof numbers for which the sign must be changed
  inline const IVect& ElementNumbering::GetNegativeDofNumber() const
  {
    return negative_dof_num;
  }
  
  
  //! sets all the dofs for which the sign must be changed
  inline void ElementNumbering::SetNegativeDofNumber(int nb, const IVect& num)
  {
    negative_dof_num.Reallocate(nb);
    for (int i = 0; i < nb; i++)
      negative_dof_num(i) = num(i);
  }
  
  
  //! informs that the element contains nf faces
  inline void ElementNumbering::ReallocateFaces(int nf)
  {
    order_face.Fill(-1);
    offset_face.Fill(0);
  }
  
  
  //! returns the order associated with a face 
  inline int ElementNumbering::GetOrderFace(int num_loc) const
  {
    return order_face(num_loc);
  }

  
  //! returns the offset associated with a face
  inline int ElementNumbering::GetOffsetFace(int num_loc) const
  {
    return offset_face(num_loc);
  }
  
  
  //! sets the order and the offset for a face
  inline void ElementNumbering::SetOrderFace(int num_loc, int r, int offset)
  {
    order_face(num_loc) = r;
    offset_face(num_loc) = offset;
  }
  
  
  //! allocates dof array 
  inline void ElementNumbering::ReallocateDof(int nb_dof_elt)
  {
    Nodle.Reallocate(nb_dof_elt);
  }
  
  
  //! sets the global number dof of a local dof
  inline void ElementNumbering::SetNumberDof(int num_loc, int num_glob)
  {
    Nodle(num_loc) = num_glob;
  }
  
  
  //! adds several dofs
  inline void ElementNumbering::AddDof(int nb_dof)
  { 
    Nodle.Resize(Nodle.GetM() + nb_dof);
  }
  
  
  //! adds a single dof
  inline void ElementNumbering::PushDof(int num_dof)
  { 
    Nodle.PushBack(num_dof);
  }
  
  
  //! erases dofs of this element
  inline void ElementNumbering::Clear()
  {
    Nodle.Clear();
    negative_dof_num.Clear();
    //phase_dof_num.Clear();
    //type_phase_dof.Clear();
  }
  
} // end namespace

#define MONTJOIE_FILE_MESH_ELEMENT_INLINE_CXX
#endif
