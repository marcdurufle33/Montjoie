#ifndef MONTJOIE_FILE_MESH_ELEMENT_CXX

namespace Montjoie
{

  //! Displays informations about an edge
  ostream& operator<<(ostream& out, const Edge_Base& e)
  {
    out<<"Edge of extremities "<<e.num_vertex(0)+1<<"  "<<e.num_vertex(1)+1
       <<"   , Reference " << e.ref << endl;
    
    return out;
  }


  //! Displays informations about 2-D edge
  ostream& operator<<(ostream& out, const Edge<Dimension2>& e)
  {
    out << static_cast<const Edge_Base&>(e);
    out<<"Edge adjacent to elements " << e.num_element[0]+1 << "  " << e.num_element[1]+1;
    out << endl;
    return out;
  }
  

  //! adds the face nf to the list of faces adjacent to the edge
  void Edge<Dimension3>::AddFace(int nf)
  {
    for (int i = 0; i < nb_faces; i++)
      if (num_elt_face(i) == nf)
	return;
    
    int n = num_elt_face.GetM();
    num_elt_face.Resize(n+1);
    for (int k = n; k > nb_faces; k--)
      num_elt_face(k) = num_elt_face(k-1);
    
    num_elt_face(nb_faces) = nf;
    nb_faces++;
  }


  //! adds the element nf to the list of elements adjacent to the edge
  void Edge<Dimension3>::AddElement(int nf)
  {
    for (int i = nb_faces; i < num_elt_face.GetM(); i++)
      if (num_elt_face(i) == nf)
        return;
    
    num_elt_face.PushBack(nf);
  }
  

  //! clears list of faces of edge e
  void Edge<Dimension3>::ClearFaces()
  {
    if (nb_faces == num_elt_face.GetM())
      {
        num_elt_face.Clear();
        nb_faces = 0;
        return;
      }
    
    IVect num(num_elt_face);
    num_elt_face.Reallocate(num.GetM()-nb_faces);
    for (int i = nb_faces; i < num.GetM(); i++)
      num_elt_face(i-nb_faces) = num(i);
    
    nb_faces = 0;
  }
  
  
  //! Displays informations about 3-D edge
  ostream& operator<<(ostream& out, const Edge<Dimension3>& e)
  {
    out << static_cast<const Edge_Base&>(e);
    out << "Faces adjacent to the edge : ";
    for (int j = 0; j < e.GetNbFaces(); j++)
      out << e.numFace(j)+1 << " ";
    
    out << endl;
    out << "Elements adjacent to the edge : "; 
    for (int j = 0; j < e.GetNbElements(); j++)
      out << e.numElement(j)+1 << " ";
    
    out << endl;
    return out;
  }
    
  
  /*************
   * Face_Base *
   *************/
  
  
  //!< modifies the local vertex num_loc
  void Face_Base::SetVertex(int num_loc, int num_glob)
  {
    num_vertex(num_loc) = num_glob;
    
    // modification of the orientation of edges
    int n = num_loc-1;
    if (n < 0)
      n += nb_vertices;
    
    if (num_glob > num_vertex(n))
      way_edge(n) = true;
    else
      way_edge(n) = false;
    
    n = num_loc+1;
    if (n >= nb_vertices)
      n = 0;
    
    if (num_glob > num_vertex(n))
      way_edge(num_loc) = false;
    else
      way_edge(num_loc) = true;
  }
  

  //! initialization of a triangular face
  void Face_Base::InitTriangular(int n1, int n2, int n3, int ref_)
  {
    nb_vertices = 3;
    num_vertex(0) = n1;
    num_vertex(1) = n2;
    num_vertex(2) = n3;
    ref = ref_;
    hybrid_type = 0;
    num_edge.Fill(-1);
    if (n2 > n1)
      way_edge(0) = true;
    else
      way_edge(0) = false;
    
    if (n3 > n2)
      way_edge(1) = true;
    else
      way_edge(1) = false;
    
    if (n1 > n3)
      way_edge(2) = true;
    else
      way_edge(2) = false;
  }
  
  
  //! initialization of a quadrilateral face
  void Face_Base::InitQuadrangular(int n1, int n2, int n3,
				   int n4, int ref_)
  {
    nb_vertices = 4;
    num_vertex(0) = n1;
    num_vertex(1) = n2;
    num_vertex(2) = n3;
    num_vertex(3) = n4;
    ref = ref_;
    hybrid_type = 1;
    num_edge.Fill(-1);
    if (n2 > n1)
      way_edge(0) = true;
    else
      way_edge(0) = false;
    
    if (n3 > n2)
      way_edge(1) = true;
    else
      way_edge(1) = false;
    
    if (n4 > n3)
      way_edge(2) = true;
    else
      way_edge(2) = false;
    
    if (n1 > n4)
      way_edge(3) = true;
    else
      way_edge(3) = false;
  }
  
  
  //! initialization of a face (triangular if num is an array of dimension 3)
  void Face_Base::Init(const IVect& num, int ref_)
  {
    if (num.GetM() == 3)
      InitTriangular(num(0), num(1), num(2), ref_);
    else
      InitQuadrangular(num(0), num(1), num(2), num(3), ref_);
  }
  
  
  //! Displays informations about a face
  ostream& operator<<(ostream& out, const Face_Base& e)
  {
    out << "Face with vertices : ";
    for (int j = 0; j < e.GetNbVertices(); j++)
      out << e.numVertex(j) + 1 << " ";
    
    out << endl << "Reference : "<<e.ref<<endl;
    out << "Edge numbers :  ";
    for (int j = 0; j < e.GetNbEdges(); j++)
      out << e.numEdge(j)+1 << " ";
    
    out << endl << "Difference of orientations for edges :  ";
    for (int j = 0; j < e.GetNbEdges(); j++)
      out << e.GetOrientationEdge(j) << " ";
    
    out << endl;
    return out;
  }
  
  
  //! returns true if the faces u and v are the same
  template<class DimensionB>
  bool operator ==(const Face<DimensionB>& u, const Face<DimensionB> & v)
  {
    int nu = u.GetNbVertices();
    int nv = v.GetNbVertices();
    if (nu == 3)
      {
        if (nv == 4)
          return false;
        else
          {
            if (u.numVertex(0) == v.numVertex(0))
              {
                if (u.numVertex(1) == v.numVertex(1))
                  {
                    if (u.numVertex(2) == v.numVertex(2))
                      return true;
                    else 
                      return false;
                  }
                else if (u.numVertex(1) == v.numVertex(2))
                  {
                    if (u.numVertex(2) == v.numVertex(1))
                      return true;
                    else 
                      return false;
                  }
                else
                  return false;
              }
            else
              return false;
          }
      }
    else
      {
        if (nv == 3)
          return false;
        else
          {
            if ((u.numVertex(0) == v.numVertex(0))&&(u.numVertex(2) == v.numVertex(2)))
              {
                if (u.numVertex(1) == v.numVertex(1))
                  {
                    if (u.numVertex(3) == v.numVertex(3))
                      return true;
                    else 
                      return false;
                  }
                else if (u.numVertex(1) == v.numVertex(3))
                  {
                    if (u.numVertex(3) == v.numVertex(1))
                      return true;
                    else 
                      return false;
                  }
                else
                  return false;
              }
            else
              return false;
          }
      }
    
    return false;
  }
   
  
  //! returns true if the vertex numbers of the face u are lower than
  //! those of the face v
  template<class DimensionB>
  bool operator <(const Face<DimensionB> & u, const Face<DimensionB>& v)
  {
    int nu = u.GetNbVertices();
    int nv = v.GetNbVertices();
    if (nu == 3)
      {
        if (nv == 4)
          return true;
        else
          {
            if (u.numVertex(0) < v.numVertex(0))
              return true;
            else if (u.numVertex(0) > v.numVertex(0))
              return false;
            else
              {
                int u1, u2, v1, v2;
                if (u.numVertex(1) < u.numVertex(2))
                  {
                    u1 = u.numVertex(1);
                    u2 = u.numVertex(2);
                  }
                else
                  {
                    u1 = u.numVertex(2);
                    u2 = u.numVertex(1);
                  }
                
                if (v.numVertex(1) < v.numVertex(2))
                  {
                    v1 = v.numVertex(1);
                    v2 = v.numVertex(2);
                  }
                else
                  {
                    v1 = v.numVertex(2);
                    v2 = v.numVertex(1);
                  }
                
                if (u1 < v1)
                  return true;
                else if (u1 > v1)
                  return false;
                else
                  {
                    if (u2 < v2)
                      return true;
                    else
                      return false;
                  }
              }
          }
      }
    else
      {
        if (nv == 3)
          return false;
        else
          {
            if (u.numVertex(0) < v.numVertex(0))
              return true;
            else if (u.numVertex(0) > v.numVertex(0))
              return false;
            else
              {
                if (u.numVertex(2) < v.numVertex(2))
                  return true;
                else if (u.numVertex(2) > v.numVertex(2))
                  return false;
                else
                  {
                    int u1, u2, v1, v2;
                    if (u.numVertex(1) < u.numVertex(3))
                      {
                        u1 = u.numVertex(1);
                        u2 = u.numVertex(3);
                      }
                    else
                      {
                        u1 = u.numVertex(3);
                        u2 = u.numVertex(1);
                      }
                    
                    if (v.numVertex(1) < v.numVertex(3))
                      {
                        v1 = v.numVertex(1);
                        v2 = v.numVertex(3);
                      }
                    else
                      {
                        v1 = v.numVertex(3);
                        v2 = v.numVertex(1);
                      }
                    
                    if (u1 < v1)
                      return true;
                    else if (u1 > v1)
                      return false;
                    else
                      {
                        if (u2 < v2)
                          return true;
                        else
                          return false;
                      }
                  }
              }
          }
      }
    
    return false;
  }
  

  //! returns true if the vertex numbers of the face u are greater than
  //! those of the face v  
  template<class DimensionB>
  bool operator >(const Face<DimensionB> & u, const Face<DimensionB>& v)
  {
    int nu = u.GetNbVertices();
    int nv = v.GetNbVertices();
    if (nu == 3)
      {
        if (nv == 4)
          return false;
        else
          {
            if (u.numVertex(0) < v.numVertex(0))
              return false;
            else if (u.numVertex(0) > v.numVertex(0))
              return true;
            else
              {
                int u1, u2, v1, v2;
                if (u.numVertex(1) < u.numVertex(2))
                  {
                    u1 = u.numVertex(1);
                    u2 = u.numVertex(2);
                  }
                else
                  {
                    u1 = u.numVertex(2);
                    u2 = u.numVertex(1);
                  }
                
                if (v.numVertex(1) < v.numVertex(2))
                  {
                    v1 = v.numVertex(1);
                    v2 = v.numVertex(2);
                  }
                else
                  {
                    v1 = v.numVertex(2);
                    v2 = v.numVertex(1);
                  }
                
                if (u1 < v1)
                  return false;
                else if (u1 > v1)
                  return true;
                else
                  {
                    if (u2 > v2)
                      return true;
                    else
                      return false;
                  }
              }
          }
      }
    else
      {
        if (nv == 3)
          return true;
        else
          {
            if (u.numVertex(0) < v.numVertex(0))
              return false;
            else if (u.numVertex(0) > v.numVertex(0))
              return true;
            else
              {
                if (u.numVertex(2) < v.numVertex(2))
                  return false;
                else if (u.numVertex(2) > v.numVertex(2))
                  return true;
                else
                  {
                    int u1, u2, v1, v2;
                    if (u.numVertex(1) < u.numVertex(3))
                      {
                        u1 = u.numVertex(1);
                        u2 = u.numVertex(3);
                      }
                    else
                      {
                        u1 = u.numVertex(3);
                        u2 = u.numVertex(1);
                      }
                    
                    if (v.numVertex(1) < v.numVertex(3))
                      {
                        v1 = v.numVertex(1);
                        v2 = v.numVertex(3);
                      }
                    else
                      {
                        v1 = v.numVertex(3);
                        v2 = v.numVertex(1);
                      }
                    
                    if (u1 < v1)
                      return false;
                    else if (u1 > v1)
                      return true;
                    else
                      {
                        if (u2 > v2)
                          return true;
                        else
                          return false;
                      }
                  }
              }
          }
      }
    
    return false;
  }
  
  
  /********************
   * Face<Dimension2> *
   ********************/  
  
  
  //! Displays informations about 2-D face
  ostream& operator<<(ostream& out, const Face<Dimension2>& e)
  {
    out << static_cast<const Face_Base&>(e);
    if (e.IsPML())
      out<<"Face inside the PML domain "<<endl;
    else if (e.IsCurved())
      out<<"Face curved "<<endl;
    
    return out;
  }
  
  
  /********************
   * Face<Dimension3> *
   ********************/  

  
  //! Displays informations about 3-D face
  ostream& operator<<(ostream& out, const Face<Dimension3>& e)
  {
    out << static_cast<const Face_Base&>(e);
    out << "Elements adjacent to the face : "<<e.num_element[0]+1<<"  "<<e.num_element[1]+1;
    out << endl;
    return out;
  }
  
  
  /**********
   * Volume *
   **********/
  
  
  //! default constructor
  Volume::Volume()
  {
    ref = 0;
    curved_elt = false; pml_num = 0; pml_type = 0;
    hybrid_type = 0;
    nb_vertices = 0; nb_edges = 0; nb_faces = 0;
  }
  

  //! initialization of vertices
  void Volume::Init(const IVect& num, int ref_)
  {
    if (num.GetM() == 4)
      InitTetrahedral(num, ref_);
    else if (num.GetM() == 5)
      InitPyramidal(num, ref_);
    else if (num.GetM() == 6)
      InitWedge(num, ref_);
    else if (num.GetM() == 8)
      InitHexahedral(num, ref_);
  }
    
  
  //! tetrahedral volume
  void Volume::InitTetrahedral(int n0, int n1, int n2, int n3, int ref_)
  {
    nb_vertices = 4; nb_edges = 6; nb_faces = 4;
    num_vertex(0) = n0; num_vertex(1) = n1; num_vertex(2) = n2; num_vertex(3) = n3;
    
    ref = ref_;
    hybrid_type = 0;
    num_edge.Fill(-1);
    num_face.Fill(-1);
    way_face.Fill(0);
    
    way_edge(0) = (n1 > n0);
    way_edge(1) = (n2 > n0);
    way_edge(2) = (n3 > n0);
    way_edge(3) = (n2 > n1);
    way_edge(4) = (n3 > n1);
    way_edge(5) = (n3 > n2);
  }
    
  
  //! pyramidal volume
  void Volume::InitPyramidal(int n0, int n1, int n2, int n3, int n4, int ref_)
  {
    nb_vertices = 5; nb_faces = 5; nb_edges = 8;
    num_vertex(0) = n0; num_vertex(1) = n1; num_vertex(2) = n2;
    num_vertex(3) = n3; num_vertex(4) = n4;
    
    ref = ref_;
    hybrid_type = 1;
    num_edge.Fill(-1);
    num_face.Fill(-1);
    way_face.Fill(0);
    
    way_edge(0) = (n1 > n0);
    way_edge(1) = (n2 > n1);
    way_edge(2) = (n2 > n3);
    way_edge(3) = (n3 > n0);
    way_edge(4) = (n4 > n0);
    way_edge(5) = (n4 > n1);
    way_edge(6) = (n4 > n2);
    way_edge(7) = (n4 > n3);
  }
  
  
  //! wedge
  void Volume::InitWedge(int n0, int n1, int n2, int n3, int n4, int n5, int ref_)
  {
    nb_vertices = 6; nb_edges = 9; nb_faces = 5;
    num_vertex(0) = n0; num_vertex(1) = n1; num_vertex(2) = n2;
    num_vertex(3) = n3; num_vertex(4) = n4; num_vertex(5) = n5;
    
    ref = ref_;
    hybrid_type = 2;
    num_edge.Fill(-1);
    num_face.Fill(-1);
    way_face.Fill(0);
    
    way_edge(0) = (n1 > n0);
    way_edge(1) = (n2 > n1);
    way_edge(2) = (n2 > n0);
    way_edge(3) = (n3 > n0);
    way_edge(4) = (n4 > n1);
    way_edge(5) = (n5 > n2);
    way_edge(6) = (n4 > n3);
    way_edge(7) = (n5 > n4);
    way_edge(8) = (n5 > n3);
  }
    
  
  //! hexahedral volume
  void Volume::InitHexahedral(int n0, int n1, int n2, int n3,
			      int n4, int n5, int n6, int n7, int ref_)
  {
    nb_vertices = 8; nb_edges = 12; nb_faces = 6;
    num_vertex(0) = n0; num_vertex(1) = n1; num_vertex(2) = n2; num_vertex(3) = n3;
    num_vertex(4) = n4; num_vertex(5) = n5; num_vertex(6) = n6; num_vertex(7) = n7;
    
    ref = ref_;
    hybrid_type = 3;
    num_edge.Fill(-1);
    num_face.Fill(-1);
    way_face.Fill(0);
    
    way_edge(0) = (n1 > n0);
    way_edge(1) = (n2 > n1);
    way_edge(2) = (n2 > n3);
    way_edge(3) = (n3 > n0);
    way_edge(4) = (n4 > n0);
    way_edge(5) = (n5 > n1);
    way_edge(6) = (n6 > n2);
    way_edge(7) = (n7 > n3);
    way_edge(8) = (n5 > n4);
    way_edge(9) = (n6 > n5);
    way_edge(10) = (n6 > n7);
    way_edge(11) = (n7 > n4);
  }
  

  //! returns true if the volume u is the same as v
  bool operator ==(const Volume& u, const Volume& v)
  {
    if (u.nb_vertices != v.nb_vertices)
      return false;

    return (u.num_vertex == v.num_vertex);
  }


  //! returns true if the volume u is different from v
  bool operator !=(const Volume& u, const Volume& v)
  {
    return !(u==v);
  }
  
  
  //! Displays information about the volume
  ostream& operator<<(ostream& out, const Volume& e)
  {
    out << "Volume with vertices :";
    for (int j = 0; j < e.GetNbVertices(); j++)
      out << e.numVertex(j)+1<<" ";
    
    out << endl << "Volume with edges : ";
    for (int j = 0; j < e.GetNbEdges(); j++)
      out << e.numEdge(j) + 1 << " ";
    
    out << endl << "Volume with faces ";
    for (int j = 0; j < e.GetNbFaces(); j++)
      out << e.numFace(j) + 1 << " ";

    out << endl << "Orientation of faces : ";
    for (int j = 0; j < e.GetNbFaces(); j++)
      out << e.GetOrientationFace(j) << " ";
    
    out << endl;
    if (e.IsPML())
      out << "Element inside the PML domain "<<endl;
    else if (e.IsCurved())
      out << "Element curved "<<endl;
    
    out << endl;
    return out;
  }


  ostream& operator<<(ostream& out, const ElementNumbering& e)
  {
    out << "Dof numbers = " << e.GetNodle() << '\n';
    return out;
  }

  
} // end namespace

#define MONTJOIE_FILE_MESH_ELEMENT_CXX
#endif
