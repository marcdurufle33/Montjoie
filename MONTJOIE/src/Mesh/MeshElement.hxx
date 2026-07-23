#ifndef MONTJOIE_FILE_MESH_ELEMENT_HXX

namespace Montjoie
{
  
  //! base class for edge
  class Edge_Base
  {
  protected:
    //! two vertices
    TinyVector<int, 2> num_vertex;
    int ref; //!< reference of the edge
    
  public :
    Edge_Base();
    
    // initialisation
    void Init(int n1, int n2);
    void Init(int n1, int n2, int ref_);
    
    int GetReference() const;
    void SetReference(int ref_);
    
    int numVertex(int ne) const;
    static int GetNbVertices();
    
    void ClearConnectivity();
        
    friend ostream& operator<<(ostream& out, const Edge_Base&);
    
    // operators to compare edges
    template<class DimensionB>
    friend bool operator==(const Edge<DimensionB>& u, const Edge<DimensionB> & v);
    
    template<class DimensionB>
    friend bool operator!=(const Edge<DimensionB> & u, const Edge<DimensionB>& v);
    
    template<class DimensionB>
    friend bool operator<(const Edge<DimensionB> & u, const Edge<DimensionB>& v);
  
    template<class DimensionB>
    friend bool operator>(const Edge<DimensionB> & u, const Edge<DimensionB>& v);

  };
  
  
  //! 1-D edge
  template<>
  class Edge<Dimension1>  : public Edge_Base
  {
  protected : 
    bool pml_elt;
    
  public :
    Edge();

    int GetMemorySize() const;
    void SetPML();
    void UnsetPML();
    bool IsPML() const;
        
    friend ostream& operator<<(ostream& out, const Edge<Dimension1>&);
    
  };
  
  
  //! 2-D edge
  template<>
  class Edge<Dimension2> : public Edge_Base
  {
  protected :
    //! array of elements that contain the edge
    int num_element[2];
    int geom_ref;
    
  public :
    Edge();

    int GetMemorySize() const;    
    int GetNbElements() const;
    int numElement(int ne) const;
    void AddElement(int nf);
    void SetElement(int nf0, int nf1);
    
    static int GetHybridType();
    
    int GetGeometryReference() const;
    void SetGeometryReference(int ref_);
    
    void ClearConnectivity();

    friend ostream& operator<<(ostream& out, const Edge<Dimension2>& );
  };
  
  
  //! 3-D edge
  template <>
  class Edge<Dimension3> : public Edge_Base
  {
  protected :
    //! list of elements and faces
    //! arrays are regrouped to save memory
    IVect num_elt_face;
    char nb_faces;
    
  public :
    Edge();
       
    int GetMemorySize() const;
    int GetNbFaces() const;
    int GetNbElements() const;
    
    int numFace(int ne) const;
    int numElement(int ne) const;
    
    void AddFace(int nf);
    void AddElement(int nf);
    
    void ClearConnectivity();
    void ClearFaces();
    
    friend ostream& operator<<(ostream& out, const Edge<Dimension3>& );
    
  };

  
  //! base class for a face
  class Face_Base
  {
  protected:
    TinyVector<int, 4> num_vertex; //!< list of vertices
    TinyVector<int, 4> num_edge; //!< list of edges
    TinyVector<bool, 4> way_edge; //!< orientation of edges
    short int ref; //!< reference of the face    
    char nb_vertices;
    char hybrid_type;
    
  public :
  
    Face_Base();

    bool IsTriangular() const;
    bool IsQuadrangular() const;
    int GetHybridType() const;
    
    bool GetOrientationEdge(int num_edge) const;
    int GetOrientationBoundary(int num_edge) const;
    
    int numVertex(int ne) const;
    int numEdge(int ne) const;
    
    int GetReference() const;
    int GetNbVertices() const;
    int GetNbEdges() const;
    
    void SetVertex(int num_loc, int num_glob);
    void SetReference(int ref_);
    
    void InitTriangular(int n1, int n2, int n3, int ref_);
    void InitQuadrangular(int n1, int n2, int n3, int n4, int ref_);
    void Init(const IVect& num, int ref_);
    
    void SetEdge(int num_loc, int ne);
    int GetPositionBoundary(int ne) const;

    void ClearConnectivity();
    void ClearEdges();
    
    int GetOrientationFace(int) const;
    
    friend ostream& operator<<(ostream& out, const Face_Base&);
    
    // operators to compare faces
    template<class DimensionB>
    friend bool operator==(const Face<DimensionB>& u, const Face<DimensionB> & v);
    
    template<class DimensionB>
    friend bool operator!=(const Face<DimensionB> & u, const Face<DimensionB>& v);
    
    template<class DimensionB>
    friend bool operator<(const Face<DimensionB> & u, const Face<DimensionB>& v);
  
    template<class DimensionB>
    friend bool operator>(const Face<DimensionB> & u, const Face<DimensionB>& v);
  };
  
  
  //! 2-D face
  template<>
  class Face<Dimension2> : public Face_Base
  {
  protected :
    bool curved_elt;
    unsigned char pml_type, pml_num;
    
  public:
    Face();

    int GetMemorySize() const;    
    void SetPML(int num, int type);
    void UnsetPML();
    bool IsPML() const;
    int GetTypePML() const;
    int GetNumberPML() const;
    
    void SetCurved();
    void UnsetCurved();
    bool IsCurved() const;
    
    int GetNbBoundary() const;
    int numBoundary(int) const;
    
    friend ostream& operator<<(ostream& out, const Face<Dimension2>& );
    
  };

  
  //! 3-D face
  template <>
  class Face<Dimension3> : public Face_Base
  {
  protected :
    //!< elements that contain the face - adjacents to the face - (one or two)
    int num_element[2];
    int geom_ref;
    
  public:
    Face(); //!< default constructor
    
    int GetMemorySize() const;
    int GetNbElements() const;
    int numElement(int ne) const;
    
    void AddElement(int ne);
    void SetElement(int n1,int n2);
    
    int GetGeometryReference() const;
    void SetGeometryReference(int ref_);
    
    void ClearConnectivity();
    
    friend ostream& operator<<(ostream& out, const Face<Dimension3>& );
    
  };
  
  
  //! class for 3-D element (tetrahedron, hexahedron ...)
  class Volume
  {
  protected :
    TinyVector<int, 8> num_vertex; //!< list of vertices
    TinyVector<int, 12> num_edge; //!< list of edges
    TinyVector<int, 6> num_face; //!< list of face
    short int ref; //!< reference of the element
    TinyVector<bool, 12> way_edge; //!< orientation of edges
    TinyVector<char, 6> way_face; //!< orientation of local faces compared to global faces
    bool curved_elt;
    unsigned char pml_type, pml_num;
    char hybrid_type;
    char nb_vertices, nb_edges, nb_faces;
    
  public:
    Volume();

    int GetMemorySize() const;    
    int GetReference() const;
    void SetReference(int ref_);
    
    int GetHybridType() const;
    
    void SetPML(int num, int type);
    void UnsetPML();
    bool IsPML() const;
    int GetTypePML() const;
    int GetNumberPML() const;
    
    void SetCurved();
    void UnsetCurved();
    bool IsCurved() const;
    
    int GetNbVertices() const;
    int GetNbEdges() const;
    int GetNbFaces() const;
    int GetNbBoundary() const;
    
    int numVertex(int ne) const;
    int numEdge(int ne) const;
    int numFace(int nf) const;
    int numBoundary(int nf) const;
    
    bool GetOrientationEdge(int num_edge) const;
    
    void SetOrientationFace(int num_face,int way);
    int GetOrientationFace(int num_face) const;
    int GetOrientationBoundary(int num_face) const;
    
    void Init(const IVect& num, int ref_);
    
    void InitTetrahedral(const IVect& num, int ref_);
    void InitTetrahedral(int n0, int n1, int n2, int n3, int ref_);
    
    void InitPyramidal(const IVect& num, int ref_);
    void InitPyramidal(int n0, int n1, int n2, int n3, int n4, int ref_);
    
    void InitWedge(const IVect& num,int ref_);
    void InitWedge(int n0, int n1, int n2, int n3, int n4, int n5, int ref_);
    
    void InitHexahedral(const IVect& num,int ref_);
    void InitHexahedral(int n0, int n1, int n2, int n3, int n4,
			int n5, int n6, int n7, int ref_);
    
    void SetEdge(int num_loc, int ne);
    
    void SetFace(int num_loc, int nf);
    void SetFace(int num_loc, int nf, int way);
    
    int GetPositionEdge(int ne) const;
    int GetPositionBoundary(int nf) const;

    void ClearConnectivity();
    void ClearEdges();
    void ClearFaces();

    friend bool operator ==(const Volume&, const Volume&);
    friend bool operator !=(const Volume&, const Volume&);
    
    friend ostream& operator<<(ostream& out, const Volume& );
    
  };
  
  class ElementNumbering
  {
  protected :
    IVect Nodle; //!< local to global dofs numbering
    IVect negative_dof_num;
    //IVect phase_dof_num;
    //IVect type_phase_dof;
    TinyVector<char, 6> order_face;
    TinyVector<int, 6> offset_face;
    
  public :
    int GetMemorySize() const;
    int GetNbDof() const;
    int GetNumberDof(int i) const;
    const IVect& GetNodle() const;
    const IVect& GetNegativeDofNumber() const;
    void SetNegativeDofNumber(int nb, const IVect& num);

    void ReallocateFaces(int nf);
    int GetOrderFace(int nf) const;
    int GetOffsetFace(int nf) const;
    void SetOrderFace(int, int rf, int);

    void ReallocateDof(int nb_dof);
    void SetNumberDof(int num_dof_loc, int num_dof_glob);
    void AddDof(int nb_dof);
    void PushDof(int num_dof);
        
    void Clear();
    
  };

  ostream& operator<<(ostream& out, const ElementNumbering&);  
  
} // end namespace

#define MONTJOIE_FILE_MESH_ELEMENT_HXX
#endif
