#ifndef MONTJOIE_FILE_MESH2D_HXX

namespace Montjoie
{
  
  //! 2-D mesh containing triangles and quadrilaterals
  template<class T>
  class Mesh<Dimension2, T> : public Mesh_Base<Dimension2, T>,
			      public MeshBoundaries<Dimension2, T>
  {
  public :
    typedef TinyVector<T, 2> R_N;
    typedef Vector<R_N> VectR_N;
    friend class MeshBoundaries<Dimension2, T>;
    friend class MeshBoundaries_Base<Dimension2, T>;
    
    // types of mesh
    enum {HYBRID_MESH, TRIANGULAR_MESH, QUADRILATERAL_MESH, RADAU_MESH};

  protected:
    TriangleGeomReference triangle_reference;
    QuadrangleGeomReference quadrangle_reference;
    
    void UpdateReferenceElement();

  public :
    
    Mesh();

    Mesh(const Mesh<Dimension2, T>& mesh);
    Mesh<Dimension2, T>& operator=(const Mesh<Dimension2, T>& mesh);
    
    Edge<Dimension2>& Boundary(int i);
    const Edge<Dimension2>& Boundary(int i) const;
    int GetNbBoundary() const;
    int GetGlobalBoundaryNumber(int i) const;
    
    void ReallocateBoundariesRef(int);
    void ResizeBoundariesRef(int);
    Edge<Dimension2>& BoundaryRef(int i);
    const Edge<Dimension2>& BoundaryRef(int i) const;
    int GetNbBoundaryRef() const;
    int GetNbEdgesRef() const;

    const TriangleGeomReference& GetTriangleReferenceElement() const;
    const QuadrangleGeomReference& GetQuadrilateralReferenceElement() const;
    
    void GetReferenceElementSurface(Vector<const ElementGeomReference<Dimension1>* >& ref) const; 
    void GetReferenceElementVolume(Vector<const ElementGeomReference<Dimension2>* >& ref) const;
    IVect GetLocalNodalNumber(int type, int num_loc);
    
    bool HybridMesh() const;
    // returns true if only triangles are in the mesh
    bool IsOnlyTriangular() const;
    // returns true if only quadrilaterals are in the mesh
    bool IsOnlyQuadrilateral() const;
    
    size_t GetMemorySize() const;
    int GetNbTriangles() const;
    int GetNbQuadrangles() const;
    int GetNbHexahedra() const;
    void SetGeometryOrder(int, bool only_quad = false);

    
    /**********************
     * Convenient methods *
     **********************/    
    
    
    int GetUpperBoundNumberOfEdges() const;
    // Destruction of informations linking edges, faces, elements
    void ClearConnectivity();
    // destruction of all the mesh
    void Clear();    
    
    void SetInputData(const string& description_field, const Vector<string>& parameters);
    int GetNodesCurvedMesh(VectR_N& PosNodes, Vector<IVect>& Nodle, int rmax,
                           bool lobatto = false) const;
    
    // Read mesh on a file
    void Read(const string & file_mesh);
    void Read(istream& file_in, const string& ext);
    
    // write mesh on a file
    void Write(const string & file_name, bool double_prec = true, bool ascii = true) const;
    void Write(ostream& file_out, const string& ext, bool double_prec = true, bool ascii = true) const;
    
    
    /******************
     * Mesh Treatment *
     ******************/
    
    
    // method, which permuts numbers of vertices of an element i (n1,n2,n3)
    // so that this element is oriented with respect to the trigonometrical way
    void ReorientElements();
    void FindConnectivity();
    // "free" vertices and "free" edges are removed
    // so that all the vertices and edges belong to the faces of the mesh
    void RemoveDuplicateVertices();
    void ForceCoherenceMesh(bool keep_edge = false);
    void ClearZeroBoundaryRef();
    // refining mesh locally around a vertex
    void ApplyLocalRefinement(int num_vertex, int nb_levels, const T& coef);    
    void ApplyRefinementEdge(const IVect& ref_cond, int ref_target, bool anisotrope, 
                             int nb_levels, const T& coef);
    
    void AddVertexConformal(const TinyVector<T, 2>& pt);
    // symmetrize the mesh
    void SymmetrizeMesh(int num_coor, T xc,
                        TinyVector<T, 2> normale, bool keep_edge = false);
    // periodize the mesh
    void PeriodizeMeshTeta(const R_N& xc);    

    // Splits each triangle and three quadrilaterals
    // and each quadrilateral in four quadrilaterals
    void SplitIntoQuadrilaterals();
    // Splits each quadrangle in two triangles
    void SplitIntoTriangles();
    // divide a general mesh
    void GetMeshSubdivision(const Vector<T>& outside_subdiv,
                            const Vector<Vector<R_N> >& points_div,
			    Vector<R_N>& PosNodes, Vector<IVect>& Nodle, Vector<IVect>& NodleSurf) const;
    
    void SubdivideMesh(const Vector<T>& step_subdiv);
    void SubdivideMesh(const Vector<T>& step_subdiv, Vector<Vector<T> >& pts_surf,
                       Vector<Vector<R_N> >& points_div, Vector<IVect>& Nodle, Vector<IVect>& NodleSurf);
    
    void SymmetrizeMeshToGetPositiveX();

    
    /*****************
     * Mesh creation *
     *****************/
    
    
    // creation of a regular mesh
    // nbPoints_x  : number of points for each direction
    //               the meshed domain is the rectangle [xmin,xmax] x [ymin,ymax]
    void CreateRegularMesh(const R_N&, const R_N&, const TinyVector<int, 2>& nbPoints,
			   int ref_domain, const TinyVector<int, 4>& ref_boundary,
			   int type_mesh, R2 ratio = R2(1, 1));

    void CreateStructuredMesh(const Vector<T>& pos_layer_x, const Vector<T>& pos_layer_y,
                              Vector<int>& nb_layers_x, Vector<int>& nb_layers_y,
                              const Matrix<int>& ref_domain_layer, const TinyVector<int, 4>& ref_boundary);
    
    // a layer is added with respect to a direction (x or y)
    void ExtrudeCoordinate(int num_coor, int nb_layers, const T& pos, const T& delta);

    // construction of a mesh 
    void ConstructMesh(int type_mesh, const VectString& parameters);
    void RearrangeElements();
    // creation of a mesh from a subset of elements
    void CreateSubmesh(Mesh<Dimension2, T>& sub_mesh, int nb_vertices_subdomain,
		       int nb_elt_subdomain, const VectBool& VertexOnSubdomain,
		       const VectBool& ElementOnSubdomain) const;
    
    // a mesh is appended to this mesh
    void AppendMesh(Mesh<Dimension2, T>& mesh, bool elimine = true, R2 u = R2());
    
        
    /*****************
     * Other methods *
     *****************/
    
    template<class T0>
    friend ostream& operator <<(ostream& out, const Mesh<Dimension2, T0>& mesh);
    
    
    void GetEdgesOnZaxis(IVect& list_vertices, IVect& list_edges,
			 VectBool& Vertex_On_Axe);    
    
    void GetTetaMinMax(const R_N& xc, T& teta_min, T& teta_max) const;
    
  protected:
    void AddElement(int& nb, int& nb_edges, const Face<Dimension2>& e);
    
  };
  
  template<class T0>
  ostream& operator <<(ostream& out, const Mesh<Dimension2, T0>& mesh);
  
  template<class T>
  void ExtrudePMLLayer(Mesh<Dimension2, T>& mesh, const TinyVector<T, 2>& axis,
		       const T& scal_max, const Vector<int>& ref_cond, int ref_target,
		       const T& thickness, int nb_layers, bool circle_layer);
  
  //! class for surfacic mesh in dimension 2
  /*!
    Used for axisymmetric integral equations
   */
  template<class T>
  class SurfacicMesh<Dimension2, T> : public SurfacicMesh_Base<Dimension2, T>
  {
  public :    
    //! default constructor
    SurfacicMesh();
        
  };
  
} // end namespace

#define MONTJOIE_FILE_MESH2D_HXX  
#endif
