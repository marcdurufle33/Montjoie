#ifndef MONTJOIE_FILE_MESH3D_HXX

namespace Montjoie
{

  //! 3-D mesh containing tetrahedra and hexahedra
  template <class T>
  class Mesh<Dimension3, T> : public Mesh_Base<Dimension3, T>,
			      public MeshBoundaries<Dimension3, T>
  {
    // main attributes
  public :
    typedef TinyVector<T, 3> R_N;
    friend class MeshBoundaries<Dimension3, T>;
    friend class MeshBoundaries_Base<Dimension3, T>;
    
    // type of mesh
    enum {HYBRID_MESH, TETRAHEDRAL_MESH, HEXAHEDRAL_MESH,
	  PYRAMID_MESH, WEDGE_MESH, CRUSH_MESH};
    
    // internal attributes
  protected :
    // faces : all the faces of the mesh, even inside faces
    // faces_ref -> referenced faces, especially on boundaries and on the interfaces
    //              between two domains
    //! faces of the mesh
    Vector<Face<Dimension3> > faces;
    //! referenced faces of the mesh
    Vector<Face<Dimension3> > faces_ref;
    
    // head_faces_minv(i) contains the first face which leaves the vertex i
    // next_face(nf) contains the next face after nf
    // if no other face leaves the vertex, -1
    //! array generated to find faces of the mesh
    IVect next_face;
    IVect head_faces_minv; //!< first face that leaves vertex i
    TetrahedronGeomReference tetrahedron_reference;
    PyramidGeomReference pyramid_reference;
    WedgeGeomReference wedge_reference;       
    HexahedronGeomReference hexahedron_reference;       
    int regular_mesh_with_slight_modification;
    enum {REGULAR_MESH_NOT_MODIFIED, REGULAR_MESH_SLIGHTLY_MODIFIED,
          REGULAR_MESH_RANDOMLY_MODIFIED};

    void UpdateReferenceElement();
    
  public :
    Mesh();

    Mesh(const Mesh<Dimension3, T>& mesh);
    Mesh<Dimension3, T>& operator=(const Mesh<Dimension3, T>& mesh);
    
    void ReallocateBoundariesRef(int);
    void ResizeBoundariesRef(int);
    void ReallocateEdgesRef(int);
    void ResizeEdgesRef(int);
    void ReallocateVertices(int);
    void ResizeVertices(int);
    void ReallocateFaces(int);
    void ResizeFaces(int);
    Edge<Dimension3>& EdgeRef(int i);
    const Edge<Dimension3>& EdgeRef(int i) const;
    Face<Dimension3>& BoundaryRef(int i);
    const Face<Dimension3>& BoundaryRef(int i) const;
    Face<Dimension3>& Boundary(int i);
    const Face<Dimension3>& Boundary(int i) const;
    int GetNbBoundary() const;
    int GetNbBoundaryRef() const;
    int GetGlobalBoundaryNumber(int i) const;
    
    const TriangleGeomReference& GetTriangleReferenceElement() const;
    const QuadrangleGeomReference& GetQuadrilateralReferenceElement() const;
    
    void GetReferenceElementSurface(Vector<const ElementGeomReference<Dimension2>* >& ref) const;
    void GetReferenceElementVolume(Vector<const ElementGeomReference<Dimension3>* >& ref) const;
    IVect GetLocalNodalNumber(int type, int num_loc);
    
    int GetNbTetrahedra() const;
    int GetNbPyramids() const;
    int GetNbWedges() const;
    int GetNbHexahedra() const;
    int GetNbTrianglesRef() const;
    int GetNbQuadranglesRef() const;
    int GetNbEdgesRef() const;
    int GetNbTriangles() const;
    int GetNbQuadrangles() const;
    int GetNbFaces() const;
    void SetGeometryOrder(int, bool only_quad = false);
    
    
    /**********************
     * Convenient methods *
     **********************/
    
    
    bool HybridMesh() const;
    bool IsOnlyTetrahedral() const;
    bool IsOnlyHexahedral() const;

    size_t GetMemorySize() const;
    
    int GetUpperBoundNumberOfFaces() const;
    int GetUpperBoundNumberOfEdges() const;
    
    void SetInputData(const string& description_field, const Vector<string>& parameters);
    int GetNodesCurvedMesh(Vector<R_N>& PosNodes, Vector<IVect>& Nodle,
			   int rmax, bool lobatto = false) const;

    int GetNodesCurvedMesh(Vector<R_N>& PosNodes, Vector<IVect>& Nodle,
			   Vector<IVect>& NodleSurf, int rmax, bool lobatto = false) const;

    // Read a mesh
    void Read(const string & file_mesh);
    void Read(istream& file_in, const string& ext);
    
    // write mesh on file    
    void Write(const string & file_name, bool double_prec = true, bool ascii = true) const;
    void Write(ostream& file_out, const string& ext,
	       bool double_prec = true, bool ascii = true) const;
    
    void ChangeLocalNumberingMesh();
    
    void ClearConnectivity();
    // destruction of all the mesh
    void Clear();
    
    
    /******************
     * Mesh Treatment *
     ******************/
    
    
    // permute numbers of vertices of element i
    // so that (n0,n1,n2,n3) is a direct triedre
    void ReorientElements();
    void ReorientFaces(int, const IVect&, const R3& normale);
    static void ReorientElementCartesian(IVect& num, VectR3& ptE);
    // find faces of the mesh
    void FindConnectivity();
    // find edges from referenced faces (instead of elements)
  protected:
    void AddEdgeOfFace(int& nb_edges, int i, int na);
    void AddSingleFace(int& nb, const Face<Dimension3>& );
    void AddFaceOfElement(int& nb, int i, int num_loc);
    void AddFacesFromElements();
    void AddEdgesFromReferencedFaces();
    void CreateConnexionEdgeToFace();

  public:
    void RemoveDuplicateVertices();
    void ForceCoherenceMesh();
    void ClearZeroBoundaryRef();

    // symmetrize the mesh
    void SymmetrizeMesh(int num_coor, T xc,
                        TinyVector<T, 3> normale, bool keep_face = false);
    
    // Cut each tetrahedron in four hexahedra
    void SplitIntoHexahedra();
    // Cut each hexahedron in six tetrahedra
    void SplitIntoTetrahedra();
    // divide a general mesh
    void GetMeshSubdivision(const Vector<T>& outside_subdiv,
                            const Vector<Vector<R_N> >& points_div,
			    const Vector<VectR2>& points_surf, Vector<R_N>& PosNodes,
			    Vector<IVect>& Nodle, Vector<IVect>& NodleSurf) const;
    
    void SubdivideMesh(const Vector<T>& step_subdiv);
    void SubdivideMesh(const Vector<T>& step_subdiv,
		       Vector<VectR2>& points_surf, Vector<VectR3>& points_div,
                       Vector<IVect>& Nodle, Vector<IVect>& NodleSurf);

    
    /*****************
     * Mesh creation *
     *****************/
    
    
    // Create a regular mesh hexahedral of tetrahedral
    void CreateRegularMesh(const R_N&, const R_N&, const TinyVector<int, 3>& nbPoints,
			   int ref_domain, const TinyVector<int, 6>& ref_boundary, int type_mesh,
                           R3 ratio = R3(1, 1, 1));
    
    void CreateStructuredMesh(const Vector<T>& pos_layer_x, const Vector<T>& pos_layer_y, const Vector<T>& pos_layer_z,
                              Vector<int>& nb_layers_x, Vector<int>& nb_layers_y, Vector<int>& nb_layers_z,
                              const Array3D<int>& ref_domain, const TinyVector<int, 6>& ref_boundary);
    
    // Create a mesh of a spherical crown (domain inside two sphere of radius a and b)
    void CreateSphericalCrown(const T& a, const T& b, const T&, bool, bool, const Vector<T>&, bool all_ref = false);
    void CreateSphericalBall(const T& a, const T& b, const T&,
                             int ref_inside, bool, int, int ref_outside = 1,
                             int nb_pts_sphere = 0, int nb_pts_line = 0,
			     int ref_extern = 3);
    
    // adding hexahedral elements on the outside boundary of the mesh
    // (assuming that the boundary is plane)
    void ExtrudeCoordinate(int num_coord, int nb_layers, const T& pos, const T& delta);
        
    // Extrude a surface from the origin
    void ExtrudeOrigin(const IVect&, const Mesh<Dimension3>&,
                       const IVect&, const IVect&, const Vector<T>&, bool all_ref = false);
    
    void ExtrudeSurfaceMesh(Mesh<Dimension2, T>& mesh_2d, const IVect& nbCells_layer,
			    const Vector<T>& step_layer, const Vector<IVect>& ref_surfacic,
                            const Vector<IVect>& ref_volumic,
			    const R3& vec_e1, const R3& vec_e2, const R3& vec_e3);
    
    void AppendMesh(const Mesh<Dimension3, T>& mesh, bool elimine = true);
    void CreateSubmesh(Mesh<Dimension3>& sub_mesh, int nb_vertices_subdomain, int nb_elt_subdomain,
		       const VectBool& VertexOnSubdomain, const VectBool& ElementOnSubdomain);

    void ConstructMesh(int type_mesh, const VectString& parameters);
    void RearrangeElements();
    
    void GetTetaMinMax(const R_N& xc, T& teta_min, T& teta_max) const;
    
    // periodize the mesh
    void PeriodizeMeshTeta(const R_N& xc);    

    void ApplyLocalRefinement(const Vector<int>& num_vertex, const Vector<int>&,
			      int nb_levels, const T& coef, int new_ref = 0);
    
    template<class T0>
    friend ostream& operator <<(ostream& out, const Mesh<Dimension3, T0>& mesh);

  public:
    IVect ref_vertices;
  };
  
  template<class T0>
  ostream& operator <<(ostream& out, const Mesh<Dimension3, T0>& mesh);
  
  class ScalingMeshVirtualFunction3D
  {
  public:
    virtual void EvaluateCoefficient(const Real_wp& z,
				     Real_wp& coef_x, Real_wp& coef_y) = 0;
  };
  
    
  class FctScalingExtrudeMesh : public ScalingMeshVirtualFunction3D
  {
  public :
    Real_wp z0, z1, max_ratio;
    
    void EvaluateCoefficient(const Real_wp& z,
                             Real_wp& coef_x, Real_wp& coef_y);
    
  };
  
  
  //! surfacic mesh in dimension 3
  template<class T>
  class SurfacicMesh<Dimension3, T> : public SurfacicMesh_Base<Dimension3, T>
  {
  public :
    //! matching array between surfacic dofs and volumic dofs
    IVect IndexDofInt_to_DofVol;
    
    //! default constructor
    SurfacicMesh();
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH3D_HXX
#endif
