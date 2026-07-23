#ifndef MONTJOIE_FILE_MESH_BASE_HXX

// declarations to use METIS as tool to split the mesh into sub-domains
extern "C"
{
  // Metis 5.1.0
#ifdef METIS_INTSIZE64
  typedef int64_t metis_int_t;
  typedef double metis_real_t;
#else
  typedef int metis_int_t;
  typedef float metis_real_t;
#endif

#ifndef MONTJOIE_WITHOUT_METIS
  int METIS_PartGraphKway(metis_int_t *nvtxs, metis_int_t *ncon, metis_int_t *xadj, metis_int_t *adjncy, 
                          metis_int_t *vwgt, metis_int_t *vsize, metis_int_t *adjwgt, metis_int_t *nparts, 
                          metis_real_t *tpwgts, metis_real_t *ubvec, metis_int_t *options, metis_int_t *objval, 
                          metis_int_t *part);
#endif
  
}

#ifdef MONTJOIE_WITH_SCOTCH
#include "scotch.h"
#endif


namespace Montjoie
{

  int GetDimensionMesh(const string& name);
  int GetDimensionMesh(const string& name, const string& extension);
  
  
  //! base class for mesh (2-D or 3-D)
  template<class Dimension, class T = Real_wp>
  class Mesh_Base
  {
    // typedef
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R2
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< R2 or R2
    typedef typename Dimension::VectR_N VectR_N; //!< vector of R2 or R3
    //! vector of edges or faces
    typedef Edge<Dimension> Edge_Rn;
    //! vector of elements (Faces in dimension 2, Volumes in dimension 3)
    typedef typename Dimension::Element_Rn Element_Rn;

    // splitting of the mesh
    enum {SUBDIV_DOMAIN_BOXES, SUBDIV_DOMAIN_METIS, SUBDIV_DOMAIN_SCOTCH};
    //! level of display
    int print_level;
    
  protected :

    //! List of vertices of the mesh
    VectR_N Vertices;
    
    //! bounding box for the mesh
    Real_wp xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;
    
    // EdgeRef are referenced edges
    // ("special" edges -> boundary edges, interface edges, curved edges)
    // Edges are all the edges of the mesh
    //! list of the edges of the mesh
    Vector<Edge_Rn> edges;
    //! list of the referenced edges of the mesh
    Vector<Edge_Rn> edges_ref;
    //! list of the elements of the mesh
    Vector<Element_Rn> elements;
    Vector<ElementGeomReference<Dimension>* > reference_element;
    
    //! path where the mesh file is searched
    string mesh_path; 
    //!< if we want to use a mesh split twice
    //! (for example, splitt in triangles, then split in quadrangles)    
    bool irregular_mesh_double_cut; 
    //! true if element i is affine
    Vector<bool> affine_element;
    int required_mesh_type; //!< type of required mesh (hexahedral, tetrahedral, etc)
    
    VectR_N ListeVertices_to_be_added; //!< list of vertices to add to the mesh
    VectR_N ListeVertices_to_be_refined; //!< list of vertices to refine around
    IVect LevelRefinement_Vertex; //!< refinement level around vertex
    Vector<T> RatioRefinement_Vertex; //!< ratio of refinement
    Vector<string> ParamRefinement_Area; //!< refinement parameters
    
    // head_minv(i) is the number of the first edge whose origin is vertex i
    // head_minv(i) is the number of the first edge whose extremity is vertex i
    // next_edge(e) is the next edge whose origin/extremity is vertex n
    //      where the vertex n is the origin/extremity of the edge e 
    IVect head_minv; //!< first edge whose origin is vertex i
    IVect head_maxv; //!< first edge whose end is vertex i
    IVect next_edge; //!< next edge after head_minv
    IVect next_edge_arr; //!< next edge after head_maxv    

    //! references given to "neighbor boundaries" initially not present in the big mesh
    IVect reference_original_neighbor;
    
    // access to leaf class
    Mesh<Dimension, T>& GetLeafClass();
    const Mesh<Dimension, T>& GetLeafClass() const;

  public :
    // public members
    
    //! 1/sqrt(J) (for Warburton trick)
    Vector<VectReal_wp> Glob_invSqrtJacobian;
    Vector<VectReal_wp> Glob_invSqrtJacobianNeighbor;
    //! gradient of jacobian (for Warburton trick)
    Vector<VectR_N> Glob_GradJacobian;
    
    //! global number for elements
    IVect GlobElementNumber_Subdomain;
    //! global number for faces (used in 3-D)
    IVect GlobFaceNumber_Subdomain;
    //! global number for edges
    IVect GlobEdgeNumber_Subdomain;
    //! global number for vertices
    IVect GlobVertexNumber_Subdomain;

    // methods        
    Mesh_Base();
    
    /*****************
     * Basic methods *
     *****************/
    
    size_t GetMemorySize() const;
    
    int GetNbElt() const;
    int GetNbVertices() const;
    int GetNbEdges() const;
    const IVect& GetOriginalNeighborReference() const;
    void SetOriginalNeighborReference(const IVect& ref);
    void SetPathName(const string&);
    const string& GetPathName() const;
    bool IsElementAffine(int i) const;
    
    void ReallocateElements(int);
    void ResizeElements(int);
    void ClearElements();
    Element_Rn& Element(int i);
    const Element_Rn& Element(int i) const;
    int GetTypeElement(int i) const;
    
    Edge_Rn& GetEdge(int i);
    const Edge_Rn& GetEdge(int i) const;
    void ReallocateEdges(int);
    void ResizeEdges(int);

    void ReallocateVertices(int);
    void ResizeVertices(int);
    void ClearVertices();
    R_N& Vertex(int i);
    const R_N& Vertex(int i) const;
    Vector<R_N>& Vertex();
    const Vector<R_N>& Vertex() const;
    
    T GetXmin() const;
    T GetXmax() const;
    T GetYmin() const;
    T GetYmax() const;
    T GetZmin() const;
    T GetZmax() const;
    int FindVertexNumber(const R_N& pt_glob) const;    
    
    const ElementGeomReference<Dimension>& GetReferenceElement(int i) const;
    const Vector<ElementGeomReference<Dimension>* >& GetReferenceElement() const;

    
    /**********************
     * Convenient methods *
     **********************/
    
    
    void GetVerticesElement(int num_elt, VectR_N& points) const;
    
    void GetDofsElement(int num_elt, VectR_N& points,
			const ElementGeomReference<Dimension>& FaceBasis) const;
    
    void SetInputData(const string& description_field, const VectString& parameters);    
    bool PointsOnSameEdge(int i, int ni, int& ne) const;
    void GetNeighboringElementsAroundVertex(int i, IVect& num) const;
    void GetNeighboringEdgesAroundVertex(int i, IVect& num) const;

    void SetVerticesToBeAdded(const VectR_N& points);
    void SetVerticesToBeRefined(const VectR_N& points, int lvl, const T& alpha);
    
    void CreateReferenceVertices(IVect& ref_vertex) const;
    T GetMeshSize() const;
    
    void ClearConnectivity();
    void ClearInput();
    void Clear();
    
  protected:
    void RemoveSingleEdge(int num_edge);
    void AddSingleEdge(int& nb_edges, const Edge<Dimension>&);
    void AddEdgeOfElement(int& nb_edges, int num_elem, int num_loc);
    void AddEdgesFromElements();

  public:
    // methods to split in several meshes
    void AddOverlapRegion(int delta_overlapping, int& nb_vertices_subdomain,
			  int& nb_elt_subdomain, Vector<bool>& VertexOnSubdomain,
			  Vector<bool>& ElementOnSubdomain) const;
    
    void SplitIntoBoxes(int nb_subdivX, int nb_subdivY, int nb_subdivZ,
			int& nb_subdomains, Vector<IVect>& NumElement_Subdomain,
			Vector<Mesh<Dimension, T> >& sub_mesh,
                        int delta_overlap, bool change_periodic_ref = true);
    
    void SplitMetis(int nb_parts, const IVect& weight_elt, Vector<IVect>& NumElement_Subdomain,
		    Vector<Mesh<Dimension, T> >& sub_mesh, int delta_overlap,
		    bool change_periodic_ref = true);
    
    void SplitScotch(int nb_parts, const IVect& weight_elt, Vector<IVect>& NumElement_Subdomain,
		     Vector<Mesh<Dimension, T> >& sub_mesh, int delta_overlap,
		     bool change_periodic_ref = true);
    
    void SplitConcentric(int nb_parts, const VectReal_wp& radius,
                         Vector<IVect>& NumElement_Subdomain,
                         Vector<Mesh<Dimension, T> >& sub_mesh,
                         int delta_overlap, bool change_periodic_ref = true);

    void SplitLayered(int nb_parts, const Vector<IVect>& ref_layer,
                      Vector<IVect>& NumElement_Subdomain,
                      Vector<Mesh<Dimension, T> >& sub_mesh,
                      int delta_overlap, bool change_periodic_ref = true);
    
    void PartMesh(int nb_parts, const IVect& epart, Vector<IVect>& NumElement_Subdomain,
                  Vector<Mesh<Dimension, T> >& sub_mesh,
                  int delta_overlapping, bool change_periodic_ref = true);

    void UpdatePeriodicReferenceSplit(const IVect& Epart, const IVect& NumLoc,
				      Vector<Mesh<Dimension, T> >& sub_mesh);
    
    void WriteOrder(const string&,
		    const MeshNumbering<Dimension, T>& mesh_num);

    void WriteColor(const string&,
		    const Vector<int>& color_elt);

    /**********************************************************************
     * Methods to use mesh for localization of points (GridInterpolation) *
     **********************************************************************/
    
    void FjElemNodal(const VectR_N&, SetPoints<Dimension>&,
		     const Mesh<Dimension>&, int) const;

    void DFjElemNodal(const VectR_N&, const SetPoints<Dimension>&,
                      SetMatrices<Dimension>&, const Mesh<Dimension>&, int) const;
    
    void Fj(const VectR_N&, const SetPoints<Dimension>&,
            const R_N&, R_N&, const Mesh<Dimension>&, int) const; 
    
    void DFj(const VectR_N&, const SetPoints<Dimension>&,
	     const R_N&, MatrixN_N&, const Mesh<Dimension>&, int) const; 
    
    int GetNbPointsNodalElt(int i) const;
    int GetNbPointsNodalBoundary(int num_elem, int num_loc) const;
    int GetNodalNumber(int i, int num_loc, int k) const;
    void GetNormale(const MatrixN_N&, R_N&, Real_wp&, int, int) const;
    
    Real_wp GetDistanceToBoundary(const R_N&, int) const;
    int ProjectPointOnBoundary(R_N&, int) const;
    void GetBoundingBox(const VectR_N& s, const Real_wp& coef,
                        VectR_N& box, TinyVector<R_N, 2>& enveloppe) const;
    
    void GetBoundingBox(int i, const VectR_N& s,
                        const SetPoints<Dimension>& Pts, TinyVector<R_N, 2>& enveloppe) const;
    
    void ComputeValuesPhiNodalRef(int i, const R_N&, VectReal_wp&) const;
    
    template<class Vect1, class Vect2>
    void ComputeNodalGradient(const SetMatrices<Dimension>&,
                              const Vect1&, Vect2&, int i) const;
    
  };
  
  
  //! base class for surfacic meshes
  template<class Dimension, class T>
  class SurfacicMesh_Base : public Mesh<Dimension, T>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    
    IVect ListeBoundaries; //!< edge/face numbers in the volumic mesh
    IVect ListeVertices; //!< vertex numbers in the volumic mesh

    IVect NumElement; //!< for each boundary, element number
    IVect NumLocalBoundary; //!< for each boundary, local position in the element
    Matrix<R_N> PointsNodal; //!< nodal points
    Matrix<R_N> NormaleNodal; //!< normale related to each point
    Matrix<T> DsjNodal; //!< surfacic integration
    Matrix<MatrixN_N> MatricesNodal; //!< jacobian matrices DF_i
    
    SurfacicMesh_Base();
    
    void SetBoundaryNodal(int i, SetPoints<Dimension>& pts,
			  SetMatrices<Dimension>& mat) const;
    
    void ReallocateNodal(int m, int n);
    
  };
  
  
  // exception raised when a method in mesh fails
  class WrongMesh : public Error
  {
  public:
    WrongMesh(string function = "", string comment = "");
    
  };
  


} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BASE_HXX
#endif
