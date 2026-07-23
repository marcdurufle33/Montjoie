#ifndef MONTJOIE_FILE_NUMBER_MESH_HXX

#include "NumberMap.hxx"

namespace Montjoie
{

  //! base class for the numbering of a mesh
  template<class T = Real_wp>
  class MeshNumbering_Base
  {
  protected :
    char order; //!< order of approximation
    int nodl; //!< number of dofs on all the mesh
    Vector<short int> nb_dof_element; //!< number of dofs for each element
    
    //! how order is computed for each element
    int type_variable_order;
    
    //! which order to use following the local mesh size
    Vector<T> mesh_size_variable_order;
    //! coefficient to apply for each reference
    Vector<T> coef_elt_variable_order;
    //! which order to use for a given reference
    Vector<char> order_for_each_ref;
    bool use_max_velocity_order;
    
    //! list of order inside each element
    Vector<char> order_inside;
    //! list of order for each edge
    Vector<char> order_edge;
    //! list of order for each face (3-D)
    Vector<char> order_face;

    //! order for each element (maximal order)
    Vector<char> order_element;
    //! order of quadrature for each face (maximal order)
    Vector<char> order_boundary;
    
    //! number of quadrature points for each edge/face (DG)
    Vector<short int> nb_points_quadrature_on_boundary;
    IVect IndexDofPML;
    int nb_dof_pml;
    
    /* Periodic dofs */
    
    IVect PeriodicityBoundary; //!< original boundaries (edge/face numbers)
    IVect type_periodicity_dof, type_periodicity_boundary;

    IVect PeriodicDofs; //!< periodic dof numbers
    IVect PeriodicDofs_Original; //!< original dof numbers
    
    //! dof numbers
    Vector<ElementNumbering> dof_element;    
    
    //! how periodicity is taken into account in the numerical scheme
    int type_formulation_periodic;    
    bool plane_wave_quasi_periodic;
    
    //! neighboring elements
    IVect order_neighbor_elt, type_neighbor_elt, num_edge_neighbor_elt, inv_num_edge_neighbor_elt;

  public:
    //! periodicity can be taken into account with the following types :
    enum {SAME_PERIODIC_DOFS, STRONG_PERIODIC, WEAK_PERIODIC};    

    enum {CONSTANT_ORDER, MAX_EDGE_ORDER,
          MEAN_EDGE_ORDER, USER_ORDER, REF_ORDER};
    
    //! used rules to number the mesh
    NumberMap number_map;
        
    //! Arrays incrementing the number of dofs associated with vertices, edges, faces, elements
    IVect OffsetDofFaceNumber;
    IVect OffsetDofEdgeNumber;
    IVect OffsetDofElementNumber;
    IVect OffsetDofVertexNumber;
    
    //! Array incrementing the number of quadrature points for each element
    IVect OffsetQuadElementNumber;
    
    //! global dof numbers (as if all subdomains were joined)
    IVect GlobDofNumber_Subdomain;

    //! global dof numbers for PML
    IVect GlobDofPML_Subdomain;
    
    bool treat_periodic_condition_during_number;
    bool compute_dof_pml;
    bool drop_interface_pml_dof;
    
  public:
    MeshNumbering_Base();
    virtual ~MeshNumbering_Base();
    
    // non-inline functions
    void CopyInputData(const MeshNumbering_Base<T>& mesh_num);
    void ReallocateNeighborElements(int n);
    void FinalizeNeighborElements(int N);
    
    size_t GetMemorySize() const;

    void SetInputData(const string& description_field, const VectString& parameters);
    
    virtual void ClearDofs();
    void Clear();

    void NumberMesh(bool surface_mesh = false);
    
    // inline functions    
    ElementNumbering& Element(int i);
    const ElementNumbering& Element(int i) const;

    /* Neighboring elements */
    
    int GetNbNeighborElt() const;
    int GetOrderNeighborElement(int i) const;
    int GetTypeNeighborElement(int i) const;
    int GetEdgeNeighborElement(int i) const;
    int GetLocalEdgeNumberNeighborElement(int i) const;

    void SetOrderNeighborElement(int i, int r);
    void SetTypeNeighborElement(int i, int type);
    void SetEdgeNeighborElement(int i, int ne);
    
    /* Periodic dofs */
    
    int GetPeriodicityTypeForDof(int i) const;
    void SetPeriodicityTypeForDof(int i, int);
    int GetPeriodicityTypeForBoundary(int i) const;
    void SetPeriodicityTypeForBoundary(int i, int);
    int GetNbPeriodicDof() const;    
    int GetPeriodicDof(int) const;
    void SetPeriodicDof(int, int);

    int GetOriginalPeriodicDof(int) const;
    void SetOriginalPeriodicDof(int, int);

    int GetPeriodicBoundary(int) const;
    void SetPeriodicBoundary(int, int);
    
    void SetSameNumberPeriodicDofs();
    void SetFormulationForPeriodicCondition(int);
    
    int GetFormulationForPeriodicCondition() const;
    bool UseSameNumberForPeriodicDofs() const;
        
    /* Dof stuff */
    
    int GetOrder() const;
    int GetNbDof() const;
    void SetNbDof(int n);
    int GetNbLocalDof(int i) const;
    
    int GetNbDofPML() const;
    int GetDofPML(int i) const;
    
    void ReallocateDofPML(int n);
    void SetDofPML(int i, int j);

    virtual void SetOrder(int);
  protected:
    virtual void SetGeometryOrder(int r, bool only_quad = false) = 0;

  public:
    bool IsOrderVariable() const;
    int GetVariableOrder() const;
    void SetVariableOrder(int);
    void SetMeshSizeVariableOrder(const Vector<T>&);
    void SetCoefficientVariableOrder(const Vector<T>&);
    
    int GetOrderElement(int num_elt) const;
    
    int GetOrderEdge(int ne) const; 
    int GetOrderFace(int ne) const; 
    int GetOrderInside(int ne) const; 

    void SetOrderEdge(int ne, int r); 
    void SetOrderFace(int ne, int r); 
    void SetOrderInside(int ne, int r); 

    int GetNbPointsQuadratureBoundary(int) const;

  protected:
    virtual void PerformMeshNumbering() = 0;
    virtual void PerformSurfaceMeshNumbering() = 0;
    
  };

  
  //! intermediate class for the numbering of a mesh
  template<class Dimension, class T = Real_wp>
  class MeshNumbering_Dim : public MeshNumbering_Base<T>
  {
  protected :
    typedef TinyVector<T, Dimension::dim_N> R_N;
    
    //! reference to the mesh
    Mesh<Dimension, T>& mesh;
    
    Vector<R_N> TranslationPeriodicBoundary;
    //!< translation vector between periodic boundary and original boundary

    Vector<R_N> TranslationPeriodicDof; 
    //!<  translation vector between periodic dof and original dof 
    
    // leaf class
    MeshNumbering<Dimension, T>& GetLeafClass();
    const MeshNumbering<Dimension, T>& GetLeafClass() const;
    
  public:    
    MeshNumbering_Dim(Mesh<Dimension, T>&);
    MeshNumbering_Dim(const MeshNumbering_Dim<Dimension, T>&);
    
    // inline methods
    int GetOrder() const;
    const Mesh<Dimension, T>& GetMesh() const;

    void InitPeriodicBoundary();
    
    const TinyVector<T, Dimension::dim_N>& GetTranslationPeriodicDof(int) const;
    void SetTranslationPeriodicDof(int, const TinyVector<T, Dimension::dim_N>&);
    const TinyVector<T, Dimension::dim_N>& GetTranslationPeriodicBoundary(int) const;
    Vector<TinyVector<T, Dimension::dim_N> >& GetTranslationPeriodicBoundary();    
    
    // other methods
    MeshNumbering_Dim<Dimension, T>& operator=(const MeshNumbering_Dim<Dimension, T>& m);
    
    void ReallocateElements(int n);
    void ReallocatePeriodicDof(int n);
    void ClearDofs();
    
    size_t GetMemorySize() const;
    
    void GetOrder(TinyVector<IVect, 4>&) const;    
    void ComputeVariableOrder(int dg_form);
    
  protected:
    void SetGeometryOrder(int r, bool only_quad = false);
    
    void PerformMeshNumbering();
    void PerformSurfaceMeshNumbering();
    
    void PerformLocalNumbering(int i, int j, int nb_dof_boundary, int offset_loc, int nf, int rf,
			       int& nb_negative_dof, IVect& NegativeDof);
    
  };
  

  template<class T>
  void InitStaticMeshData(const Dimension2& dim, const T& x);

  
  //! class to deal with mesh numbering in 2-D
  template<class T>
  class MeshNumbering<Dimension2, T> : public MeshNumbering_Dim<Dimension2, T>
  {
    friend class MeshNumbering_Dim<Dimension2, T>;
    
  public :
    // FirstExtremityEdge gives the numbers of first extremities of each edge.
    // For the quadrangle, FirstExtremityEdge(1,i) = first extremity
    // of the local edge i
    // For the triangle, it is FirstExtremityEdge(0,i)
    //! first extremity of each edge of reference elements
    static TinyMatrix<int, General, 2, 4> FirstExtremityEdge;
    //! second extremity of each edge of reference elements
    static TinyMatrix<int, General, 2, 4> SecondExtremityEdge;

  public :
    
    MeshNumbering(Mesh<Dimension2, T>&);

    // inline methods
    static int GetLocalVertexBoundary(int type, int num_loc, int i);

    int GetOrderQuadrature(int num_boundary) const;    
    void SetOrderQuadrature(int ne, int r);
    const VectReal_wp& GetReferenceQuadrature(int, int) const;

    void SetOrderElement(int num_elem, int r);
    
    static int GetRotationFace(int way_face1, int way_face2, int nv);
    static int GetOppositeOrientationFace(int rot, const Edge<Dimension2>&);

    // non-inline methods
    int GetNbDofElement(int i, const Edge<Dimension2>& elt) const;
    int GetNbDofElement(int i, const Face<Dimension2>& elt) const;
    
    void GetOrderQuadrature(TinyVector<IVect, 4>&) const;
    
    int GetBoundaryRotation(int f, int e1, int& e2, int pos1,
                            int& pos2, int& rot) const;

    void ComputeVariableOrder(int dg_form);
    void SetVariableOrder(int type);

    void ReconstructOffsetDofs();

    int GetNbDofBoundary(int ne);
    void GetDofBoundary(int i, IVect& num);
        
    static void ConstructQuadrilateralNumbering(int order, Matrix<int>& NumNodes2D,
						Matrix<int>& CoordinateNodes);
    
    static void ConstructTriangularNumbering(int order, Matrix<int>& NumNodes2D,
					     Matrix<int>& CoordinateNodes);
    
    void TreatPeriodicCondition();    

  protected:
    void NumberFaces();
    void NumberSurfaceEdges();
    
  };

  
#ifdef MONTJOIE_WITH_THREE_DIM
  template<class T>
  void InitStaticMeshData(const Dimension3& dim, const T& x);
  
  //! class to deal with mesh numbering in 3-D
  template<class T>
  class MeshNumbering<Dimension3, T> : public MeshNumbering_Dim<Dimension3, T>
  {
    friend class MeshNumbering_Dim<Dimension3, T>;
    
  public :
    
    // Extremities of faces/edges on unit tetrahedron/hexahedron
    // three extremities of a triedre are formed on each vertex N0
    // so that (N0,N1,N2,N3) is a direct triedre
    //! first extremity of edges on unit element (unit tet or unit cube)
    static TinyMatrix<int, General, 4, 12> FirstExtremityEdge;
    //! second extremity of edges on unit element
    static TinyMatrix<int, General, 4, 12> SecondExtremityEdge;
    //! first extremity of faces on unit element
    static TinyMatrix<int, General, 4, 6> FirstExtremityFace;
    //! second extremity of facs on unit element
    static TinyMatrix<int, General, 4, 6> SecondExtremityFace;
    //! third extremity of faces on unit element
    static TinyMatrix<int, General, 4, 6> ThirdExtremityFace;
    //! fourth extremity of faces on unit element
    static TinyMatrix<int, General, 4, 6> FourthExtremityFace;
    //! first edge of faces on unit element
    static TinyMatrix<int, General, 4, 6> FirstEdgeFace;
    //! second edge of faces on unit element
    static TinyMatrix<int, General, 4, 6> SecondEdgeFace;
    //! third edge of faces on unit element
    static TinyMatrix<int, General, 4, 6> ThirdEdgeFace;
    //! fourth edge of faces on unit element
    static TinyMatrix<int, General, 4, 6> FourthEdgeFace;

    static TinyMatrix<int, General, 6, 6> RotationFaceTri; //!< orientation of triangular faces
    static TinyMatrix<int, General, 8, 8> RotationFaceQuad; //!< orientation of quadrilateral faces
    
    //! true if the local face j of an elementis quadrilateral
    static TinyMatrix<bool, General, 4, 6> IsLocalFaceQuadrilateral;
    //! first vertex of triedre
    static TinyMatrix<int, General, 4, 8> FirstVertexTriedre;
    //! second vertex of triedre
    static TinyMatrix<int, General, 4, 8> SecondVertexTriedre;
    //! third vertex of triedre
    static TinyMatrix<int, General, 4, 8> ThirdVertexTriedre;
    //! first edge of triedre
    static TinyMatrix<int, General, 4, 8> FirstEdgeTriedre;
    //! second edge of triedre
    static TinyMatrix<int, General, 4, 8> SecondEdgeTriedre;
    //! third edge of triedre
    static TinyMatrix<int, General, 4, 8> ThirdEdgeTriedre;
    //! first edge of triedre
    static TinyMatrix<bool, General, 4, 8> BeginFirstEdgeTriedre;
    //! second edge of triedre
    static TinyMatrix<bool, General, 4, 8> BeginSecondEdgeTriedre;
    //! third edge of triedre
    static TinyMatrix<bool, General, 4, 8> BeginThirdEdgeTriedre;
    
  public :
    
    MeshNumbering(Mesh<Dimension3, T>&);

    // non-inline methods
    void GetOrderQuadrature(TinyVector<IVect, 4>&) const;
    void SetVariableOrder(int type);

    int GetNbDofElement(int i, const Face<Dimension3>& elt) const;
    int GetNbDofElement(int i, const Volume& elt) const;

    int GetBoundaryRotation(int f, int e1, int& e2, int pos1,
			    int& pos2, int& rot) const;
    
    void ComputeVariableOrder(int dg_form);
    void CheckMeshNumber();

  protected:
    void NumberFaces();
    void NumberSurfaceEdges();

  public:
    void ReconstructOffsetDofs();
    
    static void ConstructTetrahedralNumbering(int order, Array3D<int>& NumNodes3D,
					      Matrix<int>& CoordinateNodes);

    static void ConstructPyramidalNumbering(int order, Array3D<int>& NumNodes3D,
					  Matrix<int>& CoordinateNodes);

    static void ConstructPrismaticNumbering(int order, Array3D<int>& NumNodes3D,
					    Matrix<int>& CoordinateNodes);
    
    static void ConstructHexahedralNumbering(int order, Array3D<int>& NumNodes3D,
					     Matrix<int>& CoordinateNodes);


    static void GetRotationTriangularFace(const Matrix<int>& NumNodes2D,
					  Matrix<int>& FacesTri_Rotation);
    
    static void GetRotationTriangularFace(const VectR2& Points2D,
					  Matrix<int>& FacesTri_Rotation);
    
    static void GetRotationQuadrilateralFace(const Matrix<int>& NumNodes2D,
					     Matrix<int>& FacesQuad_Rotation);
    
    static void GetRotationQuadrilateralFace(const VectR2& Points2D,
                                             Matrix<int>& FacesQuad_Rotation);
    
    void TreatPeriodicCondition();

    // inline methods
    int GetOrderQuadrature(int num_boundary) const;        
    void SetOrderQuadrature(int ne, int r);
    const VectR2& GetReferenceQuadrature(int, int) const;

    void SetOrderElement(int num_elem, int r);
    
    static int GetRotationFace(int way_face1, int way_face2, int nv);
    
    static int GetOppositeOrientationFace(int rot, bool triangular_face);    
    static int GetOppositeOrientationFace(int rot, const Face<Dimension3>&);
        
    static int GetVertexNumberOfFace(int type, int nf, int num);
    static int GetLocalVertexBoundary(int type, int nf, int num);
    
    static int GetEdgeNumberOfFace(int type, int nf, int num);
    
  };
    
#endif
  
  
} // namespace Montjoie

#define MONTJOIE_FILE_NUMBER_MESH_HXX
#endif
