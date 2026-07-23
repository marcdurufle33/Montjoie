#ifndef MONTJOIE_FILE_MESH3D_BOUNDARIES_HXX

namespace Montjoie
{
  template<class T>
  class IntersectionCurvedEdgeFace_System;
  
  template<class T>
  class VirtualMinimizedFunction;

  template<class T>
  class SurfaceParamMinimizer;


  //! parameters specifying parameters of the PML region to add
  template<class T>
  class PmlRegionParameter<Dimension3, T> : public PmlRegionParameter_Base<Dimension3, T>
  {
  public:
    // combinations are obtained by summation
    enum {PML_Xm = 1, PML_Xp = 2, PML_Ym = 4, PML_Yp = 8, 
          PML_Zm = 16, PML_Zp = 32};
    
    T GetAverageLengthOnBoundary(const Mesh<Dimension3>&) const;

    void AddPMLElements(int nb_div, int num_pml, Mesh<Dimension3>&);
    inline int GetNbParameters() const { return 0; }
    inline void FillParameters(VectReal_wp& all_param, int& nb_param) const {}
    inline void SetRegion(const VectReal_wp& all_param, int& nb_param) {}
    
  };

  
  //! class for curved boundaries of 3-D meshes
  template<class T>
  class MeshBoundaries<Dimension3, T> : public MeshBoundaries_Base<Dimension3, T>
  {
    friend class NumberMap;
    template<class T0>
    friend class IntersectionCurvedEdgeFace_System;
    
  public :
    typedef TinyVector<T, 3> R_N;
    
    // type of 3-D curves
    enum {NO_CURVE, CURVE_CYLINDER, CURVE_SPHERE, CURVE_CONIC, CURVE_FILE};

  protected :
    Matrix<R_N> PointsEdgeRef; //!< projection of interior points of curved edges
    Vector<Vector<R_N> > PointsFaceRef; //!< projection of interior points of curved faces
    Vector<Vector<T> > value_parameter_boundary;
    Vector<TinyVector<int, 4> > integer_parameter_boundary;
    Vector<TinyVector<T, 4> > value_parameter_edge;
    Vector<TinyVector<int, 2> > integer_parameter_edge;
    
    // interpolation on curved faces
    TriangleGeomReference triangle_reference;
    QuadrangleGeomReference quadrangle_reference;
  
  public :
    NumberMap nodal_number_map;
    
    MeshBoundaries();
        
    /************************
     * Convenient functions *
     ************************/
    
    
    size_t GetMemorySize() const;    

    bool IsBoundaryCurved(int) const;
    const R_N& GetPointInsideEdge(int, int) const;
    const R_N& GetPointInsideFace(int, int) const;
    const Vector<R_N>& GetPointInsideFace(int) const;
    void SetPointInsideEdge(int, int, const R_N&);
    void SetPointInsideFace(int, int, const R_N&);
    void SetPointInsideFace(int, const Vector<R_N>&);
    int GetNbPointsInsideFace(int) const;
    void SetInputData(const string& description_field, const VectString& parameters);
    void RemoveReference(int ref);    

    const TriangleGeomReference& GetSurfaceFiniteElement() const;
    const QuadrangleGeomReference& GetSurfaceFiniteElement2() const;
    
    
    /**************************************
     * Treatment of referenced boundaries *
     **************************************/
    
    
    void SortPeriodicBoundaries();
    void SortEdgesRef();
    void SortBoundariesRef(bool keep_orientation = false);
    void GetGeometryFaceRotation(Matrix<int>& RotationTri, Matrix<int>& RotationQuad) const;
    void ConstructCrossReferenceBoundaryRef();
    void ConstructCrossReferenceEdgeRef();
    void AddBoundaryFaces();
    void AddBoundaryEdges();
    
    void ProjectVerticesOnCurves();
    void ProjectPointsOnCurves();
    void RedistributeReferences();
    void GetNeighboringBoundariesAroundVertex(int i, IVect& num) const;
    
    
    /**********************************
     * Treatment of curved boundaries *
     **********************************/

    void CheckCurveParameter(int ref);    
    void GetNodalPointsOnFace(int i, Vector<R_N>& PtFace) const;
  protected:
    void ComputeInternalPointsFaceWithCurvedEdges(int i);
    void GetParameterValueOnPoint(int ref, const R_N& pt, const VectReal_wp& param, R_N& s) const;
    void ComputeParameterValueOnCurves();
    void GetPointsOnCurvedEdge(int i, int ref, const VectReal_wp& xi, VectR3& points) const;    
    void GetPointsOnCurvedFace(int i, int ref, const VectR2& xi, VectR3& points, bool straight_face = false) const;

  public:
    void SetStraightEdge(int i);
    void SetStraightFace(int i);
    
    const R2& GetNodalPointInsideTriangle(int i) const;

    void ReadCurveType(const IVect& ref, const string& description);
    static void ProjectToSphere(R3& ptA, const R3& center, const T& radius);
    static void ProjectToConicCylinder(R3& ptM, const R3& ptC,
                                       const R3& vdirector, const T& theta);
    static void ProjectToCylinder(R3& ptM, const R3& ptC, const R3& vdirector, const T& radius);
    void ProjectToCurve(R3& ptM, int ref) const;
    
    void GetAllPointsOnReference(int ref, VectR3& Points, VectR3& Normale) const;
    void GetAllPointsOnReference(int ref, const IVect& ref_cond,
                                 VectR3& Points, VectR3& Normale) const;
    void GetDistantPointsOnReference(int ref, VectR3& Points,
                                     VectR3& Normale, int nb_points) const;
    void GetDistantPointsOnReference(int ref, const IVect& ref_cond,
                                     VectR3& Points, VectR3& Normale, int nb_points) const;
    
    void FindParametersCurve();
    void FindParametersPlane(int ref, Vector<T>& coef_plane) const;
    void FindParametersSphere(int ref, R3& center, T& radius) const;
    void FindParametersConic(int ref, R3& ptC, R3& vdirector, T& theta) const;
    void FindParametersCylinder(int ref, R3& ptC, R3& vdirector, T& radius) const;

    void FindParametersGeneric(int ref, Vector<T>& param,
			       SurfaceParamMinimizer<T>& fct) const;

    
    /*****************
     * Other methods *
     *****************/
    
    
    void ClearConnectivity();

    // create surfacic mesh associated with boundary condition ref
    void GetBoundaryMesh(int ref, SurfacicMesh<Dimension3>& mesh,
                         const IVect& ref_cond, int ref_domain = 0) const;
    
    template<class Mesh3D_>
    void GetBoundaryMesh(int ref, Mesh3D_& mesh_surf, IVect& Index_FaceSurf_to_FaceRef,
			 IVect& Index_VertexSurf_to_Vertex, IVect&, IVect&,
                         const IVect& reference_condition, int ref_domain = 0) const;
        
  };
  
  
  template<class T>
  class SurfaceParamMinimizer : public VirtualMinimizedFunction<T>,
                                public VirtualLeastSquaresFunction<T>
  {
  public :
    typedef TinyVector<T, 3> R_N;
    
    int nb_points; //!< number of points
    Vector<R_N> Points; //!< points used to find parameters
    Vector<R_N> Normales; //!< normales related to each point
    
    inline int GetM() { return this->m_; }
    
  };
  

  template<class T>
  class IntersectionCurvedEdgeFace_System
  {
  private:
    const Mesh<Dimension3>& mesh;
    int num_face1, num_edge_loc1, num_edge1;
    Real_wp pos_edge;
    int ref1; R3 center; Real_wp radius;
    bool pt_on_edge;
    
  public:
    IntersectionCurvedEdgeFace_System(const Mesh<Dimension3>& mesh0, int nf, int num_loc, const Real_wp& pos,
                                      int ref, const R3& pt0, const Real_wp& r);

    IntersectionCurvedEdgeFace_System(const Mesh<Dimension3>& mesh0, int ne,
                                      int ref, const R3& pt0, const Real_wp& r);
    
    R3 GetPointOnSurface(const Real_wp& x) const;
    
    Real_wp EvaluateFunction(const Real_wp& x) const;
    Real_wp EvaluateDerivative(const Real_wp& x) const;
    
  };
  
  
  //! class to find parameters of a cylinder (center, axis, radius)
  template<class T>
  class CylinderParam : public SurfaceParamMinimizer<T>
  {
  public :
    CylinderParam();
    void FindInitGuess(Vector<T>& param);
    void EvaluateFunction(const Vector<T>& x, T& feval);
    void EvaluateFunctionGradient(const Vector<T>& x, T& feval, Vector<T>& fjac);

    void EvaluateF(const Vector<T>& x, Vector<T>& feval);
    void EvaluateJacobian(const Vector<T>& x, Vector<T>& feval,
                          Matrix<T, General, ColMajor>& fjac);
    
  };
  
  //! class to find parameters of a conic cylinder (center, axis, angle)
  template<class T>
  class ConicParam : public SurfaceParamMinimizer<T>
  {
  public :    
    
    //! default constructor
    ConicParam();
    void FindInitGuess(Vector<T>& param);
    void EvaluateFunction(const Vector<T>& x, T& feval);
    void EvaluateFunctionGradient(const Vector<T>& x, T& feval, Vector<T>& fjac);

    void EvaluateF(const Vector<T>& x, Vector<T>& feval);
    void EvaluateJacobian(const Vector<T>& x, Vector<T>& feval, 
                          Matrix<T, General, ColMajor>& fjac);
    
  };

  
  //! intermediary class of class ConicParam
  template<class T>
  class SubConicParam : public NonLinearEquations_Base<T, Vector<T>, Matrix<T> >
  {
  public :
    typedef TinyVector<T, 3> R_N;
    
    R_N u1; //!< first vector
    R_N u2; //!< second vector
    R_N u3; //!< third vector
    
    void EvaluateFunction(const Vector<T>& param, Vector<T>& feval);
    void FindInitGuess(Vector<T>& param) const;
    void EvaluateJacobian(const Vector<T>& param, Matrix<T>& fjac);
    T Solve(Vector<T>& param);
  };

}

#define MONTJOIE_FILE_MESH3D_BOUNDARIES_HXX
#endif
