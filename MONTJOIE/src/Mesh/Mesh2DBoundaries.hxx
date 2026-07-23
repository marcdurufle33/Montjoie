#ifndef MONTJOIE_FILE_MESH2D_BOUNDARIES_HXX

namespace Montjoie
{

  //! parameters specifying parameters of the PML region to add
  template<class T>
  class PmlRegionParameter<Dimension2, T> : public PmlRegionParameter_Base<Dimension2, T>
  {
  public:
    // PML types
    // combinations are obtained by summation
    enum {PML_Xm = 1, PML_Xp = 2, PML_Ym = 4, PML_Yp = 8};
    
    T GetAverageLengthOnBoundary(const Mesh<Dimension2>&) const;
    
    void AddPMLElements(int nb_div, int num_pml, Mesh<Dimension2>&);
    
    int GetNbParameters() const;
    void FillParameters(VectReal_wp& all_param, int& nb_param) const;
    void SetRegion(const VectReal_wp& all_param, int& nb_param);
    
  };


  //! class for curved boundaries of 2-D meshes
  template<class T>
  class MeshBoundaries<Dimension2, T> : public MeshBoundaries_Base<Dimension2, T>
  {
  public :
    typedef TinyVector<T, 2> R_N;
    typedef Vector<R_N> VectR_N;
    
  protected :
    Matrix<R_N> PointsEdgeRef; //!< projection of interior points of curved edges
    //! value of parameter (for example teta for a circle) on each boundary
    Vector<R_N> value_parameter_boundary;
    //! case where parameter is an integer
    Vector<TinyVector<int, 2> > integer_parameter_boundary;
    EdgeLobatto interval_reference;

  public :
    // type of 2-D curves
    enum {NO_CURVE, CURVE_CIRCLE, CURVE_SPLINE,
	  CURVE_ELLIPSE, CURVE_PEANUT, CURVE_LOCAL_SPLINE, CURVE_FILE};
    
    MeshBoundaries();

    
    /************************
     * Convenient functions *
     ************************/
    
    
    size_t GetMemorySize() const;
    
    bool IsBoundaryCurved(int) const;
    const R_N& GetPointInsideEdge(int, int) const;
    void SetPointInsideEdge(int, int, const R_N& );
    
    int FindFollowingVertex(int n1, int n2, int ref, const IVect&, int& e) const;
    int FindFollowingVertex(int n1, int n2, int ref,
                            const IVect&, const Vector<bool>&, int& e) const;
    
    const EdgeLobatto& GetSurfaceFiniteElement() const;
    const EdgeLobatto& GetSurfaceFiniteElement2() const;

  protected:
    const Vector<R_N>& GetValueParameterOnBoundary() const;
    void ComputeParameterValueOnCurves();
    
  public:
    // is the point M on the edge num_glob_edge_ref ?
    bool IsPointOnBoundaryRef(int num_glob_edge_ref, const R_N& M) const;
    // all the references of the boundaries close to a vertex
    void GetAllReferencesAroundPoint(int num_point, IVect& ref) const;
    // get a point near num_point, and belonging 
    // to an edge with different reference from ref
    int GetPointWithOtherReference(int num_point, int ref, R_N&, int&) const;
    
    void ClearConnectivity();
    void ClearCurves();
    void Clear();
    
    /**************************************
     * Treatment of referenced boundaries *
     **************************************/
    
    void RemoveReference(int ref);
    void SortPeriodicBoundaries();
    void SortBoundariesRef();
    void ConstructCrossReferenceBoundaryRef(bool check_mesh = true);
    void AddBoundaryEdges();
    
    void ProjectPointsOnCurves();
    void RedistributeReferences();
    
    
    /**********************************
     * Treatment of curved boundaries *
     **********************************/
    
    
    void CheckCurveParameter(int ref);
    void ReadCurveType(const IVect& ref, const string& description);
    void ProjectToCurve(R_N& ptM, int ref) const;
    void GetPointsOnCurve(int i, int ref, const VectReal_wp& x, Vector<R_N>& points) const;
    
    void GetAllPointsOnReference(int ref, VectR_N& Points, VectR_N& Normale) const;
    void GetAllPointsOnReference(int ref, const IVect& ref_cond,
                                 VectR_N& Points, VectR_N& Normale) const;
    void GetDistantPointsOnReference(int ref, VectR_N& Points,
                                     VectR_N& Normale, int nb_points) const;
    void GetDistantPointsOnReference(int ref, const IVect& ref_cond,
                                     VectR_N& Points, VectR_N& Normale, int nb_points) const;
    
    void FindParametersPlane(int ref, Vector<T>& coef_plane) const;
    void FindParametersSpline(int ref, Vector<T>& coef, IVect& Point_curve);
    void FindParametersCircle(int ref, T& xc, T& yc, T& radius) const;
    void FindParametersEllipse(int ref, T& xc, T& yc, T& a, T& b) const;
    void FindParametersPeanut(int ref, T& xc, T& yc, T& a, T&b, T& c) const;
    void FindParametersCurve();

    
    /*****************
     * Other methods *
     *****************/
    

    TinyVector<T, 2> GetCenterRadialPML(T& Rmax) const;
    
    void GetBoundaryMesh(int ref, SurfacicMesh<Dimension2>& mesh,
                         const IVect& ref_cond, int ref_domain = 0) const;
    
    template<class Mesh2D_>
    void GetBoundaryMesh(int ref, Mesh2D_& mesh_surf, IVect& Index_EdgeSurf_to_EdgeRef,
			 IVect& Liste_Vertices, IVect&, IVect&,
                         const IVect& ref_cond, int ref_domain = 0) const;

#ifdef MONTJOIE_WITH_THREE_DIM
    void GenerateSurfaceOfRevolution(int ref, Mesh<Dimension3>& mesh_boundary_inside,
				     IVect& Index_FaceSurf_to_EdgeRef,
				     IVect& Index_VertexSurf_to_Vertex,
				     Vector<T>& AngleVertex, const IVect& ref_cond) const;  
    
#endif
    
  };


  //! Class used to project a point on an ellipse
  class EllipseMinimizer : public VirtualLeastSquaresFunction<Real_wp, TinyVector<Real_wp, 1>, TinyVector<Real_wp, 2>, TinyMatrix<Real_wp, General, 2, 1>, TinyMatrix<Real_wp, Symmetric, 1, 1> >
  {
    Real_wp x0, y0;
    R2 center;
    Real_wp a, b, teta0;
    
  public :
    EllipseMinimizer(const Real_wp& x, const Real_wp& y, const R2& C,
		     const Real_wp& a0, const Real_wp& b0, const Real_wp& t0);

    void SetInitialAngle(const Real_wp& t0);
    void FindInitGuess(TinyVector<Real_wp, 1>& x);
    void EvaluateF(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f);    
    void EvaluateJacobian(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f,
			  TinyMatrix<Real_wp, General, 2, 1>& df);
    
  };


  //! Class used to project a point on a peanut
  class PeanutMinimizer : public VirtualLeastSquaresFunction<Real_wp, TinyVector<Real_wp, 1>, TinyVector<Real_wp, 2>, TinyMatrix<Real_wp, General, 2, 1>, TinyMatrix<Real_wp, Symmetric, 1, 1> >
  {
    Real_wp x0, y0;
    R2 center;
    Real_wp a, b, c, teta0;
    
  public :
    PeanutMinimizer(const Real_wp& x, const Real_wp& y, const R2& C,
		    const Real_wp& a0, const Real_wp& b0, const Real_wp& c0, const Real_wp& t0);

    void SetInitialAngle(const Real_wp& t0);
    void FindInitGuess(TinyVector<Real_wp, 1>& x);
    void EvaluateF(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f);    
    void EvaluateJacobian(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f,
			  TinyMatrix<Real_wp, General, 2, 1>& df);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_MESH2D_BOUNDARIES_HXX
#endif

