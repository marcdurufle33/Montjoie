#ifndef MONTJOIE_FILE_MESH_INTERPOLATION_HXX

namespace Montjoie
{

  template<class Dimension>
  class MeshInterpolation
  {
  };
  
  //! base class to do interpolation on subdivided meshes
  template<class Dimension>
  class MeshInterpolation_Base
  {
  protected :
    typedef typename Dimension::DimensionBoundary DimensionBoundary;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::VectR_N VectR_N;
    typedef typename Dimension::R_N R_N;
    
    FiniteElementInterpolator proj;
    FiniteElementInterpolator proj_surf;
    
    //! number of subdivisions of the fine mesh
    int nb_subdivisions;
    //! references for surfacic mesh
    IVect ref_surf;
    //! numbering array
    Vector<IVect> Nodle;
    //! subdivided points on the reference elements
    Vector<typename Dimension::VectR_N > points_div;
    //! subdivided points on the surface reference elements
    Vector<typename DimensionBoundary::VectR_N > points_div_surf;
    //! weights associated with surface reference elements (if quadrature)
    TinyVector<VectReal_wp, 2> weights_surf;
    //! subdivisions on the unit edge [0,1]
    VectReal_wp step_subdiv;
    //! file where the mesh will be stored
    string file_name;
    //! type of mesh (surfacic or volumic)
    int type_mesh;
    
    //! Element where each surface is
    IVect NumElement;
    //! Local position of each surface on the element
    IVect LocalPositionBoundary;
    IVect TypeBoundary;
    //! Normales
    Vector<Vector<R_N> > NormaleBoundary, PointsBoundary, PointsNodalBoundary;
    //! w_j ds    
    Vector<VectReal_wp > WeightsBoundary, DsBoundary;
    //! jacobian matrices on quadrature points
    Vector<Vector<MatrixN_N> > DFiBoundary, invDFiBoundary;
    //! jacobian matrices on nodal points
    Vector<Vector<MatrixN_N> > DFiNodalBoundary, invDFiNodalBoundary;
    //! points, normales and weights for all processors
    Vector<R_N> AllPoints, AllNormale;
    VectReal_wp AllWeights;
    int nb_points_all_quadrature, rank_processor;
    
  public :
    IVect nb_points_quadrature_per_proc;
    IVect offset_quadrature_per_proc;
    
    //! available outputs
    enum {VOLUME_MESH, SURFACE_MESH, REFERENCE_MESH,
	  BODY_MESH, CONDITION_MESH};

    MeshInterpolation_Base();

    // Inline methods
    
    int GetNbSubdivisions() const;
    void SetRegularSubdivisions(int r);    
    void SetLobattoSubdivisions(int r);
    const VectReal_wp& GetSubdivisionStep() const;

    Vector<Vector<R_N> >& PointsReferenceVolume();
    Vector<typename DimensionBoundary::VectR_N>& PointsReferenceSurface();
    VectReal_wp& WeightsReferenceSurface(int);

    const string& GetFileName() const;
    void SetFileName(const string&);
    
    void SetType(int);
    int GetType() const;
    void SetSurfaceReference(const IVect& ref);
    const IVect& GetSurfaceReference() const;
    
    Vector<IVect>& NumberingArray();
    const Vector<IVect>& NumberingArray() const;

    void InitProjectionVolume(Mesh<Dimension>& vars);
    void InitProjectionSurface(Mesh<Dimension>& vars);
 
    template<class Vector1, class Vector2>
    void ProjectVolume(const Vector1& Uloc, Vector2& Udiv, int type_elt) const;

    template<class Vector1, class Vector2>
    void ProjectSurface(const Vector1& Uloc, Vector2& Udiv, int type_elt) const;
    
    bool IsVolumetric() const;

    bool IsSurfacic() const;

    int GetElementNumberOfSurface(int) const;
    int GetLocalPositionOfSurface(int) const;
    int GetNbBoundary() const;
    
    int GetNbPointsQuadrature() const;
    int GetNbPointsQuadrature(int i) const;

    const R_N& PointsNodal(int i, int j) const;
    const Vector<R_N>& PointsNodal(int i) const;
    
    const R_N& PointsQuadrature(int i, int j) const;
    const Vector<R_N>& PointsQuadrature(int i) const;

    const MatrixN_N& Dfjm1Quadrature(int i, int j) const;
    const Vector<MatrixN_N>& Dfjm1Quadrature(int i) const;
    
    const R_N& NormaleQuadrature(int i, int j) const;
    const Vector<R_N>& NormaleQuadrature(int i) const;
    
    Real_wp& WeightsQuadrature(int i, int j);
    const Real_wp& WeightsQuadrature(int i, int j) const;
    const Vector<Real_wp>& WeightsQuadrature(int i) const;

    const Real_wp& DsQuadrature(int i, int j) const;
    const Vector<Real_wp>& DsQuadrature(int i) const;
    
    const Vector<R_N>& GetNormale(int i) const;
    
    int GetNbAllQuadraturePoints() const;
    const R_N& GetQuadraturePoint(int k) const;
    const VectR_N& GetQuadraturePoints() const;
    const R_N& GetQuadratureNormale(int k) const;
    const Real_wp& GetQuadratureWeight(int k) const;
    
    // Other methods
    size_t GetMemorySize() const;

    void SetRegularPoints(int r);  
    void SetGaussQuadrature(int r, int type = Globatto<Real_wp>::QUADRATURE_GAUSS);
    
#ifdef SELDON_WITH_MPI
    void GatherQuadraturePoints(const MPI_Comm& comm, bool assemble = true);
#else
    void GatherQuadraturePoints(bool assemble = true);
#endif
    
  };

  
  template<>
  class MeshInterpolation<Dimension2> : public MeshInterpolation_Base<Dimension2>
  {
  public :
    void ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension2>& mesh,
			    Mesh<Dimension2>& mesh_subdiv, int ref_domain = 0);

    void StoreGeometricDataSurface(const Mesh<Dimension2>& mesh,
				   int i1, int num_elem, int num_loc);
  };


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  class MeshInterpolation<Dimension3> : public MeshInterpolation_Base<Dimension3>
  {
  public :
    void ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension3>& mesh,
			    Mesh<Dimension3>& mesh_subdiv, int ref_domain = 0);

    void StoreGeometricDataSurface(const Mesh<Dimension3>& mesh,
				   int i1, int num_elem, int num_loc);
  };
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_INTERPOLATION_HXX
#endif
