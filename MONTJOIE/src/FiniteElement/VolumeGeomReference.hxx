#ifndef MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_HXX

namespace Montjoie
{
  
  //! base class for 3-D finite element
  template<>
  class ElementGeomReference<Dimension3> : public ElementGeomReference_Base<Dimension3>
  {
    friend class ElementGeomReferenceContainer<Dimension3>;
    friend class ElementReference_Dim<Dimension3>;
    template<int t> friend class TetrahedronReference;
    template<int t> friend class PyramidReference;
    template<int t> friend class WedgeReference;
    template<int t> friend class HexahedronReference;
    
  protected :
    
    //! nodal point numbers on the faces
    Matrix<int> FacesNodal;
   
    VectR2 points_nodal2d_tri; //!< 2-D nodal points (on a face)
    VectR2 points_nodal2d_quad; //!< 2-D nodal points (on a face)
    VectR3 points_nodal3d; //!< 3-D nodal points 
    
    VectR2 points2d_tri; //!< 2-D integration points (integration on triangle)
    VectR2 points2d_quad; //!< 2-D integration points (integration on quadrangle)
    VectReal_wp weights2d_tri; //!< 2-D integration weights on triangles
    VectReal_wp weights2d_quad; //!< 2-D integration weights on quadrangles

    //! 1-D shape functions
    Globatto<Real_wp> lob_geom;
    
    VectR3 points3d; //!< 3-D integration points (integration on volume)
    VectReal_wp weights3d; //!< 3-D integration weights
    
    VectR2 points_dof2d_tri; //!< 2-D dof points (on a face)
    VectR2 points_dof2d_quad; //!< 2-D dof points (on a face)
    VectR3 points_dof3d; //!< 3-D points used to project on dof components
    
    Seldon::Array3D<int> NumNodes3D; //!< tensorial <-> scalar numbering array
    Matrix<int> CoordinateNodes; //!< tensorial <-> scalar numbering array
    
    Matrix<int> NumNodes2D_tri; //!< tensorial <-> scalar numbering array
    Matrix<int> CoordinateNodes2D_tri; //!< tensorial <-> scalar numbering array
    
    Matrix<int> NumNodes2D_quad; //!< tensorial <-> scalar numbering array
    Matrix<int> CoordinateNodes2D_quad; //!< tensorial <-> scalar numbering array

    Matrix<Real_wp> CoefCurveTri;    
    Matrix<int> FacesNodal_Rotation_Tri, FacesNodal_Rotation_Quad;
    
    Vector<bool> is_local_face_quad; // quadrilateral faces ?

    Matrix<Real_wp> RhTriGeomTrans;
    Matrix<int> NumQuad2D;
    
    VectR3 tangente_loc_x, tangente_loc_y;

  public :
    ElementGeomReference();
 
        
    /************************
     * Convenient functions *
     ************************/
    
    size_t GetMemorySize() const;
    VectR3 PointsQuadInsideND(int N) const;    
    int GetNbPointsNodalElt() const;
    int GetNbPointsDof() const;
    int GetNbPointsQuadrature() const;
        
    const Array3D<int>& GetNumNodes3D() const;
    int GetNumNodes3D(int, int, int) const;
    int GetNumNodes2D(int, int) const;
    const Matrix<int>& GetNumNodes2D_quad() const;
    const Matrix<int>& GetCoordinateNodes3D() const;
    
    int GetNodalNumber(int num_face_loc, int num_dof) const;
    const Matrix<int>& GetNodalNumber() const;
    
    const VectReal_wp& WeightsND() const;
    const VectR3& PointsQuadND() const;
    const VectR3& PointsND() const;
    const Real_wp& WeightsND(int j) const;
    const R3& PointsND(int j) const;
        
    const Vector<R3>& PointsDofND() const;
    const R3& PointsDofND(int j) const;

  protected:
    const VectR2& Points2D_tri() const;
    const VectR2& Points2D_quad() const;    
    const VectReal_wp& Weights2D_tri() const;
    const VectReal_wp& Weights2D_quad() const;    
    const VectR2& PointsNodal2D_tri() const;
    const VectR2& PointsNodal2D_quad() const;
    const VectR2& PointsDof2D_tri() const;
    const VectR2& PointsDof2D_quad() const;

    void SetPointsNodal2D_quad(const VectR2&);
    void SetPointsNodal2D_tri(const VectR2&);
    void SetPointsDof2D_quad(const VectR2&);
    void SetPointsDof2D_tri(const VectR2&);
    void SetPoints2D_quad(const VectR2&);
    void SetPoints2D_tri(const VectR2&);
    void SetWeights2D_quad(const VectReal_wp&);
    void SetWeights2D_tri(const VectReal_wp&);
    
    Matrix<int>& GetNumQuad2D();
    const Matrix<int>& GetNumQuad2D() const;
    const Matrix<int>& GetNumNodes2D_tri() const;
    const Matrix<int>& GetCoordinateNodes2D_quad() const;
    const Matrix<int>& GetCoordinateNodes2D_tri() const;

  public:
    int GetNbQuadBoundary(int) const;
    int GetNbNodalBoundary(int) const;
    bool IsLocalFaceQuadrilateral(int n) const;
    
    const R2& PointsNodalBoundary(int, int) const;
    const VectR2& PointsQuadratureBoundary(int) const;
    const VectR2& PointsDofBoundary(int) const;
    const R2& PointsQuadratureBoundary(int k, int) const;
    const Real_wp& WeightsQuadratureBoundary(int k, int) const;

    const R_N& TangenteLocX(int n) const;
    const R_N& TangenteLocY(int n) const;
    
    
    /************************
     * Overloaded functions *
     ************************/
    
    void SetMesh(Mesh<Dimension3>& mesh);

  protected:
    void FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad) const;    
    
    void GetCurvatureOnSurface(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			       SetMatrices<Dimension3>& res,
			       const Mesh<Dimension3>& mesh, int num_elem,
			       int num_face, int num_loc) const;
    
    virtual void GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
					      const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
					      R3& d_nx, R3& d_ny, R3& d_nz) const = 0;

  public:    
    virtual R3 GetTangentialVector(int num_loc, const R2& vec_u) const = 0;
    virtual R2 TransposeTangentialVector(int num_loc, const R3& vec_u) const = 0;
    virtual R2_Complex_wp TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const = 0;

    
    /*******************************************************
     * computation of geometric quantities on the boundary *
     *******************************************************/
    
    
    // we compute the enveloppe of the volume
    static void GetBoundingBox(const VectR3& s, const Real_wp coef,
			       VectR3& box, TinyVector<R3, 2>& enveloppe);

  protected:
    void InitTriangularCurvature(const TriangleGeomReference& tri);
    
    /*****************
     * Other methods *
     *****************/
  public:
    virtual void GetLocalCoordOnBoundary(int num_loc, const R2& t_loc, R3& res) const = 0;    
    virtual void GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& t_loc) const = 0;    
    
    static bool IsAffineTransformation(const VectR3& s);

    friend ostream& operator <<(ostream& out, const ElementGeomReference<Dimension3>& e);
    
  };

  ostream& operator <<(ostream& out, const ElementGeomReference<Dimension3>& e);

  template<>
  class ElementGeomReferenceContainer<Dimension3>
    : public ElementGeomReferenceContainer_Base<Dimension3>
  {
  public :
    ElementGeomReferenceContainer(ElementGeomReference<Dimension3>&);
    
    VectR3 PointsQuadInsideND(int N) const;    
    int GetNbPointsNodalElt() const;
    int GetNbPointsDof() const;
    int GetNbPointsQuadrature() const;
        
    const Array3D<int>& GetNumNodes3D() const;
    int GetNumNodes3D(int, int, int) const;
    int GetNumNodes2D(int, int) const;
    
    int GetNodalNumber(int num_face_loc, int num_dof) const;
    
    const Vector<R3>& PointsNodalND() const;
    const R3& PointsNodalND(int j) const;
    const Globatto<Real_wp>& GetNodalShapeFunctions1D() const;
    
    const VectReal_wp& WeightsND() const;
    virtual const VectR3& PointsQuadND() const;
    const VectR3& PointsND() const;
    const Real_wp& WeightsND(int j) const;
    const R3& PointsND(int j) const;
  protected:
    void SetPointsND(const VectR3& pts);
    void SetWeightsND(const VectReal_wp& w);

  public:
    const Vector<R3>& PointsDofND() const;
    const R3& PointsDofND(int j) const;

  protected:
    void SetPointsDofND(const VectR3& pts);

    const VectR2& Points2D_tri() const;
    const VectR2& Points2D_quad() const;    
    const VectReal_wp& Weights2D_tri() const;
    const VectReal_wp& Weights2D_quad() const;    
    const VectR2& PointsNodal2D_tri() const;
    const VectR2& PointsNodal2D_quad() const;
    const VectR2& PointsDof2D_tri() const;
    const VectR2& PointsDof2D_quad() const;

    void SetPointsDof2D_quad(const VectR2&);
    void SetPointsDof2D_tri(const VectR2&);
    
    const Matrix<int>& GetNumQuad2D() const;
    const Matrix<int>& GetNumNodes2D_tri() const;
    const Matrix<int>& GetNumNodes2D_quad() const;
    const Matrix<int>& GetCoordinateNodes2D_quad() const;
    const Matrix<int>& GetCoordinateNodes2D_tri() const;
    
  public:
    int GetNbQuadBoundary(int) const;
    int GetNbNodalBoundary(int) const;
    
    const R2& PointsNodalBoundary(int, int) const;
    const VectR2& PointsQuadratureBoundary(int) const;
    const VectR2& PointsDofBoundary(int) const;
    const R2& PointsQuadratureBoundary(int k, int) const;
    const Real_wp& WeightsQuadratureBoundary(int k, int) const;

    const R_N& TangenteLocX(int n) const;
    const R_N& TangenteLocY(int n) const;
    
    void GetLocalCoordOnBoundary(int num_loc, const R2& t_loc, R3& res) const;
    void GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& t_loc) const;

    R3 GetTangentialVector(int num_loc, const R2& vec_u) const;
    R2 TransposeTangentialVector(int num_loc, const R3& vec_u) const;
    R2_Complex_wp TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const;
    
  };
  
} // end namespace

#define MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_HXX
#endif
