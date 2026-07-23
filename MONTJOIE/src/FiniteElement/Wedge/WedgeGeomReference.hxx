#ifndef MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_HXX

namespace Montjoie
{
  
  //! base class for finite prismatic elements
  class WedgeGeomReference : public ElementGeomReference<Dimension3>
  {
  protected :
    //! triangular shape functions
    TriangleGeomReference function_geom_tri;
    QuadrangleGeomReference function_quad;
    Matrix<int> EdgesNodal;
    Matrix<Real_wp> LegendrePolynom;
    Vector<Real_wp> CoefLegendre;
    //! coefficients to compute quickly F_i
    Matrix<Real_wp> coefFi_curve;
    Matrix<Real_wp> coefFi, coefDFi_dx, coefDFi_dy, coefDFi_dz;
    Matrix<int> NumNodesTri;
  
  public :
    WedgeGeomReference();

    R3 GetCenterReferenceElement() const;
    
    int GetNbPointsNodalInside() const;

    const ElementGeomReference<Dimension2>& GetSurfaceFiniteElement(int n) const;
    const TriangleGeomReference& GetTriangularSurfaceFiniteElement() const;
    const QuadrangleGeomReference& GetQuadrangularSurfaceFiniteElement() const;

    FiniteElementProjector* GetNewNodalInterpolation() const;
    const Matrix<int>& GetNumNodesTri() const;

    const VectReal_wp& GetCoefLegendre() const;
    const Matrix<Real_wp>& GetLegendrePolynomial() const;
    
    /****************************
     * Initialization functions *
     ****************************/

    size_t GetMemorySize() const;
    void ConstructFiniteElement(int rgeom);
    
    void ConstructOrthogonalFunctions(int r);
    void ConstructNodalShapeFunctions(int r);

  public:
    static void ConstructRegularPoints(int r, VectReal_wp& points1d_, VectR2& points2d_tri_,
                                       const Matrix<int>& NumNodes2D_tri_,
                                       VectR2& points2d_quad_, const Matrix<int>& NumNodes2D_quad_,
                                       VectR3& points3d_, const Array3D<int>& NumNodes3D_);
    
    static void ConstructLobattoPoints(int r, VectReal_wp& points1d_, VectR2& points2d_tri_, const Matrix<int>& NumNodes2D_tri_,
				       VectR2& points2d_quad_, const Matrix<int>& NumNodes2D_quad_,
				       VectR3& points3d_, const Array3D<int>& NumNodes3D_);
    
    void ComputeCurvedTransformation();
    void ComputeCoefficientTransformation();

  protected:
    void GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
				      const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
				      R3& d_nx, R3& d_ny, R3& d_nz) const;    

  public:
    R3 GetTangentialVector(int num_loc, const R2& vec_u) const;
    R2 TransposeTangentialVector(int num_loc, const R3& vec_u) const;
    R2_Complex_wp TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const;

    
    /****************
     * Fj transform *
     ****************/
    
  public : 
    void Fj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	    const R3& point, R3& res,const Mesh<Dimension3>& mesh, int nquad) const;

    void DFj(const VectR3& s,const SetPoints<Dimension3>& PTReel,
	     const R3& point, Matrix3_3& res,
	     const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjLinear(const VectR3& s, const R3& point, R3& res) const;
    void DFjLinear(const VectR3& s, const R3& point, Matrix3_3& res) const;

  protected :
    
    void FjCurve(const SetPoints<Dimension3>& PTReel,
		 const R3& pointloc, R3& res) const;
    void DFjCurve(const SetPoints<Dimension3>& PTReel,
		  const R3& pointloc, Matrix3_3& res) const;

  public :
    Real_wp GetMinimalSize(const VectR3& s) const;
    
    bool OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc, int n) const;
    int ProjectPointOnBoundary(R3& pointloc) const;
    
    void ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const;
    
    
    /**********************
     * FjElem and DFjElem *
     **********************/
    
    
    void FjElem(const VectR3& s, SetPoints<Dimension3>& res,
		const Mesh<Dimension3>& mesh, int nquad) const;
    void FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
		     const Mesh<Dimension3>& mesh, int nquad) const;
    void FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad) const;
    void FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		 SetMatrices<Dimension3>& res,
		 const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		      SetMatrices<Dimension3>& res,
		      const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			   SetMatrices<Dimension3>& res,
			   const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res,
		    const Mesh<Dimension3>& mesh, int nquad) const;
    
  protected :
    void FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    
    void DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
        
    void FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad, const Volume& ) const;
    void FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			       const Mesh<Dimension3>& mesh, int nquad) const;
    void FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			   SetMatrices<Dimension3>& res, 
			   const Mesh<Dimension3>& mesh,int nquad) const;
    void DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
				SetMatrices<Dimension3>& res,
				const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			 SetMatrices<Dimension3>& res,
			 const Mesh<Dimension3>& mesh, int nquad) const;
    
    
    /*******************
     * Other functions *
     *******************/

    void ComputeValuesPhiFirstOrder(const R3& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiFirstOrder(const R3& pointloc, VectR3& ) const;
    
  public:
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;
    
    void ComputeValuesPhiOrthoRef(int r, const Array3D<int>&, const VectReal_wp&,
                                  const R3& pointloc, VectReal_wp& ) const;    
    
    void ComputeGradientPhiOrthoRef(int r, const Array3D<int>&, const VectReal_wp&,
                                    const R3& pointloc, VectR3& ) const;
    
    void ComputeValuesPhiNodalRef(const R3& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& ) const;
    
    void GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const;
    void GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& point_loc) const;
    
    template<int t>
    friend ostream& operator <<(ostream& out, const WedgeGeomReference& e);
    
  };

} // end namespace

#define MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_HXX
#endif
