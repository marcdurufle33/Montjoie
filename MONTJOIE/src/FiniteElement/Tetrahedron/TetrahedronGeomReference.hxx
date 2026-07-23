#ifndef MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_HXX

namespace Montjoie
{
  
  //! base class for finite tetrahedral elements
  class TetrahedronGeomReference : public ElementGeomReference<Dimension3>
  {
  public :
    enum {LOBATTO_BASIS, REGULAR_BASIS};
    int type_interpolation_nodal; //!< regular or Hesthaven (Lobatto)
    
  protected :
    int nb_points_nodal_inside;
    TriangleGeomReference function_tri;
    
    Matrix<Real_wp> coefFi_curve; //!< coefficients to compute quickly Fi
    Matrix<int> EdgesNodal;
    
    //! coefficients used to evaluate Legendre polynomials
    Matrix<Real_wp> LegendrePolynom;
    //! coefficients used to evaluate jacobi polynomials \f$ P_m^{2i+1, 0} \f$
    Vector<Matrix<Real_wp> > OddJacobiPolynom;
    //! coefficients used to evaluate jacobi polynomials \f$ P_m^{2i+2, 0} \f$
    Vector<Matrix<Real_wp> > EvenJacobiPolynom;
    Array3D<int> NumOrtho3D;
    VectReal_wp InvWeightPolynomial, CoefLegendre;
    Matrix<Real_wp> CoefOddJacobi, CoefEvenJacobi;
    //! Inverse of the VDM system
    Matrix<Real_wp> InverseVDM;
    
  public :
    TetrahedronGeomReference();

    R3 GetCenterReferenceElement() const;
    
    int GetNbPointsNodalInside() const;
    
    const ElementGeomReference<Dimension2>& GetSurfaceFiniteElement(int n) const;
    const TriangleGeomReference& GetTriangularSurfaceFiniteElement() const;
    const TriangleGeomReference& GetQuadrangularSurfaceFiniteElement() const;
    
    FiniteElementProjector* GetNewNodalInterpolation() const;
    
    const Array3D<int>& GetNumOrtho3D() const;
    const VectReal_wp& GetInvWeightPolynomial() const;
    const VectReal_wp& GetCoefLegendre() const;
    const Matrix<Real_wp>& GetCoefOddJacobi() const;
    const Matrix<Real_wp>& GetCoefEvenJacobi() const;
    const Matrix<Real_wp>& GetLegendrePolynomial() const;
    const Vector<Matrix<Real_wp> >& GetEvenJacobiPolynomial() const;
    const Vector<Matrix<Real_wp> >& GetOddJacobiPolynomial() const;
    
    /****************************
     * Initialization functions *
     ****************************/

    size_t GetMemorySize() const;
    void ConstructFiniteElement(int rgeom);

    void ConstructNodalShapeFunctions(int r);

  public:
    static void ConstructRegularPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
				       const Matrix<int>& NumNodes_tri,
				       VectR3& points_nodal3d_, const Array3D<int>& NumNodes);
    
    static void ConstructLobattoPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
				       VectR3& points_nodal3d_);

    void ComputeOrthogonalFunctions(int r);
  protected:
    void ComputeLagrangianFunctions(int r);

  public:
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
	    const R3& point, R3& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	     const R3& point, Matrix3_3& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjLinear(const VectR3& s, const R3& point, R3& res) const;
    void DFjLinear(const VectR3& s, const R3& point, Matrix3_3& res) const;

    Real_wp GetMinimalSize(const VectR3& s) const;
    
    bool OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc, int n) const;
    int ProjectPointOnBoundary(R3& pointloc) const;
    
    void ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const;

  protected :    

    void FjCurve(const SetPoints<Dimension3>& PTReel,
		 const R3& pointloc, R3& res) const;
    
    void DFjCurve(const SetPoints<Dimension3>& PTReel,
		  const R3& pointloc, Matrix3_3& res) const;


    
    /**********************
     * FjElem and DFjElem *
     **********************/
    
  public :
    void FjElem(const VectR3& s, SetPoints<Dimension3>& res,
		const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
		     const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		      SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			   SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
  protected:
    void FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    
    void DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
        
    void FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad, const Volume&) const;
    
    void FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			       const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			   SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
				SetMatrices<Dimension3>& res,
				const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			 SetMatrices<Dimension3>& res,
			 const Mesh<Dimension3>& mesh, int nquad) const;

    
    /*******************
     * Other functions *
     *******************/
    
  public :
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;
    
    void ComputeValuesPhiNodalRef(const R3& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& ) const;
    
    void ComputeValuesPhiOrthoRef(int r, const Array3D<int>&, const VectReal_wp&,
                                  const R3& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiOrthoRef(int r, const Array3D<int>&, const VectReal_wp&,
                                    const R3& pointloc, VectR3& ) const;
    
    void GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const;
    void GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& point_loc) const;
    
    template<int t>
    friend ostream& operator <<(ostream& out, const TetrahedronGeomReference& e);
    
  };


} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_HXX
#endif

