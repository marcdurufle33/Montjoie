#ifndef MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_HXX

#include "FiniteElement/Quadrangle/QuadrangleGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for finite elements on triangles
  class TriangleGeomReference : public ElementGeomReference<Dimension2>
  {    
    friend class WedgeGeomReference;
    
  public :
    //! regular interpolation or Lobatto-like interpolation
    static int type_interpolation_nodal;
    
    // type of interpolation on triangles/tetrahedras
    enum {LOBATTO_BASIS, REGULAR_BASIS, MASS_LUMPED_BASIS,
	  QUASI_LUMPED_BASIS, FEKETE_BASIS};
    
  protected :
    //! orthogonal polynomials (used for shape functions)
    Matrix<Real_wp> LegendrePolynom;
    VectReal_wp InvWeightPolynomial, CoefLegendre;
    Matrix<int> NumOrtho2D;

    Vector<Matrix<Real_wp> > OddJacobiPolynom;
    Matrix<Real_wp> CoefOddJacobi;
    
    //! Coefficients for curved triangles
    Matrix<Real_wp> coefFi_curve;
    
    //! nodal basis functions (stored as a VDM)
    Matrix<Real_wp> InverseVDM;
    
  public :
    
    TriangleGeomReference();

    size_t GetMemorySize() const;    
    int GetNbPointsNodalInside() const;
    R2 GetCenterReferenceElement() const;
    FiniteElementProjector* GetNewNodalInterpolation() const;

    const Matrix<Real_wp>& GetLegendrePolynomial() const;
    const VectReal_wp& GetInverseWeightPolynomial() const;
    const VectReal_wp& GetCoefficientLegendre() const;
    const Matrix<int>& GetNumOrtho2D() const;
    const Vector<Matrix<Real_wp> >& GetOddJacobiPolynomial() const;
    const Matrix<Real_wp>& GetCoefficientOddJacobi() const;
    const Matrix<Real_wp>& GetInverseVandermonde() const;
    
    /****************************
     * Initialization functions *
     ****************************/
        
    void ConstructFiniteElement(int rgeom);
   
    void SetNodalPoints(int r, const VectR2& nodal_points);
    
    void ComputeOrthogonalFunctions(int r);
    void ConstructNodalShapeFunctions(int r);
    
    static void ConstructRegularPoints(int r, VectReal_wp& points_nodal1d_,
				       VectR2& points_nodal2d_, const Matrix<int>& NumNodes2D_);
    
    static void ConstructLobattoPoints(int r, int type, VectReal_wp& points_nodal1d_,
                                       VectR2& points_nodal2d_);
    
  protected:
    void ComputeLagrangianFunctions(int r);

  public:
    void ComputeCurvedTransformation();
    void ComputeCoefficientTransformation();

    
    /****************
     * Fj transform *
     ****************/
    
  public :
    void Fj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	    const R2& point, R2& res, const Mesh<Dimension2>& mesh,int nquad) const;
    
    void DFj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
	     const R2& point, Matrix2_2& res,
	     const Mesh<Dimension2>& mesh, int nquad) const;
    
    Real_wp GetMinimalSize(const VectR2& s) const;
    
    bool OutsideReferenceElement(const VectR2&, const R2& Xn, const Real_wp& epsilon) const;
    Real_wp GetDistanceToBoundary(const R2& pointloc) const;
    Real_wp GetDistanceToBoundary(const R2& pointloc, int n) const;
    int ProjectPointOnBoundary(R2& pointloc) const;
    void ComputeCoefJacobian(const VectR2& s, VectReal_wp& CoefJacobian) const;
    
    void FjLinear(const VectR2& s, const R2& point, R2& res) const;
    void DFjLinear(const VectR2& s, const R2& point, Matrix2_2& res) const;

  protected:
    void FjCurve(const SetPoints<Dimension2>& PTReel, const R2& pointloc, R2& res) const;
    void DFjCurve(const SetPoints<Dimension2>& PTReel,
		  const R2& pointloc, Matrix2_2& res) const;
    void DFjLinear(const VectR2& s, Matrix2_2& res) const;
    
    
    /**********************
     * FjElem and DFjElem *
     **********************/
    
  public:
    void FjElem(const VectR2& s, SetPoints<Dimension2>& res,
		const Mesh<Dimension2>& mesh, int nquad) const;
    
    void FjElemNodal(const VectR2& s, SetPoints<Dimension2>& res,
		     const Mesh<Dimension2>& mesh, int nquad) const;
    
#ifdef MONTJOIE_WITH_THREE_DIM
    void FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
		     const Mesh<Dimension3>& mesh, int nquad, const Face<Dimension3>&) const;

    void FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad, const Face<Dimension3>&) const;
    
    void DFjElemNodal(const VectR3& pts_nodal, bool is_curved,
		      VectR3& dF_dx, VectR3& dF_dy) const;

    void DFjElemQuadrature(const VectR3& pts_nodal, bool is_curved,
			   const VectR3& dF_dx_node, const VectR3& dF_dy_node,
			   VectR3& dF_dx, VectR3& dF_dy) const;
#endif
    
    void FjElemQuadrature(const VectR2& s, SetPoints<Dimension2>& res,
			  const Mesh<Dimension2>& mesh, int nquad) const;
    
    void FjElemDof(const VectR2& s, SetPoints<Dimension2>& res,
		   const Mesh<Dimension2>& mesh, int nquad) const;
    
    void DFjElem(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		 SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
		 int nquad) const;
    
    void DFjElemNodal(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		      SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
		      int nquad) const;
    
    void DFjElemQuadrature(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			   SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			   int nquad) const;
    
    void DFjElemDof(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		    SetMatrices<Dimension2>& res,
		    const Mesh<Dimension2>& mesh, int nquad) const;
    
    template<class DimensionB, class TypeEntity>
    void FjElemNodalCurve(const typename DimensionB::VectR_N& s, SetPoints<DimensionB>& res,
			  const Mesh<DimensionB>& mesh, int nquad, const TypeEntity& elt) const;
    
  protected:
    void FjElemQuadratureCurve(const VectR2& s, SetPoints<Dimension2>& res,
			       const Mesh<Dimension2>& mesh, int nquad) const;

    void FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			       const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemDofCurve(const VectR2& s, SetPoints<Dimension2>& res,
			const Mesh<Dimension2>& mesh, int nquad) const;
    
    template<class DimensionB>
    void FjElemNodalLinear(const typename DimensionB::VectR_N& s,
                           SetPoints<DimensionB>& res) const;
    
    void FjElemQuadratureLinear(const VectR2& s, SetPoints<Dimension2>& res) const;
    void FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    
    void FjElemDofLinear(const VectR2& s, SetPoints<Dimension2>& res) const;
    
    void DFjElemNodalCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			   SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			   int nquad) const;
    
    void DFjElemQuadratureCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
				SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
				int nquad) const;
    
    void DFjElemDofCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			 SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			 int nquad) const;
    
    void DFjElemNodalLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
    void DFjElemQuadratureLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
    void DFjElemDofLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
        
  public :
    void FjSurfaceElemLinear(const VectR2& s, SetPoints<Dimension2>& res, int ne) const;
    
    
    /*******************
     * Other functions *
     *******************/
    
  public :
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;
    
    void ComputeValuesPhiOrthoRef(int r, const Matrix<int>& NumOrtho,
                                  const VectReal_wp& InvWeigh,
                                  const R2& pointloc, VectReal_wp& ) const;  

    void ComputeValuesPhiOrthoRef(const R2& pointloc, VectReal_wp& ) const;  

    void ComputeGradientPhiOrthoRef(int r, const Matrix<int>& NumOrtho,
                                    const VectReal_wp& InvWeigh,
                                    const R2& pointloc, VectR2& ) const;

    void ComputeGradientPhiOrthoRef(const R2& pointloc, VectR2& ) const;

    void ComputeValuesPhiNodalRef(const R2& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiNodalRef(const R2& pointloc, VectR2& ) const;

    void GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const;    
    void GetLocalCoordOnBoundary(int num_loc, const R2& res, Real_wp& t_loc) const;    

    friend ostream& operator <<(ostream& out, const TriangleGeomReference& e);
    
  };
    
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_HXX
#endif
