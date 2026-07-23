#ifndef MONTJOIE_FILE_PYRAMID_REFERENCE_SPLIT_HXX

namespace Montjoie
{
  
  //! base class for finite pyramidal elements
  /*!
    Case where we consider the pyramid as two tetrahedra
    This class should not be used
   */
  class PyramidReferenceSplit : public VolumeReference
  {
  public :
    typedef Dimension3 Dimension; //!< dimension3
    
    enum {LOBATTO_ELEMENT, REGULAR_ELEMENT};
    int type_interpolation; //!< Hesthaven (Lobatto)
    enum {WIENERS, LIU, WALKER};
    int type_basis;
    
  protected :
    //! nodal basis functions and derivatives
    Vector<MultivariatePolynomial<Real_wp> > LagrangePr;
    Vector<MultivariatePolynomial<Real_wp> > DxLagrangePr, DyLagrangePr, DzLagrangePr;
    Vector<Vector<Real_wp> > other_coef_lagrange;
    Vector<TinyVector<int, 4> > power_other_coef;
    
    int nb_dof_quad, nb_dof_tri;
    //! triangular basis functions
    TriangleReference function_tri;
    Matrix<int> EdgesNodal;
    //! coefficients to compute quickly F_i
    Matrix<Real_wp> coefFi_curve;
    Matrix<Real_wp> coefFi, coefDFi_dx, coefDFi_dy, coefDFi_dz,
      coefDFi_dx_curve, coefDFi_dy_curve, coefDFi_dz_curve;
      
  public :
    PyramidReferenceSplit();
    
    int GetNbDofBoundary(int) const;
    int GetNbQuadBoundary(int) const;
    int GetNbNodalBoundary(int) const;
    //int GetNbPointsNodal(TinyVector<Real_wp, 3>& x) const;
    
    const R2& PointsNodalBoundary(int, int) const;
    const R2& PointsQuadratureBoundary(int k, int) const;
    const Real_wp& WeightsQuadratureBoundary(int k, int) const;

    int GetNumNodes2D(int, int, int) const;
    
    static int GetNbVertices();
    static int GetNbBoundaries();
    
    
    /****************************
     * Initialization functions *
     ****************************/
    
    void ConstructFiniteElement(int r);
    void ConstructNumberMap(NumberMap& nmap);
    
    void ConstructQuadrature(int order_);
    void ConstructNodalPoints();
    void ConstructRegularNodalPoints();
    void ConstructLobattoPoints();
    void ComputeLagrangianFunctions();
    void ComputeCoefficientTransformation();
    

    /****************
     * Fj transform *
     ****************/
    
    
    void Fj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	    const R3& point, R3& res, const Mesh<Dimension3>& mesh, int nquad) const;

    void DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	     const R3& point, Matrix3_3& res, const Mesh<Dimension3>& mesh, int nquad) const;

    void FjLinear(const VectR3& s, const R3& point, R3& res) const;
    void DFjLinear(const VectR3& s, const R3& point, Matrix3_3& res) const;
    
    void FjCurve(const SetPoints<Dimension3>& PTReel,
		 const R3& pointloc, R3& res) const;
    void DFjCurve(const SetPoints<Dimension3>& PTReel,
		  const R3& pointloc, Matrix3_3& res) const;
    
    bool OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc) const;
    void ProjectPointOnBoundary(R3& pointloc) const;
    
    
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
		 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		      SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			   SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const;
    void DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjSurfaceElem(const VectR3& s, SetPoints<Dimension3>& res,
		       const Mesh<Dimension3>& mesh, int nquad, int ne) const;
    void DFjSurfaceElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
			int nquad, int ne) const;
    
    void FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    
    void DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
        
    void FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad) const;
    
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
    
    void ComputeValuesPhiFirstOrder(const R3& point, VectReal_wp & phi) const;
    void ComputeGradientPhiFirstOrder(const R3& point, VectR3& dphi) const;
    
    void ComputeValuesPhiNodalRef(const R3& pointloc, VectReal_wp& ) const;    
    void ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& ) const;
        
    void GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const;
    
    friend ostream& operator <<(ostream& out, const PyramidReferenceSplit& e);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_REFERENCE_SPLIT_HXX
#endif

