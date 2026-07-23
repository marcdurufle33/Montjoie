#ifndef MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_HXX

namespace Montjoie
{
  
  //! base class for hexahedral finite element
  class HexahedronGeomReference : public ElementGeomReference<Dimension3>
  {
  protected :
    //! restriction to a face
    QuadrangleGeomReference function_geom_quad;
    //! coefficients in order to compute F_i and DF_i quickly
    Matrix<Real_wp> coefFi;
    //! in the curved case too
    Vector<VectReal_wp> coefFi_curve;
    Vector<IVect> numberFi_curve; //!< in the curved case too
        
  public :
    HexahedronGeomReference();

    R3 GetCenterReferenceElement() const;
    
    size_t GetMemorySize() const;
    int GetNbPointsNodalInside() const;

    const ElementGeomReference<Dimension2>& GetSurfaceFiniteElement(int n) const;
    const QuadrangleGeomReference& GetTriangularSurfaceFiniteElement() const;
    const QuadrangleGeomReference& GetQuadrangularSurfaceFiniteElement() const;
    FiniteElementProjector* GetNewNodalInterpolation() const;
    
    
    /****************************
     * Initialization functions *
     ****************************/
    
    
    void ConstructFiniteElement(int rgeom);
    
    void ConstructNodalShapeFunctions(int r);
    
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


  public :
    /****************
     * Fj transform *
     ****************/
    
    void Fj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	    const R3& point, R3& res, const Mesh<Dimension3>& mesh, int nquad) const;
    
    void DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	     const R3& point, Matrix3_3& res,
	     const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjLinear(const VectR3& s, const R3& point, R3& res) const;
    void DFjLinear(const VectR3& s, const R3& point, Matrix3_3& res) const;

  protected:
    void FjCurve(const SetPoints<Dimension3>& PTReel,
		 const R3& pointloc, R3& res) const;
    void DFjCurve(const SetPoints<Dimension3>& PTReel,
		  const R3& pointloc, Matrix3_3& res) const;

  public :
    Real_wp GetMinimalSize(const VectR3& s) const;
    
    bool OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc, int n) const;
    Real_wp GetDistanceToBoundary(const R3& pointloc) const;
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
    
  protected:
    // in the linear case
    void FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
    void FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const;
        
    // DF_i in the linear case
    void DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    void DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const;
    
    // in the curvilinear case
    void FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
			  const Mesh<Dimension3>& mesh, int nquad, const Volume& ) const;
    
    void FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			       const Mesh<Dimension3>& mesh, int nquad) const;
    
    void FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const;
    
    // DF_i in the curvilinear case
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

    template <class T0>
    void ComputeNodalGradientRefT(const Vector<T0>&, Vector<TinyVector<T0, 3> >&) const;
    
  public :
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;
    
    void ComputeValuesPhiNodalRef(const R3& pointloc, VectReal_wp&) const;    
    void ComputeGradientPhiNodalRef(const R3& pointloc, VectR3&) const;    
    
    void GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const;
    void GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& point_loc) const;
    
    void ComputeNodalGradientRef(const Vector<Real_wp>& Uloc, VectR3& gradU) const;
    void ComputeNodalGradientRef(const Vector<Complex_wp>& Uloc, Vector<R3_Complex_wp>& gradU) const;

    friend ostream& operator <<(ostream& out, const HexahedronGeomReference& e);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_HXX
#endif

