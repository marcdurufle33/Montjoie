#ifndef MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_HXX

namespace Montjoie
{
  
  //! Base class for quadrilateral finite element
  class QuadrangleGeomReference : public ElementGeomReference<Dimension2>
  {
    friend class MeshBoundaries<Dimension3>;
    
  protected :
    IVect SymmetricNodes;
    
    //! coefficients in order to compute F_i and DF_i quickly
    Matrix<Real_wp> coefFi;
    //! in the curved case too
    Vector<VectReal_wp> coefFi_curve;
    Vector<IVect> numberFi_curve; //!< in the curved case too
    
  public :
    QuadrangleGeomReference();
    
    size_t GetMemorySize() const;
    int GetNbPointsNodalInside() const;
    R2 GetCenterReferenceElement() const;    
    FiniteElementProjector* GetNewNodalInterpolation() const;
    
    /****************************
     * Initialization functions *
     ****************************/
    
    // sets map in order to number dofs of the mesh
    void ConstructFiniteElement(int rgeom);
    
    void SetNodalPoints(int r, const Vector<R2>& pts);
    void ConstructNodalShapeFunctions(int r);
    
    void ComputeCurvedTransformation();
    void ComputeCoefficientTransformation();
    
    
    /****************
     * Fj transform *
     ****************/
    
  public :
    // res = Fj(point)
    void Fj(const VectR2& s, const SetPoints<Dimension2>& PTReel, const R2& point,
	    R2& res, const Mesh<Dimension2>& mesh, int nquad) const;
    
    // res = DFj(point)
    void DFj(const VectR2& s, const SetPoints<Dimension2>& PTReel, const R2& point,
	     Matrix2_2& res, const Mesh<Dimension2>& mesh, int nquad) const;
    
    Real_wp GetMinimalSize(const VectR2& s) const;

    // returns true if Xn is outside the unit cube
    bool OutsideReferenceElement(const VectR2& s, const R2& Xn, const Real_wp& epsilon) const;
    // returns distance of Xn to the boundaries of the unit cube
    Real_wp GetDistanceToBoundary(const R2& pointloc) const;
    Real_wp GetDistanceToBoundary(const R2& pointloc, int n) const;
    // projects pointloc on boundaries of the unit cube, if pointloc is outside
    int ProjectPointOnBoundary(R2& pointloc) const;
    void ComputeCoefJacobian(const VectR2& s, VectReal_wp& CoefJacobian) const;
    
    void FjLinear(const VectR2& s, const R2& point, R2& res) const;
    void DFjLinear(const VectR2& s, const R2& point, Matrix2_2& res) const;
    
  protected:
    void FjCurve(const SetPoints<Dimension2>& PTReel, const R2& pointloc, R2& res) const;
    
    void DFjCurve(const SetPoints<Dimension2>& PTReel,
		  const R2& pointloc, Matrix2_2& res) const;
    
    
    /**********************
     * FjElem and DFjElem *
     **********************/
    
  public :
    void FjElem(const VectR2& s, SetPoints<Dimension2>& res,
		const Mesh<Dimension2>& mesh, int nquad) const;
    
    void FjElemNodal(const VectR2& s, SetPoints<Dimension2>& res,
		     const Mesh<Dimension2>& mesh,int nquad) const;
    
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
			  const Mesh<Dimension2>& mesh,int nquad) const;
    
    void FjElemDof(const VectR2& s, SetPoints<Dimension2>& res,
		   const Mesh<Dimension2>& mesh, int nquad) const;
    
    void DFjElem(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		 SetMatrices<Dimension2>& res, const Mesh<Dimension2>&, int nquad) const;
    
    void DFjElemNodal(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		      SetMatrices<Dimension2>&, const Mesh<Dimension2>&, int) const;
    
    void DFjElemQuadrature(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			   SetMatrices<Dimension2>&, const Mesh<Dimension2>&, int) const;
    
    void DFjElemDof(const VectR2& s, const SetPoints<Dimension2>& PTReel,
		    SetMatrices<Dimension2>&, const Mesh<Dimension2>& mesh, int) const;
    
  protected:
    template<class DimensionB, class TypeEntity>
    void FjElemNodalCurve(const typename DimensionB::VectR_N& s, SetPoints<DimensionB>& res,
			  const Mesh<DimensionB>& mesh, int nquad, const TypeEntity& elt) const;
    
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
				SetMatrices<Dimension2>& res, const Mesh<Dimension2>&,
				int nquad) const;
    void DFjElemDofCurve(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			 SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			 int nquad) const;
    
    void DFjElemNodalLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
    void DFjElemQuadratureLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
    void DFjElemDofLinear(const VectR2& s, SetMatrices<Dimension2>& res) const;
    
    void FjSurfaceElemLinear(const VectR2& s, SetPoints<Dimension2>& res, int ne) const;
    
    
    /*******************
     * Other functions *
     *******************/

    template <class T0>
    void ComputeNodalGradientRefT(const Vector<T0>&, Vector<TinyVector<T0, 2> >&) const;
    
  public:
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;
    
    void ComputeValuesPhiNodalRef(const R2& pointloc, VectReal_wp&) const;    
    void ComputeGradientPhiNodalRef(const R2& pointloc, VectR2&) const;    
    
    void ComputeNodalGradientRef(const Vector<Real_wp>& Uloc, VectR2& gradU) const;
    void ComputeNodalGradientRef(const Vector<Complex_wp>& Uloc, Vector<R2_Complex_wp>& gradU) const;
  
    // computes 2-D point res located at the local coordinate t_loc on the edge num_loc
    void GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const;    
    void GetLocalCoordOnBoundary(int num_loc, const R2& res, Real_wp& t_loc) const;    

    template<int t>
    friend ostream& operator <<(ostream& out,const QuadrangleReference<t>& e);
    
  };

}

#define MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_HXX
#endif
