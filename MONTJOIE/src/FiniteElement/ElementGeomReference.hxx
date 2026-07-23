#ifndef MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_HXX

namespace Montjoie
{

  template<class Dimension>
  class ElementGeomReferenceContainer {};

  template<class Dimension>
  class ElementGeomReferenceContainer_Base;
    
  //! Generic class for finite element (only shape functions)
  template<class Dimension>
  class ElementGeomReference_Base
  {
    template<int t> friend class FaceReference;
#ifdef MONTJOIE_WITH_THREE_DIM
    template<int t> friend class VolumeReference;
#endif
    friend class ElementGeomReferenceContainer_Base<Dimension>;
    
  public :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::VectR_N VectR_N;

    //! equality between quadrature and nodal points ?
    bool quadrature_equal_nodal;
    bool dof_equal_nodal; //!< equality between dof and nodal points ?
    bool dof_equal_quadrature; //!< equality between dof and quadrature points ?
    
  protected :

    //! order of approximation for geometry
    int order_geom;
    
    //! 1-D nodal points
    VectReal_wp points_nodal1d;

   //! 1-D Weights of integration
    VectReal_wp weights1d;
    //! 1-D Points of integration
    VectReal_wp points1d;
    
    //! 1-D dof points
    VectReal_wp points_dof1d;
    
    //! matrix used for interpolation from mesh to finite element
    Matrix<Real_wp> CoefCurve1D;

    //! normale to each edge (or face in 3-D)
    VectR_N normale;

    //! 2-D or 3-D nodal points
    VectR_N points_nodal_nd;
    VectReal_wp weights_nodal_nd;
    
    //! Values of gradient of shape functions on nodal points
    Matrix<R_N> GradientPhi_Nodal; 

    //! hybrid type of the element (0, 1, 2, 3)
    int hybrid_type_elt;
    
    //! number of vertices, edges, boundaries for this element
    int nb_vertices_elt, nb_edges_elt, nb_boundaries_elt;
    
  public :
    
    ElementGeomReference_Base();
    virtual ~ElementGeomReference_Base();

    /*****************
     * Basic methods *
     *****************/
    
    
    int GetGeometryOrder() const;    
    
    const R_N& NormaleLoc(int n) const;
    
    const VectReal_wp& PointsNodal1D() const;
    const Real_wp& PointsNodal1D(int i) const;
    
    const Vector<R_N>& PointsNodalND() const;
    const R_N& PointsNodalND(int j) const;

    const Vector<Real_wp>& WeightsNodalND() const;
    const Real_wp& WeightsNodalND(int j) const;

    const VectReal_wp& PointsDof1D() const;
    const Real_wp& PointsDof1D(int i) const;
    
    const VectReal_wp& Points1D() const;
    const Real_wp& Points1D(int i) const;
    const VectReal_wp& Weights1D() const;
    const Real_wp& Weights1D(int i) const;
    
    int GetHybridType() const;
    int GetNbVertices() const;
    int GetNbEdges() const;
    int GetNbBoundaries() const;
    void SetAxisymmetricGeometry();

    virtual const ElementGeomReference<typename Dimension::DimensionBoundary>&
    GetSurfaceFiniteElement(int n) const = 0;
    
    
    /******************
     * Common methods *
     ******************/
    
    
    template<class DimensionB>
    void SetMesh(Mesh<DimensionB>& mesh);

    // computes normale and surface element ds
    void GetNormale(const MatrixN_N& dfjm1, R_N& normale_fj,
                    Real_wp& dsj, int num_edge) const;
    
    void FjCurvePhi(const SetPoints<Dimension>& PTReel,
		    const R_N& point, R_N& res, VectReal_wp& phi) const;
    
    void DFjCurvePhi(const SetMatrices<Dimension>&, const R_N& point, MatrixN_N& res, 
		     const VectReal_wp& phi) const;
    
    /*******************
     * Virtual methods *
     *******************/

    virtual size_t GetMemorySize() const = 0;    
    virtual R_N GetCenterReferenceElement() const = 0;
    virtual void ConstructFiniteElement(int rgeom) = 0;
    virtual FiniteElementProjector* GetNewNodalInterpolation() const = 0;
    
    virtual void Fj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		    const R_N& point, R_N& res, const Mesh<Dimension>& mesh,int nquad) const = 0;
    
    virtual void FjLinear(const VectR_N& s, const R_N& point, R_N& res) const = 0;
    virtual void DFjLinear(const VectR_N& s, const R_N& point, MatrixN_N& res) const = 0;

    virtual void DFj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		     const R_N& point, MatrixN_N& res,
		     const Mesh<Dimension>& mesh, int nquad) const = 0;
    
    virtual Real_wp GetMinimalSize(const VectR_N& s) const = 0;
    
    virtual bool OutsideReferenceElement(const VectR_N& s, const R_N& Xn, const Real_wp& epsilon) const = 0;
    
    virtual Real_wp GetDistanceToBoundary(const R_N& pointloc) const = 0;
    virtual Real_wp GetDistanceToBoundary(const R_N& pointloc, int n) const = 0;
    
    virtual int ProjectPointOnBoundary(R_N& pointloc) const = 0;
    
    virtual void ComputeCoefJacobian(const VectR_N& s, VectReal_wp& CoefJacobian) const = 0;

    virtual void FjElem(const VectR_N& s, SetPoints<Dimension>& res,
			const Mesh<Dimension>& mesh, int nquad) const = 0;
    
    virtual void FjElemNodal(const VectR_N& s, SetPoints<Dimension>& res,
			     const Mesh<Dimension>& mesh, int nquad) const = 0;    
    
    virtual void FjElemQuadrature(const VectR_N& s, SetPoints<Dimension>& res,
				  const Mesh<Dimension>& mesh, int nquad) const = 0;
    
    virtual void FjElemDof(const VectR_N& s, SetPoints<Dimension>& res,
			   const Mesh<Dimension>& mesh, int nquad) const = 0;
    
    virtual void DFjElem(const VectR_N& s, const SetPoints<Dimension>& PTReel,
			 SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
			 int nquad) const = 0;
    
    virtual void DFjElemNodal(const VectR_N& s, const SetPoints<Dimension>& PTReel,
			      SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
			      int nquad) const = 0;
    
    virtual void DFjElemQuadrature(const VectR_N& s, const SetPoints<Dimension>& PTReel,
				   SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
				   int nquad) const = 0;
    
    virtual void DFjElemDof(const VectR_N& s, const SetPoints<Dimension>& PTReel,
			    SetMatrices<Dimension>& res,
			    const Mesh<Dimension>& mesh, int nquad) const = 0;
    
    virtual void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const = 0;

    virtual void ComputeValuesPhiNodalRef(const R_N& pointloc, VectReal_wp&) const = 0;    
    virtual void ComputeGradientPhiNodalRef(const R_N& pointloc, VectR_N& ) const = 0;

    virtual void ComputeCoefficientTransformation();
    
    /*****************
     * Other methods *
     *****************/
    
    virtual void ComputeNodalGradientRef(const Vector<Real_wp>& Un, VectR_N& Unode) const;
    virtual void ComputeNodalGradientRef(const Vector<Complex_wp>& Un,
					 Vector<TinyVector<Complex_wp, Dimension::dim_N> >& Unode) const;
    
    // not virtual since it will not be overloaded
    template<class T1>
    void ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
			      const Vector<T1>& Uloc, Vector<TinyVector<T1, Dimension::dim_N> >& grad_Uloc) const;
    
    template<class T>
    void ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
			      const Vector<Vector<T> >& U, Vector<Vector<T> >& gradU) const;
    
  };


  template<class Dimension>
  class ElementGeomReferenceContainer_Base
  {
  protected : 
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::VectR_N VectR_N;

    ElementGeomReference<Dimension>& elt_geom;
    
  public :    
    ElementGeomReferenceContainer_Base(ElementGeomReference<Dimension>& elt);
    
    bool QuadratureEqualNodal() const;
    bool DofEqualNodal() const;
    bool DofEqualQuadrature() const;
    
    int GetGeometryOrder() const;    
    void SetGeometryOrder(int r);

    ElementGeomReference<Dimension>& GetGeometricElement();
    const ElementGeomReference<Dimension>& GetGeometricElement() const;
    
    const R_N& NormaleLoc(int n) const;
    
    const VectReal_wp& PointsNodal1D() const;
    const Real_wp& PointsNodal1D(int i) const;
    
    const VectReal_wp& PointsDof1D() const;
    const Real_wp& PointsDof1D(int i) const;
  protected:
    void SetPointsDof1D(const VectReal_wp&);

  public:
    const VectReal_wp& Points1D() const;
    const Real_wp& Points1D(int i) const;
    const VectReal_wp& Weights1D() const;
    const Real_wp& Weights1D(int i) const;

  protected:
    void SetPoints1D(const VectReal_wp&);
    void SetWeights1D(const VectReal_wp&);

  public:
    int GetHybridType() const;
    int GetNbVertices() const;
    int GetNbEdges() const;
    int GetNbBoundaries() const;

    FiniteElementProjector* GetNewNodalInterpolation() const;
    void SetMesh(Mesh<Dimension>& mesh);
    
    void GetNormale(const MatrixN_N& dfjm1, R_N& normale_fj,
		    Real_wp& dsj, int num_edge) const;

    void Fj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
	    const R_N& point, R_N& res, const Mesh<Dimension>& mesh,int nquad) const;
    
    void FjLinear(const VectR_N& s, const R_N& point, R_N& res) const ;
    void DFjLinear(const VectR_N& s, const R_N& point, MatrixN_N& res) const ;
    
    void DFj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
	     const R_N& point, MatrixN_N& res,
	     const Mesh<Dimension>& mesh, int nquad) const;
    
    Real_wp GetMinimalSize(const VectR_N& s) const;
    
    bool OutsideReferenceElement(const VectR_N& s, const R_N& Xn, const Real_wp& epsilon) const;
    
    Real_wp GetDistanceToBoundary(const R_N& pointloc) const;
    
    int ProjectPointOnBoundary(R_N& pointloc) const;
    
    void ComputeCoefJacobian(const VectR_N& s, VectReal_wp& CoefJacobian) const;
    
    void FjElem(const VectR_N& s, SetPoints<Dimension>& res,
		const Mesh<Dimension>& mesh, int nquad) const;
    
    void FjElemNodal(const VectR_N& s, SetPoints<Dimension>& res,
		     const Mesh<Dimension>& mesh, int nquad) const;    
    
    void FjElemQuadrature(const VectR_N& s, SetPoints<Dimension>& res,
			  const Mesh<Dimension>& mesh, int nquad) const;
    
    void FjElemDof(const VectR_N& s, SetPoints<Dimension>& res,
		   const Mesh<Dimension>& mesh, int nquad) const;
    
    void DFjElem(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		 SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
		 int nquad) const;
    
    void DFjElemNodal(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		      SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
		      int nquad) const;
    
    void DFjElemQuadrature(const VectR_N& s, const SetPoints<Dimension>& PTReel,
			   SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
			   int nquad) const;
    
    void DFjElemDof(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		    SetMatrices<Dimension>& res,
		    const Mesh<Dimension>& mesh, int nquad) const;
    
    void ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const;

    void ComputeValuesPhiNodalRef(const R_N& pointloc, VectReal_wp&) const;    
    void ComputeGradientPhiNodalRef(const R_N& pointloc, VectR_N& ) const;
    
    void ComputeNodalGradientRef(const Vector<Real_wp>& Un, VectR_N& Unode) const;
    void ComputeNodalGradientRef(const Vector<Complex_wp>& Un,
				 Vector<TinyVector<Complex_wp, Dimension::dim_N> >& Unode) const;

    template<class T>
    void ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
			      const Vector<Vector<T> >& U, Vector<Vector<T> >& gradU) const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_HXX
#endif
