#ifndef MONTJOIE_FILE_FACE_GEOM_REFERENCE_HXX

namespace Montjoie
{
  //! Base class for 2-D finite element (shape functions only)
  template<>
  class ElementGeomReference<Dimension2> : public ElementGeomReference_Base<Dimension2>
  {
    friend class ElementGeomReferenceContainer<Dimension2>;
    friend class ElementReference_Dim<Dimension2>;
    
  public :
    typedef Dimension2 Dimension; //!< dimension
    
  protected :    
    //! 1-D nodal shape functions
    Globatto<Real_wp> lob_geom;    
    //! 2-D Weights of integration
    VectReal_wp weights2d;
    //! 2-D Points of integration
    VectR2 points2d;
        
    //! 2-D dof points
    VectR2 points_dof2d;
    
    //! numbers of the nodal points on the edges of the element
    Matrix<int> EdgesNode;

    //! NumNodes2D(i,j) is the number of the dof located at \f$ (\xi_i, \xi_j) \f$
    Matrix<int> NumNodes2D;
    //!  i = CoordinatesNodes(num_dof,0) is linked with \f$ \xi_i \f$
    //!  j = CoordinatesNodes(num_dof,1) is linked with \f$ \xi_j \f$
    Matrix<int> CoordinateNodes;
    
    //! if true, the geometry is considered as axisymmetric (for computing the curvatures)
    bool axisym_geometry;

    EdgeLobatto edge_geom;
    
  public :

    ElementGeomReference();

    
    /**********************
     * Convenient methods *
     **********************/
    
    size_t GetMemorySize() const;
    VectR2 PointsQuadInsideND(int N) const;    
    int GetNbPointsNodalElt() const;
    
    int GetNbNodalBoundary(int) const;
        
    const Matrix<int>& GetNumNodes2D() const;
    const Matrix<int>& GetCoordinateNodes2D() const;
    int GetNumNodes2D(int, int) const;
    
    const Real_wp& PointsNodalBoundary(int j, int k) const;
    
    int GetNodalNumber(int j, int k) const;
    
    const VectReal_wp& WeightsND() const;
    const VectR2& PointsQuadND() const;
    const VectR2& PointsND() const;
    const Real_wp& WeightsND(int j) const;
    const R2& PointsND(int j) const;

    const VectReal_wp& PointsDofBoundary(int num_loc) const;
    const VectReal_wp& PointsQuadratureBoundary(int num_loc) const;
    const Real_wp& PointsQuadratureBoundary(int k, int num_loc) const;
    const Real_wp& WeightsQuadratureBoundary(int k, int num_loc) const;
    
    const Vector<R2>& PointsDofND() const;
    const R2& PointsDofND(int j) const;
    
    int GetNbPointsDof() const;
    int GetNbPointsQuadrature() const;
    int GetNbQuadBoundary(int) const;
    void SetAxisymmetricGeometry();
    
    const ElementGeomReference<Dimension1>& GetSurfaceFiniteElement(int n) const;
    
    
    /**********************
     * Overloaded methods *
     **********************/
    

#ifdef MONTJOIE_WITH_THREE_DIM
    virtual void FjElemNodal(const VectR2& s, SetPoints<Dimension2>& res,
			     const Mesh<Dimension2>& mesh, int nquad) const = 0;    

    virtual void FjElemQuadrature(const VectR2& s, SetPoints<Dimension2>& res,
				  const Mesh<Dimension2>& mesh, int nquad) const = 0;
    
    virtual void DFjElemNodal(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			      SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			      int nquad) const = 0;
    
    virtual void DFjElemQuadrature(const VectR2& s, const SetPoints<Dimension2>& PTReel,
				   SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
				   int nquad) const = 0;

    virtual void FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
			     const Mesh<Dimension3>& mesh, int nquad, const Face<Dimension3>&) const = 0;

    virtual void FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
				  const Mesh<Dimension3>& mesh, int nquad, const Face<Dimension3>&) const = 0;
    
    virtual void DFjElemNodal(const VectR3& pts_nodal, bool is_curved,
			      VectR3& dF_dx, VectR3& dF_dy) const = 0;

    virtual void DFjElemQuadrature(const VectR3& pts_nodal, bool is_curved,
				   const VectR3& dF_dx_node, const VectR3& dF_dy_node,
				   VectR3& dF_dx, VectR3& dF_dy) const = 0;
#endif
    
    virtual void GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const = 0;    
    virtual void GetLocalCoordOnBoundary(int num_loc, const R2& res, Real_wp& t_loc) const = 0;    
    
    
    /*******************************************************
     * computation of geometric quantities on the boundary *
     *******************************************************/

    static void GetCurvatureNodal(const Globatto<Real_wp>& lob,
                                  const VectR2& normale, const VectReal_wp& ds,
				  VectReal_wp& curvature);
    
  protected:        
    // computes curvature
    void GetCurvatureOnSurface(const VectR2& s, const SetPoints<Dimension2>& PTReel,
				SetMatrices<Dimension2>& res,
				const Mesh<Dimension2>& mesh, int num_elem,
				int num_edge, int num_loc_edge) const;
    
    // computes the curvature, when no curved element is used
    void GetCurvatureOnLinearSurface(const VectR2& s, const SetPoints<Dimension2>&,
				     SetMatrices<Dimension2>& res,
				     const Mesh<Dimension2>& mesh, int num_elem,
				     int num_edge, int num_loc_edge) const;
    
    // computes the curvature of edge, when curved element is used
    void GetCurvatureOnCurvedSurface(const VectR2& s, const SetPoints<Dimension2>&,
				     SetMatrices<Dimension2>& res,
				     const Mesh<Dimension2>& mesh, int num_elem,
				     int num_edge, int num_loc_edge) const;
    
    // for axisymmetric code, we compute the two curvatures k1 and k2
    void GetCurvatureOnAxisymmetricSurface(const VectR2& s,
					   const SetPoints<Dimension2>& PTReel,
					   SetMatrices<Dimension2>& res,
					   const Mesh<Dimension2>& mesh, int num_elem,
					   int ne, int num_loc) const;
    
    // in the linear case
    void GetCurvatureOnAxisymmetricLinearSurface(const VectR2& s,
						 const SetPoints<Dimension2>& PTReel,
						 SetMatrices<Dimension2>& res) const;
    
    /* methods used to compute F_i^{-1} */
  public:
    // returns the bounding box of an element defined by the vertices s
    static void GetBoundingBox(const VectR2& s, const Real_wp& coef,
			       VectR2& box, TinyVector<R2, 2>& enveloppe);
    
    // intersection of the edge [pointA, pointB] with the polygon defined by box
    static bool ComputeIntersection(int nb_vertices, const VectR2& box, const R2& pointA,
                                    const R2& pointB, VectR2& res);
    
    
    /*****************
     * Other methods *
     *****************/

    
    static bool IsAffineTransformation(const VectR2& s);        
        
    friend ostream& operator <<(ostream& out, const ElementGeomReference<Dimension2>& e);
    
  };

  ostream& operator <<(ostream& out, const ElementGeomReference<Dimension2>& e);
  
  template<>
  class ElementGeomReferenceContainer<Dimension2>
    : public ElementGeomReferenceContainer_Base<Dimension2>
  {
  public :
    ElementGeomReferenceContainer(ElementGeomReference<Dimension2>&);

    VectR2 PointsQuadInsideND(int N) const;    
    int GetNbPointsNodalElt() const;
    
    int GetNbNodalBoundary(int) const;
        
    const Matrix<int>& GetNumNodes2D() const;
    const Matrix<int>& GetCoordinateNodes2D() const;
    int GetNumNodes2D(int, int) const;
    
    const Vector<R2>& PointsNodalND() const;
    const R2& PointsNodalND(int j) const;
    const Real_wp& PointsNodalBoundary(int j, int k) const;
    
    int GetNodalNumber(int j, int k) const;
    const Matrix<int>& GetNodalNumber() const;
    const Globatto<Real_wp>& GetNodalShapeFunctions1D() const;

    const VectReal_wp& WeightsND() const;
    virtual const VectR2& PointsQuadND() const;
    const VectR2& PointsND() const;
    const Real_wp& WeightsND(int j) const;
    const R2& PointsND(int j) const;
  protected:
    void SetPointsND(const VectR2&);
    void SetWeightsND(const VectReal_wp&);

  public:
    const VectReal_wp& PointsDofBoundary(int num_loc) const;
    const VectReal_wp& PointsQuadratureBoundary(int num_loc) const;
    const Real_wp& PointsQuadratureBoundary(int k, int num_loc) const;
    const Real_wp& WeightsQuadratureBoundary(int k, int num_loc) const;
    
    const Vector<R2>& PointsDofND() const;
    const R2& PointsDofND(int j) const;
  protected:
    void SetPointsDofND(const VectR2&);

  public:
    int GetNbPointsDof() const;
    int GetNbPointsQuadrature() const;
    int GetNbQuadBoundary(int) const;
    
    void GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const;    
    void GetLocalCoordOnBoundary(int num_loc, const R2& res, Real_wp& t_loc) const;    

  };

} // end namespace


#define MONTJOIE_FILE_FACE_GEOM_REFERENCE_HXX
#endif
