#ifndef MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
    
  /****************************
   * Initialization functions *
   ****************************/
  
  
  //! return the center of the unit prism
  inline R3 WedgeGeomReference::GetCenterReferenceElement() const
  {
    return R3(1.0/3, 1.0/3, 0.5);
  }


  //! returns the number of nodal points strictly inside the element
  inline int WedgeGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - (4*this->order_geom*this->order_geom+2);
  }


  //! returns the triangular finite element that is the restriction of 
  //! current finite element on a triangular face
  inline const ElementGeomReference<Dimension2>& WedgeGeomReference::GetSurfaceFiniteElement(int n) const
  {
    if (this->IsLocalFaceQuadrilateral(n))
      return function_quad;
    
    return function_geom_tri;
  }

  
  //! returns the restriction of the finite element on a triangular face  
  inline const TriangleGeomReference& WedgeGeomReference::GetTriangularSurfaceFiniteElement() const
  {
    return function_geom_tri;
  }
  
  
  //! returns the restriction of the finite element on a quadrilateral face
  inline const QuadrangleGeomReference&
  WedgeGeomReference::GetQuadrangularSurfaceFiniteElement() const
  {
    return function_quad;
  }


  inline FiniteElementProjector* WedgeGeomReference::GetNewNodalInterpolation() const
  {
    return new DenseProjector<Dimension3>();
  }


  inline const Matrix<int>& WedgeGeomReference::GetNumNodesTri() const
  {
    return NumNodesTri;
  }
  
  
  inline const VectReal_wp& WedgeGeomReference::GetCoefLegendre() const
  {
    return CoefLegendre;
  }
  
  
  inline const Matrix<Real_wp>& WedgeGeomReference::GetLegendrePolynomial() const
  {
    return LegendrePolynom;
  }
  

  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the Wedge
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit Wedge
    \param[out] res the result of transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void WedgeGeomReference::
  Fj(const VectR3& s,const SetPoints<Dimension3>& PTReel,
     const R3& point, R3& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel,point,res);
    else
      FjLinear(s,point,res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the Wedge
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit Wedge
    \param[out] res the jacobian matrix
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void WedgeGeomReference::
  DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
      const R3& point, Matrix3_3& res,
      const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel, point, res);
    else
      DFjLinear(s, point, res);
  }
    
} // end namespace

#define MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_INLINE_CXX
#endif
