#ifndef MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  //! returns the number of nodal points strictly inside the element
  inline int TriangleGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - 3*this->order_geom;
  }
  
  
  //! return the center of the unit triangle
  inline R2 TriangleGeomReference::GetCenterReferenceElement() const
  {
    return R2(1.0/3, 1.0/3);
  }

  
  //! allocates a new projector
  inline FiniteElementProjector* TriangleGeomReference::GetNewNodalInterpolation() const
  {
    return new DenseProjector<Dimension2>();
  }

  
  inline const Matrix<Real_wp>& TriangleGeomReference::GetLegendrePolynomial() const
  {
    return LegendrePolynom;
  }
  
  
  inline const VectReal_wp& TriangleGeomReference::GetInverseWeightPolynomial() const
  {
    return InvWeightPolynomial;
  }
  
  
  inline const VectReal_wp& TriangleGeomReference::GetCoefficientLegendre() const
  {
    return CoefLegendre;
  }
  
  
  inline const Matrix<int>& TriangleGeomReference::GetNumOrtho2D() const
  {
    return NumOrtho2D;
  }
  
  
  inline const Vector<Matrix<Real_wp> >& TriangleGeomReference::GetOddJacobiPolynomial() const
  {
    return OddJacobiPolynom;
  }
  
  
  inline const Matrix<Real_wp>& TriangleGeomReference::GetCoefficientOddJacobi() const
  {
    return CoefOddJacobi;
  }
  
  
  inline const Matrix<Real_wp>& TriangleGeomReference::GetInverseVandermonde() const
  {
    return InverseVDM;
  }


  /************************
   * Convenient functions *
   ************************/
  
  
  //! evaluating the 1-D nodal functions at x
  inline void TriangleGeomReference
  ::ComputeValuesNodalPhi1D(const Real_wp& x, VectReal_wp& phi) const
  {
    lob_geom.ComputeValuesPhiRef(x, phi);
  }
  

  inline void TriangleGeomReference
  ::ComputeValuesPhiOrthoRef(const R2& pointloc, VectReal_wp& phi) const
  {
    ComputeValuesPhiOrthoRef(this->order_geom, NumOrtho2D, InvWeightPolynomial, pointloc, phi);
  }


  inline void TriangleGeomReference
  ::ComputeGradientPhiOrthoRef(const R2& pointloc, VectR2& phi) const
  {
    ComputeGradientPhiOrthoRef(this->order_geom, NumOrtho2D, InvWeightPolynomial, pointloc, phi);
  }
  

  //! returns 2-D point of a point on a edge
  /*!
    \param[in] num_loc local edge number in the triangle
    \param[in] t_loc local coordinate on the edge
    \param[out] res local coordinates of the point in the unit triangle
  */
  inline void TriangleGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const
  {
    if (num_loc == 0)
      res.Init(t_loc, 0);
    else if (num_loc == 1)
      res.Init(Real_wp(1)-t_loc, t_loc);
    else
      res.Init(0, Real_wp(1)-t_loc);
  }
  
  
  //! computes local abscissa from a 2-D point
  /*!
    \param[in] num_loc local edge number in the triangle
    \param[out] t_loc local coordinate on the edge
    \param[res] pt local coordinates of the point in the unit triangle
  */
  inline void TriangleGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R2& pt, Real_wp& t_loc) const
  {
    if (num_loc == 0)
      t_loc = pt(0);
    else if (num_loc == 1)
      t_loc = pt(1);
    else
      t_loc = Real_wp(1)-pt(1);
  }
  
  
  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the triangle
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit triangle
    \param[out] res the result of transformation Fi
    \param[in] mesh mesh considered
    \param[in] nquad triangular element number
  */
  inline void TriangleGeomReference
  ::Fj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
       const R2& point, R2& res,
       const Mesh<Dimension2>& mesh,int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel, point, res);
    else
      FjLinear(s, point, res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the triangle
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit triangle
    \param[out] res the jacobian matrix
    \param[in] mesh triangular mesh considered
    \param[in] nquad triangular element number
  */
  inline void TriangleGeomReference
  ::DFj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
        const R2& point, Matrix2_2& res,
        const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel,point,res);
    else
      DFjLinear(s,res);
  }
  
  
  //! not implemented
  inline void TriangleGeomReference
  ::ComputeCoefJacobian(const VectR2& s, VectReal_wp& CoefJacobian) const
  {
    abort();
  }

  
  //! transformation Fi in the case of straight triangle
  inline void TriangleGeomReference
  ::FjLinear(const VectR2& s, const R2& point, R2& res) const
  {
    Real_wp x, y, l1;
    x = point(0); y = point(1);
    l1 = Real_wp(1) - x - y;
    res.Init(s(0)(0)*l1 + s(1)(0)*x + s(2)(0)*y,
	     s(0)(1)*l1 + s(1)(1)*x + s(2)(1)*y);
  }
  
  
  //! transformation DFi in the case of straight triangle
  inline void TriangleGeomReference
  ::DFjLinear(const VectR2& s, Matrix2_2& res) const
  {
    res(0,0) = s(1)(0)-s(0)(0);
    res(0,1) = s(2)(0)-s(0)(0);
    res(1,0) = s(1)(1)-s(0)(1);
    res(1,1) = s(2)(1)-s(0)(1);
  }
  
  
  //! computes DF_i(point) for straight elements 
  inline void TriangleGeomReference
  ::DFjLinear(const VectR2& s, const R2& point, Matrix2_2& res) const
  {
    DFjLinear(s, res);
  }
   
}

#define MONTJOIE_FILE_TRIANGLE_GEOM_REFERENCE_INLINE_CXX
#endif
