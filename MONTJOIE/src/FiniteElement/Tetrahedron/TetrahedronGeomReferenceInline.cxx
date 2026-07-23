#ifndef MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  //! return the center of the unit tetrahedron
  inline R3 TetrahedronGeomReference::GetCenterReferenceElement() const
  {
    return R3(0.25, 0.25, 0.25);
  }


  //! returns the number of nodal points strictly inside the element
  inline int TetrahedronGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - (2*this->order_geom*this->order_geom+2);
  }
  
  
  //! returns the restriction of the finite element on a triangular face
  inline const ElementGeomReference<Dimension2>&
  TetrahedronGeomReference::GetSurfaceFiniteElement(int n) const
  {
    return function_tri;
  }
  
     
  //! returns the restriction of the finite element on a triangular face
  inline const TriangleGeomReference&
  TetrahedronGeomReference::GetTriangularSurfaceFiniteElement() const
  {
    return function_tri;
  }


  //! returns the restriction of the finite element on a triangular face
  inline const TriangleGeomReference&
  TetrahedronGeomReference::GetQuadrangularSurfaceFiniteElement() const
  {
    return function_tri;
  }
    

  inline FiniteElementProjector* TetrahedronGeomReference::GetNewNodalInterpolation() const
  {
    return new DenseProjector<Dimension3>();
  }


  inline const Array3D<int>& TetrahedronGeomReference::GetNumOrtho3D() const
  {
    return NumOrtho3D;
  }
  
  
  inline const VectReal_wp& TetrahedronGeomReference::GetInvWeightPolynomial() const
  {
    return InvWeightPolynomial;
  }
  
  
  inline const VectReal_wp& TetrahedronGeomReference::GetCoefLegendre() const
  {
    return CoefLegendre;
  }
  
  
  inline const Matrix<Real_wp>& TetrahedronGeomReference::GetCoefOddJacobi() const
  {
    return CoefOddJacobi;
  }
  
  
  inline const Matrix<Real_wp>& TetrahedronGeomReference::GetCoefEvenJacobi() const
  {
    return CoefEvenJacobi;
  }
  
  
  inline const Matrix<Real_wp>& TetrahedronGeomReference::GetLegendrePolynomial() const
  {
    return LegendrePolynom;
  }
  
  
  inline const Vector<Matrix<Real_wp> >& TetrahedronGeomReference::GetEvenJacobiPolynomial() const
  {
    return EvenJacobiPolynom;
  }
  
  
  inline const Vector<Matrix<Real_wp> >& TetrahedronGeomReference::GetOddJacobiPolynomial() const
  {
    return OddJacobiPolynom;
  }
  
  
  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the tetrahedron
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit tetrahedron
    \param[out] res the result of transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void TetrahedronGeomReference::
  Fj(const VectR3& s, const SetPoints<Dimension3>& PTReel, const R3& point,
     R3& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel,point,res);
    else
      FjLinear(s,point,res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the tetrahedron
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit tetrahedron
    \param[out] res the jacobian matrix
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void TetrahedronGeomReference::
  DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel, const R3& point,
      Matrix3_3& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel, point, res);
    else
      DFjLinear(s, point, res);
  }
  
  
  //! transformation Fi in the case of straight tetrahedron
  inline void TetrahedronGeomReference::
  FjLinear(const VectR3& s, const R3& point, R3& res) const
  {
    // Fi = (1-x-y-z)*S0 + x*S1 + y*S2 + z*S3 
    Real_wp x = point(0), y = point(1), z = point(2),
      l1 = Real_wp(1)-x-y-z;
    
    for (int j = 0; j < 3; j++)
      res(j) = l1*s(0)(j) + x*s(1)(j) + y*s(2)(j) + z*s(3)(j);
  }
  
  
  //! transformation DFi in the case of straight tetrahedron
  inline void TetrahedronGeomReference::
  DFjLinear(const VectR3& s, const R3& pt, Matrix3_3& res) const
  {
    for (int j = 0; j < 3; j++)
      {
	// dFi/dx = S1 - S0
	// dFi/dy = S2 - S0
	// dFi/dz = S3 - S0
	res(j,0) = s(1)(j)-s(0)(j);
	res(j,1) = s(2)(j)-s(0)(j);
	res(j,2) = s(3)(j)-s(0)(j);
      }
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_GEOM_REFERENCE_INLINE_CXX
#endif
