#ifndef MONTJOIE_FILE_PYRAMID_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  /****************************
   * Initialization functions *
   ****************************/
  
  
  //! return the center of the symmetric pyramid
  inline R3 PyramidGeomReference::GetCenterReferenceElement() const
  {
    return R3(0.0, 0.0, 0.25);
  }

  
  //! returns the number of nodal points strictly inside the element
  inline int PyramidGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - (3*this->order_geom*this->order_geom+2);
  }

  
  //! returns the triangular finite element that is the restriction of 
  //! current finite element on a triangular face
  inline const ElementGeomReference<Dimension2>& PyramidGeomReference::GetSurfaceFiniteElement(int n) const
  {
    if (this->IsLocalFaceQuadrilateral(n))
      return function_quad;
    
    return function_tri;
  }
  

  //! returns the triangular finite element that is the restriction of 
  //! current finite element on a triangular face    
  inline const TriangleGeomReference& PyramidGeomReference
  ::GetTriangularSurfaceFiniteElement() const
  {
    return function_tri;
  }


  //! returns the triangular finite element that is the restriction of 
  //! current finite element on a triangular face    
  inline TriangleGeomReference& PyramidGeomReference::GetTriangularSurfaceFiniteElement()
  {
    return function_tri;
  }
  
  
  //! returns the finite element that is the restriction of 
  //! current finite element on a quadrangular face
  inline const QuadrangleGeomReference& PyramidGeomReference
  ::GetQuadrangularSurfaceFiniteElement() const
  {
    return function_quad;
  }


  inline FiniteElementProjector* PyramidGeomReference::GetNewNodalInterpolation() const
  {
    return new DenseProjector<Dimension3>();
  }


  inline const VectReal_wp& PyramidGeomReference::GetCoefLegendre() const
  {
    return CoefLegendre;
  }
  
  
  inline const Matrix<Real_wp>& PyramidGeomReference::GetLegendrePolynomial() const
  {
    return LegendrePolynom;
  }
  
  
  inline const Matrix<Real_wp>& PyramidGeomReference::GetCoefEvenJacobi() const
  {
    return CoefEvenJacobi;
  }
  
  
  inline const Vector<Matrix<Real_wp> >& PyramidGeomReference::GetEvenJacobiPolynomial() const
  {
    return EvenJacobiPolynom;
  }
  
    
  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the pyramid
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit pyramid
    \param[out] res the result of transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
  */
  inline void PyramidGeomReference::
  Fj(const VectR3& s,const SetPoints<Dimension3>& PTReel,
     const R3& point, R3& res,
     const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel, point, res);
    else
      FjLinear(s, point, res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the pyramid
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit pyramid
    \param[out] res the jacobian matrix
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void PyramidGeomReference::
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

#define MONTJOIE_FILE_PYRAMID_GEOM_REFERENCE_INLINE_CXX
#endif
