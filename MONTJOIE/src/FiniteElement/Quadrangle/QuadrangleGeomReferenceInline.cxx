#ifndef MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  //! returns the number of points strictly inside the element
  inline int QuadrangleGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - 4*this->order_geom;
  }


  //! return the center of the unit square
  inline R2 QuadrangleGeomReference::GetCenterReferenceElement() const
  {
    return R2(0.5, 0.5);
  }


  inline FiniteElementProjector* QuadrangleGeomReference::GetNewNodalInterpolation() const
  {
    return new TensorizedProjector<Dimension2>();
  }

  
  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the quad
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit square
    \param[out] res the result of transformation Fi
    \param[in] mesh quadrilateral mesh considered
    \param[in] nquad quadrilateral element number
   */  
  inline void QuadrangleGeomReference
  ::Fj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
       const R2& point, R2& res, const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel,point,res);
    else
      FjLinear(s,point,res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the quad
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit square
    \param[out] res the jacobian matrix
    \param[in] mesh quadrilateral mesh considered
    \param[in] nquad quadrilateral element number
   */  
  inline void QuadrangleGeomReference
  ::DFj(const VectR2& s, const SetPoints<Dimension2>& PTReel,
        const R2& point, Matrix2_2& res,
        const Mesh<Dimension2>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel,point,res);
    else
      DFjLinear(s,point,res);
  }
  
  
  //! not implemented
  inline void QuadrangleGeomReference
  ::ComputeCoefJacobian(const VectR2& s, VectReal_wp& CoefJacobian) const
  {
    abort();
  }

  
  //! evaluation of 1-D nodal basis functions at a given point x
  inline void QuadrangleGeomReference
  ::ComputeValuesNodalPhi1D(const Real_wp& x, VectReal_wp& phi) const
  {
    lob_geom.ComputeValuesPhiRef(x, phi);
  }
    

  inline void QuadrangleGeomReference
  ::ComputeNodalGradientRef(const Vector<Real_wp>& Uloc, VectR2& gradU) const
  {
    ComputeNodalGradientRefT(Uloc, gradU);
  }
  
  
  inline void QuadrangleGeomReference
  ::ComputeNodalGradientRef(const Vector<Complex_wp>& Uloc, Vector<R2_Complex_wp>& gradU) const
  {
    ComputeNodalGradientRefT(Uloc, gradU);
  }
  
}
  
#define MONTJOIE_FILE_QUADRANGLE_GEOM_REFERENCE_INLINE_CXX
#endif
