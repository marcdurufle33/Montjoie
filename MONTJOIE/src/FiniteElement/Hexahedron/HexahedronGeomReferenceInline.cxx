#ifndef MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{

  //! return the center of the unit cube
  inline R3 HexahedronGeomReference::GetCenterReferenceElement() const
  {
    return R3(0.5, 0.5, 0.5);
  }


  //! returns the number of nodal points strictly inside the element
  inline int HexahedronGeomReference::GetNbPointsNodalInside() const
  {
    return this->points_nodal_nd.GetM() - (6*this->order_geom*this->order_geom+2);
  }


  //! returns the finite element that is the restriction of the current finite element to a face
  inline const ElementGeomReference<Dimension2>&
  HexahedronGeomReference::GetSurfaceFiniteElement(int n) const
  {
    return function_geom_quad;
  }
  

  //! returns the finite element that is the restriction of the current finite element to a face  
  inline const QuadrangleGeomReference&
  HexahedronGeomReference::GetTriangularSurfaceFiniteElement() const
  {
    return function_geom_quad;
  }
  

  //! returns the finite element that is the restriction of the current finite element to a face  
  inline const QuadrangleGeomReference&
  HexahedronGeomReference::GetQuadrangularSurfaceFiniteElement() const
  {
    return function_geom_quad;
  }

  
  inline FiniteElementProjector* HexahedronGeomReference::GetNewNodalInterpolation() const
  {
    return new TensorizedProjector<Dimension3>();
  }
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the hexahedron
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit cube
    \param[out] res the result of transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void HexahedronGeomReference::
  Fj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
     const R3& point, R3& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel, point, res);
    else
      FjLinear(s, point, res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the hexahedron
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit cube
    \param[out] res the jacobian matrix
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void HexahedronGeomReference::
  DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
      const R3& point, Matrix3_3& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel, point, res);
    else
      DFjLinear(s, point, res);
  }
    
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_INLINE_CXX
#endif
  


