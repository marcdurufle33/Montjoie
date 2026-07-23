#ifndef MONTJOIE_FILE_FACE_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  inline ElementGeomReference<Dimension2>::ElementGeomReference()
    : ElementGeomReference_Base<Dimension2>()
  {
    axisym_geometry = false;
  }
  
  
  /**********************
   * Convenient methods *
   **********************/
  
  
  //! returns the number of nodal points in the face
  inline int ElementGeomReference<Dimension2>::GetNbPointsNodalElt() const
  {
    return this->points_nodal_nd.GetM();
  }
  
  
  //! returns the number of nodal points along an edge
  inline int ElementGeomReference<Dimension2>::GetNbNodalBoundary(int) const
  {
    return this->points_nodal1d.GetM();
  }
  
  
  //! node numbers
  inline const Matrix<int>& ElementGeomReference<Dimension2>::GetNumNodes2D() const
  {
    return NumNodes2D;
  }


  //! node numbers
  inline const Matrix<int>& ElementGeomReference<Dimension2>::GetCoordinateNodes2D() const
  {
    return CoordinateNodes;
  }
  

  //! node number
  inline int ElementGeomReference<Dimension2>::GetNumNodes2D(int i, int j) const
  {
    return NumNodes2D(i, j);
  }

  
  //! nodal point on the boundary
  /*!
    \param[in] j local node number on the edge
    \param[in] ne edge number
   */
  inline const Real_wp& ElementGeomReference<Dimension2>::PointsNodalBoundary(int j, int ne) const
  {
    return this->points_nodal1d(j);
  }


  //! returns number of the node k of the edge j
  /*!
    \param[in] k local node number on the edge
    \param[in] j edge number
   */
  inline int ElementGeomReference<Dimension2>::GetNodalNumber(int j, int k) const
  {
    return EdgesNode(k, j);
  }
  
  
  //! 2-D Weight of integration
  inline const VectReal_wp& ElementGeomReference<Dimension2>::WeightsND() const
  {
    return weights2d;
  }


  //! 2-D integration points
  inline const VectR2& ElementGeomReference<Dimension2>::PointsQuadND() const
  {
    return points2d;
  }


  //! 2-D integration points
  inline const VectR2& ElementGeomReference<Dimension2>::PointsND() const
  {
    return points2d;
  }
  
  
  //! 2-D Weight of integration
  inline const Real_wp& ElementGeomReference<Dimension2>::WeightsND(int j) const
  {
    return weights2d(j);
  }
  
  
  //! 2-D Point of integration
  inline const R2& ElementGeomReference<Dimension2>::PointsND(int j) const
  {
    return points2d(j);
  }
  

  //! returns quadrature points on the boundary ne
  /*!
    \param[in] ne edge number
   */
  inline const VectReal_wp& ElementGeomReference<Dimension2>::PointsDofBoundary(int ne) const
  {
    return this->points_dof1d;
  }


  //! returns quadrature points on the boundary ne
  /*!
    \param[in] ne edge number
   */
  inline const VectReal_wp& ElementGeomReference<Dimension2>::PointsQuadratureBoundary(int ne) const
  {
    return this->points1d;
  }

  
  //! quadrature point on the boundary
  /*!
    \param[in] k local qudrature point number on the edge
    \param[in] ne edge number
   */
  inline const Real_wp& ElementGeomReference<Dimension2>::PointsQuadratureBoundary(int k, int ne) const
  {
    return this->points1d(k);
  }
  
  
  //! weight on the boundary
  /*!
    \param[in] k local qudrature weight number on the edge
    \param[in] ne edge number
   */
  inline const Real_wp& ElementGeomReference<Dimension2>::WeightsQuadratureBoundary(int k, int ne) const
  {
    return this->weights1d(k);
  }
  
  
  //! 2-D dof points
  inline const Vector<R2>& ElementGeomReference<Dimension2>::PointsDofND() const
  {
    return points_dof2d;
  }


  //! 2-D dof point
  inline const R2& ElementGeomReference<Dimension2>::PointsDofND(int j) const
  {
    return points_dof2d(j);
  }

  
  //! returns the number of quadrature points on an edge
  inline int ElementGeomReference<Dimension2>::GetNbQuadBoundary(int) const
  {
    return this->points1d.GetM();
  }
  

  //! return the number of dof points
  inline int ElementGeomReference<Dimension2>::GetNbPointsDof() const
  {
    return points_dof2d.GetM();
  }
  
  
  //! returns the number of quadrature points (interior + boundary)
  inline int ElementGeomReference<Dimension2>::GetNbPointsQuadrature() const
  {
    return points2d.GetM();
  }
  
  
  inline void ElementGeomReference<Dimension2>::SetAxisymmetricGeometry()
  {
    axisym_geometry = true;
  }


  inline const ElementGeomReference<Dimension1>& 
  ElementGeomReference<Dimension2>::GetSurfaceFiniteElement(int n) const
  {
    return edge_geom;
  }
  
  
  /**************************************
   * Container for ElementGeomReference *
   **************************************/
  

  inline ElementGeomReferenceContainer<Dimension2>
  ::ElementGeomReferenceContainer(ElementGeomReference<Dimension2>& elt)
    : ElementGeomReferenceContainer_Base<Dimension2>(elt)
  {
  }
  
  inline VectR2 ElementGeomReferenceContainer<Dimension2>::PointsQuadInsideND(int N) const
  {
    return this->elt_geom.PointsQuadInsideND(N);
  }
    
  inline int ElementGeomReferenceContainer<Dimension2>::GetNbPointsNodalElt() const
  {
    return this->elt_geom.GetNbPointsNodalElt();
  }
    
  inline int ElementGeomReferenceContainer<Dimension2>::GetNbNodalBoundary(int n) const
  {
    return this->elt_geom.GetNbNodalBoundary(n);
  }
  
  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension2>::GetNumNodes2D() const
  {
    return this->elt_geom.GetNumNodes2D();
  }
  
  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension2>::GetCoordinateNodes2D() const
  {
    return this->elt_geom.GetCoordinateNodes2D();
  }
  
  inline int ElementGeomReferenceContainer<Dimension2>::GetNumNodes2D(int i, int j) const
  {
    return this->elt_geom.GetNumNodes2D(i, j);
  }
    
  inline const Vector<R2>& ElementGeomReferenceContainer<Dimension2>::PointsNodalND() const
  {
    return this->elt_geom.PointsNodalND();
  }
  
  inline const R2& ElementGeomReferenceContainer<Dimension2>::PointsNodalND(int j) const
  {
    return this->elt_geom.PointsNodalND(j);
  }
  
  inline const Real_wp& ElementGeomReferenceContainer<Dimension2>::PointsNodalBoundary(int j, int k) const
  {
    return this->elt_geom.PointsNodalBoundary(j, k);
  }
  
  inline int ElementGeomReferenceContainer<Dimension2>::GetNodalNumber(int j, int k) const
  {
    return this->elt_geom.GetNodalNumber(j, k);
  }
  
  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension2>::GetNodalNumber() const
  {
    return this->elt_geom.EdgesNode;
  }
  
  inline const Globatto<Real_wp>& ElementGeomReferenceContainer<Dimension2>
  ::GetNodalShapeFunctions1D() const
  {
    return this->elt_geom.lob_geom;
  }
  
  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension2>::WeightsND() const
  {
    return this->elt_geom.WeightsND();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension2>::PointsQuadND() const
  {
    return this->elt_geom.PointsQuadND();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension2>::PointsND() const
  {
    return this->elt_geom.PointsND();
  }
  
  inline const Real_wp& ElementGeomReferenceContainer<Dimension2>::WeightsND(int j) const
  {
    return this->elt_geom.WeightsND(j);
  }
  
  inline const R2& ElementGeomReferenceContainer<Dimension2>::PointsND(int j) const
  {
    return this->elt_geom.PointsND(j);
  }

  inline void ElementGeomReferenceContainer<Dimension2>::SetPointsND(const VectR2& pts)
  {
    this->elt_geom.points2d = pts;
  }
  
  inline void ElementGeomReferenceContainer<Dimension2>::SetWeightsND(const VectReal_wp& w)
  {
    this->elt_geom.weights2d = w;
  }

  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension2>
  ::PointsQuadratureBoundary(int num_loc) const
  {
    return this->elt_geom.PointsQuadratureBoundary(num_loc);
  }

  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension2>
  ::PointsDofBoundary(int num_loc) const
  {
    return this->elt_geom.PointsDofBoundary(num_loc);
  }
    
  inline const Real_wp& ElementGeomReferenceContainer<Dimension2>
  ::PointsQuadratureBoundary(int k, int num_loc) const
  {
    return this->elt_geom.PointsQuadratureBoundary(k, num_loc);
  }
  
  inline const Real_wp& ElementGeomReferenceContainer<Dimension2>
  ::WeightsQuadratureBoundary(int k, int num_loc) const
  {
    return this->elt_geom.WeightsQuadratureBoundary(k, num_loc);
  }
    
  inline const Vector<R2>& ElementGeomReferenceContainer<Dimension2>::PointsDofND() const
  {
    return this->elt_geom.PointsDofND();
  }
  
  inline const R2& ElementGeomReferenceContainer<Dimension2>::PointsDofND(int j) const
  {
    return this->elt_geom.PointsDofND(j);
  }

  inline void ElementGeomReferenceContainer<Dimension2>::SetPointsDofND(const VectR2& pts)
  {
    this->elt_geom.points_dof2d = pts;
  }
    
  inline int ElementGeomReferenceContainer<Dimension2>::GetNbPointsDof() const
  {
    return this->elt_geom.GetNbPointsDof();
  }
  
  inline int ElementGeomReferenceContainer<Dimension2>::GetNbPointsQuadrature() const
  {
    return this->elt_geom.GetNbPointsQuadrature();
  }
  
  inline int ElementGeomReferenceContainer<Dimension2>::GetNbQuadBoundary(int n) const
  {
    return this->elt_geom.GetNbQuadBoundary(n);
  }
    
  inline void ElementGeomReferenceContainer<Dimension2>
  ::GetLocalCoordOnBoundary(int num_loc, const Real_wp& t_loc, R2& res) const
  {
    return this->elt_geom.GetLocalCoordOnBoundary(num_loc, t_loc, res);
  }
  
  inline void ElementGeomReferenceContainer<Dimension2>
  ::GetLocalCoordOnBoundary(int num_loc, const R2& res, Real_wp& t_loc) const
  {
    return this->elt_geom.GetLocalCoordOnBoundary(num_loc, res, t_loc);
  }
  
}
  
#define MONTJOIE_FILE_FACE_GEOM_REFERENCE_INLINE_CXX
#endif
