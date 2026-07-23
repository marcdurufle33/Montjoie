#ifndef MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  inline ElementGeomReference<Dimension3>::ElementGeomReference()
    : ElementGeomReference_Base<Dimension3>()
  {
  }
  
  
  /************************
   * Convenient functions *
   ************************/
  
  
  //! returns the number of nodal points in the element
  inline int ElementGeomReference<Dimension3>::GetNbPointsNodalElt() const
  {
    return this->points_nodal_nd.GetM();
  }
  
  
  //! return the number of dof points
  inline int ElementGeomReference<Dimension3>::GetNbPointsDof() const
  {
    return points_dof3d.GetM();
  }
  
  
  //! returns the number of quadrature points (inside + boundary)
  inline int ElementGeomReference<Dimension3>::GetNbPointsQuadrature() const
  {
    return points3d.GetM();
  }


  //! returns numbering array for nodal points
  inline const Array3D<int>& ElementGeomReference<Dimension3>::GetNumNodes3D() const
  {
    return NumNodes3D;
  }
  
  
  //! returns node number of the node (i, j, k)
  inline int ElementGeomReference<Dimension3>::GetNumNodes3D(int i, int j, int k) const
  {
    return NumNodes3D(i, j, k);
  }

  
  //! returns node numbers of quadrilateral face
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetNumNodes2D_quad() const
  {
    return NumNodes2D_quad;
  }
  

  //! returns scalar to tensorial matching for 3-D nodes
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetCoordinateNodes3D() const
  {
    return CoordinateNodes;
  }
  
  
  //! returns node number of the node num_dof on face num_face_loc
  inline int ElementGeomReference<Dimension3>::GetNodalNumber(int num_face_loc, int num_dof) const
  {
    return FacesNodal(num_dof,num_face_loc);
  }
  
  
  //! return the array containing node numbers on the faces
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetNodalNumber() const
  {
    return FacesNodal;
  }
  
  
  //! returns 3-D integration weights
  inline const VectReal_wp& ElementGeomReference<Dimension3>::WeightsND() const
  {
    return weights3d;
  }
  

  //! returns 3-D integration points
  inline const VectR3& ElementGeomReference<Dimension3>::PointsQuadND() const
  {
    return points3d;
  }


  //! returns 3-D integration points
  inline const VectR3& ElementGeomReference<Dimension3>::PointsND() const
  {
    return points3d;
  }
  
  
  //! returns 3-D integration weight
  inline const Real_wp& ElementGeomReference<Dimension3>::WeightsND(int j) const
  {
    return weights3d(j);
  }
  

  //! returns 3-D integration point
  inline const R3& ElementGeomReference<Dimension3>::PointsND(int j) const
  {
    return points3d(j);
  }  
  
 
  //! returns 3-D dof points
  inline const Vector<R3>& ElementGeomReference<Dimension3>::PointsDofND() const
  {
    return points_dof3d;
  }
  
  
  //! returns 3-D dof point
  inline const R3& ElementGeomReference<Dimension3>::PointsDofND(int j) const
  {
    return points_dof3d(j);
  }


  inline const VectR2& ElementGeomReference<Dimension3>::Points2D_tri() const
  {
    return points2d_tri;
  }


  inline const VectR2& ElementGeomReference<Dimension3>::Points2D_quad() const
  {
    return points2d_quad;
  }


  inline const VectReal_wp& ElementGeomReference<Dimension3>::Weights2D_tri() const
  {
    return weights2d_tri;
  }


  inline const VectReal_wp& ElementGeomReference<Dimension3>::Weights2D_quad() const
  {
    return weights2d_quad;
  }


  inline const VectR2& ElementGeomReference<Dimension3>::PointsNodal2D_tri() const
  {
    return points_nodal2d_tri;
  }


  inline const VectR2& ElementGeomReference<Dimension3>::PointsNodal2D_quad() const
  {
    return points_nodal2d_quad;
  }


  inline const VectR2& ElementGeomReference<Dimension3>::PointsDof2D_tri() const
  {
    return points_dof2d_tri;
  }


  inline const VectR2& ElementGeomReference<Dimension3>::PointsDof2D_quad() const
  {
    return points_dof2d_quad;
  }
  
  
  inline void ElementGeomReference<Dimension3>::SetPointsNodal2D_tri(const VectR2& pts)
  {
    points_nodal2d_tri = pts;
  }


  inline void ElementGeomReference<Dimension3>::SetPointsNodal2D_quad(const VectR2& pts)
  {
    points_nodal2d_quad = pts;
  }


  inline void ElementGeomReference<Dimension3>::SetPointsDof2D_tri(const VectR2& pts)
  {
    points_dof2d_tri = pts;
  }


  inline void ElementGeomReference<Dimension3>::SetPointsDof2D_quad(const VectR2& pts)
  {
    points_dof2d_quad = pts;
  }


  inline void ElementGeomReference<Dimension3>::SetPoints2D_tri(const VectR2& pts)
  {
    points2d_tri = pts;
  }

  
  inline void ElementGeomReference<Dimension3>::SetPoints2D_quad(const VectR2& pts)
  {
    points2d_quad = pts;
  }


  inline void ElementGeomReference<Dimension3>::SetWeights2D_tri(const VectReal_wp& w)
  {
    weights2d_tri = w;
  }

  
  inline void ElementGeomReference<Dimension3>::SetWeights2D_quad(const VectReal_wp& w)
  {
    weights2d_quad = w;
  }


  inline Matrix<int>& ElementGeomReference<Dimension3>::GetNumQuad2D()
  {
    return NumQuad2D; 
  }
  
  
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetNumQuad2D() const
  {
    return NumQuad2D; 
  }

  
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetNumNodes2D_tri() const
  {
    return NumNodes2D_tri;
  }
  
  
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetCoordinateNodes2D_quad() const
  {
    return CoordinateNodes2D_quad;
  }
  
  
  inline const Matrix<int>& ElementGeomReference<Dimension3>::GetCoordinateNodes2D_tri() const
  {
    return CoordinateNodes2D_tri;
  }


  //! returns the number of quadrature points on a face
  inline int ElementGeomReference<Dimension3>::GetNbQuadBoundary(int n) const
  {
    if (is_local_face_quad(n))
      return this->points2d_quad.GetM();

    return this->points2d_tri.GetM();
  }


  //! returns the number of nodal points on a face
  inline int ElementGeomReference<Dimension3>::GetNbNodalBoundary(int n) const
  {
    if (is_local_face_quad(n))
      return this->points_nodal2d_quad.GetM();
    
    return this->points_nodal2d_tri.GetM();
  }

  
  //! returns true if the face n is a quadrilateral
  inline bool ElementGeomReference<Dimension3>::IsLocalFaceQuadrilateral(int n) const
  {
    return is_local_face_quad(n);
  }
  
    
  //! return a nodal point of a face
  /*!
    \param[in] k local nodal point number
    \param[in] nf face number
   */
  inline const R2& ElementGeomReference<Dimension3>::PointsNodalBoundary(int k, int nf) const
  {
    if (is_local_face_quad(nf))
      return this->points_nodal2d_quad(k);
    
    return this->points_nodal2d_tri(k);
  }


  //! returns quadrature points of a face
  /*!
    \param[in] nf face number
   */
  inline const VectR2& ElementGeomReference<Dimension3>::PointsQuadratureBoundary(int nf) const
  {
    if (is_local_face_quad(nf))      
      return this->points2d_quad;
    
    return this->points2d_tri;
  }


  //! returns quadrature points of a face
  /*!
    \param[in] nf face number
   */
  inline const VectR2& ElementGeomReference<Dimension3>::PointsDofBoundary(int nf) const
  {
    if (is_local_face_quad(nf))      
      return this->points_dof2d_quad;
    
    return this->points_dof2d_tri;
  }

  
  //! returns a quadrature point of a face
  /*!
    \param[in] k local quadrature point number
    \param[in] nf face number
   */
  inline const R2& ElementGeomReference<Dimension3>::PointsQuadratureBoundary(int k, int nf) const
  {
    if (is_local_face_quad(nf))      
      return this->points2d_quad(k);
    
    return this->points2d_tri(k);
  }
  
  
  //! returns a quadrature weight of a face
  /*!
    \param[in] k local quadrature point number
    \param[in] nf face number
   */
  inline const Real_wp& ElementGeomReference<Dimension3>::WeightsQuadratureBoundary(int k, int nf) const
  {
    if (is_local_face_quad(nf))      
      return this->weights2d_quad(k);
    
    return this->weights2d_tri(k);
  }
      

  //! returns node numbering of a quadrilateral face
  inline int ElementGeomReference<Dimension3>::GetNumNodes2D(int i, int j) const
  {
    return this->NumNodes2D_quad(i, j);
  }


  inline const R3& ElementGeomReference<Dimension3>::TangenteLocX(int n) const
  {
    return tangente_loc_x(n);
  }
   
  
  inline const R3& ElementGeomReference<Dimension3>::TangenteLocY(int n) const
  {
    return tangente_loc_y(n);
  }

  
  /**************************************
   * Container for ElementGeomReference *
   **************************************/
  

  inline ElementGeomReferenceContainer<Dimension3>
  ::ElementGeomReferenceContainer(ElementGeomReference<Dimension3>& elt)
    : ElementGeomReferenceContainer_Base<Dimension3>(elt)
  {
  }
  
  inline VectR3 ElementGeomReferenceContainer<Dimension3>::PointsQuadInsideND(int N) const
  {
    return this->elt_geom.PointsQuadInsideND(N);
  }
  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNbPointsNodalElt() const
  {
    return this->elt_geom.GetNbPointsNodalElt();
  }
  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNbPointsDof() const
  {
    return this->elt_geom.GetNbPointsDof();
  }
  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNbPointsQuadrature() const
  {
    return this->elt_geom.GetNbPointsQuadrature();
  }
        
  inline const Array3D<int>& ElementGeomReferenceContainer<Dimension3>::GetNumNodes3D() const
  {
    return this->elt_geom.GetNumNodes3D();
  }

  inline int ElementGeomReferenceContainer<Dimension3>::GetNumNodes3D(int i, int j, int k) const
  {
    return this->elt_geom.GetNumNodes3D(i, j, k);
  }
  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNumNodes2D(int i, int j) const
  {
    return this->elt_geom.GetNumNodes2D(i, j);
  }
    
  inline int ElementGeomReferenceContainer<Dimension3>::GetNodalNumber(int num_face_loc, int num_dof) const
  {
    return this->elt_geom.GetNodalNumber(num_face_loc, num_dof);
  }
    
  inline const Vector<R3>& ElementGeomReferenceContainer<Dimension3>::PointsNodalND() const
  {
    return this->elt_geom.PointsNodalND();
  }
  
  inline const R3& ElementGeomReferenceContainer<Dimension3>::PointsNodalND(int j) const
  {
    return this->elt_geom.PointsNodalND(j);
  }
        
  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension3>::WeightsND() const
  {
    return this->elt_geom.WeightsND();
  }
  
  inline const VectR3& ElementGeomReferenceContainer<Dimension3>::PointsQuadND() const
  {
    return this->elt_geom.PointsQuadND();
  }
  
  inline const VectR3& ElementGeomReferenceContainer<Dimension3>::PointsND() const
  {
    return this->elt_geom.PointsND();
  }
  
  inline const Real_wp& ElementGeomReferenceContainer<Dimension3>::WeightsND(int j) const
  {
    return this->elt_geom.WeightsND(j);
  }
  
  inline const R3& ElementGeomReferenceContainer<Dimension3>::PointsND(int j) const
  {
    return this->elt_geom.PointsND(j);
  }
        
  inline void ElementGeomReferenceContainer<Dimension3>::SetPointsND(const VectR3& pts)
  {
    this->elt_geom.points3d = pts;
  }
  
  inline void ElementGeomReferenceContainer<Dimension3>::SetWeightsND(const VectReal_wp& w)
  {
    this->elt_geom.weights3d = w;
  }

  inline const Vector<R3>& ElementGeomReferenceContainer<Dimension3>::PointsDofND() const
  {
    return this->elt_geom.PointsDofND();
  }

  inline void ElementGeomReferenceContainer<Dimension3>::SetPointsDofND(const VectR3& pts)
  {
    this->elt_geom.points_dof3d = pts;
  }
  
  inline const R3& ElementGeomReferenceContainer<Dimension3>::PointsDofND(int j) const
  {
    return this->elt_geom.PointsDofND(j);
  }

  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::Points2D_tri() const
  {
    return this->elt_geom.Points2D_tri();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::Points2D_quad() const
  {
    return this->elt_geom.Points2D_quad();
  }
  
  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension3>::Weights2D_tri() const
  {
    return this->elt_geom.Weights2D_tri();
  }
  
  inline const VectReal_wp& ElementGeomReferenceContainer<Dimension3>::Weights2D_quad() const
  {
    return this->elt_geom.Weights2D_quad();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsNodal2D_tri() const
  {
    return this->elt_geom.PointsNodal2D_tri();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsNodal2D_quad() const
  {
    return this->elt_geom.PointsNodal2D_quad();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsDof2D_tri() const
  {
    return this->elt_geom.PointsDof2D_tri();
  }
  
  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsDof2D_quad() const
  {
    return this->elt_geom.PointsDof2D_quad();
  }

  inline void ElementGeomReferenceContainer<Dimension3>::SetPointsDof2D_quad(const VectR2& pts)
  {
    this->elt_geom.SetPointsDof2D_quad(pts);
  }
  
  inline void ElementGeomReferenceContainer<Dimension3>::SetPointsDof2D_tri(const VectR2& pts)
  {
    this->elt_geom.SetPointsDof2D_tri(pts);
  }

  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension3>::GetNumQuad2D() const
  {
    return this->elt_geom.GetNumQuad2D();
  }

  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension3>::GetNumNodes2D_tri() const
  {
    return this->elt_geom.GetNumNodes2D_tri();
  }
  
  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension3>::GetNumNodes2D_quad() const
  {
    return this->elt_geom.GetNumNodes2D_quad();
  }

  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension3>::GetCoordinateNodes2D_quad() const
  {
    return this->elt_geom.GetCoordinateNodes2D_quad();
  }
  
  inline const Matrix<int>& ElementGeomReferenceContainer<Dimension3>::GetCoordinateNodes2D_tri() const
  {
    return this->elt_geom.GetCoordinateNodes2D_tri();
  }

  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNbQuadBoundary(int n) const
  {
    return this->elt_geom.GetNbQuadBoundary(n);
  }
  
  inline int ElementGeomReferenceContainer<Dimension3>::GetNbNodalBoundary(int n) const
  {
    return this->elt_geom.GetNbNodalBoundary(n);
  }
      
  inline const R2& ElementGeomReferenceContainer<Dimension3>::PointsNodalBoundary(int k, int n) const
  {
    return this->elt_geom.PointsNodalBoundary(k, n);
  }

  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsQuadratureBoundary(int n) const
  {
    return this->elt_geom.PointsQuadratureBoundary(n);
  }

  inline const VectR2& ElementGeomReferenceContainer<Dimension3>::PointsDofBoundary(int n) const
  {
    return this->elt_geom.PointsDofBoundary(n);
  }
  
  inline const R2& ElementGeomReferenceContainer<Dimension3>::PointsQuadratureBoundary(int k, int n) const
  {
    return this->elt_geom.PointsQuadratureBoundary(k, n);
  }
  
  inline const Globatto<Real_wp>& ElementGeomReferenceContainer<Dimension3>
  ::GetNodalShapeFunctions1D() const
  {
    return this->elt_geom.lob_geom;
  }

  inline const Real_wp& ElementGeomReferenceContainer<Dimension3>
  ::WeightsQuadratureBoundary(int k, int n) const
  {
    return this->elt_geom.WeightsQuadratureBoundary(k, n);
  }

  inline const R3& ElementGeomReferenceContainer<Dimension3>::TangenteLocX(int n) const
  {
    return this->elt_geom.TangenteLocX(n);
  }
     
  inline const R3& ElementGeomReferenceContainer<Dimension3>::TangenteLocY(int n) const
  {
    return this->elt_geom.TangenteLocY(n);
  }
    
  inline void ElementGeomReferenceContainer<Dimension3>
  ::GetLocalCoordOnBoundary(int num_loc, const R2& t_loc, R3& res) const
  {
    return this->elt_geom.GetLocalCoordOnBoundary(num_loc, t_loc, res);
  }

  inline void ElementGeomReferenceContainer<Dimension3>
  ::GetLocalCoordOnBoundary(int num_loc, const R3& res, R2& t_loc) const
  {
    return this->elt_geom.GetLocalCoordOnBoundary(num_loc, res, t_loc);
  }

  inline R3 ElementGeomReferenceContainer<Dimension3>::GetTangentialVector(int num_loc, const R2& vec_u) const
  {
    return this->elt_geom.GetTangentialVector(num_loc, vec_u);
  }

  inline R2 ElementGeomReferenceContainer<Dimension3>::TransposeTangentialVector(int num_loc, const R3& vec_u) const
  {
    return this->elt_geom.TransposeTangentialVector(num_loc, vec_u);
  }

  inline R2_Complex_wp ElementGeomReferenceContainer<Dimension3>::TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const
  {
    return this->elt_geom.TransposeTangentialVector(num_loc, vec_u);
  }
  
} // end namespace

#define MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_INLINE_CXX
#endif
