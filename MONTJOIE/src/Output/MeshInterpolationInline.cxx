#ifndef MONTJOIE_FILE_MESH_INTERPOLATION_INLINE_CXX

namespace Montjoie
{
    
  //! returns the number of subdivisions
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetNbSubdivisions() const
  {
    return nb_subdivisions;
  }
  
  
  //! returns the subdivision used by the interpolation mesh
  template<class Dimension>
  inline const VectReal_wp& MeshInterpolation_Base<Dimension>::GetSubdivisionStep() const
  {
    return step_subdiv;
  }
  
  
  //! returns the subdivision on each reference element (unit tetrahedron, unit cube, etc)
  template<class Dimension>
  inline Vector<Vector<typename Dimension::R_N> >&
  MeshInterpolation_Base<Dimension>::PointsReferenceVolume()
  {
    return points_div;
  }

  
  //! returns the subdivision on each reference face (unit triangle, unit square)
  template<class Dimension>
  inline Vector<typename Dimension::DimensionBoundary::VectR_N>&
  MeshInterpolation_Base<Dimension>::PointsReferenceSurface()
  {
    return points_div_surf;
  }


  //! returns the weights for a given shape
  template<class Dimension>
  inline VectReal_wp& MeshInterpolation_Base<Dimension>::WeightsReferenceSurface(int n)
  {
    return weights_surf(n);
  }
  
  
  //! returns the file name for the output mesh
  template<class Dimension>
  inline const string& MeshInterpolation_Base<Dimension>::GetFileName() const
  {
    return file_name;
  }
  
    
  //! sets the file name for the output mesh
  template<class Dimension>
  inline void MeshInterpolation_Base<Dimension>::SetFileName(const string& name)
  {
    file_name = name;
  }
  
  
  //! sets the type of output mesh (VOLUME_MESH, SURFACE_MESH, etc)
  template<class Dimension>
  inline void MeshInterpolation_Base<Dimension>::SetType(int type)
  {
    type_mesh = type;
  }
  
  
  //! returns the type of output mesh (VOLUME_MESH, SURFACE_MESH, etc)
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetType() const
  {
    return type_mesh;
  }
  
  
  //! sets the references of surfaces where an output is required
  template<class Dimension>
  inline void MeshInterpolation_Base<Dimension>::SetSurfaceReference(const IVect& ref)
  {
    ref_surf = ref;
  }
  
  
  //! returns the references of surfaces where an output is required
  template<class Dimension>
  inline const IVect& MeshInterpolation_Base<Dimension>::GetSurfaceReference() const
  {
    return ref_surf;
  }
  
    
  //! returns the node numbers for each element of the mesh
  template<class Dimension>
  inline Vector<IVect>& MeshInterpolation_Base<Dimension>::NumberingArray()
  {
    return Nodle;
  }
  
  
  //! returns the node numbers for each element of the mesh
  template<class Dimension>
  inline const Vector<IVect>& MeshInterpolation_Base<Dimension>::NumberingArray() const
  {
    return Nodle;
  }
  
  
  //! computation of volumetric projectors
  template<class Dimension>
  inline void MeshInterpolation_Base<Dimension>::InitProjectionVolume(Mesh<Dimension>& mesh)
  {
    Vector<const ElementGeomReference<Dimension>* > reference_element;
    mesh.GetReferenceElementVolume(reference_element);
    proj.InitProjection(reference_element, step_subdiv, points_div);
  }
  
  
  //! computation of surfacic projectors
  template<class Dimension>
  inline void MeshInterpolation_Base<Dimension>::InitProjectionSurface(Mesh<Dimension>& mesh)
  {
    Vector<const ElementGeomReference<DimensionBoundary>* > reference_element;
    mesh.GetReferenceElementSurface(reference_element);
    proj_surf.InitProjection(reference_element, step_subdiv, points_div_surf);
  }
  
  
  //! volumetric projection from nodal points to interpolation points
  /*!
    \param[in] Uloc nodal values of U on element i
    \param[out] Udiv values of U on subdivisions (usually regular)
    \param[in] vars given problem
    \param[in] type_elt element type
  */
  template<class Dimension> template<class Vector1, class Vector2>
  inline void MeshInterpolation_Base<Dimension>
  ::ProjectVolume(const Vector1& Uloc, Vector2& Udiv, int type_elt) const
  {
    proj.ProjectScalar(Uloc, Udiv, type_elt);
  }
  
  
  //! volumetric projection from nodal points to interpolation points
  /*!
    \param[in] Uloc nodal values of U on element i
    \param[out] Udiv values of U on subdivisions (usually regular)
    \param[in] vars given problem
    \param[in] type_elt element type
  */
  template<class Dimension> template<class Vector1, class Vector2>
  inline void MeshInterpolation_Base<Dimension>
  ::ProjectSurface(const Vector1& Uloc, Vector2& Udiv, int type_elt) const
  {
    proj_surf.Project(Uloc, Udiv, type_elt);
  }
  
  
  //! volumetric mesh interpolator ?
  template<class Dimension>
  inline bool MeshInterpolation_Base<Dimension>::IsVolumetric() const
  {
    return (type_mesh == VOLUME_MESH);
  }
  
  
  //! surfacic mesh interpolator ?
  template<class Dimension>
  inline bool MeshInterpolation_Base<Dimension>::IsSurfacic() const
  {
    return (type_mesh != VOLUME_MESH);
  }
  
  
  //! returns the element number associated with face i
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetElementNumberOfSurface(int i) const
  {
    return NumElement(i);
  }
  
  
  //! returns the local face number in the element
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetLocalPositionOfSurface(int i) const
  {
    return LocalPositionBoundary(i);
  }
  
  
  //! returns the number of faces for which an output is required
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetNbBoundary() const
  {
    return NumElement.GetM();
  }
  
  
  //! returns the number of stored quadrature points (on the current processor)
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetNbPointsQuadrature() const
  {
    return nb_points_quadrature_per_proc(rank_processor);
  }
  
  
  //! returns the number of quadrature points on edge/face i
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetNbPointsQuadrature(int i) const
  {
    return PointsBoundary(i).GetM();
  }

  
  //! returns the nodal point j of edge/face i
  template<class Dimension>
  inline const typename Dimension::R_N&
  MeshInterpolation_Base<Dimension>::PointsNodal(int i, int j) const
  {
    return PointsNodalBoundary(i)(j);
  }
  
  
  //! returns the nodal points of edge/face i
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>&
  MeshInterpolation_Base<Dimension>::PointsNodal(int i) const
  {
    return PointsNodalBoundary(i);
  }
    
  
  //! returns the quadrature point j of edge/face i
  template<class Dimension>
  inline const typename Dimension::R_N&
  MeshInterpolation_Base<Dimension>::PointsQuadrature(int i, int j) const
  {
    return PointsBoundary(i)(j);
  }
  

  //! returns the quadrature points of edge/face i  
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>&
  MeshInterpolation_Base<Dimension>::PointsQuadrature(int i) const
  {
    return PointsBoundary(i);
  }

  
  //! returns the inverse of jacobian matrix (i.e. DF_i^{-1})
  //! for the quadrature point j of edge/face i
  template<class Dimension>
  inline const typename Dimension::MatrixN_N&
  MeshInterpolation_Base<Dimension>::Dfjm1Quadrature(int i, int j) const
  {
    return invDFiBoundary(i)(j);
  }
  

  //! returns the inverse of jacobian matrix (i.e. DF_i^{-1})
  //! for the quadrature points of edge/face i  
  template<class Dimension>
  inline const Vector<typename Dimension::MatrixN_N>&
  MeshInterpolation_Base<Dimension>::Dfjm1Quadrature(int i) const
  {
    return invDFiBoundary(i);
  }
    

  //! returns the normale at the quadrature point j of edge/face i
  template<class Dimension>
  inline const typename Dimension::R_N&
  MeshInterpolation_Base<Dimension>::NormaleQuadrature(int i, int j) const
  {
    return NormaleBoundary(i)(j);
  }
  
  
  //! returns the normale at the quadrature points of edge/face i
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>&
  MeshInterpolation_Base<Dimension>::NormaleQuadrature(int i) const
  {
    return NormaleBoundary(i);
  }
    
  
  //! return the weight of quadrature point j of edge/face i
  template<class Dimension>
  inline Real_wp& MeshInterpolation_Base<Dimension>::WeightsQuadrature(int i, int j)
  {
    return WeightsBoundary(i)(j);
  }
  

  //! return the weight of quadrature point j of edge/face i
  template<class Dimension>
  inline const Real_wp& MeshInterpolation_Base<Dimension>::WeightsQuadrature(int i, int j) const
  {
    return WeightsBoundary(i)(j);
  }
  
  
  //! return the weight of quadrature points of edge/face i
  template<class Dimension>
  inline const Vector<Real_wp>& MeshInterpolation_Base<Dimension>::WeightsQuadrature(int i) const
  {
    return WeightsBoundary(i);
  }

  
  //! return the surfacic element ds of quadrature point j of edge/face i
  template<class Dimension>
  inline const Real_wp& MeshInterpolation_Base<Dimension>::DsQuadrature(int i, int j) const
  {
    return DsBoundary(i)(j);
  }
  
 
  //! return the surfacic element ds of quadrature points of edge/face i
  template<class Dimension>
  inline const Vector<Real_wp>& MeshInterpolation_Base<Dimension>::DsQuadrature(int i) const
  {
    return DsBoundary(i);
  }
    

  //! returns the normales at nodal ponts of edge/face i
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>&
  MeshInterpolation_Base<Dimension>::GetNormale(int i) const
  {
    return NormaleBoundary(i);
  }

  
  //! returns all the quadrature points (gathered between processors)
  template<class Dimension>
  inline int MeshInterpolation_Base<Dimension>::GetNbAllQuadraturePoints() const
  {
    return AllPoints.GetM();
  }
  
  
  //! returns the global quadrature point k
  template<class Dimension>
  inline const typename Dimension::R_N&
  MeshInterpolation_Base<Dimension>::GetQuadraturePoint(int k) const
  {
    return AllPoints(k);
  }


  //! returns the global quadrature points
  template<class Dimension>
  inline const typename Dimension::VectR_N&
  MeshInterpolation_Base<Dimension>::GetQuadraturePoints() const
  {
    return AllPoints;
  }
  
  
  //! returns the normale at the global quadrature point k
  template<class Dimension>
  inline const typename Dimension::R_N&
  MeshInterpolation_Base<Dimension>::GetQuadratureNormale(int k) const
  {
    return AllNormale(k);
  }
  
  
  //! returns the weight of the global quadrature point k
  template<class Dimension>
  inline const Real_wp& MeshInterpolation_Base<Dimension>::GetQuadratureWeight(int k) const
  {
    return AllWeights(k);
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_INTERPOLATION_INLINE_CXX
#endif
