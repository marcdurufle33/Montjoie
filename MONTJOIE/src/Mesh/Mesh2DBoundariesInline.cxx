#ifndef MONTJOIE_FILE_MESH2D_BOUNDARIES_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor
  template<class T>
  inline MeshBoundaries<Dimension2, T>::MeshBoundaries()
    : MeshBoundaries_Base<Dimension2, T>()
  {
  }

  
  /************************
   * Convenient functions *
   ************************/
  
  
  template<class T>
  inline const TinyVector<T, 2>& MeshBoundaries<Dimension2, T>::
  GetPointInsideEdge(int ne, int k) const
  {
    return PointsEdgeRef(ne, k);
  }
  
  
  template<class T>
  inline void MeshBoundaries<Dimension2, T>::
  SetPointInsideEdge(int ne, int k, const TinyVector<T, 2>& point)
  {
    PointsEdgeRef(ne, k) = point;
  }
  
  
  template<class T>
  inline bool MeshBoundaries<Dimension2, T>::
  IsBoundaryCurved(int ne) const
  {
    if (ne >= this->GetLeafClass().GetNbBoundaryRef())
      return false;
    
    if (this->Type_curve(this->GetLeafClass().edges_ref(ne).GetReference()) >= 1)
      return true;
    
    return false;
  }

  
  //! returns finite element over the unit interval
  template<class T>
  inline const EdgeLobatto& MeshBoundaries<Dimension2, T>::GetSurfaceFiniteElement() const
  {
    return interval_reference;
  }
  
  
  //! returns finite element over the unit interval
  template<class T>
  inline const EdgeLobatto& MeshBoundaries<Dimension2, T>::GetSurfaceFiniteElement2() const
  {
    return interval_reference;
  }


  //! returns parameters associated with curved boundaries
  template<class T>
  inline const Vector<TinyVector<T, 2> >& MeshBoundaries<Dimension2, T>::GetValueParameterOnBoundary() const
  {
    return value_parameter_boundary;
  }
  
  
  template<class T>
  inline void MeshBoundaries<Dimension2, T>::ClearConnectivity()
  {
    MeshBoundaries_Base<Dimension2, T>::ClearConnectivity();
  }
  
  
  template<class T>
  inline void MeshBoundaries<Dimension2, T>::ClearCurves()
  {
    MeshBoundaries_Base<Dimension2, T>::ClearCurves();
    PointsEdgeRef.Clear(); value_parameter_boundary.Clear(); integer_parameter_boundary.Clear();
    this->lob_curve.ConstructQuadrature(1, this->lob_curve.QUADRATURE_LOBATTO);
  }
  
  
  template<class T>
  inline void MeshBoundaries<Dimension2, T>::Clear()
  {
    MeshBoundaries_Base<Dimension2, T>::Clear();
    PointsEdgeRef.Clear();
    integer_parameter_boundary.Clear();
    value_parameter_boundary.Clear();
  }
    
}

#define MONTJOIE_FILE_MESH2D_BOUNDARIES_INLINE_CXX
#endif
