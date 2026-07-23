#ifndef MONTJOIE_FILE_MESH_BOUNDARIES_INLINE_CXX

namespace Montjoie
{

  /**********************
   * PmlRegionParameter *
   **********************/
  
  
  //! returns the thickness of PML layers
  template <class Dimension, class T>
  inline const T& PmlRegionParameter_Base<Dimension, T>::GetThicknessPML() const
  {
    return thicknessPML;
  }
  

  //! returns the origin associated with PML layers
  template <class Dimension, class T>
  inline const TinyVector<T, Dimension::dim_N>&
  PmlRegionParameter_Base<Dimension, T>::GetOriginRadialPML() const
  {
    return originPML;
  }

  
  //! returns the radius associated with PML layers
  template <class Dimension, class T>
  inline const T& PmlRegionParameter_Base<Dimension, T>::GetRadiusPML() const
  {
    return radiusPML;
  }


  //! sets the thickness of PML layers
  template <class Dimension, class T>
  inline void PmlRegionParameter_Base<Dimension, T>::SetThicknessPML(const T& delta)
  {
    thicknessPML = delta;
  }
  
  
  //! returns the type of addition for PML layers along the coordinate num_coor
  template <class Dimension, class T>
  inline int PmlRegionParameter_Base<Dimension, T>::GetTypeAdditionPML(int num_coor) const
  {
    switch (num_coor)
      {
      case 0:
	return pml_add_axis_x;
      case 1:
	return pml_add_axis_y;
      case 2:
	return pml_add_axis_z;
      }
    
    return -1;
  }
  
  
  //! returns the number of PML layers
  template <class Dimension, class T>
  inline int PmlRegionParameter_Base<Dimension, T>::GetNbLayersPML() const
  {
    return nb_layers_pml;
  }
  
  
  //! sets the radius associated with PML layers
  template <class Dimension, class T>
  inline void PmlRegionParameter_Base<Dimension, T>::SetRadiusPML(const T& r)
  {
    radiusPML = r;
  }


  //! sets the origin associated with PML layers
  template <class Dimension, class T>
  inline void
  PmlRegionParameter_Base<Dimension, T>::SetOriginRadialPML(const TinyVector<T, Dimension::dim_N>& pt)
  {
    originPML = pt;
  }

  
  //! returns the leaf class
  template <class Dimension, class T>
  inline PmlRegionParameter<Dimension, T>& PmlRegionParameter_Base<Dimension, T>::GetLeafClass()
  {
    return static_cast<PmlRegionParameter<Dimension, T>& >(*this);
  }
  
  
  /***********************
   * MeshBoundaries_Base *
   ***********************/
  
  
  //! returns leaf class
  template <class Dimension, class T>
  inline Mesh<Dimension, T>& MeshBoundaries_Base<Dimension, T>::GetLeafClass()
  {
    return static_cast<Mesh<Dimension, T>& >(*this);
  }
  
  
  //! returns leaf class
  template <class Dimension, class T>
  inline const Mesh<Dimension, T>& MeshBoundaries_Base<Dimension, T>::GetLeafClass() const
  {
    return static_cast<const Mesh<Dimension, T>& >(*this);
  }
  
  
  //! returns the number of references in the mesh
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetNbReferences() const
  {
    return this->Type_curve.GetM()-1;
  }

  
  //! returns the type of the curved boundary i
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetCurveType(int i) const
  {
    if (i < 1)
      {
        // reference 0 is affected for internal faces -> no curve
        
        // reference -1 is affected to "curved faces" found in a mesh file 
        // without physical reference    
        
        // reference -2 is affected to "curved faces" found for predefined curves        
        if (i == -1)
          return this->GetLeafClass().CURVE_FILE;
        else 
          return this->GetLeafClass().NO_CURVE;
      }
    
    return Type_curve(i);
  }
  
  
  //! returns the type of the curved boundaries
  template <class Dimension, class T>
  inline const IVect& MeshBoundaries_Base<Dimension, T>::GetCurveType() const
  {
    return Type_curve;
  }
  
  
  //! returns the boundary conditions associated with faces of reference i
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetBoundaryCondition(int i) const
  {
    if (i < 0)
      return 0;
    
    return Cond_curve(i);
  }

  
  //! returns the boundary conditions associated with references of faces
  template <class Dimension, class T>
  inline const IVect& MeshBoundaries_Base<Dimension, T>::GetBoundaryCondition() const
  {
    return Cond_curve;
  }
  
  
  //! returns the body number associated with faces of reference i
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetBodyNumber(int i) const
  {
    if (i < 0)
      return 0;

    return Body_curve(i);
  }
  
  
  //! returns the body numbers associated with references of faces
  template <class Dimension, class T>
  inline const IVect& MeshBoundaries_Base<Dimension, T>::GetBodyNumber() const
  {
    return Body_curve;
  }
  
  
  //! returns the parameters of the reference i (for example radius and center of sphere)
  template <class Dimension, class T>
  inline void MeshBoundaries_Base<Dimension, T>::GetCurveParameter(int i, Vector<T>& param) const
  {
    if (i < 0)
      param.Clear();
    else
      param = Param_curve(i);
  }

  
  //! returns the parameters of the reference i (for example radius and center of sphere)
  template <class Dimension, class T>
  inline const Vector<T>& MeshBoundaries_Base<Dimension, T>::GetCurveParameter(int i) const
  {
    if (i < 0)
      return Param_curve(0);
    
    return Param_curve(i);
  }
  
  
  //! returns the parameters of all the references (for example radius and center of spheres)
  template <class Dimension, class T>
  inline const Vector<Vector<T> >& MeshBoundaries_Base<Dimension, T>::GetCurveParameter() const
  {
    return Param_curve;
  }
  
  
  //! returns 1-D basis functions
  template <class Dimension, class T>
  inline const Globatto<Real_wp>& MeshBoundaries_Base<Dimension, T>::GetCurveFunctions1D() const
  {
    return lob_curve;
  }
  
  
  //! returns the number of PML areas
  template <class Dimension, class T>
  int MeshBoundaries_Base<Dimension, T>::GetNbPmlAreas() const
  {
    return pml_areas.GetM();
  }
  
  
  //! changes the number of PML areas
  template <class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::ReallocatePmlAreas(int n)
  {
    pml_areas.Reallocate(n);
  }
  

  //! returns parameters defining the PML i
  template <class Dimension, class T>
  PmlRegionParameter<Dimension, T>& MeshBoundaries_Base<Dimension, T>::GetPmlArea(int i)
  {
    return pml_areas(i);
  }


  //! returns parameters defining the PML i
  template <class Dimension, class T>
  const PmlRegionParameter<Dimension, T>& MeshBoundaries_Base<Dimension, T>::GetPmlArea(int i) const
  {
    return pml_areas(i);
  }


  //! returns the number of periodic references
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetNbPeriodicReferences() const
  {
    return periodicity_references.GetM();
  }
  

  //! returns the periodic references
  template <class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetPeriodicReference(int i, int j) const
  {
    return periodicity_references(i)(j);
  }
  
  
  //! adds a pair of periodic references and type of periodicity
  template<class Dimension, class T>
  inline void MeshBoundaries_Base<Dimension, T>::
  AddPeriodicCondition(const TinyVector<int, 2>& ref, int type_coor)
  {
    periodicity_references.PushBack(ref);
    type_coord_periodicity.PushBack(type_coor);
  }
  
  
  //! returns the type of periodicity
  template<class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetPeriodicityTypeReference(int i) const
  {
    return type_coord_periodicity(i);;
  }
  
  
  //! returns the cyclicity angle alpha
  template <class Dimension, class T>
  inline Real_wp MeshBoundaries_Base<Dimension, T>::GetPeriodicAlpha() const
  {
    return periodic_alpha_polar;
  }
  
  
  //! sets the cyclicity angle alpha
  template <class Dimension, class T>
  inline void MeshBoundaries_Base<Dimension, T>::SetPeriodicAlpha(const Real_wp& alpha)
  {
    periodic_alpha_polar = alpha;
  }
  
  
  /* Methods for curved boundaries */
  
  
  //! returns safety coefficient used to bound an element by a rectangular box
  template <class Dimension, class T>
  inline T MeshBoundaries_Base<Dimension, T>::GetSafetyCoef(int i) const
  {
    const Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    if (mesh_leaf.Element(i).IsCurved() )
      return 0.2;
    else
      return TinyVector<T, Dimension::dim_N>::threshold;
  }

  
  //! returns the order of approximation for the geometry
  template<class Dimension, class T>
  inline int MeshBoundaries_Base<Dimension, T>::GetGeometryOrder() const
  {
    return lob_curve.GetOrder();
  }
  
  
  //! evaluates 1-D shape functions (Lagrange function)
  template<class Dimension, class T>
  inline Real_wp MeshBoundaries_Base<Dimension, T>::GetInterpolate1D(int i, const Real_wp& x)
  {
    return lob_curve.EvaluatePhi(i, x);
  }
  
  
  //! returns an internal nodal point on interval [0, 1]
  template<class Dimension, class T>
  inline const Real_wp& MeshBoundaries_Base<Dimension, T>::GetNodalPointInsideEdge(int i) const
  {
    return lob_curve.Points(i+1);
  }
  
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BOUNDARIES_INLINE_CXX
#endif
