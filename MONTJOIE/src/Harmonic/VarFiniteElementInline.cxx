#ifndef MONTJOIE_FILE_VAR_FINITE_ELEMENT_INLINE_CXX

namespace Montjoie
{

  /************************
   * VarFiniteElementEnum *
   ************************/
  
  
  //! returns identity number associated with element name_elt (t = hybrid type)
  template<class Dimension, int type>
  inline int VarFiniteElementEnum_Base<Dimension, type>
  ::GetIdentityNumber(const string& name_elt, int t)
  {
    return list_element(t)[name_elt];
  }

  
  //! initializes list_element
  template<class Dimension, int type>
  inline void VarFiniteElementEnum<Dimension, type>::InitStaticData()
  {
  }


  //! returns a new pointer to the appropriate ElementReference class
  template<class Dimension, int type>
  inline ElementReference<Dimension, type>* VarFiniteElementEnum<Dimension, type>
  ::GetNewReferenceElement(int id)
  {
    return NULL;
  }


  /********************
   * VarFiniteElement *
   ********************/
  

  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline const ElementReference_Dim<Dimension>& VarFiniteElement<Dimension>
  ::GetReferenceElement(int i, int n) const
  {
    if (n == 0)
      return *reference_element(i);

    return *reference_element(i + n*this->nb_elt_);
  }
  
  
  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline ElementReference_Dim<Dimension>& VarFiniteElement<Dimension>
  ::GetReferenceElement(int i, int n)
  {
    if (n == 0)
      return *reference_element(i);

    return *reference_element(i + n*this->nb_elt_);
  }



  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline const ElementReference<Dimension, 1>& VarFiniteElement<Dimension>
  ::GetReferenceElementH1(int i, int n) const
  {
    if (n == 0)
      return dynamic_cast<const ElementReference<Dimension, 1>& >(*reference_element(i));

    return dynamic_cast<const ElementReference<Dimension, 1>& >(*reference_element(i + n*this->nb_elt_));
  }
  
  
  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline ElementReference<Dimension, 1>& VarFiniteElement<Dimension>
  ::GetReferenceElementH1(int i, int n)
  {
    if (n == 0)
      return dynamic_cast<ElementReference<Dimension, 1>& >(*reference_element(i));

    return dynamic_cast<ElementReference<Dimension, 1>& >(*reference_element(i + n*this->nb_elt_));
  }


    //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline const ElementReference<Dimension, 2>& VarFiniteElement<Dimension>
  ::GetReferenceElementHcurl(int i, int n) const
  {
    if (n == 0)
      return dynamic_cast<const ElementReference<Dimension, 2>& >(*reference_element(i));

    return dynamic_cast<const ElementReference<Dimension, 2>& >(*reference_element(i + n*this->nb_elt_));
  }
  
  
  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline ElementReference<Dimension, 2>& VarFiniteElement<Dimension>
  ::GetReferenceElementHcurl(int i, int n)
  {
    if (n == 0)
      return dynamic_cast<ElementReference<Dimension, 2>& >(*reference_element(i));

    return dynamic_cast<ElementReference<Dimension, 2>& >(*reference_element(i + n*this->nb_elt_));
  }


    //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline const ElementReference<Dimension, 3>& VarFiniteElement<Dimension>
  ::GetReferenceElementHdiv(int i, int n) const
  {
    if (n == 0)
      return dynamic_cast<const ElementReference<Dimension, 3>& >(*reference_element(i));

    return dynamic_cast<const ElementReference<Dimension, 3>& >(*reference_element(i + n*this->nb_elt_));
  }
  
  
  //! returns a reference to the reference element associated with element i
  template<class Dimension>
  inline ElementReference<Dimension, 3>& VarFiniteElement<Dimension>
  ::GetReferenceElementHdiv(int i, int n)
  {
    if (n == 0)
      return dynamic_cast<ElementReference<Dimension, 3>& >(*reference_element(i));

    return dynamic_cast<ElementReference<Dimension, 3>& >(*reference_element(i + n*this->nb_elt_));
  }

  
  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline const ElementReference_Dim<typename Dimension::DimensionBoundary>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElement(int i, int n) const
  {
    if (n == 0)
      return *surface_element(i);

    return *surface_element(i + n*this->nb_surf_);
  }

  
  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline ElementReference_Dim<typename Dimension::DimensionBoundary>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElement(int i, int n)
  {
    if (n == 0)
      return *surface_element(i);

    return *surface_element(i + n*this->nb_surf_);
  }


  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline const ElementReference<typename Dimension::DimensionBoundary, 1>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElementH1(int i, int n) const
  {
    if (n == 0)
      return dynamic_cast<const ElementReference<typename Dimension::DimensionBoundary, 1>& >(*surface_element(i));

    return dynamic_cast<const ElementReference<typename Dimension::DimensionBoundary, 1>& >(*surface_element(i + n*this->nb_surf_));
  }

  
  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline ElementReference<typename Dimension::DimensionBoundary, 1>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElementH1(int i, int n)
  {
    if (n == 0)
      return dynamic_cast<ElementReference<typename Dimension::DimensionBoundary, 1>& >(*surface_element(i));

    return dynamic_cast<ElementReference<typename Dimension::DimensionBoundary, 1>& >(*surface_element(i + n*this->nb_surf_));
  }


  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline const ElementReference<typename Dimension::DimensionBoundary, 2>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElementHcurl(int i, int n) const
  {
    if (n == 0)
      return dynamic_cast<const ElementReference<typename Dimension::DimensionBoundary, 2>& >(*surface_element(i));

    return dynamic_cast<const ElementReference<typename Dimension::DimensionBoundary, 2>& >(*surface_element(i + n*this->nb_surf_));
  }

  
  //! returns a reference to the reference element associated with face i
  template<class Dimension>
  inline ElementReference<typename Dimension::DimensionBoundary, 2>&
  VarFiniteElement<Dimension>::GetSurfaceFiniteElementHcurl(int i, int n)
  {
    if (n == 0)
      return dynamic_cast<ElementReference<typename Dimension::DimensionBoundary, 2>& >(*surface_element(i));

    return dynamic_cast<ElementReference<typename Dimension::DimensionBoundary, 2>& >(*surface_element(i + n*this->nb_surf_));
  }
  

  //! returns the reference element associated with the i-th neighbor element (parallel)
  template<class Dimension>
  inline const ElementReference_Dim<Dimension>& 
  VarFiniteElement<Dimension>::GetNeighborReferenceElement(int i, int n) const
  {
    if (n == 0)
      return *reference_neighbor_element(i);
    
    return *reference_neighbor_element(i + n*this->nb_surf_neigh_);
  }
  

  //! returns the reference element associated with the i-th neighbor element (parallel)
  template<class Dimension>
  inline ElementReference_Dim<Dimension>&
  VarFiniteElement<Dimension>::GetNeighborReferenceElement(int i, int n)
  {
    if (n == 0)
      return *reference_neighbor_element(i);
    
    return *reference_neighbor_element(i + n*this->nb_surf_neigh_);
  }
    
}

#define MONTJOIE_FILE_VAR_FINITE_ELEMENT_INLINE_CXX
#endif

