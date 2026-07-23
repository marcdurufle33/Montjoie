#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_INLINE_CXX

namespace Montjoie
{

  //! returns weights of integration
  template<class Dim>
  inline const VectReal_wp& ElementReference_Dim<Dim>::WeightsND() const
  {
    return ElementGeomReferenceContainer<Dim>::WeightsND();
  }


  //! returns j-weight of integration
  template<class Dim>
  inline const Real_wp& ElementReference_Dim<Dim>::WeightsND(int j) const
  {
    return ElementGeomReferenceContainer<Dim>::WeightsND(j);
  }


#ifdef MONTJOIE_WITH_THREE_DIM
#ifdef MONTJOIE_WITH_NODAL_HCURL
  inline int ElementReferenceType<Dimension3, 2>::GetCoordinateDof(int i) const
  {
    return -1;
  }

  inline int ElementReferenceType<Dimension3, 2>::GetPointNode_FromDofNumber(int i) const
  {
    return -1;
  }
#endif
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_REFERENCE_INLINE_CXX
#endif
