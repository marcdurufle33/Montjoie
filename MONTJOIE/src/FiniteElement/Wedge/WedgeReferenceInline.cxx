#ifndef MONTJOIE_FILE_WEDGE_REFERENCE_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  template<int type>
  inline WedgeReference<type>::WedgeReference() : ElementReference<Dimension3, type>(Fb_geom)
  {
  }
    
} // end namespace

#define MONTJOIE_FILE_WEDGE_REFERENCE_INLINE_CXX
#endif
