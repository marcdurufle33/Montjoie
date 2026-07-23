#ifndef MONTJOIE_FILE_HEXAHEDRON_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int type>
  inline HexahedronReference<type>::HexahedronReference() : ElementReference<Dimension3, type>(Fb_geom)
  {
    // serendip_element = true;
    serendip_element = false;
  }
    
}

#define MONTJOIE_FILE_HEXAHEDRON_REFERENCE_INLINE_CXX
#endif

    
