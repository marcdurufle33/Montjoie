#ifndef MONTJOIE_FILE_PYRAMID_REFERENCE_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<int type>
  inline PyramidReference<type>::PyramidReference() : ElementReference<Dimension3, type>(Fb_geom)
  {
    type_interpolation = PyramidGeomReference::LOBATTO_BASIS;
  }
    
}

#define MONTJOIE_FILE_PYRAMID_REFERENCE_INLINE_CXX
#endif

  
