#ifndef MONTJOIE_FILE_TETRAHEDRON_REFERENCE_INLINE_CXX

namespace Montjoie
{
  //! default constructor
  template<int type>
  inline TetrahedronReference<type>::TetrahedronReference()
    : ElementReference<Dimension3, type>(Fb_geom)
  {
    type_interpolation = TetrahedronGeomReference::LOBATTO_BASIS;
  }
  
  
  //! modification of class with a line of the data file
  template<int type>
  inline void TetrahedronReference<type>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (!keyword.compare("TypeInterpolation"))
      {
	if (!param(0).compare("REGULAR"))
	  type_interpolation = TetrahedronGeomReference::REGULAR_BASIS;
	else
	  type_interpolation = TetrahedronGeomReference::LOBATTO_BASIS;
      }
  }
  
}

#define MONTJOIE_FILE_TETRAHEDRON_REFERENCE_INLINE_CXX
#endif
