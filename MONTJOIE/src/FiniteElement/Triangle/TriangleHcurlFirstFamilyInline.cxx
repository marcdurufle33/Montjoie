#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{
  
  inline void TriangleHcurlFirstFamily
  ::ComputeValuesPhiHRef(const R2& point_loc, VectReal_wp& res) const 
  {
    abort(); 
  }
  

  //! default constructor
  inline TriangleHcurlOptimalFirstFamily::TriangleHcurlOptimalFirstFamily()
    : TriangleHcurlFirstFamily()
  {
    type_nodal_basis = NODAL_LOBATTO;
  }
  
}

#define MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_INLINE_CXX
#endif
