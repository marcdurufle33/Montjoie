#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_INLINE_CXX

namespace Montjoie
{
  
  inline void TriangleHcurlSecondFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void TriangleHcurlSecondFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_INLINE_CXX
#endif

  
