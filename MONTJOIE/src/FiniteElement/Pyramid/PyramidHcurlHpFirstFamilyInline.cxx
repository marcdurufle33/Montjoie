#ifndef MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void PyramidHcurlHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void PyramidHcurlHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_INLINE_CXX
#endif
