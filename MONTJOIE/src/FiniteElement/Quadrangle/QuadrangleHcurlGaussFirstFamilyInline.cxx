#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void QuadrangleHcurlGaussFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlGaussFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_INLINE_CXX
#endif
