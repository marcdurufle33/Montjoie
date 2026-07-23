#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline const Matrix<int>& TriangleHcurlOptimalHpFirstFamily::GetNumDofsX() const
  {
    return NumDofsX;
  }
  
  inline const Matrix<int>& TriangleHcurlOptimalHpFirstFamily::GetNumDofsY() const
  {
    return NumDofsY;
  }

  inline void TriangleHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void TriangleHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }

  
}

#define MONTJOIE_FILE_TRIANGLE_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX
#endif

