#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline const Matrix<int>& QuadrangleHcurlOptimalHpFirstFamily::GetNumDofsX() const
  {
    return NumDofsX;
  }
  
  inline const Matrix<int>& QuadrangleHcurlOptimalHpFirstFamily::GetNumDofsY() const
  {
    return NumDofsY;
  }

  inline void QuadrangleHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX
#endif
