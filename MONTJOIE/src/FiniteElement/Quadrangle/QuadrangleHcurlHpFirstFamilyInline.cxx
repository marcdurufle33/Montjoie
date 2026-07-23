#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline const Matrix<int>& QuadrangleHcurlHpFirstFamily::GetNumDofsX() const
  {
    return NumDofsX;
  }
  
  inline const Matrix<int>& QuadrangleHcurlHpFirstFamily::GetNumDofsY() const
  {
    return NumDofsY; 
  }
  
  inline void QuadrangleHcurlHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_HP_FIRST_FAMILY_INLINE_CXX
#endif
