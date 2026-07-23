#ifndef MONTJOIE_FILE_TRIANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{
  
  inline void TriangleHdivOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
    
  inline void TriangleHdivOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
}

#define MONTJOIE_FILE_TRIANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX  
#endif
