#ifndef MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void WedgeHdivHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void WedgeHdivHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
    
  inline void WedgeHdivHpFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void WedgeHdivHpFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

}

#define MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_INLINE_CXX
#endif
