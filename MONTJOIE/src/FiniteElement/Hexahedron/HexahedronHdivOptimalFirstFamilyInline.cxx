#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void HexahedronHdivOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHdivOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
    
  inline void HexahedronHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void HexahedronHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX
#endif
  
