#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void HexahedronHdivFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHdivFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
    
  inline void HexahedronHdivFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void HexahedronHdivFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

}

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_FIRST_FAMILY_INLINE_CXX
#endif
