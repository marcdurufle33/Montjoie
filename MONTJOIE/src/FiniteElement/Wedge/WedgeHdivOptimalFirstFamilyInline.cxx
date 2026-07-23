#ifndef MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void WedgeHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void WedgeHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
}

#define MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX
#endif
