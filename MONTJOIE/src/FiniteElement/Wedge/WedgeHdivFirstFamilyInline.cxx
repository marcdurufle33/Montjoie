#ifndef MONTJOIE_FILE_WEDGE_HDIV_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void WedgeHdivFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void WedgeHdivFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
}

#define MONTJOIE_FILE_WEDGE_HDIV_FIRST_FAMILY_INLINE_CXX
#endif
