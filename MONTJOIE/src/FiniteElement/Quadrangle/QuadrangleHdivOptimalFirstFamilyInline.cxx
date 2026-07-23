#ifndef MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void QuadrangleHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void QuadrangleHdivOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
}

#define MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_FIRST_FAMILY_INLINE_CXX
#endif
