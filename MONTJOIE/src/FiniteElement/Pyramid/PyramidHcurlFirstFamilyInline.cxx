#ifndef MONTJOIE_FILE_PYRAMID_HCURL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void PyramidHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void PyramidHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HCURL_FIRST_FAMILY_INLINE_CXX
#endif
