#ifndef MONTJOIE_FILE_TRIANGLE_CLASSICAL_INLINE_CXX

namespace Montjoie
{

  inline void TriangleClassical
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleClassical
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleClassical
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void TriangleClassical
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void TriangleClassical
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleClassical
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleClassical::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				       VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void TriangleClassical
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  
  /*******************
   * TriangleLobatto *
   *******************/


  inline void TriangleLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void TriangleLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
}

#define MONTJOIE_FILE_TRIANGLE_CLASSICAL_INLINE_CXX
#endif
