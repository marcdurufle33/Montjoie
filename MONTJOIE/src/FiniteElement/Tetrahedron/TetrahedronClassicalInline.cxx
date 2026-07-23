#ifndef MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_INLINE_CXX

namespace Montjoie
{

  inline void TetrahedronClassical
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TetrahedronClassical
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TetrahedronClassical
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void TetrahedronClassical
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void TetrahedronClassical
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TetrahedronClassical
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TetrahedronClassical::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
					  VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void TetrahedronClassical
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

}

#define MONTJOIE_FILE_TETRAHEDRON_CLASSICAL_INLINE_CXX
#endif
