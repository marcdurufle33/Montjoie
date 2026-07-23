#ifndef MONTJOIE_FILE_WEDGE_HIERARCHIC_INLINE_CXX

namespace Montjoie
{

  inline void WedgeHierarchic::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void WedgeHierarchic::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void WedgeHierarchic::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void WedgeHierarchic::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void WedgeHierarchic::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void WedgeHierarchic::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void WedgeHierarchic::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void WedgeHierarchic::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void WedgeHierarchic::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void WedgeHierarchic::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void WedgeHierarchic::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }

  inline void WedgeHierarchic::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeHierarchic::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }

  inline void WedgeHierarchic
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeHierarchic::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
					 VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void WedgeHierarchic::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void WedgeHierarchic::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
				       VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

}

#define MONTJOIE_FILE_WEDGE_HIERARCHIC_INLINE_CXX
#endif
