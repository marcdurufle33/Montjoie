#ifndef MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_INLINE_CXX

namespace Montjoie
{

  inline const Matrix<int>& QuadrangleHierarchic::GetNumDofs2D() const
  {
    return NumDofsQuad;
  }
  
  inline void QuadrangleHierarchic
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHierarchic
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  inline void QuadrangleHierarchic::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void QuadrangleHierarchic::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void QuadrangleHierarchic::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void QuadrangleHierarchic::ApplyRhQuadrature(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhQuadratureGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyRhQuadrature(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhQuadratureGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleHierarchic::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }

  inline void QuadrangleHierarchic
  ::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHierarchic
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHierarchic
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHierarchic
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHierarchic
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
}

#define MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_INLINE_CXX
#endif
