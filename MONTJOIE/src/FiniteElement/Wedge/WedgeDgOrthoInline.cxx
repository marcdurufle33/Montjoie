#ifndef MONTJOIE_FILE_WEDGE_DG_ORTHO_INLINE_CXX

namespace Montjoie
{

  inline void WedgeDgOrtho::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void WedgeDgOrtho::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void WedgeDgOrtho::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void WedgeDgOrtho::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void WedgeDgOrtho::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void WedgeDgOrtho::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void WedgeDgOrtho::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void WedgeDgOrtho::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void WedgeDgOrtho::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void WedgeDgOrtho::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void WedgeDgOrtho::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void WedgeDgOrtho::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void WedgeDgOrtho
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgOrtho::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
						VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgOrtho::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
					      VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeDgOrtho
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
						     VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void WedgeDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
						     VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void WedgeDgOrtho::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
					    VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgOrtho
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgOrtho::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void WedgeDgOrtho::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_WEDGE_DG_ORTHO_INLINE_CXX
#endif
