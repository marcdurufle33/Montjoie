#ifndef MONTJOIE_FILE_PYRAMID_DG_LEGENDRE_INLINE_CXX

namespace Montjoie
{

  //! x is overwritten by M^-1 x where M is the mass matrix  
  template<class Vector1>
  inline void PyramidDgLegendre::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity, we don't change x
  }
  

  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector1>
  inline void PyramidDgLegendre::MltMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity
  }
  

  inline void PyramidDgLegendre::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void PyramidDgLegendre::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void PyramidDgLegendre::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void PyramidDgLegendre::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void PyramidDgLegendre::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void PyramidDgLegendre::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void PyramidDgLegendre::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void PyramidDgLegendre::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void PyramidDgLegendre::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void PyramidDgLegendre::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void PyramidDgLegendre::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void PyramidDgLegendre
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void PyramidDgLegendre
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgLegendre
  ::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
			       VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgLegendre
  ::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void PyramidDgLegendre
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void PyramidDgLegendre::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
								 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void PyramidDgLegendre::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
								 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void PyramidDgLegendre
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgLegendre
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgLegendre
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void PyramidDgLegendre
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_PYRAMID_DG_LEGENDRE_INLINE_CXX
#endif

