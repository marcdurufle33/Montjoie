#ifndef MONTJOIE_FILE_WEDGE_DG_LEGENDRE_INLINE_CXX

namespace Montjoie
{

  //! x is overwritten by M^-1 x where M is the mass matrix    
  template<class Vector1>
  inline void WedgeDgLegendre::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity, we don't change x
  }
  
  
  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector1>
  inline void WedgeDgLegendre::MltMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity
  }
  

  inline void WedgeDgLegendre::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void WedgeDgLegendre::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void WedgeDgLegendre::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void WedgeDgLegendre::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void WedgeDgLegendre::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void WedgeDgLegendre::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void WedgeDgLegendre::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void WedgeDgLegendre::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void WedgeDgLegendre::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void WedgeDgLegendre::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void WedgeDgLegendre::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void WedgeDgLegendre::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void WedgeDgLegendre
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgLegendre::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							  VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgLegendre
  ::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeDgLegendre
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeDgLegendre::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
						     VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void WedgeDgLegendre::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
						     VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void WedgeDgLegendre
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgLegendre
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeDgLegendre
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void WedgeDgLegendre
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_WEDGE_DG_LEGENDRE_INLINE_CXX
#endif
