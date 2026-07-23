#ifndef MONTJOIE_FILE_TRIANGLE_DG_ORTHO_INLINE_CXX

namespace Montjoie
{

  inline const VectReal_wp& TriangleDgOrtho::GetInverseWeightFunction() const
  {
    return InvWeightFct; 
  }
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  inline void TriangleDgOrtho::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct2D, InvWeightFct, pointloc, phi);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  inline void TriangleDgOrtho::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct2D, InvWeightFct, pointloc, grad_phi);
  }

  inline void TriangleDgOrtho::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void TriangleDgOrtho::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void TriangleDgOrtho::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void TriangleDgOrtho::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void TriangleDgOrtho::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void TriangleDgOrtho::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void TriangleDgOrtho
  ::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void TriangleDgOrtho
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void TriangleDgOrtho::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void TriangleDgOrtho::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void TriangleDgOrtho::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void TriangleDgOrtho::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void TriangleDgOrtho
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void TriangleDgOrtho
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void TriangleDgOrtho
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleDgOrtho::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							  VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void TriangleDgOrtho
  ::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void TriangleDgOrtho
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void TriangleDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
							       VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void TriangleDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							       VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void TriangleDgOrtho
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void TriangleDgOrtho
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

}

#define MONTJOIE_FILE_TRIANGLE_DG_ORTHO_INLINE_CXX
#endif
