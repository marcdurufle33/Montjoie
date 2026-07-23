#ifndef MONTJOIE_FILE_PYRAMID_DG_ORTHO_INLINE_CXX

namespace Montjoie
{

  //! x is overwritten by M^-1 x where M is the mass matrix  
  template<class Vector1>
  inline void PyramidDgOrtho::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity, we don't change x
  }
  

  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector1>
  inline void PyramidDgOrtho::MltMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on point
  */  
  inline void PyramidDgOrtho::ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& phi) const
  {
    Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct3D, InvWeightFct, pointloc, phi);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point
  */
  inline void PyramidDgOrtho::ComputeGradientPhiRef(const R3& pointloc, VectR3& grad_phi) const
  {
    Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct3D, InvWeightFct, pointloc, grad_phi);
  }
  

  inline void PyramidDgOrtho::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void PyramidDgOrtho::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void PyramidDgOrtho::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void PyramidDgOrtho::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void PyramidDgOrtho::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void PyramidDgOrtho::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void PyramidDgOrtho::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void PyramidDgOrtho::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void PyramidDgOrtho::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void PyramidDgOrtho::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void PyramidDgOrtho::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void PyramidDgOrtho::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void PyramidDgOrtho
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgOrtho::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							 VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgOrtho::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
						VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void PyramidDgOrtho
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void PyramidDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
						       VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void PyramidDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							      VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void PyramidDgOrtho
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgOrtho
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void PyramidDgOrtho
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void PyramidDgOrtho
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

}

#define MONTJOIE_FILE_PYRAMID_DG_ORTHO_INLINE_CXX
#endif

