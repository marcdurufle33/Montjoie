#ifndef MONTJOIE_FILE_TETRAHEDRON_DG_ORTHO_INLINE_CXX

namespace Montjoie
{

  //! x is overwritten by M^-1 x where M is the mass matrix    
  template<class Vector1>
  inline void TetrahedronDgOrtho::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity, we don't change x
  }
  
  
  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector1>
  inline void TetrahedronDgOrtho::MltMassMatrixGen(Vector1& y) const
  {
    // mass matrix equal to identity
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on point
  */  
  inline void TetrahedronDgOrtho::ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& phi) const
  {
    Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct3D, InvWeightFct, pointloc, phi);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point
  */
  inline void TetrahedronDgOrtho::ComputeGradientPhiRef(const R3& pointloc, VectR3& grad_phi) const
  {
    Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct3D, InvWeightFct, pointloc, grad_phi);
  }
  
  
  inline void TetrahedronDgOrtho
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void TetrahedronDgOrtho
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  inline void TetrahedronDgOrtho::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void TetrahedronDgOrtho::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }

  inline void TetrahedronDgOrtho::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void TetrahedronDgOrtho::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void TetrahedronDgOrtho::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void TetrahedronDgOrtho::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void TetrahedronDgOrtho::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void TetrahedronDgOrtho::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void TetrahedronDgOrtho::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }

  inline void TetrahedronDgOrtho::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void TetrahedronDgOrtho::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void TetrahedronDgOrtho
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void TetrahedronDgOrtho
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void TetrahedronDgOrtho
  ::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
			       VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void TetrahedronDgOrtho
  ::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void TetrahedronDgOrtho
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void TetrahedronDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
							   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void TetrahedronDgOrtho::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
							   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void TetrahedronDgOrtho
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void TetrahedronDgOrtho
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

}

#define MONTJOIE_FILE_TETRAHEDRON_DG_ORTHO_INLINE_CXX
#endif
