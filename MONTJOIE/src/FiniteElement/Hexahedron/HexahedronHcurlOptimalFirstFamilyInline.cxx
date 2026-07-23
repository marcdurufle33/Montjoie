#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{
  
  inline const Array3D<int>& HexahedronHcurlOptimalFirstFamily::GetNumDofsX() const
  {
    return NumDofs_X;
  }
  
  inline const Array3D<int>& HexahedronHcurlOptimalFirstFamily::GetNumDofsY() const
  {
    return NumDofs_Y;
  }
  
  inline const Array3D<int>& HexahedronHcurlOptimalFirstFamily::GetNumDofsZ() const
  {
    return NumDofs_Z;
  }
  
  inline const Matrix<Real_wp>& HexahedronHcurlOptimalFirstFamily::GetGradGL_GL() const
  {
    return dGL_GL;
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlOptimalFirstFamily::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
    
  inline void HexahedronHcurlOptimalFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronHcurlOptimalFirstFamily::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlOptimalFirstFamily::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalFirstFamily::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalFirstFamily::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlOptimalFirstFamily::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r ) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
		     VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
			  VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
			  VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void HexahedronHcurlOptimalFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
    
}

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_FIRST_FAMILY_INLINE_CXX
#endif
