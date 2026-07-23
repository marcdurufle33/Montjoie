#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
    
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
		     VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
			  VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
			  VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void HexahedronHcurlOptimalHpFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
}

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_INLINE_CXX
#endif
