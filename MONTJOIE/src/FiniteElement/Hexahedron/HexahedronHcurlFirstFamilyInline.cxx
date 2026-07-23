#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{
  inline const Matrix<int>& HexahedronHcurlFirstFamily::GetCoordinateDofs() const
  {
    return CoordinateDofs;
  }
  
  inline const Array3D<int>& HexahedronHcurlFirstFamily::GetNumDofsX() const
  {
    return NumDofs_X;
  }
  
  inline const Array3D<int>& HexahedronHcurlFirstFamily::GetNumDofsY() const
  {
    return NumDofs_Y;
  }
  
  inline const Array3D<int>& HexahedronHcurlFirstFamily::GetNumDofsZ() const
  {
    return NumDofs_Z;
  }
  
  inline const Matrix<Real_wp>& HexahedronHcurlFirstFamily::GetGradGL_GL() const
  {
    return dGL_GL;
  }
  
  inline const Matrix<Real_wp>& HexahedronHcurlFirstFamily::GetGradGL_G() const
  {
    return dGL_G;
  }
  
  inline const Matrix<Real_wp>& HexahedronHcurlFirstFamily::GetValGauss_GL() const
  {
    return G_GL;
  }

  inline int HexahedronHcurlFirstFamily::GetCoordinateDof(int i) const
  {
    return this->CoordinateDofs(i, 3);
  }

  inline int HexahedronHcurlFirstFamily::GetPointNode_FromDofNumber(int i) const
  {
    return i;
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlFirstFamily
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void HexahedronHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlFirstFamily
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlFirstFamily::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlFirstFamily::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
    
  inline void HexahedronHcurlFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlFirstFamily::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlFirstFamily::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronHcurlFirstFamily::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void HexahedronHcurlFirstFamily
  ::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlFirstFamily::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlFirstFamily::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlFirstFamily::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlFirstFamily::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlFirstFamily::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
		     VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
			  VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
			  VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }

  inline void HexahedronHcurlFirstFamily
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void HexahedronHcurlFirstFamily
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }

  inline void HexahedronHcurlFirstFamily
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			       VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableStiffnessMatrixGen(off_row, off_col, C, mat);
  }
    
  inline void HexahedronHcurlFirstFamily
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			       VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableStiffnessMatrixGen(off_row, off_col, C, mat);
  }
      
}

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_INLINE_CXX
#endif
