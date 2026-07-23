#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_INLINE_CXX

namespace Montjoie
{

  inline const Real_wp& HexahedronHcurlLobatto::dGL_GL(int i, int j) const
  {
    return lob_basis.GradPhi(i, j);
  }    
  
  inline const Matrix<Real_wp>& HexahedronHcurlLobatto::GradPhi() const
  {
    return lob_basis.GradPhi();
  }

  inline int HexahedronHcurlLobatto::GetCoordinateDof(int i) const
  {
    return this->CoordinateDofs(i, 3);
  }

  inline int HexahedronHcurlLobatto::GetPointNode_FromDofNumber(int i) const
  {
    return Dof_to_Node(i);
  }

  inline void HexahedronHcurlLobatto::GetDofNumber_FromPointNode(int i, int& ix, int& iy, int& iz) const
  {
    ix = Node_to_DofX(i); iy = Node_to_DofY(i); iz = Node_to_DofZ(i);
  }

  inline void HexahedronHcurlLobatto
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  inline void HexahedronHcurlLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlLobatto
  ::ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronHcurlLobatto
  ::ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeCurlBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronHcurlLobatto::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlLobatto::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
    
  inline void HexahedronHcurlLobatto
  ::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlLobatto
  ::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void HexahedronHcurlLobatto::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronHcurlLobatto::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronHcurlLobatto::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlLobatto::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronHcurlLobatto::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronHcurlLobatto
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void HexahedronHcurlLobatto::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlLobatto::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlLobatto
  ::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronHcurlLobatto
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r ) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronHcurlLobatto
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
		     VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
			  VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
			  VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void HexahedronHcurlLobatto
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
			       VirtualMatrix<Real_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlLobatto
  ::AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
			       VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantStiffnessMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronHcurlLobatto
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void HexahedronHcurlLobatto
  ::AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, C, A);
  }
  
  inline void HexahedronHcurlLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void HexahedronHcurlLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
}

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_INLINE_CXX
#endif
