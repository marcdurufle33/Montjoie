#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_INLINE_CXX

namespace Montjoie
{

#ifdef MONTJOIE_FAST_HEXA 
  inline int HexahedronDgLegendre::GetNbPointsUsedForSh() const
  {
    return 3*(order+1)*(order+2);
  }
  
  inline int HexahedronDgLegendre::GetNbPointsUsedForFaceSh(int num_loc) const
  {
    return (order+1)*(order+2)/2;
  }
  
  inline int HexahedronDgLegendre::GetOffsetSh(int num_loc) const
  {
    return num_loc*(order+1)*(order+2)/2; 
  }
#endif

  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector1>
  inline void HexahedronDgLegendre::SolveMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity
  }
  
  
  //! x is overwritten by M x where M is the mass matrix
  template<class Vector1>
  inline void HexahedronDgLegendre::MltMassMatrixGen(Vector1& x) const
  {
    // mass matrix equal to identity
  }
  

  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix
  template<class Vector1>
  inline void HexahedronDgLegendre::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
  }
  
  
  inline void HexahedronDgLegendre
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void HexahedronDgLegendre::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }

  inline void HexahedronDgLegendre
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronDgLegendre
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void HexahedronDgLegendre
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
   
  inline void HexahedronDgLegendre
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void HexahedronDgLegendre
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronDgLegendre
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronDgLegendre
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void HexahedronDgLegendre
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  inline void HexahedronDgLegendre::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronDgLegendre::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }

  inline void HexahedronDgLegendre::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void HexahedronDgLegendre::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void HexahedronDgLegendre::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronDgLegendre::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronDgLegendre::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronDgLegendre::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronDgLegendre::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronDgLegendre::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void HexahedronDgLegendre::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronDgLegendre::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronDgLegendre
  ::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void HexahedronDgLegendre
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }

  inline void HexahedronDgLegendre::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronDgLegendre::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronDgLegendre::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void HexahedronDgLegendre::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
    
  inline void HexahedronDgLegendre
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void HexahedronDgLegendre
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			 VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
    
  inline void HexahedronDgLegendre::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronDgLegendre::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
    
  inline void HexahedronDgLegendre
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronDgLegendre
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void HexahedronDgLegendre
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronDgLegendre
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronDgLegendre
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronDgLegendre
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void HexahedronDgLegendre::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				    VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronDgLegendre
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronDgLegendre
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronDgLegendre
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronDgLegendre
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 3, 3>& C,
			  const R3& D, const R3& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }  
  
  inline void HexahedronDgLegendre
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  const R3_Complex_wp& D, const R3_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
    
  inline void HexahedronDgLegendre
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void HexahedronDgLegendre
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void HexahedronDgLegendre
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  const Vector<R3>& D, const Vector<R3>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
    
  inline void HexahedronDgLegendre
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }

}

#define MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_INLINE_CXX
#endif
