#ifndef MONTJOIE_FILE_HEXAHEDRON_LOBATTO_INLINE_CXX

namespace Montjoie
{

  //! returns 3-D quadrature points
  inline const Vector<R3>& HexahedronLobatto::PointsQuadND() const
  {
    return points3d_all;
  }

  inline void HexahedronLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void HexahedronLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void HexahedronLobatto
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void HexahedronLobatto
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void HexahedronLobatto::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void HexahedronLobatto::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }

  inline void HexahedronLobatto::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void HexahedronLobatto::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void HexahedronLobatto::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void HexahedronLobatto::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void HexahedronLobatto::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void HexahedronLobatto::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void HexahedronLobatto::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void HexahedronLobatto::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void HexahedronLobatto::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronLobatto::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void HexahedronLobatto
  ::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void HexahedronLobatto
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }

  inline void HexahedronLobatto::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronLobatto::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void HexahedronLobatto::ApplyRhBoundary(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhBoundaryGen(Vh, Uh);
  }

  inline void HexahedronLobatto::ApplyRhBoundary(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhBoundaryGen(Vh, Uh);
  }

  inline void HexahedronLobatto::ApplyRhBoundaryTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhBoundaryTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronLobatto::ApplyRhBoundaryTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhBoundaryTransposeGen(Uh, Vh);
  }

  inline void HexahedronLobatto
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronLobatto
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void HexahedronLobatto
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronLobatto
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void HexahedronLobatto::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				    VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void HexahedronLobatto::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
				       VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void HexahedronLobatto::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void HexahedronLobatto::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
    
  inline void HexahedronLobatto
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void HexahedronLobatto
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			 VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
    
  inline void HexahedronLobatto::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
  
  inline void HexahedronLobatto::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
    
  inline void HexahedronLobatto
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronLobatto
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void HexahedronLobatto
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 3, 3>& C,
			  const R3& D, const R3& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }  
  
  inline void HexahedronLobatto
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  const R3_Complex_wp& D, const R3_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
    
  inline void HexahedronLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void HexahedronLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void HexahedronLobatto
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  const Vector<R3>& D, const Vector<R3>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
    
  inline void HexahedronLobatto
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
}

#define MONTJOIE_FILE_HEXAHEDRON_LOBATTO_INLINE_CXX
#endif

