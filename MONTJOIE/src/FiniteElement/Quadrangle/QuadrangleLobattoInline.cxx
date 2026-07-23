#ifndef MONTJOIE_FILE_QUADRANGLE_LOBATTO_INLINE_CXX

namespace Montjoie
{
  
  //! returns 2-D quadrature points
  inline const Vector<R2>& QuadrangleLobatto::PointsQuadND() const
  {
    return points2d_all;
  }

  inline void QuadrangleLobatto::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }

  inline void QuadrangleLobatto::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }

  inline void QuadrangleLobatto::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void QuadrangleLobatto::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
      
  inline void QuadrangleLobatto::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void QuadrangleLobatto::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void QuadrangleLobatto::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void QuadrangleLobatto::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleLobatto::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void QuadrangleLobatto::ApplyRhBoundary(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhBoundaryGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyRhBoundary(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhBoundaryGen(Vh, Uh);
  }

  inline void QuadrangleLobatto::ApplyRhBoundaryTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhBoundaryTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleLobatto::ApplyRhBoundaryTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhBoundaryTransposeGen(Uh, Vh);
  }

  inline void QuadrangleLobatto::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void QuadrangleLobatto::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }

  inline void QuadrangleLobatto
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy);
  }

  inline void QuadrangleLobatto
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy);
  }

  inline void QuadrangleLobatto::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }

  inline void QuadrangleLobatto::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }

  inline void QuadrangleLobatto
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void QuadrangleLobatto
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void QuadrangleLobatto
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleLobatto
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }

  inline void QuadrangleLobatto
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void QuadrangleLobatto
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleLobatto
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  inline void QuadrangleLobatto
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void QuadrangleLobatto
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
    
  inline void QuadrangleLobatto
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  const R2_Complex_wp& D, const R2_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }

  inline void QuadrangleLobatto
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 2, 2>& C,
			  const R2& D, const R2& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleLobatto
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  const Vector<R2>& D, const Vector<R2>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleLobatto
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
    
  inline void QuadrangleLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
  
  inline void QuadrangleLobatto
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }  
    
  inline void QuadrangleLobatto
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void QuadrangleLobatto
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
    
}
  
#define MONTJOIE_FILE_QUADRANGLE_LOBATTO_INLINE_CXX
#endif
