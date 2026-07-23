#ifndef MONTJOIE_FILE_QUADRANGLE_GAUSS_INLINE_CXX

namespace Montjoie
{

  inline void QuadrangleGauss::SetRadauInterpolationPoints()
  {
  }
  
  inline int QuadrangleGauss::GetNumDofs(int i, int j) const
  {
    return NumFct2D(i, j);
  }

  inline void QuadrangleGauss
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void QuadrangleGauss
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void QuadrangleGauss::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleGauss::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleGauss
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void QuadrangleGauss
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void QuadrangleGauss
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void QuadrangleGauss
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void QuadrangleGauss::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void QuadrangleGauss::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
    
  inline void QuadrangleGauss::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void QuadrangleGauss::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void QuadrangleGauss::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void QuadrangleGauss::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void QuadrangleGauss::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleGauss::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleGauss::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void QuadrangleGauss::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void QuadrangleGauss::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleGauss::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleGauss::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleGauss::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleGauss::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy);
  }
  
  inline void QuadrangleGauss
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy);
  }

  inline void QuadrangleGauss::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
  
  inline void QuadrangleGauss::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureGen(Uh, Vh);
  }
    
  inline void QuadrangleGauss::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleGauss::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleGauss
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void QuadrangleGauss
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleGauss::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
					    VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleGauss
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void QuadrangleGauss
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
    
  inline void QuadrangleGauss
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  const R2_Complex_wp& D, const R2_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }

  inline void QuadrangleGauss
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 2, 2>& C,
			  const R2& D, const R2& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleGauss
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  const Vector<R2>& D, const Vector<R2>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleGauss
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
    
  inline void QuadrangleGauss
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
  
  inline void QuadrangleGauss
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }  

}

#define MONTJOIE_FILE_QUADRANGLE_GAUSS_INLINE_CXX
#endif
