#ifndef MONTJOIE_FILE_WEDGE_CLASSICAL_INLINE_CXX

namespace Montjoie
{

  inline void WedgeClassical
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void WedgeClassical
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void WedgeClassical
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void WedgeClassical::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }

  inline void WedgeClassical
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
   
  inline void WedgeClassical
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void WedgeClassical
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void WedgeClassical
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void WedgeClassical::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void WedgeClassical::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }

  inline void WedgeClassical::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void WedgeClassical::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void WedgeClassical
  ::ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void WedgeClassical
  ::ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyRhSplitGen(Uh, Vx, Vy, Vz);
  }

  inline void WedgeClassical::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeClassical::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeClassical::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void WedgeClassical::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  inline void WedgeClassical::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void WedgeClassical::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void WedgeClassical
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeClassical
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void WedgeClassical
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeClassical
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeClassical
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void WedgeClassical
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void WedgeClassical::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				    VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void WedgeClassical::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
				    VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void WedgeClassical::ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
  
  inline void WedgeClassical::ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhGen(Uh, Vh);
  }
    
  inline void WedgeClassical
  ::ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
  
  inline void WedgeClassical
  ::ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			 VectComplex_wp& Vy, VectComplex_wp& Vz) const
  {
    ApplyConstantRhSplitGen(Uh, Vx, Vy, Vz);
  }
    
  inline void WedgeClassical::ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
  
  inline void WedgeClassical::ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyConstantRhTransposeGen(Uh, Vh);
  }
    
  inline void WedgeClassical
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void WedgeClassical
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void WedgeClassical
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 3, 3>& C,
			  const R3& D, const R3& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }  
  
  inline void WedgeClassical
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 3, 3>& C,
			  const R3_Complex_wp& D, const R3_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
    
  inline void WedgeClassical
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A,
			  VirtualMatrix<Real_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void WedgeClassical
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
    
  inline void WedgeClassical
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			  const Vector<R3>& D, const Vector<R3>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
    
  inline void WedgeClassical
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			  const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
}

#define MONTJOIE_FILE_WEDGE_CLASSICAL_INLINE_CXX
#endif
