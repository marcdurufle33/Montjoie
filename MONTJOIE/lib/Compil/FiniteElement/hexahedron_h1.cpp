#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/Hexahedron/HexahedronGauss.cxx"
#include "FiniteElement/Hexahedron/HexahedronLobatto.cxx"
#include "FiniteElement/Hexahedron/HexahedronHierarchic.cxx"
#include "FiniteElement/Hexahedron/HexahedronDgGauss.cxx"
#include "FiniteElement/Hexahedron/HexahedronDgLegendre.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template void HexahedronLobatto::ComputeGaussIntegralSurfaceGen(const VectReal_wp & feval,
                                                             VectReal_wp & res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ComputeGradientBoundaryGen(const VectReal_wp& feval, VectReal_wp& res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronLobatto::SolveMassMatrixGen(VectReal_wp&) const;
  
  SELDON_EXTERN template void HexahedronLobatto::SolveCholeskyGen(const SeldonTranspose& TransA, VectReal_wp&) const;
    
  SELDON_EXTERN template void HexahedronLobatto::MltMassMatrixGen(VectReal_wp&) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyChGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyChTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyRhBoundaryGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhBoundaryTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyRhSplitGen(const VectReal_wp& Uh, VectReal_wp& Vh, VectReal_wp&, VectReal_wp&) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyShTransposeGen(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyShGen(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyNablaShTransposeGen(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyNablaShGen(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
                                              VectReal_wp& Vh, int r) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhGen(const VectReal_wp&, VectReal_wp& ) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhTransposeGen(const VectReal_wp&, VectReal_wp& ) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhSplitGen(const VectReal_wp& Uh, VectReal_wp& Vh, VectReal_wp&, VectReal_wp&) const;
    
  SELDON_EXTERN template void HexahedronLobatto::AddConstantMassMatrixGen(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const;
    
  SELDON_EXTERN template void HexahedronLobatto::AddConstantElemMatrixGen(int m, int n, const Real_wp& mass, const TinyMatrix<Real_wp, General, 3, 3>& C,
                                                       const TinyVector<Real_wp, 3>& D, const TinyVector<Real_wp, 3>& E, 
                                                       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;

  SELDON_EXTERN template void HexahedronLobatto::AddVariableMassMatrixGen(int off_row, int off_col,
                                                       const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const;
  
  SELDON_EXTERN template void HexahedronLobatto::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<Real_wp>& A,
                                                       const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
                                                       const Vector<TinyVector<Real_wp, 3> >& D,
                                                       const Vector<TinyVector<Real_wp, 3> >& E,
                                                       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const;


  SELDON_EXTERN template void HexahedronLobatto::ComputeGaussIntegralSurfaceGen(const VectComplex_wp & feval,
                                                             VectComplex_wp & res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ComputeGradientBoundaryGen(const VectComplex_wp& feval, VectComplex_wp& res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronLobatto::SolveMassMatrixGen(VectComplex_wp&) const;
  
  SELDON_EXTERN template void HexahedronLobatto::SolveCholeskyGen(const SeldonTranspose& TransA, VectComplex_wp&) const;
    
  SELDON_EXTERN template void HexahedronLobatto::MltMassMatrixGen(VectComplex_wp&) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyChGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyChTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyRhBoundaryGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyRhBoundaryTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyRhSplitGen(const VectComplex_wp& Uh, VectComplex_wp& Vh, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyShTransposeGen(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyShGen(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;

  SELDON_EXTERN template void HexahedronLobatto::ApplyNablaShTransposeGen(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyNablaShGen(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
                                              VectComplex_wp& Vh, int r) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhGen(const VectComplex_wp&, VectComplex_wp& ) const;
  
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhTransposeGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void HexahedronLobatto::ApplyConstantRhSplitGen(const VectComplex_wp& Uh, VectComplex_wp& Vh, VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void HexahedronLobatto::AddConstantMassMatrixGen(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const;
    
  SELDON_EXTERN template void HexahedronLobatto::AddConstantElemMatrixGen(int m, int n, const Complex_wp& mass, const TinyMatrix<Complex_wp, General, 3, 3>& C,
                                                       const TinyVector<Complex_wp, 3>& D, const TinyVector<Complex_wp, 3>& E, 
                                                       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;

  SELDON_EXTERN template void HexahedronLobatto::AddVariableMassMatrixGen(int off_row, int off_col,
                                                       const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const;
  
  SELDON_EXTERN template void HexahedronLobatto::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<Complex_wp>& A,
                                                       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
                                                       const Vector<TinyVector<Complex_wp, 3> >& D,
                                                       const Vector<TinyVector<Complex_wp, 3> >& E,
                                                       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const;
 

  SELDON_EXTERN template void HexahedronGauss::ComputeGaussIntegralSurfaceGen(const VectReal_wp & feval,
                                                                              VectReal_wp& res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronGauss::ComputeNodalValuesGen(const VectReal_wp& Un, VectReal_wp & Unode) const;    
  
  SELDON_EXTERN template void HexahedronGauss::ComputeValueBoundaryGen(const VectReal_wp& Un, VectReal_wp & Unode, int num_loc) const;
    
  SELDON_EXTERN template void HexahedronGauss::ComputeGradientBoundaryGen(const VectReal_wp& Un, VectReal_wp & Unode, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronGauss::SolveMassMatrixGen(VectReal_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::SolveCholeskyGen(const SeldonTranspose& TransA, VectReal_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::MltMassMatrixGen(VectReal_wp& Vh) const;    
  
  SELDON_EXTERN template void HexahedronGauss::ApplyChGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
  
  SELDON_EXTERN template void HexahedronGauss::ApplyChTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyRhGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyRhTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhSplitGen(const VectReal_wp& Uh, VectReal_wp& Vh, VectReal_wp&, VectReal_wp&) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhGen(const VectReal_wp&, VectReal_wp& ) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhTransposeGen(const VectReal_wp&, VectReal_wp& ) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhSplitGen(const VectReal_wp& Uh, VectReal_wp& Vh, VectReal_wp&, VectReal_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhQuadratureGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhQuadratureTransposeGen(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyShTransposeGen(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyShGen(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyNablaShTransposeGen(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyNablaShGen(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
                       VectReal_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::AddConstantMassMatrixGen(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const;
    
  SELDON_EXTERN template void HexahedronGauss::AddConstantElemMatrixGen(int m, int n, const Real_wp& mass, const TinyMatrix<Real_wp, General, 3, 3>& C,
                                const TinyVector<Real_wp, 3>& D, const TinyVector<Real_wp, 3>& E, 
                                const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
  
  SELDON_EXTERN template void HexahedronGauss::AddVariableMassMatrixGen(int off_row, int off_col,
                                const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const;

  SELDON_EXTERN template void HexahedronGauss::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<Real_wp>& A,
                                const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
                                const Vector<TinyVector<Real_wp, 3> >& D,
                                const Vector<TinyVector<Real_wp, 3> >& E,
                                const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const;

  SELDON_EXTERN template void HexahedronGauss::ComputeGaussIntegralSurfaceGen(const VectComplex_wp & feval,
                                                                              VectComplex_wp& res, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronGauss::ComputeNodalValuesGen(const VectComplex_wp& Un, VectComplex_wp & Unode) const;    
  
  SELDON_EXTERN template void HexahedronGauss::ComputeValueBoundaryGen(const VectComplex_wp& Un, VectComplex_wp & Unode, int num_loc) const;
    
  SELDON_EXTERN template void HexahedronGauss::ComputeGradientBoundaryGen(const VectComplex_wp& Un, VectComplex_wp & Unode, int num_loc) const;
  
  SELDON_EXTERN template void HexahedronGauss::SolveMassMatrixGen(VectComplex_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::SolveCholeskyGen(const SeldonTranspose& TransA, VectComplex_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::MltMassMatrixGen(VectComplex_wp& Vh) const;    
  
  SELDON_EXTERN template void HexahedronGauss::ApplyChGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
  
  SELDON_EXTERN template void HexahedronGauss::ApplyChTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyRhGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyRhTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhSplitGen(const VectComplex_wp& Uh, VectComplex_wp& Vh, VectComplex_wp&, VectComplex_wp&) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhTransposeGen(const VectComplex_wp&, VectComplex_wp& ) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyConstantRhSplitGen(const VectComplex_wp& Uh, VectComplex_wp& Vh, VectComplex_wp&, VectComplex_wp&) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhQuadratureGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyRhQuadratureTransposeGen(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyShTransposeGen(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyShGen(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;

  SELDON_EXTERN template void HexahedronGauss::ApplyNablaShTransposeGen(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::ApplyNablaShGen(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
                       VectComplex_wp& Vh, int r) const;
    
  SELDON_EXTERN template void HexahedronGauss::AddConstantMassMatrixGen(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const;
    
  SELDON_EXTERN template void HexahedronGauss::AddConstantElemMatrixGen(int m, int n, const Complex_wp& mass, const TinyMatrix<Complex_wp, General, 3, 3>& C,
                                const TinyVector<Complex_wp, 3>& D, const TinyVector<Complex_wp, 3>& E, 
                                const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
  
  SELDON_EXTERN template void HexahedronGauss::AddVariableMassMatrixGen(int off_row, int off_col,
                                const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const;

  SELDON_EXTERN template void HexahedronGauss::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<Complex_wp>& A,
                                const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
                                const Vector<TinyVector<Complex_wp, 3> >& D,
                                const Vector<TinyVector<Complex_wp, 3> >& E,
                                const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const;
  
}
