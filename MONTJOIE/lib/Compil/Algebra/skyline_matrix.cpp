#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/GeneralSkyLineMatrix.cxx"
#include "Algebra/SymmetricSkyLineMatrix.cxx"
#endif

namespace Seldon
{
  /* non-symmetric skyline matrices */
  
  SELDON_EXTERN template class Matrix<Real_wp, General, RowSkyLine>;
  SELDON_EXTERN template class GeneralSkyLineMatrix<Real_wp, General, RowSkyLine>;

  SELDON_EXTERN template void GeneralSkyLineMatrix<Real_wp, General, RowSkyLine>::Fill(const int&);
  SELDON_EXTERN template void GeneralSkyLineMatrix<Real_wp, General, RowSkyLine>::Fill(const Real_wp&);
  
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, RowSkyLine>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, RowSkyLine>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, RowSkyLine>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, RowSkyLine>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, General, RowSkyLine>&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, General, RowSkyLine>&, Vector<Real_wp>&);
  SELDON_EXTERN template void SolveLuVector(const Matrix<Real_wp, General, RowSkyLine>&, Vector<Real_wp>&);
  
  /* symmetric skyline matrices */
  
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, SymColSkyLine>;
  SELDON_EXTERN template class SymmetricSkyLineMatrix<Real_wp, Symmetric, SymColSkyLine>;

  SELDON_EXTERN template void SymmetricSkyLineMatrix<Real_wp, Symmetric, SymColSkyLine>::Fill(const int&);
  SELDON_EXTERN template void SymmetricSkyLineMatrix<Real_wp, Symmetric, SymColSkyLine>::Fill(const Real_wp&);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, SymColSkyLine>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, SymColSkyLine>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, SymColSkyLine>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, SymColSkyLine>&);
  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, SymColSkyLine>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, SymColSkyLine>&, Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(Matrix<Real_wp, Symmetric, SymColSkyLine>&);
  SELDON_EXTERN template void SolveLuVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, SymColSkyLine>&, Vector<Real_wp>&);
  
}
