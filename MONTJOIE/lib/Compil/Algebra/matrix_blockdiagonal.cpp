#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/MatrixBlockDiagonal.cxx"
#endif

namespace Seldon
{
  /* BlockDiagRow */
  
  SELDON_EXTERN template class Matrix<Real_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class Matrix<Complex_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>;

  SELDON_EXTERN template void Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>::Fill(const Real_wp&);
  SELDON_EXTERN template void Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>::Fill(const int&);
  SELDON_EXTERN template Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>& Matrix_BlockDiagonal<Real_wp, General, BlockDiagRow>::operator*=(const Real_wp&);
  SELDON_EXTERN template void Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>::Fill(const int&);
  SELDON_EXTERN template void Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>::Fill(const Complex_wp&);
  SELDON_EXTERN template Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>& Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>::operator*=(const Complex_wp&);
  SELDON_EXTERN template Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>& Matrix_BlockDiagonal<Complex_wp, General, BlockDiagRow>::operator*=(const Real_wp&);
  
  /* BlockDiagRowSym */
  
  SELDON_EXTERN template class Matrix<Real_wp, Symmetric, BlockDiagRowSym>;
  SELDON_EXTERN template class MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>;
  SELDON_EXTERN template class Matrix<Complex_wp, Symmetric, BlockDiagRowSym>;
  SELDON_EXTERN template class MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>;

  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>::Fill(const int&);
  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>::Fill(const Real_wp&);
  SELDON_EXTERN template MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>& MatrixSym_BlockDiagonal<Real_wp, Symmetric, BlockDiagRowSym>::operator*=(const Real_wp&);
  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>::Fill(const int&);
  SELDON_EXTERN template void MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>::Fill(const Complex_wp&);
  SELDON_EXTERN template MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>& MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>::operator*=(const Complex_wp&);
  SELDON_EXTERN template MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>& MatrixSym_BlockDiagonal<Complex_wp, Symmetric, BlockDiagRowSym>::operator*=(const Real_wp&);
  
  /* Functions */
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, General, BlockDiagRow>&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void ConvertToSparse(const Matrix<Real_wp, General, BlockDiagRow>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void ConvertToSparse(const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, General, BlockDiagRow>&, Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, RowMajor>&);
  
  
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void AddMatrix(const Real_wp&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void ConvertToSparse(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void ConvertToSparse(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  
  /* Conversion between sparse matrices and block-diagonal matrices */
  
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&, int, int);
    
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void FillBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, General, BlockDiagRow>&, int, int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);
  SELDON_EXTERN template void ConvertToBlockDiagonal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, General, BlockDiagRow>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, BlockDiagRow>&, Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);
  
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void CopyMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
  /* matrix-vector products */
  
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, General, BlockDiagRow>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, General, BlockDiagRow>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, General, BlockDiagRow>&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, General, BlockDiagRow>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, General, BlockDiagRow>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, General, BlockDiagRow>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, General, BlockDiagRow>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  SELDON_EXTERN template void Mlt(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template void Mlt(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltVector(const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp >&, Vector<Complex_wp >&);  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const Vector<Complex_wp >&, const Complex_wp&, Vector<Complex_wp >&);
  
  /* matrix functions */
  
  SELDON_EXTERN template void MltAddMatrix(const Real_wp&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const Real_wp&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void MltAddMatrix(const Complex_wp&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const Complex_wp&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void GetInverse(Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void GetInverse(Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
  SELDON_EXTERN template void GetCholesky(Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Vector<Real_wp>&);

  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void MltScalar(const Real_wp&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
}
