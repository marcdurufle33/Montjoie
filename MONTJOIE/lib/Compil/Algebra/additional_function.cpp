#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/Eigenvalue.cxx"
#include "Algebra/FactorisationLU.cxx"
#include "Algebra/FunctionMatrixExtraction.cxx"
#include "matrix_sparse/Functions_MatrixArray.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&, const Matrix<Real_wp, Symmetric, RowSymPacked>&, Vector<Complex_wp>&, LapackInfo&);
    
#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void ExtractDistributedInfo(const DistributedMatrix_Base<Real_wp>&, const Vector<int>&, const Vector<int>&, DistributedMatrixIntegerArray&);
  SELDON_EXTERN template void ExtractDistributedInfo(const DistributedMatrix_Base<Complex_wp>&, const Vector<int>&, const Vector<int>&, DistributedMatrixIntegerArray&);
#endif
  
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, ArrayRowSparse>&, const Vector<int>&, Vector<Vector<Real_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, Vector<Vector<Real_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, ArrayRowSparse>&, const Vector<int>&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const Vector<int>&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, const Vector<int>&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const Vector<int>&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);

  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, General, BlockDiagRow>&, const IVect&, Vector<Vector<Real_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, General, BlockDiagRow>&, const IVect&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, BlockDiagRowSym>&, const IVect&, Vector<Vector<Real_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&, const IVect&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Real_wp, Symmetric, DiagonalRow>&, const IVect&, Vector<Vector<Real_wp, VectSparse>, VectSparse>&);
  SELDON_EXTERN template void GetCol(const Matrix<Complex_wp, Symmetric, DiagonalRow>&, const IVect&, Vector<Vector<Complex_wp, VectSparse>, VectSparse>&);

  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void EraseRow(const IVect&, Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void EraseCol(const IVect&, Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrowCol>&);

  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void GetColSum(Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrowCol>&);

  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, General, ArrowCol>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, ArrowCol>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Real_wp, General, BandedCol>&);
  SELDON_EXTERN template void GetRowColSum(Vector<Real_wp>&, Vector<Real_wp>&, const Matrix<Complex_wp, General, BandedCol>&);
  
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Real_wp, General, ArrowCol>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void ScaleLeftMatrix(Matrix<Complex_wp, General, ArrowCol>&, const Vector<Real_wp>&);

  SELDON_EXTERN template void GetSubMatrix(const Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, int, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, int, int, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, int, int, int, int, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, int, Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, int, int, Matrix<Complex_wp, General, ArrayRowSparse>&);

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, int, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, int, int, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, int, int, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, int, int, int, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, int, int, int, int, DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, int, int, int, DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void GetSubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, int, int, int, int, DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);

#endif
  
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, Symmetric, RowSymComplexSparse>&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, Symmetric, RowSymSparse>&, Matrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, General, ArrayRowComplexSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, General, ArrayRowSparse>&, Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, General, RowComplexSparse>&, Matrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void CopyReal(const Matrix<Complex_wp, General, RowSparse>&, Matrix<Real_wp, General, RowSparse>&);

  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Real_wp, General, RowSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, General, RowSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&, const IVect&);

  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, const IVect&);

  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, General, RowComplexSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&, const IVect&);

  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, const IVect&);
  SELDON_EXTERN template void CompressMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, const IVect&);

  SELDON_EXTERN template void GetExponential(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetExponential(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void GetExponential(Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void GetExponential(Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void GetExponential(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetExponential(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void GetExponential(Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void GetExponential(Matrix<Complex_wp, Symmetric, RowSym>&);

  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Real_wp, General, RowMajor>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Real_wp, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Real_wp, Symmetric, RowSymPacked>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Real_wp, Symmetric, RowSym>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Complex_wp, General, RowMajor>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Complex_wp, General, ColMajor>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Complex_wp, Symmetric, RowSymPacked>&, bool);
  SELDON_EXTERN template void GetAbsoluteValue(Matrix<Complex_wp, Symmetric, RowSym>&, bool);

  SELDON_EXTERN template void GetSquareRoot(Matrix<Real_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Real_wp, General, ColMajor>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Real_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Real_wp, Symmetric, RowSym>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Complex_wp, General, RowMajor>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Complex_wp, General, ColMajor>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Complex_wp, Symmetric, RowSymPacked>&);
  SELDON_EXTERN template void GetSquareRoot(Matrix<Complex_wp, Symmetric, RowSym>&);

  SELDON_EXTERN template void ExtractSubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, const IVect&, const IVect&, Matrix<Real_wp, General, ArrayRowSparse>& Asub);
  SELDON_EXTERN template void ExtractSubMatrix(const Matrix<Complex_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, const IVect&, const IVect&, Matrix<Complex_wp, General, ArrayRowSparse>& Asub);
  SELDON_EXTERN template void ExtractSubMatrix(const Matrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, const IVect&, const IVect&, Matrix<Real_wp, Symmetric, ArrayRowSymSparse>& Asub);

  SELDON_EXTERN template void ExtractSubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, DistributedMatrixIntegerArray&, DistributedMatrix<Real_wp, General, ArrayRowSparse>& Asub);
  SELDON_EXTERN template void ExtractSubMatrix(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, DistributedMatrixIntegerArray&, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Asub);
  SELDON_EXTERN template void ExtractSubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&, const IVect&, const IVect&, DistributedMatrixIntegerArray&, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Asub);

}
