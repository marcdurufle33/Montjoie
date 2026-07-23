#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#ifdef SELDON_WITH_MPI

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/DistributedBlockDiagonalMatrix.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class DistributedMatrix_BlockDiag<Real_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class DistributedMatrix_BlockDiag<Real_wp, Symmetric, BlockDiagRowSym>;
  SELDON_EXTERN template class DistributedMatrix_BlockDiag<Complex_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class DistributedMatrix_BlockDiag<Complex_wp, Symmetric, BlockDiagRowSym>;

  SELDON_EXTERN template class DistributedMatrix<Real_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>;
  SELDON_EXTERN template class DistributedMatrix<Complex_wp, General, BlockDiagRow>;
  SELDON_EXTERN template class DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>;
  
  SELDON_EXTERN template void DistributedMatrix_BlockDiag<Real_wp, General, BlockDiagRow>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void DistributedMatrix_BlockDiag<Real_wp, Symmetric, BlockDiagRowSym>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void DistributedMatrix_BlockDiag<Complex_wp, General, BlockDiagRow>::SetPattern(const Vector<IVect>&);
  SELDON_EXTERN template void DistributedMatrix_BlockDiag<Complex_wp, Symmetric, BlockDiagRowSym>::SetPattern(const Vector<IVect>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					DistributedMatrix<Real_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					DistributedMatrix<Complex_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					DistributedMatrix<Complex_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void AddMatrix(const Real_wp&, const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void AddMatrix(const Complex_wp&, const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					      DistributedMatrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					      DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					      DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					      DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					      DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void ConvertToSparse(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					      DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&,
						     DistributedMatrix<Real_wp, General, BlockDiagRow>&, int , int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&,
						     DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&, int , int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&,
						     DistributedMatrix<Complex_wp, General, BlockDiagRow>&, int, int);
  
  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&,
						     DistributedMatrix<Complex_wp, General, BlockDiagRow>&, int, int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
						     DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
						     DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);

  SELDON_EXTERN template void ConvertToBlockDiagonal(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&,
						     DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
				   DistributedMatrix<Real_wp, General, ArrayRowSparse>&);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
				   DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
				   DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
				   DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
				   DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
				   DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);


  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>&,
				   DistributedMatrix<Real_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowSparse>&,
				   DistributedMatrix<Complex_wp, General, BlockDiagRow>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&,
				   DistributedMatrix<Complex_wp, General, BlockDiagRow>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&,
				   DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&,
				   DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);
  
  SELDON_EXTERN template void Copy(const DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&,
				   DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void Mlt(const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
				  Vector<Real_wp>& X, bool);

  SELDON_EXTERN template void Mlt(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
				  Vector<Complex_wp>& X, bool);
  
  SELDON_EXTERN template void Mlt(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& A,
				  Vector<Real_wp>& X, bool);

  SELDON_EXTERN template void Mlt(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>& A,
				  Vector<Complex_wp>& X, bool);

  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& A);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Real_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>& A);
  SELDON_EXTERN template void MltScalar(const Real_wp&, DistributedMatrix<Complex_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>& A);
  SELDON_EXTERN template void MltScalar(const Complex_wp&, DistributedMatrix<Complex_wp, General, BlockDiagRow>& A);


  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Real_wp>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					const Vector<Real_wp>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltVector(const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Real_wp>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					const Vector<Real_wp>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltVector(const SeldonTranspose&,
					const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					const Vector<Complex_wp>&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Real_wp&,
					   const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&,
					   const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Real_wp&,
					   const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&,
					   const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  
  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Real_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					   const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Real_wp, General, BlockDiagRow>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&,
					   const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&,
					   const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  
  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);

  SELDON_EXTERN template void MltAddVector(const Complex_wp&, const SeldonTranspose&,
					   const DistributedMatrix<Complex_wp, General, BlockDiagRow>&,
					   const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&, bool);
  
  SELDON_EXTERN template void GetInverse(DistributedMatrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void GetInverse(DistributedMatrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void GetInverse(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void GetInverse(DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void GetCholesky(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void SolveCholesky(const SeldonTranspose&,
					    const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					    Vector<Real_wp>&);

  SELDON_EXTERN template void MltCholesky(const SeldonTranspose&,
					  const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&,
					  Vector<Real_wp>&);

  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& A);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void EraseRow(const IVect&, DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>& A);

  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& A);
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, General, BlockDiagRow>& A);
  SELDON_EXTERN template void EraseCol(const IVect&, DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>& A);
  
}

#endif
