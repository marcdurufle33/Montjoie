#ifndef SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_HXX

namespace Seldon
{

#ifdef SELDON_WITH_MPI
  void ExtractDistributedInfo(const DistributedMatrixIntegerArray&,
			      const Vector<int>& num, const Vector<int>& index,
			      DistributedMatrixIntegerArray&);

  template<class T>
  void ExtractDistributedInfo(const DistributedMatrix_Base<T>&,
			      const Vector<int>& num, const Vector<int>& index,
			      DistributedMatrixIntegerArray&);
#endif
  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSparse, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowComplexSparse, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, BlockDiagRowSym, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);  
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, BlockDiagRow, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, DiagonalRow, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
 
  // complex -> real : forbidden case
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T1>, Prop, ArrayRowSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<complex<T1>, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T1>, Prop, ArrayRowSymComplexSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T1>, Prop, ArrayRowComplexSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T>, Prop, BlockDiagRowSym, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T>, Prop, BlockDiagRow, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V);  
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T>, Prop, DiagonalRow, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V);
  
  // non-implemented functions
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T>, Prop, RowSymComplexSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, RowSymSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<complex<T>, Prop, RowComplexSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, RowSparse, Allocator>& A,
	      const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V);

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, DiagonalRow, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, DiagonalRow, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BandedCol, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, BandedCol, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRow, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRow, Allocator>& A);

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRowSym, Allocator>& A);
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRowSym, Allocator>& A);

  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 ArrowCol, Allocator> & mat_direct);

  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General,
                 ArrowCol, Allocator> & mat_direct);
  
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
		    BandedCol, Allocator> & A);
  
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row,
                    Vector<T>& sum_col,
                    const Matrix<Complexe, General,
		    ArrowCol, Allocator> & A);
  
  template<class T, class Allocator, class T0>
  void ScaleLeftMatrix(Matrix<T, General, ArrowCol, Allocator>& A,
                       const Vector<T0>& coef_row);

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSymComplexSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSymSparse, Allocator2>& B);
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSymSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSymSparse, Allocator2>& B);

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSymComplexSparse, Allocator1>& A,
		Matrix<T, Prop, RowSymSparse, Allocator2>& B);
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSymSparse, Allocator1>& A,
		Matrix<T, Prop, RowSymSparse, Allocator2>& B);

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowComplexSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSparse, Allocator2>& B);
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSparse, Allocator2>& B);

  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowComplexSparse, Allocator1>& A,
		Matrix<T, Prop, RowSparse, Allocator2>& B);
  
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSparse, Allocator1>& A,
		Matrix<T, Prop, RowSparse, Allocator2>& B);

  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void CopyReal(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
		Matrix<T2, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m, int n, Matrix<T, General, ArrayRowSparse, Allocator>& B);

  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, General, ArrayRowSparse, Allocator>& B);

  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B);

  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
                    int m, int n, Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B);
  
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B);

  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
                    int m, int n, Matrix<T, Prop, Storage, Allocator>& B);

  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Prop, Storage, Allocator>& B);

#ifdef SELDON_WITH_MPI
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m, int n,
		    DistributedMatrix<T, Prop, Storage, Allocator>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    DistributedMatrix<T, Prop, Storage, Allocator>& B);
#endif

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSymSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
		      const IVect& IndexRow);

    template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowComplexSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSymComplexSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
		      const IVect& IndexRow);

  template<class T>
  void ExtractSubMatrix(const Matrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			const IVect& col_num, const IVect& index_col,
			Matrix<T, General, ArrayRowSparse>& Asub);

  template<class T>
  void ExtractSubMatrix(const Matrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			const IVect& col_num, const IVect& index_col,
			Matrix<T, Symmetric, ArrayRowSymSparse>& Asub);

  template<class T>
  void ExtractSubMatrix(const DistributedMatrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			DistributedMatrixIntegerArray& Asub_info,
			DistributedMatrix<T, General, ArrayRowSparse>& Asub);

  template<class T>
  void ExtractSubMatrix(const DistributedMatrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			DistributedMatrixIntegerArray& Asub_info,
			DistributedMatrix<T, Symmetric, ArrayRowSymSparse>& Asub);

}

#define SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_HXX
#endif

