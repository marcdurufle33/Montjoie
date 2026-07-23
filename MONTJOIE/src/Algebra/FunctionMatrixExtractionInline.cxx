#ifndef SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_INLINE_CXX

namespace Seldon
{

  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T1>, Prop, ArrayRowSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<complex<T1>, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T1>, Prop, ArrayRowSymComplexSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  template<class T1, class Prop, class Allocator,
	   class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T1>, Prop, ArrayRowComplexSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T1, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T>, Prop, BlockDiagRowSym, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V) 
  { abort(); }
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T>, Prop, BlockDiagRow, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V) 
  { abort(); }
  
  template<class T, class Prop, class Allocator,
           class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T>, Prop, DiagonalRow, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  // non-implemented functions
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T>, Prop, RowSymComplexSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<T, Prop, RowSymSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }

  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<complex<T>, Prop, RowComplexSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }
  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  inline void GetCol(const Matrix<T, Prop, RowSparse, Allocator>& A,
		     const IVect& col_number,
		     Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  { abort(); }

}

#define SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_INLINE_CXX
#endif
