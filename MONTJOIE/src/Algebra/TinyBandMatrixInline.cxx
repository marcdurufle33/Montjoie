#ifndef SELDON_FILE_TINY_BAND_MATRIX_INLINE_CXX

#include "TinyBandMatrix.hxx"

namespace Seldon
{
  
  /******************
   * TinyBandMatrix *
   ******************/
  
  
  //! returns the number of elements stored
  template<class T, int d>
  inline int TinyBandMatrix<T, d>::GetDataSize() const
  {
    return (2*d+1)*diag_.GetM();
  }
  

  //! returns the memory used by the matrix in bytes
  template<class T, int d>
  inline size_t TinyBandMatrix<T, d>::GetMemorySize() const
  {
    return diag_.GetMemorySize();
  }
  
  
  //! adding a term val in a_{i,j}
  template<class T, int d>    
  inline void TinyBandMatrix<T, d>::AddInteraction(int i, int j, const T& val)
  {
    int k = j-i+d;
    if ((k >= 0)&&(k <= 2*d))
      diag_(j)(k) += val;
    else
      {
	cout << "Matrix not compatible " << endl;
	abort();
      }
  }

  
  //! returns entry (i, j) of the matrix
  template<class T, int d>    
  inline T TinyBandMatrix<T, d>::operator()(int i, int j) const
  {
    T zero; SetComplexZero(zero);
    int k = j-i+d;
    if ((k >= 0)&&(k <= 2*d))
      return diag_(j)(k);
    
    return zero;
  }
  
  
  //! modification of entry (i, j)
  template<class T, int d>    
  inline T& TinyBandMatrix<T, d>::Get(int i, int j)
  {
    int k = j-i+d;
    if ((k >= 0)&&(k <= 2*d))
      return diag_(j)(k);
    else
      {
	cout << "Element not accessible " << endl;
	abort();
      }
  }


  //! returns entry (i, j)
  template<class T, int d>    
  inline const T& TinyBandMatrix<T, d>::Get(int i, int j) const
  {
    int k = j-i+d;
    if ((k >= 0)&&(k <= 2*d))
      return diag_(j)(k);
    else
      {
	cout << "Element not accessible " << endl;
	abort();
      }
  }


  //! modification of entry (i, j)
  template<class T, int d>    
  inline T& TinyBandMatrix<T, d>::Val(int i, int j)
  {
    return Get(i, j);
  }


  //! returns entry (i, j)
  template<class T, int d>    
  inline const T& TinyBandMatrix<T, d>::Val(int i, int j) const
  {
    return Get(i, j);
  }
  

  //! modification of entry (i, j)
  template<class T, int d>    
  inline void TinyBandMatrix<T, d>::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! modification of entry (i, j)
  template<class T, int d>    
  inline void TinyBandMatrix<T, d>::SetEntry(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }

  
  //! Gaussian elimination of column j
  template<int p> template<class T, int d>
  inline void TinyBandMatrixLoop<p>::EliminateColumn(TinyBandMatrix<T, d>& A, int j)
  {
    T pivot = A.diag_(j)(p-1)*A.diag_(j)(d);
    TinyBandMatrixDoubleLoop<p, d>::CombineRow(A, j, pivot); 
    A.diag_(j)(p-1) = pivot;
    TinyBandMatrixLoop<p-1>::EliminateColumn(A, j);
  }
  
  
  //! solving by L
  template<int p> template<class T, int d, class T0>
  inline void TinyBandMatrixLoop<p>::
  SolveLower(const TinyBandMatrix<T, d>& A, int j, Vector<T0>& x)
  {
    x(j) -= A.diag_(j-p)(d-p)*x(j-p);
    TinyBandMatrixLoop<p-1>::SolveLower(A, j, x);
  }
  
  
  //! solving by U
  template<int p> template<class T, int d, class T0>
  inline void TinyBandMatrixLoop<p>::SolveUpper(const TinyBandMatrix<T, d>& A,
                                                int j, Vector<T0>& x)
  {
    x(j) -= A.diag_(j+p)(d+p)*x(j+p);
    TinyBandMatrixLoop<p-1>::SolveUpper(A, j, x);
  }
  

  //! the i-th row of matrix A is multiplied by coefficient coef
  template<int p> template<class T, int d, class T0>
  inline void TinyBandMatrixLoop<p>::MltRow(int i, const T0& coef, TinyBandMatrix<T, d>& A)
  {
    A.diag_(p+i-d-1)(p-1) *= coef;
    TinyBandMatrixLoop<p-1>::MltRow(i, coef, A);
  }
  
  
  //! matrix-vector product y = y + A x
  template<int p> template<class T, int size, class T1>
  inline void TinyBandMatrixLoop<p>::MltAdd(int offset, const TinyVector<T, size>& diag,
                                            const T1& xj, Vector<T1>& y)
  {
    y(offset-p) += diag(p-1)*xj;
    TinyBandMatrixLoop<p-1>::MltAdd(offset, diag, xj, y);
  }
  

  //! matrix-vector product y = y + A^T x
  template<int p> template<class T, int size, class T1>
  inline void TinyBandMatrixLoop<p>
  ::MltTransAdd(int offset, const TinyVector<T, size>& diag, const Vector<T1>& x, T1& yj)
  {
    yj += diag(p-1)*x(offset-p);
    TinyBandMatrixLoop<p-1>::MltTransAdd(offset, diag, x, yj);
  }


  //! matrix-vector product y = y + A^H x
  template<int p> template<class T, int size, class T1>
  inline void TinyBandMatrixLoop<p>
  ::MltConjTransAdd(int offset, const TinyVector<T, size>& diag, const Vector<T1>& x, T1& yj)
  {
    yj += conjugate(diag(p-1))*x(offset-p);
    TinyBandMatrixLoop<p-1>::MltConjTransAdd(offset, diag, x, yj);
  }

  
  //! combination of row to perform Gaussian elimination
  template<int p, int q> template<class T, int d>
  inline void TinyBandMatrixDoubleLoop<p, q>::
  CombineRow(TinyBandMatrix<T, d>& A, int j, const T& pivot)
  {
    A.diag_(j+q)(p+q-1) -= pivot*A.diag_(j+q)(d+q);
    TinyBandMatrixDoubleLoop<p, q-1>::CombineRow(A, j, pivot); 
  }
  
  
  /*******************
   * TinyArrowMatrix *
   *******************/

  
  //! returns the number of rows
  template<class T, int d, int m>    
  inline int TinyArrowMatrix<T, d, m>::GetM() const
  {
    if (this->diag_.GetM() <= 0)
      return 0;
    
    return this->diag_.GetM() + m;
  }
  
  
  //! returns the number of columns
  template<class T, int d, int m>    
  inline int TinyArrowMatrix<T, d, m>::GetN() const
  {
    if (this->diag_.GetM() <= 0)
      return 0;
    
    return this->diag_.GetM() + m;
  }
  

  //! returns the number of elements stored
  template<class T, int d, int m>    
  inline int TinyArrowMatrix<T, d, m>::GetDataSize() const
  {
    int size = TinyBandMatrix<T, d>::GetDataSize();
    size += 2*last_row_.GetM()*m + last_block_.GetSize();
    return size;
  }
  
  
  //! returns the memory used by the matrix in bytes
  template<class T, int d, int m>
  inline size_t TinyArrowMatrix<T, d, m>::GetMemorySize() const
  {
    return this->diag_.GetMemorySize() + last_row_.GetMemorySize()
      + last_col_.GetMemorySize() + sizeof(last_block_);
  }
  
  
  //! eliminates last columns in Gaussian elimination
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::
  EliminateLastColumn(TinyArrowMatrix<T, d, m>& A, int j)
  {
    T pivot = A.diag_(j)(d-p);
    TinyArrowMatrixDoubleLoop<p, m>::CombineLastColumn(A, j, pivot);
    TinyArrowMatrixLoop<p-1>::EliminateLastColumn(A, j);
  }
  
  
  //! eliminates last rows in Gaussian elimination
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::EliminateLastRow(TinyArrowMatrix<T, d, m>& A, int j)
  {
    T pivot = A.last_row_(j)(p-1)*A.diag_(j)(d);
    A.last_row_(j)(p-1) = pivot;
    TinyArrowMatrixDoubleLoop<p, d>::CombineLastRow(A, j, pivot);    
    TinyArrowMatrixDoubleLoop<p, m>::CombineLastBlock(A, j, pivot);
        
    TinyArrowMatrixLoop<p-1>::EliminateLastRow(A, j);
  }
  
  
  //! solving by L
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::
  SolveLower(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x)
  {
    x(j) -= A.diag_(j-p)(d-p)*x(j-p);
    TinyArrowMatrixLoop<p-1>::SolveLower(A, j, x);
  }
  
  
  //! function used in SolveLU
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::
  SolveRow(const TinyArrowMatrix<T, d, m>& A, int j, int n, Vector<T>& x)
  {
    x(n+p-1) -= x(j)*A.last_row_(j)(p-1);
    TinyArrowMatrixLoop<p-1>::SolveRow(A, j, n, x);
  }
  
  
  //! function used in SolveLU
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::
  SolveColumn(const TinyArrowMatrix<T, d, m>& A, int j, int n, Vector<T>& x)
  {
    x(j) -= A.last_col_(j)(p-1)*x(n+p-1);
    TinyArrowMatrixLoop<p-1>::SolveColumn(A, j, n, x);
  }
  
  
  //! function used in SolveLU
  template<int p> template<class T, int d, int m>
  inline void TinyArrowMatrixLoop<p>::
  SolveUpper(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x)
  {
    x(j) -= A.diag_(j+p)(d+p)*x(j+p);
    TinyArrowMatrixLoop<p-1>::SolveUpper(A, j, x);
  }

  
  //! function used for matrix-vector product
  template<int p> template<class T, int m, class T1>
  inline void TinyArrowMatrixLoop<p>
  ::MltAdd(int n, const TinyVector<T, m>& coef, const Vector<T1>& x, T1& yj)
  {
    yj += coef(p-1)*x(n+p-1);
    TinyArrowMatrixLoop<p-1>::MltAdd(n, coef, x, yj);
  }


  //! function used for matrix-vector product
  template<int p> template<class T, int m, class T1>
  inline void TinyArrowMatrixLoop<p>
  ::MltConjAdd(int n, const TinyVector<T, m>& coef, const Vector<T1>& x, T1& yj)
  {
    yj += conjugate(coef(p-1))*x(n+p-1);
    TinyArrowMatrixLoop<p-1>::MltConjAdd(n, coef, x, yj);
  }

  
  //! function used for matrix-vector product
  template<int p> template<class T, int m, class T1>
  inline void TinyArrowMatrixLoop<p>
  ::MltTransAdd(int n, const TinyVector<T, m>& coef, const T1& xj, Vector<T1>& y)
  {
    y(n+p-1) += coef(p-1)*xj;
    TinyArrowMatrixLoop<p-1>::MltTransAdd(n, coef, xj, y);
  }


  //! function used for matrix-vector product
  template<int p> template<class T, int m, class T1>
  inline void TinyArrowMatrixLoop<p>
  ::MltConjTransAdd(int n, const TinyVector<T, m>& coef, const T1& xj, Vector<T1>& y)
  {
    y(n+p-1) += conjugate(coef(p-1))*xj;
    TinyArrowMatrixLoop<p-1>::MltConjTransAdd(n, coef, xj, y);
  }

  
  //! function used in GetLU
  template<int p, int q> template<class T, int d, int m>
  inline void TinyArrowMatrixDoubleLoop<p, q>::
  CombineLastColumn(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot)
  {
    A.last_col_(j+p)(q-1) -= pivot*A.last_col_(j)(q-1);
    TinyArrowMatrixDoubleLoop<p, q-1>::CombineLastColumn(A, j, pivot);
  }


  //! function used in GetLU
  template<int p, int q> template<class T, int d, int m>
  inline void TinyArrowMatrixDoubleLoop<p, q>::
  CombineLastRow(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot)
  {
    A.last_row_(j+q)(p-1) -= pivot*A.diag_(j+q)(d+q);
    TinyArrowMatrixDoubleLoop<p, q-1>::CombineLastRow(A, j, pivot);
  }
  
  
  //! function used in GetLU
  template<int p, int q> template<class T, int d, int m>
  inline void TinyArrowMatrixDoubleLoop<p, q>::
  CombineLastBlock(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot)
  {
    A.last_block_(p-1, q-1) -= pivot*A.last_col_(j)(q-1);  
    TinyArrowMatrixDoubleLoop<p, q-1>::CombineLastBlock(A, j, pivot);
  }
  
}

#define SELDON_FILE_TINY_BAND_MATRIX_INLINE_CXX
#endif

