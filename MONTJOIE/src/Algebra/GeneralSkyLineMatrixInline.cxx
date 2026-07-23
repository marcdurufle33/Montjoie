#ifndef MONTJOIE_FILE_GENERAL_SKYLINE_MATRIX_INLINE_CXX

#include "GeneralSkyLineMatrix.hxx"

namespace Seldon
{
  
  //! Default constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::GeneralSkyLineMatrix()
  {
    this->m_ = 0;
    this->n_ = 0;
  }
  
  
  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::Clear()
  {
    val.Clear();
    first_col.Clear();
    this->n_ = 0;
    this->m_ = 0;
  }
  
  
  //! allocation of a row
  template <class T, class Prop, class Storage, class Allocator>
  inline void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  ReallocateRow(int i, int size_row, int first_col_n)
  {
    val(i).Reallocate(size_row);
    first_col(i) = first_col_n;
  }
	
  
  //! returns size of row
  template <class T, class Prop, class Storage, class Allocator>
  inline int GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  GetRowSize(int i) const
  {
    return val(i).GetM();
  }


  //! returns column numbers of first non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  inline int GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  GetFirstColNumber(int i) const
  {
    return this->first_col(i);
  }
  
  
  //! returns j-th non-zero entry of row i
  template <class T, class Prop, class Storage, class Allocator>
  inline T& GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::Value(int i, int j)
  {
    return val(i)(j);
  }
  
  
  //! returns j-th non-zero entry of row i
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  Value(int i, int j) const
  {
    return val(i)(j);
  }
  
}
  
#define MONTJOIE_FILE_GENERAL_SKYLINE_MATRIX_INLINE_CXX
#endif
