#ifndef SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_INLINE_CXX

#include "SymmetricSkyLineMatrix.hxx"

namespace Seldon
{
  
  //! Default constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::SymmetricSkyLineMatrix()
  {
  }
  
  
  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::Clear()
  {
    ind.Clear();
    val.Clear();
    this->m_ = 0; this->n_ = 0;
  }
  
  
  //! nothing
  template<class T, class Prop, class Storage, class Allocator>
  inline void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // matrix can only be constructed with SetData or Copy
  }
  
  
  //! returns non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  inline T* SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::GetData() const
  {
    return this->val.GetData();
  }
	
  
  //! returns indices to access beginning of columns
  template <class T, class Prop, class Storage, class Allocator>
  inline int* SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::GetInd() const
  {
    return this->ind.GetData();
  }
  

  //! Returns the size used by the object in bytes.
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    return val.GetMemorySize() + ind.GetMemorySize();
  }

  
  //! returns number of non-zero entries stored
  template <class T, class Prop, class Storage, class Allocator>
  inline int SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return val.GetM();
  }
  
}

#define SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_INLINE_CXX
#endif
