#ifndef SELDON_FILE_MATRIX_BLOCK_DIAGONAL_INLINE_CXX

#include "MatrixBlockDiagonal.hxx"

namespace Seldon
{
  
  template <class T, class Prop, class Storage, class Allocator>
  inline size_t Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + sizeof(int)*(3*nblok_+1+2*this->m_);
    taille += this->nz_*sizeof(T);
    return taille;
  }    
  

  //! returns the number of non-zero entries on row i
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowSize(int i) const
  {
    int n = num_blok_[i];
    return ptr_[n+1] - ptr_[n];
  }
  
  
  //! returns the maximal size for all the blocks contained in the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetSizeMaxBlock() const
  {
    return size_max_blok_;
  }
  
  
  //! returns the block number of the row i
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockNumber(int i) const
  {
    return num_blok_[i];
  }
  
  
  //! returns the number of non-zero entries stored
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->nz_;
  }
  
  
  //! return the number of blocks
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbBlocks() const
  {
    return nblok_;
  }
  
  
  //! return the size of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockSize(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->nblok_)
      throw WrongRow("Matrix_BlockDiagonal::GetBlockSize(int)",
		     string("Index should be in [0, ")
		     + to_str(this->nblok_-1) + "], but is equal to "
		     + to_str(i) + ".");
#endif
    
    return ptr_[i+1] - ptr_[i];
  }
  
  
  //! returns the number of row numbers (including distant rows in parallel)
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbRowNumbers() const
  {
    return size_ind_;
  }
  
  
  //! returns the pointer of the array storing row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowNumData() const
  {
    return ind_;
  }
  
  
  //! returns the pointer of the array "ptr" (beginning of each block)
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowPtr() const
  {
    return ptr_;
  }
  
  
  //! returns the pointer to the array giving local position of rows
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetReverseNum() const
  {
    return rev_;
  }
  
  
  //! returns global row number of local row j of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Index(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->nblok_, this->ptr_[i+1]-this->ptr_[i],
		"Matrix_BlockDiagonal");
#endif

    return ind_[ptr_[i] + j];
  }
 
  
  //! returns non-zero entry (j, k) of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  Value(int i, int j, int k) const
  {
    int size_row = ptr_[i+1] - ptr_[i];
    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, this->nblok_, size_row, size_row, "Matrix_BlockDiagonal");
#endif

    return this->data_[offset_[i] + size_row*j + k];
  }
  
  
  //! returns non-zero entry (j, k) of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline T& Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  Value(int i, int j, int k)
  {
    int size_row = ptr_[i+1] - ptr_[i];
    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, this->nblok_, size_row, size_row, "Matrix_BlockDiagonal");
#endif

    return this->data_[offset_[i] + size_row*j + k];
  }
  
  
  //! returns first index for each block
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }
  
  
  //! returns row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }
  
  
  //! returns local row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetRev() const
  {
    return rev_;
  }
  
  
  //! returns block numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetNumBlock() const
  {
    return num_blok_;
  }
  
  
  //! returns offset for accessing non-zero entries for each block
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::GetOffset() const
  {
    return offset_;
  }
  
  
  //! Adds a to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  AddInteraction(int i, int j, const T& a)
  {
    Vector<T> val(1);
    IVect col(1); val(0) = a; col(0) = j;
    AddInteractionRow(i, 1, col, val);
  }
  
  
  /***************************
   * MatrixSym_BlockDiagonal *
   ***************************/


  template <class T, class Prop, class Storage, class Allocator>
  inline size_t MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + sizeof(int)*(3*nblok_+1+2*this->m_);
    taille += this->nz_*sizeof(T);
    return taille;
  }    
  
  
  //! returns the number of non-zero entries on row i
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowSize(int i) const
  {
    int n = num_blok_[i];
    return ptr_[n+1] - ptr_[n];
  }
  
  
  //! returns the maximal size for all the blocks contained in the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetSizeMaxBlock() const
  {
    return size_max_blok_;
  }
  

  //! returns the block number of the row i
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockNumber(int i) const
  {
    return num_blok_[i];
  }

  
  //! returns the number of non-zero entries stored
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->nz_;
  }
  
  
  //! return the number of blocks
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbBlocks() const
  {
    return nblok_;
  }
  
  
  //! return the size of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockSize(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->nblok_)
      throw WrongRow("MatrixSym_BlockDiagonal::GetBlockSize(int)",
		     string("Index should be in [0, ")
		     + to_str(this->nblok_-1) + "], but is equal to "
		     + to_str(i) + ".");
#endif
    
    return ptr_[i+1] - ptr_[i];
  }
  
  
  //! returns the number of row numbers (including distant rows in parallel)
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbRowNumbers() const
  {
    return size_ind_;
  }
  
  
  //! returns the pointer of the array storing row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowNumData() const
  {
    return ind_;
  }
  
  
  //! returns the pointer of the array "ptr" (beginning of each block)
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowPtr() const
  {
    return ptr_;
  }
  
  
  //! returns the pointer to the array giving local position of rows
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetReverseNum() const
  {
    return rev_;
  }
  
  
  //! returns global row number of local row j of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Index(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->nblok_, this->ptr_[i+1]-this->ptr_[i],
		"MatrixSym_BlockDiagonal");
#endif

    return ind_[ptr_[i] + j];
  }
 
  
  //! returns non-zero entry (j, k) of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  Value(int i, int j, int k) const
  {
    int size_row = ptr_[i+1] - ptr_[i];
    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, this->nblok_, size_row, size_row, "MatrixSym_BlockDiagonal");
#endif

    if (j < k)
      return this->data_[offset_[i]+ j*size_row - j*(j+1)/2 + k];
    else
      return this->data_[offset_[i]+ k*size_row - k*(k+1)/2 + j];
  }
  
  
  //! returns non-zero entry (j, k) of block i
  template <class T, class Prop, class Storage, class Allocator>
  inline T& MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  Value(int i, int j, int k)
  {
    int size_row = ptr_[i+1] - ptr_[i];
    
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, k, this->nblok_, size_row, size_row, "MatrixSym_BlockDiagonal");
#endif

    if (j < k)
      return this->data_[offset_[i]+ j*size_row - j*(j+1)/2 + k];
    else
      return this->data_[offset_[i]+ k*size_row - k*(k+1)/2 + j];
  }
  
  
  //! returns first index for each block
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }
  
  
  //! returns row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }
  
  
  //! returns local row numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRev() const
  {
    return rev_;
  }
  
  
  //! returns block numbers
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetNumBlock() const
  {
    return num_blok_;
  }
  
  
  //! returns offset for accessing non-zero entries for each block
  template <class T, class Prop, class Storage, class Allocator>
  inline int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetOffset() const
  {
    return offset_;
  }
  
  
  //! adds value a to entry (i,j) of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixSym_BlockDiagonal<T,Prop,Storage,Allocator>::
  AddInteraction(int i, int j, const T& a)
  {
    Vector<T> val(1);
    IVect col(1); val(0) = a; col(0) = j;
    AddInteractionRow(i, 1, col, val);
  }


  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  inline void CopySubMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
			    const IVect& row, const IVect& col,
			    Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {}
  
  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  inline void CopySubMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
			    const IVect& row, const IVect& col,
			    Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {}

#ifndef SELDON_WITH_MPI
  template<class T0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, General, BlockDiagRow, Allocator0>& A,
		  Vector<T0>& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<T0, General, BlockDiagRow, Allocator0>& >(A), Y);    
  }

  template<class T0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
                  Vector<T0>& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& >(A), Y); 
  }
#endif
    
}

#define SELDON_FILE_MATRIX_BLOCK_DIAGONAL_INLINE_CXX
#endif
