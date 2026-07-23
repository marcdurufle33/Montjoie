#ifndef SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_INLINE_CXX

namespace Seldon
{
  
  //! returns MPI communicator (contains processors that will share the matrix)
  template<class T, class Prop, class Storage, class Allocator>
  inline MPI_Comm& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::GetCommunicator()
  {
    return comm_;
  }
  
  
  //! returns MPI communicator (contains processors that will share the matrix)
  template<class T, class Prop, class Storage, class Allocator>
  inline const MPI_Comm& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::GetCommunicator() const
  {
    return comm_;
  }
  
  
  //! returns the global number of rows
  /*!
    the global number of rows is equal to the number of rows in a sequential execution
   */
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetGlobalM() const
  {
    return nglob_;
  }
  
  
  //! returns the number of scalar unknowns
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetNodlScalar() const
  {
    return nodl_scalar_;
  }
  

  //! returns the number of scalar unknowns
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetNbScalarUnknowns() const
  {
    return nb_unknowns_scal_;
  }
  
  
  //! returns pointer to the array containing processor number for each row
  template<class T, class Prop, class Storage, class Allocator>
  inline int* DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetProcNumber() const
  {
    return proc_row.GetData();
  }


  //! returns the array containing processor number for each row
  template<class T, class Prop, class Storage, class Allocator>
  inline const IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetProcRowNumber() const
  {
    return proc_row;
  }


  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  inline void CopySubMatrix(const DistributedMatrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
			    const IVect& row, const IVect& col,
			    DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {}

  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  inline void CopySubMatrix(const DistributedMatrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
			    const IVect& row, const IVect& col,
			    DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {}
  
}

#define SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_INLINE_CXX
#endif
