#ifndef SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_HXX

namespace Seldon
{
  template<class T, class Prop, class Storage, class Allocator
	   = typename SeldonDefaultAllocator<Storage, T>::allocator >
  class DistributedMatrix_BlockDiag : public Matrix<T, Prop, Storage, Allocator>
  {    
  protected :
    //! row numbers shared with other processors
    /*!
      It is the same array required in a distributed vector.
      Some rows can be shared by several processors.
    */
    IVect* OverlapRowNumbers;
    
    //! processor where each shared row should be assembled
    /*!
      For rows which have already been counted in other processors,
      you specify on which processor these rows are assembled
    */
    IVect* OverlapProcNumbers;
    
    //! global row numbers 
    IVect* GlobalRowNumbers;

    //! list of processors sharing rows with the current one
    IVect* ProcSharingRows;
    
    //! for each processor sharing rows, list of "local numbers" shared
    /*!
      It is assumed that these local numbers are "sorted", such that
      local numbers shared with processor j on processor i are corresponding
      with local numbers shared with processor i on processor j
     */
    Vector<IVect>* SharingRowNumbers;
    
    //! number of "scalar" unknowns
    int nodl_scalar_, nb_unknowns_scal_;
    
    //! total number of rows (on all processors)
    int nglob_;
    
    //! MPI communicator
    MPI_Comm comm_;
    //! for each row, processor associated with
    IVect proc_row;

    //! global row/col numbers (needed for MltAdd)
    IVect global_col_to_recv;
    IVect ptr_global_col_to_recv;
    
    //! local row/col numbers (needed for MltAdd)
    Vector<IVect> local_col_to_send;
    
    //! processor numbers (needed for MltAdd)
    IVect proc_col_to_recv, proc_col_to_send;
    
    //! if true local numbers are present in row_num instead of global numbers
    bool local_number_distant_values;
    
    //! number of distant non-zero entries
    int size_max_distant_col;
    
    //! true if the matrix is symmetric
    bool sym_matrix;
    
    // internal functions
    void EraseArrayForMltAdd();
    void SwitchToGlobalNumbers();
    void CompleteWithAssociatedNumbers(IVect& num);
    
    void SortAndAssembleDistantInteractions(const IVect& proc_num, IVect& dist_num,
                                            IVect& glob_num, IVect& ptr_glob_num, IVect& proc_glob,
                                            Vector<IVect>& local_num, IVect& proc_local);
    template<class T2>
    void ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                       const IVect&, const IVect& proc_recv,
                       const Vector<IVect>& num_send,
                       const IVect& proc_send, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                        const IVect&, const IVect& proc_recv,
                        const Vector<IVect>& num_send,
                        const IVect& proc_send, Vector<T2>& X) const;
    
    void AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                           const IVect& num_recv, const IVect& ptr_num_recv,
                           const IVect& proc_recv,
                           const Vector<IVect>& num_send, const IVect& proc_send,
                           IVect& Y, IVect& Yproc) const;
    
    void AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const;
    
  public :
    DistributedMatrix_BlockDiag();
    
    void Clear();
    
    // initialisation of pointers
    void Init(int n, IVect*, IVect*, IVect*,
              int, int, IVect*, Vector<IVect>*, const MPI_Comm&);
    
    template<class T0, class Prop0, class Storage0, class Allocator0>
    void Init(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A);
    
    // parallel stuff (inline)
    MPI_Comm& GetCommunicator();
    const MPI_Comm& GetCommunicator() const;
    
    int GetGlobalM() const;
    int GetNodlScalar() const;
    int GetNbScalarUnknowns() const;
    int* GetProcNumber() const;
    const IVect& GetProcRowNumber() const;
    
    // non inline
    IVect& GetGlobalRowNumber();
    const IVect& GetGlobalRowNumber() const;
    IVect& GetOverlapRowNumber();
    const IVect& GetOverlapRowNumber() const;
    IVect& GetOverlapProcNumber();
    const IVect& GetOverlapProcNumber() const;
    IVect& GetProcessorSharingRows();
    Vector<IVect>& GetSharingRowNumbers();
    
    // construction of matrix
    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
                    const Vector<IVect>& ProcDofs);
    
    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks);

    // modifies matrix
    void AddDistantInteraction(int i, int jglob, int proc, const T& val);
    void AddRowDistantInteraction(int iglob, int j, int proc, const T& val);
    
    int IndexGlobal(int i, int j) const;

    // conversion
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;
    
    template<class T0, class Prop0, class Storage0>
    void ConvertToSparse(DistributedMatrix<T0, Prop0, Storage0>& A) const;

    // functions for matrix-vector product
    bool IsReadyForMltAdd() const;
    
    void PrepareMltAdd();
    
    template<class T2>
    void ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const;
    
    template<class T2>
    void AssembleVec(Vector<T2>&) const;
    
    // Blas/Lapack functions
    template<class T2, class Allocator2>
    void Mlt_(Vector<T2, VectFull, Allocator2>&, bool assemble) const;

    template<class T0, class T1, class Allocator1, class T3, class T2, class Allocator2>
    void MltAdd_(const T0& alpha, const SeldonTranspose&,
                 const Vector<T1, VectFull, Allocator1>& X, const T3& beta,
                 Vector<T2, VectFull, Allocator2>& Y, bool assemble) const;

    template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
    void Add_(const T0& alpha, const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);
    
    void Assemble();

    template<class T1, class Prop1, class Allocator1>
    friend void EraseRow(const IVect& num,
			 DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A);

    template<class T1, class Prop1, class Allocator1>
    friend void EraseRow(const IVect& num,
			 DistributedMatrix<T1, Prop1, BlockDiagRowSym, Allocator1>& A);

    template<class T1, class Prop1, class Allocator1>
    friend void EraseCol(const IVect& num,
			 DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A);

  };

  template<class T, class Prop, class Allocator>
  class DistributedMatrix<T, Prop, BlockDiagRow, Allocator>
    : public DistributedMatrix_BlockDiag<T, Prop, BlockDiagRow, Allocator>
  {
  };
  
  
  template<class T, class Prop, class Allocator>
  class DistributedMatrix<T, Prop, BlockDiagRowSym, Allocator>
    : public DistributedMatrix_BlockDiag<T, Prop, BlockDiagRowSym, Allocator>
  {
  };


  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		 DistributedMatrix<T2, General, BlockDiagRow, Allocator2>& B);
  
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, General, BlockDiagRow, Alloc>& A,
                       DistributedMatrix<T, General, ArrayRowSparse, Alloc2>& B);

  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, General, BlockDiagRow, Alloc>& A,
                       DistributedMatrix<T, General, ArrayRowComplexSparse, Alloc2>& B);

  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		 DistributedMatrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& B);

  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       DistributedMatrix<T, Symmetric, ArrayRowSymSparse, Alloc2>& B);

  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       DistributedMatrix<T, Symmetric, ArrayRowSymComplexSparse, Alloc2>& B);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
			      DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B,
			      int m1 = 0, int m2 = 0);
  
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRow, Allocator1>& B);

  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRow, Allocator1>& B);
  
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRowSym, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRow, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSparse, Allocator1>& B);

  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRowSym, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSymComplexSparse, Allocator1>& B);

  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRow, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowComplexSparse, Allocator1>& B);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const DistributedMatrix<T1, Prop, BlockDiagRow, Allocator1>& A,
		 Vector<T2, VectFull, Allocator2>& X, bool assemble = true);
  
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const DistributedMatrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X, bool assemble = true);

  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble = true);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble = true);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble = true);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble = true);
  
  template<class T, class Allocator>
  void GetInverse(DistributedMatrix<T, General, BlockDiagRow, Allocator>& A);
    
  template<class T, class Allocator>
  void GetInverse(DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, class Allocator>
  void GetCholesky(DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, class Allocator, class Allocator2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                     Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator, class Allocator2>
  void MltCholesky(const SeldonTranspose& trans,
                   const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                   Vector<T, VectFull, Allocator2>& x);

  // inline functions
  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const DistributedMatrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
		     const IVect& row, const IVect& col,
		     DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const DistributedMatrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
		     const IVect& row, const IVect& col,
		     DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);
  
  template<class T1, class Prop1, class Allocator1>
  void EraseRow(const IVect& num,
		DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A);

  template<class T1, class Prop1, class Allocator1>
  void EraseRow(const IVect& num,
		DistributedMatrix<T1, Prop1, BlockDiagRowSym, Allocator1>& A);

  template<class T1, class Prop1, class Allocator1>
  void EraseCol(const IVect& num,
		DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A);

  template<class T1, class Prop1, class Allocator1>
  void EraseCol(const IVect& num,
		DistributedMatrix<T1, Prop1, BlockDiagRowSym, Allocator1>& A);
    
  // last functions (non-inline)
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const DistributedMatrix<T, Prop, BlockDiagRow, Allocator>& A);

  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out,
                      const DistributedMatrix<T, Prop, BlockDiagRowSym, Allocator>& A);
  
}

#define SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_HXX
#endif
