#ifndef SELDON_FILE_MATRIX_BLOCK_DIAGONAL_HXX

namespace Seldon
{
  
  class BlockDiagRow : public RowSparse
  {
  };
  
  class BlockDiagRowSym : public RowSymSparse
  {
  };
  
  //! Block-diagonal matrix (the pattern of the matrix is symmetric)
  /*!
    Block diagonal matrices are defined by: (1) the number of rows and columns;
    (2) the size for each block; (3) the indices of rows for each block
    (4) the values for each block
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator >
  class Matrix_BlockDiagonal : public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename SeldonDefaultAllocator<VectFull, int>::allocator AllocatorInt;
    
    // Attributes.
  protected:
    
    //! Index of first index for each block
    /*!
      size of each block is ptr_[i+1] - ptr-[i]
    */
    int* ptr_;    
    //! row numbers
    int* ind_;
    //! local row number inside a block
    /*!
      For a global row i, rev_[i] gives the local position 
      of this row in the block
     */
    int* rev_;
    //! offset indices to access each block in the array data
    int* offset_;
    //! for each row, block to which it belongs
    int* num_blok_;
    //! number of non-zero elements
    int nz_;
    //! number of blocks
    int nblok_;
    //! maximal block size
    int size_max_blok_;
    //! size of ind_
    int size_ind_;
    
    // Methods.
  public:
    // Constructors.
    Matrix_BlockDiagonal();
    Matrix_BlockDiagonal(const Matrix_BlockDiagonal<T, Prop,
			 Storage, Allocator>& A);
    
    // Destructor.
    ~Matrix_BlockDiagonal();

    void Clear();
    void ClearRow(int i);
    void ClearColumn(int i);
    
    size_t GetMemorySize() const;
    int GetRowSize(int i) const;
    int GetSizeMaxBlock() const;
    int GetBlockNumber(int i) const;
    
    // Memory management.
    void SetPattern(const IVect& num_dof_blocks, const IVect& size_blocks,
                    const Vector<bool>& ignore_row, IVect& permut);

    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
                    const Vector<Vector<bool> >& ignore_row, IVect& permut);
    
    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks);
    
    void Copy(const Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& B);
    void Reallocate(int i, int j);
    void Resize(int i, int j);
    
    // Basic methods.
    int GetDataSize() const;
    int GetNbBlocks() const;
    int GetBlockSize(int i) const;

    int GetNbRowNumbers() const;
    int* GetRowNumData() const;
    int* GetRowPtr() const;
    int* GetReverseNum() const;
    
    int Index(int i, int j) const;
    
    const T& Value(int i, int j, int k) const;
    T& Value(int i, int j, int k);
    
    int* GetPtr() const;
    int* GetInd() const;
    int* GetRev() const;
    int* GetNumBlock() const;
    int* GetOffset() const;
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Zero();
    
    // Element acess and affectation.
    const value_type operator() (int i, int j) const;
    
    reference Get(int i, int j);
    const_reference Get(int i, int j) const;
    
    void Set(int, int j, const T& val);
    void SetEntry(int, int j, const T& val);
    
    Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& 
    operator= (const Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& A);
    
    template<class T0>
    Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& operator *= (const T0& alpha);
    
    void AddInteraction(int i, int j, const T& a);
    
    void AddInteractionRow(int i,int nb_interac, const IVect& col_interac,
			   const Vector<T>& val_interac, bool sorted = false);
    
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;

    void RemoveSmallEntry(typename ClassComplexType<T>::Treal eps){} 
    
  };
  
  //! block-diagonal matrix with symmetric profile
  template<class T, class Allocator>
  class Matrix<T, General, BlockDiagRow, Allocator>
    : public Matrix_BlockDiagonal<T, General, BlockDiagRow, Allocator>
  {
  };
  

  //! Symmetric block-diagonal matrix
  /*!
    Block diagonal matrices are defined by: (1) the number of rows and columns;
    (2) the size for each block; (3) the indices of rows for each block
    (4) the values for each block
  */
  template <class T, class Prop, class Storage,
            class Allocator = typename SeldonDefaultAllocator<Storage, T>::allocator >
  class MatrixSym_BlockDiagonal : public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename SeldonDefaultAllocator<VectFull, int>::allocator AllocatorInt;
    
    // Attributes.
  protected:

    //! Index of first index for each block
    /*!
      size of each block is ptr_[i+1] - ptr-[i]
    */
    int* ptr_;
    //! numbers of rows
    int* ind_;
    //! local row number inside a block
    /*!
      For a global row i, rev_[i] gives the local position 
      of this row in the block
     */
    int* rev_;
    //! offset indices to access each block in the array data
    int* offset_;
    //! for each row, block to which it belongs
    int* num_blok_;
    //! number of non-zero elements
    int nz_;
    //! number of blocks
    int nblok_;
    //! maximal block size
    int size_max_blok_;
    //! size of ind_
    int size_ind_;
    
    // Methods.
  public:
    // Constructors.
    MatrixSym_BlockDiagonal();
    MatrixSym_BlockDiagonal(const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& A);
    
    // Destructor.
    ~MatrixSym_BlockDiagonal();
    
    void Clear();
    void ClearRow(int i);
    void ClearColumn(int i);

    size_t GetMemorySize() const;
    int GetRowSize(int i) const;
    int GetSizeMaxBlock() const;
    int GetBlockNumber(int i) const;
    
    // Memory management.
    void SetPattern(const IVect& num_dof_blocks, const IVect& size_blocks,
                    const Vector<bool>& ignore_row, IVect& permut);

    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
                    const Vector<Vector<bool> >& ignore_row, IVect& permut);
    
    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks);
    
    void Copy(const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& B);
    void Reallocate(int i, int j);
    void Resize(int i, int j);
    
    // Basic methods.
    int GetDataSize() const;
    int GetNbBlocks() const;
    int GetBlockSize(int i) const;
    
    int GetNbRowNumbers() const;
    int* GetRowNumData() const;
    int* GetRowPtr() const;
    int* GetReverseNum() const;
    
    int Index(int i, int j) const;
    
    const T& Value(int i, int j, int k) const;
    T& Value(int i, int j, int k);
    
    int* GetPtr() const;
    int* GetInd() const;
    int* GetRev() const;
    int* GetNumBlock() const;
    int* GetOffset() const;
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Zero();
    
    // Element access and affectation.
    value_type operator() (int i, int j) const;

    reference Get(int i, int j);
    const_reference Get(int i, int j) const;
    
    void Set(int, int j, const T& val);
    void SetEntry(int, int j, const T& val);
        
    MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>&
    operator= (const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& A);
    
    template<class T0>
    MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& operator *= (const T0& alpha);

    void AddInteraction(int i, int j, const T& a);
    
    void AddInteractionRow(int i,int nb_interac, const IVect& col_interac,
			   const Vector<T>& val_interac, bool sorted = false);
  
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;

    void RemoveSmallEntry(typename ClassComplexType<T>::Treal eps){}
    
  };


  //! Symmetric block-diagonal matrix
  template<class T, class Allocator>
  class Matrix<T, Symmetric, BlockDiagRowSym, Allocator>
    : public MatrixSym_BlockDiagonal<T, Symmetric, BlockDiagRowSym, Allocator>
  {
  };


  template<class T0, class T1, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		 Matrix<T1, General, BlockDiagRow, Allocator2>& B);
  
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowSparse, Alloc2>& B);

  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowComplexSparse, Alloc2>& B);
  
  template<class T, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, General, BlockDiagRow, Allocator1>& A,
		  Matrix<T, General, RowMajor, Allocator2>& B);
  
  template<class T0, class T1, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		 Matrix<T1, Symmetric, BlockDiagRowSym, Allocator2>& B);
  
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       Matrix<T, Symmetric, ArrayRowSymSparse, Alloc2>& B);

  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       Matrix<T, Symmetric, ArrayRowSymComplexSparse, Alloc2>& B);
  
  template<class T, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Symmetric, BlockDiagRowSym, Allocator1>& A,
		  Matrix<T, Symmetric, RowSymPacked, Allocator2>& B);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2);
    
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B, int m1, int m2);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B, int m1, int m2);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X);


  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const SeldonTranspose&,
		 const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const SeldonTranspose&,
		 const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y);

  template<class T0, class Allocator0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4>
  void MltAddMatrix(const T3& alpha, const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>&,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>&,
		    const T4& beta, Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>&);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Prop, class Allocator1,
	   class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, BlockDiagRow, Allocator>& A);
    
  template<class T, class Allocator>
  void GetInverse(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, class Allocator>
  void GetCholesky(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, class Allocator, class Allocator2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                     Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator, class Allocator2>
  void MltCholesky(const SeldonTranspose& trans,
                   const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, General, BlockDiagRow, Allocator>& A);

  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  // inline functions
  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Matrix<T1, Prop1, Storage1, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void CopySubMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
                     const IVect& row, const IVect& col,
                     Matrix<T1, Prop1, Storage1, Allocator1>& B);

#ifndef SELDON_WITH_MPI
  template<class T0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, General, BlockDiagRow, Allocator0>& A,
	   Vector<T0>& Y, bool assemble = true);

  template<class T0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
	   Vector<T0>& Y, bool assemble = true);
#endif
  
}

#define SELDON_FILE_MATRIX_BLOCK_DIAGONAL_HXX
#endif
