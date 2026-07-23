#ifndef SELDON_FILE_MATRIX_BLOCK_DIAGONAL_HXX

namespace Seldon
{

  class BlockDiagRow : public RowMajor
  {
  public :
    typedef RowMajor StorageDense;
  };
  
  class BlockDiagRowSym : public RowSymPacked
  {
  public :
    typedef RowSymPacked StorageDense;
  };
  
  //! Symmetric block-diagonal matrix
  /*!
    Block diagonal matrices are defined by: (1) the number of rows and columns;
    (2) the size for each block; (3) the indices of rows for each block
    (4) the values for each block
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator >
  class MatrixSym_BlockDiagonal
  {
  public :
    typedef T entry_type;
    
  protected:
    //! number of rows
    int m_, n_;
    //! row numbers
    Vector<int> row_num_;
    //! for each row, local number of the row in the block
    Vector<int> rev_num_;
    //! offset for each blok
    Vector<int> ptr_, size_block_, offset_blok_;
    //! values for diagonal part of the matrix
    Vector<T> diag_value_;
    //! for each row, block number where it belongs
    Vector<int> blok_num_;
    //! maximal size for blocks
    int size_max_blok_;
    //! 2x2 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 2, 2> > mat22_;
    //! 3x3 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 3, 3> > mat33_;
    //! 4x4 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 4, 4> > mat44_;
    //! 5x5 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 5, 5> > mat55_;
    //! 6x6 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 6, 6> > mat66_;
    //! 7x7 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 7, 7> > mat77_;
    //! 8x8 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 8, 8> > mat88_;
    //! 9x9 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 9, 9> > mat99_;
    //! 10x10 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 10, 10> > mat1010_;
    //! 11x11 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 11, 11> > mat1111_;
    //! 12x12 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 12, 12> > mat1212_;
    //! 13x13 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 13, 13> > mat1313_;
    //! 14x14 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 14, 14> > mat1414_;
    //! 15x15 block-diagonal matrices
    /*Vector<TinyMatrix<T, Prop, 15, 15> > mat1515_;
    //! 16x16 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 16, 16> > mat1616_;
    //! 17x17 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 17, 17> > mat1717_;
    //! 15x15 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 18, 18> > mat1818_;
    //! 16x16 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 19, 19> > mat1919_;
    //! 17x17 block-diagonal matrices
    Vector<TinyMatrix<T, Prop, 20, 20> > mat2020_; */
    enum{max_size_small = 15};
    //! other block-diagonal matrices of size > 9
    Vector<Matrix<T, Prop, typename Storage::StorageDense, Allocator> > all_mat_;
    
    // Methods.
  public:
    // Constructors.
    MatrixSym_BlockDiagonal();
    
    void Clear();
    
    int GetRowSize(int i) const;
    int GetSizeMaxBlock() const;
    int GetBlockNumber(int i) const;
    
    // Memory management.
    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks);

    template<class Allocator1>
    void SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
                    const Vector<Vector<bool> >& ignore_row, IVect& permut);

    void Reallocate(int i, int j);
    void Resize(int i, int j);
        
    // Basic methods.
    int GetDataSize() const;
    int GetM() const;
    int GetN() const;
    int GetNbBlocks() const;
    int GetBlockSize(int i) const;

    int GetNbRowNumbers() const;
    int* GetRowNumData() const;
    int* GetRowPtr() const;
    int* GetReverseNum() const;
    
    int Index(int i, int j) const;
    
    const T& Value(int i, int j, int k) const;
    T& Value(int i, int j, int k);
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Zero();
    
    // Element access
    const T operator() (int i, int j) const;
    
    MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>&
    operator *= (const T& alpha);
    
    // input/output
    template<class T0, class Prop0, class Storage0>
    void ConvertToSparse(Matrix<T0, Prop0, Storage0>& A) const;
    
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;
    
    // methods used for blas/lapack functions
    template<class T2, class Allocator2>
    void Mlt_(Vector<T2, VectFull, Allocator2>& ) const;

    template<class T0, class Prop0, class Storage0, class Allocator0,
             class T1, class Prop1, class Storage1, class Allocator1>
    void Mlt_(const MatrixSym_BlockDiagonal<T0, Prop0, Storage0, Allocator0>&,
              MatrixSym_BlockDiagonal<T1, Prop1, Storage1, Allocator1>&) const;
    
    template<class T0, class T1, class Allocator1, class T2, class Allocator2>
    void MltAdd_(const T0& alpha, const class_SeldonNoTrans&,
                 const Vector<T1, VectFull, Allocator1>& X,
                 Vector<T2, VectFull, Allocator2>& Y) const;

    template<class T0, class T1, class Allocator1, class T2, class Allocator2>
    void MltAdd_(const T0& alpha, const class_SeldonTrans&,
                 const Vector<T1, VectFull, Allocator1>& X,
                 Vector<T2, VectFull, Allocator2>& Y) const;

    template<class T0, class T1, class Allocator1, class T2, class Allocator2>
    void MltAdd_(const T0& alpha, const class_SeldonConjTrans&,
                 const Vector<T1, VectFull, Allocator1>& X,
                 Vector<T2, VectFull, Allocator2>& Y) const;
    
    void GetInverse_();
    
    template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
    void Add_(const T0& alpha, const Matrix<T1, Prop1, Storage1, Allocator1>& A);
    
  };
  
  
  //! Symmetric block-diagonal matrix stored by rows
  template<class T, class Allocator>
  class Matrix<T, Symmetric, BlockDiagRowSym, Allocator>
    : public MatrixSym_BlockDiagonal<T, Symmetric, BlockDiagRowSym, Allocator>
  {
  public :
    
    void AddInteraction(int i, int j, const T& a);
    
    void AddInteractionRow(int i, int nb_interac, const IVect& col_interac,
			   const Vector<T>& val_interac);
    
    template<class T0, class Allocator0>
    friend void GetCholesky(Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A);
    
    template<class T0, class Allocator0, class T2, class Allocator2>
    friend void SolveCholesky(const SeldonTranspose& trans,
                              const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
                              Vector<T2, VectFull, Allocator2>& x);
    
    template<class T0, class Allocator0, class T2, class Allocator2>
    friend void MltCholesky(const SeldonTranspose& trans,
                            const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
                            Vector<T2, VectFull, Allocator2>& x);
    
  };


  //! General block-diagonal matrix stored by rows
  template<class T, class Allocator>
  class Matrix<T, General, BlockDiagRow, Allocator>
    : public MatrixSym_BlockDiagonal<T, General, BlockDiagRow, Allocator>
  {
  public :
    
    void AddInteraction(int i, int j, const T& a);
    
    void AddInteractionRow(int i, int nb_interac, const IVect& col_interac,
			   const Vector<T>& val_interac);

  };

  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		 Matrix<T2, General, BlockDiagRow, Allocator2>& B);
  
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowSparse, Alloc2>& B);

  template<class T, class Allocator1, class Allocator2>
  void Copy(const Matrix<T, General, BlockDiagRow, Allocator1>& A,
            Matrix<T, General, RowMajor, Allocator2>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		 Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& B);
  
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       Matrix<T, Symmetric, ArrayRowSymSparse, Alloc2>& B);
  
  template<class T, class Allocator1, class Allocator2>
  void Copy(const Matrix<T, Symmetric, BlockDiagRowSym, Allocator1>& A,
            Matrix<T, Symmetric, RowSymPacked, Allocator2>& B);

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B,
                              int m1, int m2);

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
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B,
                         int m1, int m2);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B,
                         int m1, int m2);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B);
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B);

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X);

  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X);

  template<class T0, class Allocator0, class T1, class Allocator1,
           class T2, class Allocator2>
  void MltMatrix(const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>&,
		 const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>&,
		 Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>&);
  
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonTrans& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonTrans& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonConjTrans& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonConjTrans& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B);
  
  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, BlockDiagRow, Allocator>& A);
    
  template<class T, class Allocator>
  void GetInverse(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, General, BlockDiagRow, Allocator>& A);

  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A);

  template<class T, int m>
  void MltMatrix(const TinyMatrix<T, Symmetric, m, m>& A, const TinyMatrix<T, Symmetric, m, m>& B,
                 TinyMatrix<T, Symmetric, m, m>& C);

  template<class T, int m>
  void GetInverseMatrix(TinyMatrix<T, General, m, m>& A);
  
  template<class T, int m>
  void GetInverseMatrix(TinyMatrix<T, Symmetric, m, m>& A);
  
  template<class T, int m>
  void GetCholeskyMatrix(TinyMatrix<T, Symmetric, m, m>& A);

  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, BlockDiagRow, Allocator>& A);

  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, BlockDiagRowSym, Allocator>& A);

}

#define SELDON_FILE_MATRIX_BLOCK_DIAGONAL_HXX
#endif
