#ifndef MONTJOIE_FILE_DIAGONAL_MATRIX_INLINE_CXX

namespace Seldon
{
  //! default constructor
  template<class T, class Prop, class Allocator>
  inline Matrix<T, Prop, DiagonalRow, Allocator>::Matrix()
  {
  }
  
  
  //! matrix with i rows
  template<class T, class Prop, class Allocator>
  inline Matrix<T, Prop, DiagonalRow, Allocator>::Matrix(int i, int j)
    : VirtualMatrix<T>(i, j)
  {    
    T zero;
    SetComplexZero(zero);
    
    diag.Reallocate(i);
    diag.Fill(zero);
  }
  
  
  //! returns the memory in bytes used by the object
  template<class T, class Prop, class Allocator>
  inline size_t Matrix<T, Prop, DiagonalRow, Allocator>::GetMemorySize() const 
  {
    return diag.GetMemorySize();
  }
  
  
  //! matrix with i rows
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::Reallocate(int i, int j)
  {
    diag.Reallocate(i);
    T zero;
    SetComplexZero(zero);

    diag.Fill(zero);
    this->m_ = i;
    this->n_ = i;
  }

 
  
  //! changing the size of the matrix
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::Resize(int i, int j)
  {
    diag.Resize(i);
    this->m_ = i;
    this->n_ = i;
  }
  
  
  //! returns A(i, j)
  template<class T, class Prop, class Allocator>
  inline T Matrix<T, Prop, DiagonalRow, Allocator>::operator() (int i, int j) const
  {
    if (i == j) 
      return this->diag(i); 
    
    T zero;
    SetComplexZero(zero);
    
    return zero;
  }
  
  
  //! returns access to A(i, j)
  template<class T, class Prop, class Allocator>
  inline T& Matrix<T, Prop, DiagonalRow, Allocator>::Get(int i, int j)
  {
    if (i == j)
      return this->diag(i);
    
    throw WrongArgument("Matrix::Get(i, j)", "only diagonal value can be modified");
  }
  
  
  //! adds val, if i == j
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    if (i == j)
      this->diag(i) += val;      
  }


  //! sets val, if i == j
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::SetEntry(int i, int j, const T& val)
  {
    if (i == j)
      this->diag(i) = val;      
  }

  
  //! not in the pattern of the matrix => nothing to do
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::AddDistantInteraction(int i, int j, int proc, const T& val)
  {
  }
  
  
  //! not in the pattern of the matrix => nothing to do
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::AddRowDistantInteraction(int i, int j, int proc, const T& val)
  {
  }
  
  
  //! adds val, if i == col(j)
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::AddInteractionRow(int i, int n, const IVect& col,
		      const Vector<T>& val, bool sorted)
  {    
    for (int k = 0; k < n; k++)
      if (i == col(k))
	this->diag(i) += val(k);      
  }
  

  //! extracts a row of the matrix, row = A(i, :)
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>
  ::GetSparseRow(int i, Vector<T, VectSparse>& row) const
  {
    row.Reallocate(1);
    row.Index(0) = i;
    row.Value(0) = this->diag(i);
  }
  

  //! extracts a row of the matrix, row = A(i, :)
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::GetDenseRow(int i, Vector<T>& row) const
  {
    row.Reallocate(this->m_); row.Zero(); row(i) = this->diag(i);
  }


  //! extracts a column of the matrix, col = A(:, i)
  template<class T, class Prop, class Allocator>
  void Matrix<T, Prop, DiagonalRow, Allocator>::GetDenseCol(int i, Vector<T>& col) const
  {
    col.Reallocate(this->m_); col.Zero(); col(i) = this->diag(i);
  }


  //! sets a row of the matrix, A(i, :) = row
  template<class T, class Prop, class Allocator>
  void Matrix<T, Prop, DiagonalRow, Allocator>::SetDenseRow(int i, const Vector<T>& row)
  {
    this->diag(i) = row(i);
  }

  
  //! sets a column of the matrix, A(i, :) = col
  template<class T, class Prop, class Allocator>
  void Matrix<T, Prop, DiagonalRow, Allocator>::SetDenseCol(int i, const Vector<T>& col)
  {
    this->diag(i) = col(i);
  }

  
  //! returns 1
  template<class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, DiagonalRow, Allocator>::GetRowSize(int i) const
  {
    return 1;
  }
  

  //! matrix with i rows
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::Zero()
  {
    T zero;
    SetComplexZero(zero);

    diag.Fill(zero);
  }
  

  //! Releases memory used by the matrix
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::Clear()
  {
    diag.Clear();
    this->m_ = 0;
    this->n_ = 0;
  }


  //! Clears a row of the matrix
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::ClearRow(int i)
  {
    SetComplexZero(diag(i));
  }
  
  
  //! sets the pointer storing entries of the matrix
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::SetData(int n, T* data)
  {
    diag.SetData(n, data);
    this->m_ = n;
    this->n_ = n;
  }
  
  
  //! nullifies the pointer without releasing memory
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::Nullify()
  {
    diag.Nullify();
    this->m_ = 0;
    this->n_ = 0;
  }


  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::WriteText(const string& FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  template<class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, DiagonalRow, Allocator>::WriteText(ostream& FileStream) const
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    // the matrix is written in the same format as sparse matrices
    for (int i = 0; i < this->GetM(); i++)
      FileStream << i+1 << " " << i+1 << " " << (*this)(i, i) << '\n';
    
#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif
  }
  

  //! matrix vector product
  template<class T, class Prop, class Allocator, class AllocatorVect>
  inline void MltVector(const Matrix<T, Prop, DiagonalRow, Allocator>  M,
			const Vector<T, VectFull, AllocatorVect> A, Vector<T, VectFull, AllocatorVect> & B)
  {
    for (int i = 0; i < M.GetM(); i++)
      B(i) = M(i, i)*A(i);    
  }
  
  
  //! matrix vector product y = beta y + alpha A x
  template<class T0, class T1, class T2, class T3, class T4,
           class Prop, class Allocator, class Allocator1,
           class Allocator2, class Allocator3>
  inline void MltAddVector(const T0& alpha, const Matrix<T1, Prop, DiagonalRow, Allocator1>  A,
			   const Vector<T2, VectFull, Allocator2> x, 
			   const T3& beta, Vector<T4, VectFull, Allocator3> & y)
  {
    for (int i = 0; i < A.GetM(); i++)
      y(i) = beta*y(i) + alpha*A(i, i)*x(i);    
  }
  
}

#define MONTJOIE_FILE_DIAGONAL_MATRIX_INLINE_CXX
#endif
  
