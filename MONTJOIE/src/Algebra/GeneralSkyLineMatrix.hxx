#ifndef SELDON_FILE_GENERAL_SKYLINE_MATRIX_HXX

namespace Seldon
{
  
  class RowSkyLine : public RowSparse
  {
  };
  
  //! Skyline matrix
  /*!
    SkyLine matrices stored by rows
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class GeneralSkyLineMatrix : public VirtualMatrix<T>
  {
    
    // Attributes.
  protected:
  
    //! Vector of values
    /*!
      all the non-zero entries stored by rows
     */
    Vector<Vector<T>, VectFull, NewAlloc<Vector<T> > > val;
    //! Beginning of colum numbers for each row
    /*!
      A(i, j) = 0 for all j < first_col(i)
      A(i, j) <> 0 for j = first_col(i)
    */
    Vector<int> first_col;
    
    // Methods.
  public:
    // Constructors.
    GeneralSkyLineMatrix();
    
    // Memory management
    void Clear();
    void Reallocate(int i, int j);

    // Basic methods.    
    void ReallocateRow(int, int, int);
    void AddInteraction(int i, int j, const T& val);
    void AddInteractionRow(int i, int n, const IVect& j,
			   const Vector<T>& val, bool sorted = false);
    
    int GetRowSize(int i) const;
    int GetFirstColNumber(int i) const;
    size_t GetMemorySize() const;
    
    T& Value(int i, int j);
    const T& Value(int i, int j) const;
    
    int GetDataSize() const;
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Zero();
    
    // Element access
    const T operator() (int i, int j) const;
    
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;
    
  };
  
  
  //! General skyline matrix stored by rows
  template<class T, class Prop, class Allocator> 
  class Matrix<T, Prop, RowSkyLine, Allocator>
    : public GeneralSkyLineMatrix <T, Prop, RowSkyLine, Allocator>
  {
  };


  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, RowSkyLine, Allocator>& A);

  
  // Conversion from any matrix to skyline matrix
  template<class T, class Prop, class Storage,
           class Allocator, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
		  Matrix<T, General, RowSkyLine, Allocator2>& B);

  // LU factorisation without pivot
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, RowSkyLine, Allocator>& A);
  
  template<class T, class Allocator, class T2>
  void SolveLuVector(const SeldonTranspose& trans,
		     const Matrix<T, General, RowSkyLine, Allocator>&, Vector<T2>&);
  
  template<class T, class Allocator, class T2>
  void SolveLuVector(const Matrix<T, General, RowSkyLine, Allocator>&,
		     Vector<T2>&);

  template<class T1, class Allocator1, class T2,
           class Allocator2, class T4, class Allocator4>
  void MltVector(const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& x,
		 Vector<T4, VectFull, Allocator4>& y);
  
  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y);

  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& trans,
		    const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y);
    
}

#define SELDON_FILE_GENERAL_SKYLINE_MATRIX_HXX
#endif
