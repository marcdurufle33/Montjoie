#ifndef SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_HXX

namespace Seldon
{
  
  class SymColSkyLine : public ColSparse
  {
  };
  
  //! Skyline matrix
  /*!
    SkyLine matrices are defined by: (1) the vector with the skyline values;
    (2) the vector with the indices of the first non-zero value for each row in (1)
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class SymmetricSkyLineMatrix : public VirtualMatrix<T>
  {
    
    // Attributes.
  protected:
  
    //! Vector of values
    /*!
      Upper part of the matrix is stored by columns
      (i.e. Lower part of the matrix is stored by rows)
     */
    Vector<T> val;
    //! Vector of indices
    /*!
      Column j is comprised between val(ind(j)) and val(ind(j+1)-1)
     */
    Vector<int> ind;
    
    // Methods.
  public:
    // Constructors.
    SymmetricSkyLineMatrix();
 
    void Clear();
    
    // Memory management.
    void Reallocate(int i, int j);
    
    // Basic methods.
    T* GetData() const;
    int* GetInd() const;
    
    void SetData(IVect&, Vector<T>&);

    size_t GetMemorySize() const;
    int GetDataSize() const;
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Zero();
    
    // Element access
    T operator() (int i, int j) const;
    
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;
    
  };
  
  //! Symmetric skyline matrix stored by columns
  template<class T, class Prop, class Allocator> 
  class Matrix<T, Prop, SymColSkyLine, Allocator>
    : public SymmetricSkyLineMatrix <T, Prop, SymColSkyLine, Allocator>
  {
  };

  
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, SymColSkyLine, Allocator>& A);
  
  // Conversion from any matrix to skyline matrix
  template<class T, class Prop, class Storage,
           class Allocator, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
		  Matrix<T, Symmetric, SymColSkyLine, Allocator2>& B);

  // Cholesky factorisation
  template<class T, class Allocator>
  void GetCholesky(Matrix<T, Symmetric, SymColSkyLine, Allocator>& A);
  
  template<class T, class Allocator, class Alloc2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const Matrix<T, Symmetric, SymColSkyLine, Allocator>&,
                     Vector<T, VectFull, Alloc2>&);
    
  template<class T, class Allocator, class Alloc2>
  void MltCholesky(const SeldonTranspose& trans,
                   const Matrix<T, Symmetric, SymColSkyLine, Allocator>&,
                   Vector<T, VectFull, Alloc2>&);
  
  template <class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, SymColSkyLine, Allocator>& A);
  
  template<class T, class Allocator, class T2>
  void SolveLuVector(const SeldonTranspose& trans,
		     const Matrix<T, Symmetric, SymColSkyLine, Allocator>& A,
		     Vector<T2>& x);

  template<class T1, class Allocator1, class T2,
           class Allocator2, class T4, class Allocator4>
  void MltVector(const Matrix<T1, Symmetric, SymColSkyLine, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& x,
		 Vector<T4, VectFull, Allocator4>& y);
  
  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Symmetric, SymColSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y);
    
}

#define SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_HXX
#endif
