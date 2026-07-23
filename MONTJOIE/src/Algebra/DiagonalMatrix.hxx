#ifndef MONTJOIE_FILE_DIAGONAL_MATRIX_HXX

namespace Seldon
{
  
  class DiagonalRow
  {
  };
  
  //! class implementing diagonal matrices
  template<class T, class Prop, class Allocator>
  class Matrix<T, Prop, DiagonalRow, Allocator> : public VirtualMatrix<T>
  {
  protected:
    Vector<T> diag;

  public :
    typedef T entry_type;
    
    Matrix();
    explicit Matrix(int i, int j);

    size_t GetMemorySize() const ;
    
    void Reallocate(int i, int j);
    void Resize(int i, int j);

    T operator() (int i, int j) const;
    T& Get(int i, int j);

    void SetEntry(int i, int j, const T& val);
    void AddInteraction(int i, int j, const T& val);
    void AddDistantInteraction(int i, int j, int proc, const T& val);
    void AddRowDistantInteraction(int i, int j, int proc, const T& val);
    void AddInteractionRow(int i, int n, const IVect& col, const Vector<T>& val,
			   bool sorted = false);
    
    void GetSparseRow(int i, Vector<T, VectSparse>& row) const;
    void GetDenseRow(int i, Vector<T>& row) const;
    void GetDenseCol(int i, Vector<T>& col) const;
    void SetDenseRow(int i, const Vector<T>& row);
    void SetDenseCol(int i, const Vector<T>& col);

    int GetRowSize(int i) const;

    void Zero();
    void Clear();
    void ClearRow(int i);
    void SetData(int n, T* data);
    void Nullify();
    
    void WriteText(const string& file_name) const;
    void WriteText(ostream& FileStream) const;
    
    inline void RemoveSmallEntry(typename ClassComplexType<T>::Treal eps){}
    inline bool IsSymmetric() const { return true; }
    inline bool IsSparse() const { return true; }

  };


  template<class T, class Prop, class Allocator, class AllocatorVect>
  void MltVector(const Matrix<T, Prop, DiagonalRow, Allocator>  M,
		 const Vector<T, VectFull, AllocatorVect> A, Vector<T, VectFull, AllocatorVect> & B);

  template<class T0, class T1, class T2, class T3, class T4,
           class Prop, class Allocator, class Allocator1,
           class Allocator2, class Allocator3>
  void MltAddVector(const T0& alpha, const Matrix<T1, Prop, DiagonalRow, Allocator1>  A,
		    const Vector<T2, VectFull, Allocator2> x, 
		    const T3& beta, Vector<T4, VectFull, Allocator3> & y);
  
}

#define MONTJOIE_FILE_DIAGONAL_MATRIX_HXX
#endif
  
