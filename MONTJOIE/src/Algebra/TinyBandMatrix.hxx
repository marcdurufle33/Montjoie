#ifndef SELDON_FILE_TINY_BAND_MATRIX_HXX

namespace Seldon
{

  /*******************
   * Banded matrices *
   *******************/
  
  
  template<int p>
  class TinyBandMatrixLoop;

  template<int p, int q>
  class TinyBandMatrixDoubleLoop;
  
  //! class for a band matrix whose band width is known at the compilation
  template<class T, int d>
  class TinyBandMatrix : public VirtualMatrix<T>
  {
    template<int p>
    friend class TinyBandMatrixLoop;
    
    template<int p, int q>
    friend class TinyBandMatrixDoubleLoop;

  protected :
    //! storing all the diagonals from lower-diagonal to upper-diagonal 
    //! diag_(:)(d) is then the main diagonal
    Vector<TinyVector<T, 2*d+1> > diag_;
    
  public :
    typedef T entry_type;
    
    int GetDataSize() const;
    size_t GetMemorySize() const;
    
    void Clear();
    void Zero();
    void HideMessages();
    void Reallocate(int m, int n);
    
    void AddInteraction(int i, int j, const T& val);
    void AddInteractionRow(int i, int n, const Vector<int>&,
			   const Vector<T>& val, bool sorted = false);
        
    T operator()(int i, int j) const;

    TinyBandMatrix<T, d>& operator *=(const T& alpha);
    
    T& Get(int i, int j);
    const T& Get(int i, int j) const;

    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    
    void Set(int, int, const T&);
    void SetEntry(int, int, const T&);
    
    void ClearRow(int i);
    
    void SetIdentity();
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Copy(const Matrix<T, General, ArrayRowSparse>& A);
    
    void Factorize();
    
    template<class T0>
    void Add_(const T0& alpha, const TinyBandMatrix<T, d>& A);
    
    template<class T0, class T1>
    void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		const Vector<T1>& x, Vector<T1>& y) const;
    
    template<class T0>
    void Solve(Vector<T0>& x);
    
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    
  };
  
  
  //! class for unrolling simple loops for band matrices
  template<int p>
  class TinyBandMatrixLoop
  {
  public :
    template<class T, int d>
    static void EliminateColumn(TinyBandMatrix<T, d>& A, int j);

    template<class T, int d, class T0>
    static void SolveLower(const TinyBandMatrix<T, d>& A, int j, Vector<T0>& x);
    
    template<class T, int d, class T0>
    static void SolveUpper(const TinyBandMatrix<T, d>& A, int j, Vector<T0>& x);
    
    template<class T, int d, class T0>
    static void MltRow(int i, const T0& coef, TinyBandMatrix<T, d>& A);
    
    template<class T, int size, class T1>
    static void MltAdd(int offset, const TinyVector<T, size>& diag, const T1& xj, Vector<T1>& y);

    template<class T, int size, class T1>
    static void MltTransAdd(int offset, const TinyVector<T, size>& diag,
                            const Vector<T1>& xj, T1& val);

    template<class T, int size, class T1>
    static void MltConjTransAdd(int offset, const TinyVector<T, size>& diag,
                                const Vector<T1>& xj, T1& val);

  };
  
  
  //! class terminating simple loops for band matrices
  template<>
  class TinyBandMatrixLoop<0>
  {
  public :
    template<class T, int d>
    static inline void EliminateColumn(TinyBandMatrix<T, d>& A, int j) {}

    template<class T, int d, class T0>
    static inline void SolveLower(const TinyBandMatrix<T, d>& A, int j, Vector<T0>& x) {}
    
    template<class T, int d, class T0>
    static inline void SolveUpper(const TinyBandMatrix<T, d>& A, int j, Vector<T0>& x) {}

    template<class T, int d, class T0>
    static inline void MltRow(int i, const T0& coef, TinyBandMatrix<T, d>& A) {}

    template<class T, int size, class T1>
    static inline void MltAdd(int offset, const TinyVector<T, size>& diag, const T1& xj, Vector<T1>& y) {}

    template<class T, int size, class T1>
    static inline void MltTransAdd(int offset, const TinyVector<T, size>& diag,
				   const Vector<T1>& xj, T1& val) {}

    template<class T, int size, class T1>
    static inline void MltConjTransAdd(int offset, const TinyVector<T, size>& diag,
				       const Vector<T1>& xj, T1& val) {}
    
  };
  
  
  //! class for unrolling double loops for band matrices
  template<int p, int q>
  class TinyBandMatrixDoubleLoop
  {
  public :
    template<class T, int d>
    static void CombineRow(TinyBandMatrix<T, d>& A, int j, const T& pivot);
    
  };
  
  
  //! class terminating double loops for band matrices
  template<int p>
  class TinyBandMatrixDoubleLoop<p, 0>
  {
  public :
    template<class T, int d>
    static inline void CombineRow(TinyBandMatrix<T, d>& A, int j, const T& pivot) {}
    
  };
  
  // functions  for banded matrices

  template<class T, int d>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& row_sum, const TinyBandMatrix<T, d>& A);
      
  template<class T, int d>
  void ScaleLeftMatrix(TinyBandMatrix<T, d>& A, const Vector<typename ClassComplexType<T>::Treal>& scale);
  
  template<class T, int d>
  void GetLU(TinyBandMatrix<T, d>& A,
             TinyBandMatrix<T, d>& mat_lu, bool keep_matrix = false);

  template<class T, int d>
  void GetLU(TinyBandMatrix<T, d>& A);

  template<class T, int d>
  void Copy(const Matrix<T, General, ArrayRowSparse>& A, TinyBandMatrix<T, d>& B);

  template<class T, int d>
  void SolveLU(TinyBandMatrix<T, d>& mat_lu, Vector<T>& x);

  template<class T0, class T1, int d>
  void Add(const T0& alpha, const TinyBandMatrix<T1, d>& A, TinyBandMatrix<T1, d>& B);

  template<class T0, class T, int d, class T2>
  void MltAdd(const T0& alpha, const TinyBandMatrix<T, d>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y);
  
  template<class T0, class T, int d, class T2>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans, const TinyBandMatrix<T, d>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y);
  
  template<class T, int d, class T1>
  void Mlt(const SeldonTranspose& trans, const TinyBandMatrix<T, d>& A,
	   const Vector<T1>& x, Vector<T1>& y);
  
  template<class T, int d, class T1>
  void Mlt(const TinyBandMatrix<T, d>& A,
	   const Vector<T1>& x, Vector<T1>& y);

  template<class T, int d>
  void Mlt(const T& alpha, TinyBandMatrix<T, d>& A);
  
  template<class T, int d>
  ostream& operator<<(ostream& out, const TinyBandMatrix<T, d>& A);
  
  
  /******************
   * Arrow matrices *
   ******************/
  
  
  template<int p>
  class TinyArrowMatrixLoop;
  
  template<int p, int q>
  class TinyArrowMatrixDoubleLoop;

  //! class for an arrow matrix whose band width and size of last columns/rows
  //! are known at the compilation
  template<class T, int d, int m>
  class TinyArrowMatrix : public TinyBandMatrix<T, d>
  {
    template<int p>
    friend class TinyArrowMatrixLoop;

    template<int p, int q>
    friend class TinyArrowMatrixDoubleLoop;

    template<class T0, int p, int q>
    friend class TinyArrowMatrix;
    
  protected :
    //! last columns/rows of the matrix
    Vector<TinyVector<T, m> > last_row_, last_col_;
    //! last tiny block 
    TinyMatrix<T, General, m, m> last_block_;
    
  public :
    int GetM() const;
    int GetN() const;
    int GetDataSize() const;
    size_t GetMemorySize() const;
    
    void Clear();
    void Zero();
    void Reallocate(int, int);
    
    void AddInteraction(int i, int j, const T& val);
    void AddInteractionRow(int i, int nb, const IVect& j,
			   const Vector<T>& val, bool sorted = false);
    
    T operator()(int i, int j) const;    

    TinyArrowMatrix<T, d, m>& operator *=(const T& alpha);
    
    T& Get(int i, int j);
    const T& Get(int i, int j) const;

    T& Val(int i, int j);
    const T& Val(int i, int j) const;
    
    void Set(int, int, const T&);
    void SetEntry(int, int, const T&);
    
    template<int d2, int m2>
    const TinyArrowMatrix<T, d, m>& operator=(const TinyArrowMatrix<T, d2, m2>& A);
    
    template<class T0>
    const TinyArrowMatrix<T, d, m>& operator *=(const T0& x);
    
    void ClearRow(int i);
    
    void SetIdentity();
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    void Copy(const Matrix<T, General, ArrayRowSparse>& A);
    
    void Factorize();
    
    template<class T0, int d2, int m2>
    void Add_(const T0& alpha, const TinyArrowMatrix<T, d2, m2>& A);
    
    template<class T0>
    void Add_(const T0& alpha, const TinyArrowMatrix<T, d, m>& A);
    
    template<class T0, class T1>
    void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		const Vector<T1>& x, Vector<T1>& y) const;
    
    void Solve(Vector<T>& x);
    
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    
  };
  
    
  //! class for unrolling simple loops for arrow matrices
  template<int p>
  class TinyArrowMatrixLoop
  {
  public :
    template<class T, int d, int m>
    static void EliminateLastColumn(TinyArrowMatrix<T, d, m>& A, int j);

    template<class T, int d, int m>
    static void EliminateLastRow(TinyArrowMatrix<T, d, m>& A, int j);
    
    template<class T, int d, int m>
    static void SolveLower(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x);
    
    template<class T, int d, int m>
    static void SolveRow(const TinyArrowMatrix<T, d, m>& A, int j, int n, Vector<T>& x);
    
    template<class T, int d, int m>
    static void SolveColumn(const TinyArrowMatrix<T, d, m>& A, int j, int n, Vector<T>& x);
    
    template<class T, int d, int m>
    static void SolveUpper(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x);
    
    template<class T, int m, class T1>
    static void MltAdd(int n, const TinyVector<T, m>& coef, const Vector<T1>& x, T1& yj);
    
    template<class T, int m, class T1>
    static void MltConjAdd(int n, const TinyVector<T, m>& coef, const Vector<T1>& x, T1& yj);

    template<class T, int m, class T1>
    static void MltTransAdd(int n, const TinyVector<T, m>& coef, const T1& xj, Vector<T1>& y);

    template<class T, int m, class T1>
    static void MltConjTransAdd(int n, const TinyVector<T, m>& coef, const T1& xj, Vector<T1>& y);
    
  };


  //! class terminating simple loops for band matrices  
  template<>
  class TinyArrowMatrixLoop<0>
  {
  public :
    template<class T, int d, int m>
    static inline void EliminateLastColumn(TinyArrowMatrix<T, d, m>& A, int j) {}
    
    template<class T, int d, int m>
    static inline void EliminateLastRow(TinyArrowMatrix<T, d, m>& A, int j) {}
    
    template<class T, int d, int m>
    static inline void SolveLower(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x) {}
    
    template<class T, int d, int m>
    static inline void SolveRow(const TinyArrowMatrix<T, d, m>& A,
                                int j, int n, Vector<T>& x) {}
    
    template<class T, int d, int m>
    static inline void SolveColumn(const TinyArrowMatrix<T, d, m>& A,
                                   int j, int n, Vector<T>& x) {}

    template<class T, int d, int m>
    static inline void SolveUpper(const TinyArrowMatrix<T, d, m>& A, int j, Vector<T>& x) {}

    template<class T, int m, class T1>
    static inline void MltAdd(int n, const TinyVector<T, m>& coef, const Vector<T1>& x, T1& yj) {}

    template<class T, int m, class T1>
    static inline void MltConjAdd(int n, const TinyVector<T, m>& coef,
                                  const Vector<T1>& x, T1& yj) {}

    template<class T, int m, class T1>
    static inline void MltTransAdd(int n, const TinyVector<T, m>& coef,
                                   const T1& xj, Vector<T1>& y) {}

    template<class T, int m, class T1>
    static inline void MltConjTransAdd(int n, const TinyVector<T, m>& coef,
                                       const T1& xj, Vector<T1>& y) {}
    
  };
  
  
  //! class for unrolling double loops for arrow matrices
  template<int p, int q>
  class TinyArrowMatrixDoubleLoop
  {
  public :
    template<class T, int d, int m>
    static void CombineLastColumn(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot);
    
    template<class T, int d, int m>
    static void CombineLastRow(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot);

    template<class T, int d, int m>
    static void CombineLastBlock(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot);
    
  };
  
  
  //! class terminating double loops for arrow matrices
  template<int p>
  class TinyArrowMatrixDoubleLoop<p, 0>
  {
  public :
    template<class T, int d, int m>
    static inline void CombineLastColumn(TinyArrowMatrix<T, d, m>& A,
                                         int j, const T& pivot) {}
    
    template<class T, int d, int m>
    static inline void CombineLastRow(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot) {}

    template<class T, int d, int m>
    static inline void CombineLastBlock(TinyArrowMatrix<T, d, m>& A, int j, const T& pivot) {}
    
  };

  template<class T, int d, int m>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& row_sum, const TinyArrowMatrix<T, d, m>& A);
      
  template<class T, int d, int m>
  void ScaleLeftMatrix(TinyArrowMatrix<T, d, m>& A, const Vector<typename ClassComplexType<T>::Treal>& scale);
 
  template<class T, int d, int m>
  void GetLU(TinyArrowMatrix<T, d, m>& A,
             TinyArrowMatrix<T, d, m>& mat_lu, bool keep_matrix = false);

  template<class T, int d, int m>
  void GetLU(TinyArrowMatrix<T, d, m>& A);

  template<class T, int d, int m>
  void SolveLU(TinyArrowMatrix<T, d, m>& mat_lu, Vector<T>& x);

  template<class T, int d, int m>
  void Copy(const Matrix<T, General, ArrayRowSparse>& A, TinyArrowMatrix<T, d, m>& B);

  template<class T0, class T1, int d, int m, int d2, int m2>
  void Add(const T0& alpha, const TinyArrowMatrix<T1, d, m>& A,
	   TinyArrowMatrix<T1, d2, m2>& B);
  
  template<class T0, class T, int d, int m, class T2>
  void MltAdd(const T0& alpha, const TinyArrowMatrix<T, d, m>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y);
  
  template<class T0, class T, int d, int m, class T2>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans, const TinyArrowMatrix<T, d, m>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y);
  
  template<class T, int d, int m, class T2>
  void Mlt(const SeldonTranspose& trans, const TinyArrowMatrix<T, d, m>& A,
	   const Vector<T2>& x, Vector<T2>& y);
  
  template<class T, int d, int m, class T2>
  void Mlt(const TinyArrowMatrix<T, d, m>& A,
           const Vector<T2>& x, Vector<T2>& y);
  
  template<class T, int d, int m>
  void Mlt(const T& alpha, TinyArrowMatrix<T, d, m>& A);
  
  template<class T, int d, int m>
  ostream& operator<<(ostream& out, const TinyArrowMatrix<T, d, m>& A);
  
}

#define SELDON_FILE_TINY_BAND_MATRIX_HXX
#endif

