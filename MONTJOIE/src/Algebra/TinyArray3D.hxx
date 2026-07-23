#ifndef SELDON_FILE_TINY_ARRAY3D_HXX

#include "TinyArray3DExpression.hxx"

namespace Seldon
{
  template<int p>
  class TinyArray3DLoop;

  template<int p, int q>
  class TinyArray3DDoubleLoop;

  template<int p, int q, int r>
  class TinyArray3DTripleLoop;
  
  //! class for tiny 3-D arrays, whose dimensions are known at compilation time
  template <class T, int m, int n, int p>
  class TinyArray3D : public TinyArray3DExpression<T, m, n, p, TinyArray3D<T, m, n, p> >
  {
    template<int q>
    friend class TinyArray3DLoop;

    template<int pp, int q>
    friend class TinyArray3DDoubleLoop;

    template<int pp, int q, int r>
    friend class TinyArray3DTripleLoop;
    
  protected: 
    //! all the elements of the tiny 3-D array
    T data_[m*n*p];
    
  public :
    enum{ size_ = m*n*p };
    
    TinyArray3D();
    
    // Basic methods.
    static int GetLength1();
    static int GetLength2();
    static int GetLength3();
    static int GetSize();
    
    T* GetData();
    void FillRand();
    
    // access operator
    T& Val(int i);
    const T& Val(int i) const;
    
    T& operator()(int i, int j, int k);
    const T& operator()(int i, int j, int k) const;

    // operation this = u
    template<class T1, class E>
    TinyArray3D<T, m, n, p>& operator =(const TinyArray3DExpression<T1, m, n, p, E> & u);

    // this += u
    template<class T1, class E>
    TinyArray3D<T, m, n, p>& operator +=(const TinyArray3DExpression<T1, m, n, p, E> & u);

    // this -= u
    template<class T1, class E>
    TinyArray3D<T, m, n, p>& operator -=(const TinyArray3DExpression<T1, m, n, p, E> & u);

    template<class T0>
    const TinyArray3D<T, m, n, p>& operator *=(const T0& alpha);    
    
    // Convenient functions
    void Zero();         // Fills the 3-D array with 0
    void Fill(const T& a); // Fills the 3-D array with a
        
  };
  

  //! class for unrolling simple loops in TinyArray3D functions
  template<int q>
  class TinyArray3DLoop
  {
  public :
    template<int m, int n, int p, class T, class T1>
    static void Fill(TinyArray3D<T, m, n, p>& A, const T1& alpha);
    
    template<int m, int n, int p, class T>
    static void FillRand(TinyArray3D<T, m, n, p>& A);
    
    template<int m, int n, int p, class T, class T0>
    static void Mlt(TinyArray3D<T, m, n, p>& A, const T0& alpha);
    
    template<class T0, int p, int r, int s>
    static void MltVect(const TinyArray3D<T0, p, r, s>& A,
                        const TinyVector<T0, r>& x, TinyMatrix<T0, General, p, s>& y);
	  
    template<class T0, int p, int r, int s>
    static void ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
                              int theta, TinyMatrix<T0, General, p, s>& y);
    
    template<class T0, int p, int r, int s>
    static void ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
                              TinyVector<int,s>& theta, TinyMatrix<T0, General, p, s>& y);

    template<class T0, int m, int n, int s, class E, class T1>
    static void Copy(const TinyArray3DExpression<T0, m, n, s, E>& u,
		     TinyArray3D<T1, m, n, s>& v);
    
    template<class T0, int m, int n, int s, class E, class T1>
    static void AddCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
			TinyArray3D<T1, m, n, s>& v);

    template<class T0, int m, int n, int s, class E, class T1>
    static void DiffCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
			 TinyArray3D<T1, m, n, s>& v);

  };


  //! class terminating simple loops in TinyArray3D functions
  template<>
  class TinyArray3DLoop<0>
  {
  public :
    template<int m, int n, int p, class T, class T1>
    static inline void Fill(TinyArray3D<T, m, n, p>& A, const T1& alpha) {}
	  
    template<int m, int n, int p, class T>
    static inline void FillRand(TinyArray3D<T, m, n, p>& A){}
    
    template<int m, int n, int p, class T, class T0>
    static inline void Mlt(TinyArray3D<T, m, n, p>& A, const T0& alpha) {}
    
    template<class T0, int p, int r, int s>
    static inline void MltVect(const TinyArray3D<T0, p, r, s>& A,
                               const TinyVector<T0, r>& x,
                               TinyMatrix<T0, General, p, s>& y) {}
    
    template<class T0, int p, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
				     int theta, TinyMatrix<T0, General, p, s>& y){}
    
    template<class T0, int p, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
				     TinyVector<int,s>& theta, TinyMatrix<T0, General, p, s>& y){}

    template<class T0, int m, int n, int s, class E, class T1>
    static inline void Copy(const TinyArray3DExpression<T0, m, n, s, E>& u,
			    TinyArray3D<T1, m, n, s>& v) {}
    
    template<class T0, int m, int n, int s, class E, class T1>
    static inline void AddCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
			       TinyArray3D<T1, m, n, s>& v) {}

    template<class T0, int m, int n, int s, class E, class T1>
    static inline void DiffCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
				TinyArray3D<T1, m, n, s>& v) {}
    
  };
  
  
  //! class for unrolling double loops in TinyArray3D functions
  template<int p, int q>
  class TinyArray3DDoubleLoop
  {
  public :
    template<class T0, int m, int n, int s>
    static void MltVect(const TinyArray3D<T0, m, n, s>& A,
                        const TinyVector<T0, n>& x, TinyMatrix<T0, General, m, s>& y);
    
    template<class T0, int t, int r, int s>
    static void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                              int theta, TinyMatrix<T0, General, t, s>& y);
    
    template<class T0, int t, int r, int s>
    static void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                              TinyVector<int,s>& theta, TinyMatrix<T0, General, t, s>& y);
    
  };
  
  //! class terminating double loops in TinyArray3D functions
  template<int p>
  class TinyArray3DDoubleLoop<p, 0>
  {
  public :
    template<class T0, int m, int n, int s>
    static inline void MltVect(const TinyArray3D<T0, m, n, s>& A,
                               const TinyVector<T0, n>& x,
                               TinyMatrix<T0, General, m, s>& y) {}
    
    template<class T0, int t, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A, int theta,
                              TinyMatrix<T0, General, t, s>& y){}
	  
    template<class T0, int t, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                                     TinyVector<int,s>& theta,
                                     TinyMatrix<T0, General, t, s>& y){}
    
  };
  
  
  //! class terminating double loops in TinyArray3D functions
  template<int q>
  class TinyArray3DDoubleLoop<0, q>
  {
  public :
    
    template<class T0, int t, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                                     int theta, TinyMatrix<T0, General, t, s>& y){}
    
    template<class T0, int t, int r, int s>
    static inline void ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                                     TinyVector<int,s>& theta,
                                     TinyMatrix<T0, General, t, s>& y){}
    
  };
  
  
  //! class for unrolling triple loops in TinyArray3D functions
  template<int p, int q, int r>
  class TinyArray3DTripleLoop
  {
  public :
    template<class T0, int m, int n, int s>
    static void MltVect(const TinyArray3D<T0, m, n, s>& A,
                        const TinyVector<T0, n>& x, TinyMatrix<T0, General, m, s>& y);

  };

  
  //! class terminating triple loops in TinyArray3D functions
  template<int p, int q>
  class TinyArray3DTripleLoop<p, q, 0>
  {
  public :
    template<class T0, int m, int n, int s>
    static inline void MltVect(const TinyArray3D<T0, m, n, s>& A,
                               const TinyVector<T0, n>& x,
                               TinyMatrix<T0, General, m, s>& y) {}

  };

  
  // inline functions
  template<class T0, int m, int n, int p>
  TinyMatrix<T0, General, m, p>
  dot(const TinyArray3D<T0, m, n, p>& A, const TinyVector<T0, n>& x);
    
  template<int m, int n, int q, class T>
  void ExtractMatrix(const TinyArray3D<T, n, q, m> & A, int theta,
		     TinyMatrix<T, General, n, m>& y);

  template<int m, int n, int q, class T>
  void ExtractMatrix(const TinyArray3D<T, n, q, m>& A, TinyVector<int, m> & theta,
		     TinyMatrix<T, General, n, m>& y);
  
  // non-inline functions
  template<class T, int m, int n, int p> inline
  ostream & operator << (ostream & out , const TinyArray3D<T, m, n, p> & A);

} // end namespace


#define SELDON_FILE_TINY_ARRAY3D_HXX
#endif

