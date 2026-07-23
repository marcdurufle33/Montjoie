#ifndef SELDON_FILE_TINY_ARRAY3D_INLINE_CXX

#include "TinyArray3D.hxx"

#include "TinyArray3DExpressionInline.cxx"

namespace Seldon
{  

  /***************
   * TinyArray3D *
   ***************/

  
  //! sets all elements to 0
  template <class T, int m, int n, int p>
  inline TinyArray3D<T, m, n, p>::TinyArray3D()
  {
    Zero();
  }
  
  
  //! returns first dimension of 3-D array
  template <class T, int m, int n, int p>
  inline int TinyArray3D<T, m, n, p>::GetLength1()
  {
    return m;
  }
  
  
  //! returns second dimension of 3-D array
  template <class T, int m, int n, int p>
  inline int TinyArray3D<T, m, n, p>::GetLength2()
  {
    return n;
  }
  
  
  //! returns third dimension of 3-D array
  template <class T, int m, int n, int p>
  inline int TinyArray3D<T, m, n, p>::GetLength3()
  {
    return p;
  }
  
  
  //! returns the number of element stored in the array
  template <class T, int m, int n, int p>
  inline int TinyArray3D<T, m, n, p>::GetSize()
  {
    return size_;
  }
  
  
  //! returns pointer to the array containing all the elements
  template <class T, int m, int n, int p>
  inline T* TinyArray3D<T, m, n, p>::GetData()
  {
    return data_;
  }
  
  
  //! sets randomly all elements of the 3-D array
  template <class T, int m, int n, int p>
  inline void TinyArray3D<T, m, n, p>::FillRand()
  {
    TinyArray3DLoop<size_>::FillRand(*this);
  }


  //! returns i-th element stored
  template <class T, int m, int n, int p>
  inline T& TinyArray3D<T, m, n, p>::Val(int i)
  {
    return data_[i];
  }


  //! returns i-th element stored
  template <class T, int m, int n, int p>
  inline const T& TinyArray3D<T, m, n, p>::Val(int i) const
  {
    return data_[i];
  }

  
  //! returns A(i, j, k)
  template <class T, int m, int n, int p>
  inline T& TinyArray3D<T, m, n, p>::operator()(int i, int j, int k)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= m)
      throw WrongRow("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(m-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >=  n)
      throw WrongCol("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(n-1) + "], but is equal to "
		     + to_str(j) + ".");

    if (k < 0 || k >=  p)
      throw WrongCol("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(p-1) + "], but is equal to "
		     + to_str(k) + ".");
#endif

    return data_[p*(n*i + j) + k]; 
  }
  
  
  //! returns A(i, j, k)
  template <class T, int m, int n, int p> inline
  const T& TinyArray3D<T, m, n, p>::operator()(int i, int j, int k) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= m)
      throw WrongRow("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(m-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >=  n)
      throw WrongCol("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(n-1) + "], but is equal to "
		     + to_str(j) + ".");

    if (k < 0 || k >=  p)
      throw WrongCol("TinyArray3D::operator()",
		     string("Index should be in [0, ")
		     + to_str(p-1) + "], but is equal to "
		     + to_str(k) + ".");
#endif
    
    return data_[p*(n*i + j) + k]; 
  }


  //! operation this = u
  template <class T, int m, int n, int p> template<class T1, class E>
  TinyArray3D<T, m, n, p>& TinyArray3D<T, m, n, p>
  ::operator =(const TinyArray3DExpression<T1, m, n, p, E> & u)
  {
    TinyArray3DLoop<size_>::Copy(u, *this);
    return *this;
  }
  
  
  //! this += u
  template <class T, int m, int n, int p> template<class T1, class E>
  TinyArray3D<T, m, n, p>& TinyArray3D<T, m, n, p>
  ::operator +=(const TinyArray3DExpression<T1, m, n, p, E> & u)
  {
    TinyArray3DLoop<size_>::AddCopy(u, *this);
    return *this;
  }
  

  //! this -= u
  template <class T, int m, int n, int p> template<class T1, class E>
  TinyArray3D<T, m, n, p>& TinyArray3D<T, m, n, p>
  ::operator -=(const TinyArray3DExpression<T1, m, n, p, E> & u)
  {
    TinyArray3DLoop<size_>::DiffCopy(u, *this);
    return *this;
  }

  
  //! multiplies all the elements by alpha
  template <class T, int m, int n, int p> template<class T0>
  inline const TinyArray3D<T, m, n, p>& TinyArray3D<T, m, n, p>::operator*=(const T0& alpha)
  {
    TinyArray3DLoop<size_>::Mlt(*this, alpha);
    return *this;
  }

  
  //! sets all elements to 0
  template <class T, int m, int n, int p>
  inline void TinyArray3D<T, m, n, p>::Zero()
  {
    T zero; SetComplexZero(zero);
    TinyArray3DLoop<size_>::Fill(*this, zero);
  }
  
  
  //! sets all elements to a given value
  template <class T, int m, int n, int p>
  inline void TinyArray3D<T, m, n, p>::Fill(const T& a)
  {
    TinyArray3DLoop<size_>::Fill(*this, a);
  }
  
  
  /*******************
   * TinyArray3DLoop *
   *******************/

  
  //! sets all elements to a given value
  template<int q> template<int m, int n, int p, class T, class T1>
  inline void TinyArray3DLoop<q>::Fill(TinyArray3D<T, m, n, p>& A, const T1& alpha)
  {
    SetComplexReal(alpha, A.data_[q-1]);
    TinyArray3DLoop<q-1>::Fill(A, alpha);
  }
  

  //! sets randomly all elements of the array
  template<int q> template<int m, int n, int p, class T>
  inline void TinyArray3DLoop<q>::FillRand(TinyArray3D<T, m, n, p>& A)
  {
    A.data_[q-1] = rand();
    TinyArray3DLoop<q-1>::FillRand(A);
  }
  
  
  //! multiplies A by alpha
  template<int q> template<int m, int n, int p, class T, class T0>
  inline void TinyArray3DLoop<q>::Mlt(TinyArray3D<T, m, n, p>& A, const T0& alpha)
  {
    A.data_[q-1] *= alpha;
    TinyArray3DLoop<q-1>::Mlt(A, alpha);
  }
  
  
  //! computes y with y(i, j) = \sum_k A(i, k, j) x(k)
  template<int q> template<class T0, int p, int r, int s>
  inline void TinyArray3DLoop<q>::
  MltVect(const TinyArray3D<T0, p, r, s>& A,
          const TinyVector<T0, r>& x, TinyMatrix<T0, General, p, s>& y)
  {
    TinyArray3DDoubleLoop<q, s>::MltVect(A, x, y);
    TinyArray3DLoop<q-1>::MltVect(A, x, y);
  }
  
  
  //! extracts matrix y with y(:,:) = A(:, theta, :)
  template<int q> template<class T0, int p, int r, int s>
  inline void TinyArray3DLoop<q>::
  ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
                int theta, TinyMatrix<T0, General, p, s>& y)
  {
    TinyArray3DDoubleLoop<p, q>::ExtractMatrix(A, theta, y);
    TinyArray3DLoop<q-1>::ExtractMatrix(A, theta, y);
  }
  

  //! extracts matrix y with y(:,q) = A(:, theta(q), q)
  template<int q> template<class T0, int p, int r, int s>
  inline void TinyArray3DLoop<q>::
  ExtractMatrix(const TinyArray3D<T0, p, r, s>& A,
                TinyVector<int, s>& theta, TinyMatrix<T0, General, p, s>& y)
  {
    TinyArray3DDoubleLoop<p, q>::ExtractMatrix(A, theta, y);
    TinyArray3DLoop<q-1>::ExtractMatrix(A, theta, y);
  }


  //! v = u
  template<int q> template<class T0, int m, int n, int s, class E, class T1>
  inline void TinyArray3DLoop<q>::Copy(const TinyArray3DExpression<T0, m, n, s, E>& u,
				       TinyArray3D<T1, m, n, s>& v)
  {
    v.data_[q-1] = u.Val(q-1);
    TinyArray3DLoop<q-1>::Copy(u, v);
  }
  
  
  //! v = v + u
  template<int q> template<class T0, int m, int n, int s, class E, class T1>
  inline void TinyArray3DLoop<q>::AddCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
					  TinyArray3D<T1, m, n, s>& v)
  {
    v.data_[q-1] += u.Val(q-1);
    TinyArray3DLoop<q-1>::AddCopy(u, v);
  }


  //! v = v - u
  template<int q> template<class T0, int m, int n, int s, class E, class T1>
  inline void TinyArray3DLoop<q>::DiffCopy(const TinyArray3DExpression<T0, m, n, s, E>& u,
					   TinyArray3D<T1, m, n, s>& v)
  {
    v.data_[q-1] -= u.Val(q-1);
    TinyArray3DLoop<q-1>::DiffCopy(u, v);
  }
  
  
  /*************************
   * TinyArray3DDoubleLoop *
   *************************/
  
  
  //! computes y with y(i, j) = \sum_k A(i, k, j) x(k)
  template<int p, int q> template<class T0, int m, int n, int s>
  inline void TinyArray3DDoubleLoop<p, q>::
  MltVect(const TinyArray3D<T0, m, n, s>& A, const TinyVector<T0, n>& x,
          TinyMatrix<T0, General, m, s>& y)
  {
    TinyArray3DTripleLoop<p, q, n>::MltVect(A, x, y);
    TinyArray3DDoubleLoop<p, q-1>::MltVect(A, x, y);
  }


  //! computes y with y(i, j) = \sum_k A(i, k, j) x(k)
  template<int p, int q, int r> template<class T0, int m, int n, int s>
  inline void TinyArray3DTripleLoop<p, q, r>::
  MltVect(const TinyArray3D<T0, m, n, s>& A, const TinyVector<T0, n>& x,
          TinyMatrix<T0, General, m, s>& y)
  {
    y(p-1, q-1) += A(p-1, r-1, q-1)*x(r-1);
    TinyArray3DTripleLoop<p, q, r-1>::MltVect(A, x, y);
  }
  
  
  //! extracts matrix y with y(:,:) = A(:, theta, :)
  template<int p, int q> template<class T0, int t, int r, int s>
  inline void TinyArray3DDoubleLoop<p, q>::
  ExtractMatrix(const TinyArray3D<T0, t, r, s>& A,
                int theta, TinyMatrix<T0, General, t, s>& y)
  {
    y(p-1, q-1) = A(p-1, theta, q-1);
    TinyArray3DDoubleLoop<p-1, q>::ExtractMatrix(A, theta, y);
  }
  

  //! extracts matrix y with y(:,q) = A(:, theta(q), q)  
  template<int p, int q> template<class T0, int t, int r, int s>
  inline void TinyArray3DDoubleLoop<p, q>::
  ExtractMatrix(const TinyArray3D<T0, t, r, s>& A, TinyVector<int, s>& theta,
                TinyMatrix<T0, General, t, s>& y)
  {
    y(p-1, q-1) = A(p-1, theta(q-1), q-1);
    TinyArray3DDoubleLoop<p-1, q>::ExtractMatrix(A, theta, y);
  }


  /*************
   * Functions *
   *************/
  

  //! returns matrix B with B(i, j) = \sum_k A(i, k, j) x(k)
  template<class T0, int m, int n, int p>
  inline TinyMatrix<T0, General, m, p>
  dot(const TinyArray3D<T0, m, n, p>& A, const TinyVector<T0, n>& x)
  {
    TinyMatrix<T0, General, m, p> y;
    TinyArray3DLoop<m>::MltVect(A, x, y);
    return y;
  }


  //! extracts matrix y with y(:,:) = A(:, theta, :)
  template<int m, int n, int q, class T>
  inline void ExtractMatrix(const TinyArray3D<T, n, q, m> & A, int theta,
			    TinyMatrix<T, General, n, m>& y)
  {
    TinyArray3DLoop<m>::ExtractMatrix(A, theta, y);
  }
  

  //! extracts matrix y with y(:,q) = A(:, theta(q), q)  
  template<int m, int n, int q, class T>
  inline void ExtractMatrix(const TinyArray3D<T, n, q, m>& A, TinyVector<int, m> & theta,
			    TinyMatrix<T, General, n, m>& y)
  {
    TinyArray3DLoop<m>::ExtractMatrix(A, theta, y);
  }

} // end namespace

#define SELDON_FILE_TINY_ARRAY3D_INLINE_CXX
#endif
