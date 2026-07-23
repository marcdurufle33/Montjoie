#ifndef SELDON_FILE_TINY_ARRAY3D_EXPRESSION_INLINE_CXX

namespace Seldon
{
  
  //! returns the element i of expression
  template<class T, int m, int n, int p, class E>
  inline const T TinyArray3DExpression<T, m, n, p, E>::operator()(int i, int j, int k) const
  {
    return static_cast<const E&>(*this)(i, j, k);
  }


  //! returns the element i of expression
  template<class T, int m, int n, int p, class E>
  inline const T TinyArray3DExpression<T, m, n, p, E>::Val(int i) const
  {
    return static_cast<const E&>(*this).Val(i);
  }

  
  //! Constructor u-v with two expressions u and v
  template<class T, int m, int n, int p, class E1, class E2> template<class T1, class T2>
  inline TinyArray3DDifference<T, m, n, p, E1, E2>::
  TinyArray3DDifference(const TinyArray3DExpression<T1, m, n, p, E1>& u,
			const TinyArray3DExpression<T2, m, n, p, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the difference
  template<class T, int m, int n, int p, class E1, class E2> 
  inline const T TinyArray3DDifference<T, m, n, p, E1, E2>::operator()(int i, int j, int k) const
  {
    return u_(i, j, k) - v_(i, j, k);
  }


  //! returns the i-th element of the difference
  template<class T, int m, int n, int p, class E1, class E2> 
  inline const T TinyArray3DDifference<T, m, n, p, E1, E2>::Val(int i) const
  {
    return u_.Val(i) - v_.Val(i);
  }

  
  //! Constructor u+v with two expressions u and v
  template<class T, int m, int n, int p, class E1, class E2> template<class T1, class T2>
  inline TinyArray3DSum<T, m, n, p, E1, E2>
  ::TinyArray3DSum(const TinyArray3DExpression<T1, m, n, p, E1>& u,
		   const TinyArray3DExpression<T2, m, n, p, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the sum  
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DSum<T, m, n, p, E1, E2>::operator()(int i, int j, int k) const
  {
    return u_(i, j, k) + v_(i, j, k);
  }


  //! returns the i-th element of the sum  
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DSum<T, m, n, p, E1, E2>::Val(int i) const
  {
    return u_.Val(i) + v_.Val(i);
  }


  //! Constructor u*v with two expressions u and v
  template<class T, int m, int n, int p, class E1, class E2> template<class T1, class T2>
  inline TinyArray3DProduct<T, m, n, p, E1, E2>
  ::TinyArray3DProduct(const TinyArray3DExpression<T1, m, n, p, E1>& u,
		       const TinyArray3DExpression<T2, m, n, p, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise product 
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DProduct<T, m, n, p, E1, E2>::operator()(int i, int j, int k) const
  {
    return u_(i, j, k) * v_(i, j, k);
  }


  //! returns the i-th element of the element-wise product 
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DProduct<T, m, n, p, E1, E2>::Val(int i) const
  {
    return u_.Val(i) * v_.Val(i);
  }


  //! Constructor u / v with two expressions u and v
  template<class T, int m, int n, int p, class E1, class E2> template<class T1, class T2>
  inline TinyArray3DDivision<T, m, n, p, E1, E2>
  ::TinyArray3DDivision(const TinyArray3DExpression<T1, m, n, p, E1>& u,
			const TinyArray3DExpression<T2, m, n, p, E2>& v)
    : u_(u), v_(v)
  {
  }
  
  
  //! returns the i-th element of the element-wise division
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DDivision<T, m, n, p, E1, E2>::operator()(int i, int j, int k) const
  {
    return u_(i, j, k) / v_(i, j, k);
  }


  //! returns the i-th element of the element-wise division
  template<class T, int m, int n, int p, class E1, class E2>
  inline const T TinyArray3DDivision<T, m, n, p, E1, E2>::Val(int i) const
  {
    return u_.Val(i) / v_.Val(i);
  }


  //! Constructor alpha * u with a scalar alpha and an expression u
  template<class T, int m, int n, int p, class T0, class E> template<class T1>
  inline TinyArray3DScaled<T, m, n, p, T0, E>
  ::TinyArray3DScaled(const T0& alpha,
		      const TinyArray3DExpression<T1, m, n, p, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha*u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScaled<T, m, n, p, T0, E>::operator()(int i, int j, int k) const
  {
    return alpha_*u_(i, j, k);
  }


  //! returns the i-th element of alpha*u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScaled<T, m, n, p, T0, E>::Val(int i) const
  {
    return alpha_*u_.Val(i);
  }


  //! Constructor alpha / u with a scalar alpha and an expression u
  template<class T, int m, int n, int p, class T0, class E> template<class T1>
  inline TinyArray3DScalDiv<T, m, n, p, T0, E>
  ::TinyArray3DScalDiv(const T0& alpha,
		       const TinyArray3DExpression<T1, m, n, p, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha / u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalDiv<T, m, n, p, T0, E>::operator()(int i, int j, int k) const
  {
    return alpha_ / u_(i, j, k);
  }


  //! returns the i-th element of alpha / u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalDiv<T, m, n, p, T0, E>::Val(int i) const
  {
    return alpha_ / u_.Val(i);
  }


  //! Constructor alpha + u with a scalar alpha and an expression u
  template<class T, int m, int n, int p, class T0, class E> template<class T1>
  inline TinyArray3DScalSum<T, m, n, p, T0, E>
  ::TinyArray3DScalSum(const T0& alpha,
		       const TinyArray3DExpression<T1, m, n, p, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha + u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalSum<T, m, n, p, T0, E>::operator()(int i, int j, int k) const
  {
    return alpha_ + u_(i, j, k);
  }


  //! returns the i-th element of alpha + u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalSum<T, m, n, p, T0, E>::Val(int i) const
  {
    return alpha_ + u_.Val(i);
  }


  //! Constructor alpha - u with a scalar alpha and an expression u
  template<class T, int m, int n, int p, class T0, class E> template<class T1>
  inline TinyArray3DScalDiff<T, m, n, p, T0, E>
  ::TinyArray3DScalDiff(const T0& alpha,
			const TinyArray3DExpression<T1, m, n, p, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the i-th element of alpha - u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalDiff<T, m, n, p, T0, E>::operator()(int i, int j, int k) const
  {
    return alpha_ - u_(i, j, k);
  }


  //! returns the i-th element of alpha - u
  template<class T, int m, int n, int p, class T0, class E>
  inline const T TinyArray3DScalDiff<T, m, n, p, T0, E>::Val(int i) const
  {
    return alpha_ - u_.Val(i);
  }
  

  //! Constructor -u with an expression u
  template<class T, int m, int n, int p, class E>
  inline TinyArray3DOpposite<T, m, n, p, E>::TinyArray3DOpposite(const TinyArray3DExpression<T, m, n, p, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the i-th element of -u
  template<class T, int m, int n, int p, class E>
  inline const T TinyArray3DOpposite<T, m, n, p, E>::operator()(int i, int j, int k) const
  {
    return -u_(i, j, k);
  }


  //! returns the i-th element of -u
  template<class T, int m, int n, int p, class E>
  inline const T TinyArray3DOpposite<T, m, n, p, E>::Val(int i) const
  {
    return -u_.Val(i);
  }

  
  /*************
   * Operators *
   *************/
  
  
  //! returns u+v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DSum<T, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DSum<T, m, n, p, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DSum<complex<T>, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DSum<complex<T>, m, n, p, E1, E2>(u, v);
  }
  

  //! returns u+v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DSum<complex<T>, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v)
  {
    return TinyArray3DSum<complex<T>, m, n, p, E1, E2>(u, v);
  }
 

  //! returns u-v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDifference<T, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DDifference<T, m, n, p, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDifference<complex<T>, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DDifference<complex<T>, m, n, p, E1, E2>(u, v);
  }


  //! returns u-v
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDifference<complex<T>, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v)
  {
    return TinyArray3DDifference<complex<T>, m, n, p, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DProduct<T, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DProduct<T, m, n, p, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DProduct<complex<T>, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DProduct<complex<T>, m, n, p, E1, E2>(u, v);
  }


  //! returns u*v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DProduct<complex<T>, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v)
  {
    return TinyArray3DProduct<complex<T>, m, n, p, E1, E2>(u, v);
  }
  

  //! returns u/v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDivision<T, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DDivision<T, m, n, p, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDivision<complex<T>, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v)
  {
    return TinyArray3DDivision<complex<T>, m, n, p, E1, E2>(u, v);
  }


  //! returns u/v element-wise
  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDivision<complex<T>, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v)
  {
    return TinyArray3DDivision<complex<T>, m, n, p, E1, E2>(u, v);
  }
  
  
  //! returns alpha*u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<T, m, n, p, T, E>
  operator *(const T& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScaled<T, m, n, p, T, E>(alpha, u);
  }


  //! returns alpha*u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }

  
  //! returns alpha*u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, T, E>
  operator *(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u)
  {
    return TinyArray3DScaled<complex<T>, m, n, p, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<T, m, n, p, T, E> operator *(const TinyArray3DExpression<T, m, n, p, E>& u,
						       const T& alpha)
  {
    return TinyArray3DScaled<T, m, n, p, T, E>(alpha, u);
  }


  //! returns u*alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator *(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha)
  {
    return TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }

  
  //! returns u*alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, T, E>
  operator *(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha)
  {
    return TinyArray3DScaled<complex<T>, m, n, p, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalDiv<T, m, n, p, T, E> operator /(const T& alpha,
							const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalDiv<T, m, n, p, T, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalDiv<complex<T>, m, n, p, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalDiv<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }


  //! returns alpha / u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalDiv<complex<T>, m, n, p, T, E>
  operator /(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u)
  {
    return TinyArray3DScalDiv<complex<T>, m, n, p, T, E>(alpha, u);
  }

  
  //! returns u / alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<T, m, n, p, T, E> operator /(const TinyArray3DExpression<T, m, n, p, E>& u,
						       const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyArray3DScaled<T, m, n, p, T, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator /(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha)
  {
    T one; SetComplexOne(one);
    return TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>(one/alpha, u);
  }


  //! returns u / alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScaled<complex<T>, m, n, p, T, E> 
  operator /(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha)
  {
    T one; SetComplexOne(one);
    return TinyArray3DScaled<complex<T>, m, n, p, T, E>(one/alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<T, m, n, p, T, E> operator +(const T& alpha,
							const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalSum<T, m, n, p, T, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }


  //! returns alpha + u 
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator +(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, T, E>(alpha, u);
  }

  
  //! returns u + alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<T, m, n, p, T, E> operator +(const TinyArray3DExpression<T, m, n, p, E>& u,
							const T& alpha)
  {
    return TinyArray3DScalSum<T, m, n, p, T, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator +(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }


  //! returns u + alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator +(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, T, E>(alpha, u);
  }

    
  //! returns u - alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<T, m, n, p, T, E> operator -(const TinyArray3DExpression<T, m, n, p, E>& u,
							const T& alpha)
  {
    return TinyArray3DScalSum<T, m, n, p, T, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator -(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>(-alpha, u);
  }


  //! returns u - alpha
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator -(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha)
  {
    return TinyArray3DScalSum<complex<T>, m, n, p, T, E>(-alpha, u);
  }

    
  //! returns alpha - u
  template<class T, int m, int n, int p, class E> 
  inline const TinyArray3DScalDiff<T, m, n, p, T, E> operator -(const T& alpha,
						      const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalDiff<T, m, n, p, T, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, int n, int p, class E> 
  inline const TinyArray3DScalDiff<complex<T>, m, n, p, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DScalDiff<complex<T>, m, n, p, complex<T>, E>(alpha, u);
  }


  //! returns alpha - u
  template<class T, int m, int n, int p, class E> 
  inline const TinyArray3DScalDiff<complex<T>, m, n, p, T, E>
  operator -(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u)
  {
    return TinyArray3DScalDiff<complex<T>, m, n, p, T, E>(alpha, u);
  }
  
  
  //! returns -u
  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DOpposite<T, m, n, p, E> operator-(const TinyArray3DExpression<T, m, n, p, E>& u)
  {
    return TinyArray3DOpposite<T, m, n, p, E>(u);
  }
  
}

#define SELDON_FILE_TINY_ARRAY3D_EXPRESSION_INLINE_CXX
#endif
