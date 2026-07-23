#ifndef MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor, null polynom
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>::UnivariatePolynomial()
    : Vector<T, Vect_Full, Allocator>(1)
  {
    n_ = 0; 
    this->Fill(T(0));
  }
  
  
  //! polynome of order n
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>::UnivariatePolynomial(int n)
    : Vector<T, Vect_Full, Allocator>(n+1)
  {
    n_ = n;
    this->Fill(T(0));
  }
  
  
  //! default destructor
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>::~UnivariatePolynomial()
  {
  }
    
  
  /************************
   * Convenient Functions *
   ************************/
  
  
  //! returns order of polynome
  template<class T, class Allocator>
  inline int UnivariatePolynomial<T, Allocator>::GetOrder() const
  {
    return n_;
  }
  
  
  //! Evaluation of the polynome for a given x
  /*! we compute \f$ P(x) = a_0 + a_1 x + a_2 x^2 + a_3 x^3 + ... + a_n x^n \f$
    with the Horner algorithm
    \f$ P(x) = a_0 + x*( a_1 + x*( a_2 + x*(a_3 ...+x*a_n) ) ) \f$
  */
  template<class T, class Allocator> template<class T0>
  inline T0 UnivariatePolynomial<T, Allocator>::Evaluate(const T0& x) const
  {
    T0 res = this->data_[n_];
    int i = n_;
    while (i > 0)
      res = res*x + this->data_[--i];
    
    return res;
  }
  
  
  /*************
   * Operators *
   *************/
  
  
  //! addition of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator +=(const UnivariatePolynomial<T, Allocator>& P)
  {
    Add(T(1), P, *this);
    return *this;
  }
  
  
  //! difference of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator -=(const UnivariatePolynomial<T, Allocator>& P)
  {
    Add(T(-1), P, *this);
    return *this;
  }
  
  
  //! multiplication of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator *=(const UnivariatePolynomial<T, Allocator>& P)
  {
    UnivariatePolynomial<T, Allocator> Q(*this);
    MltAdd(T(1), Q, P, T(0), *this);
    return *this;
  }
  
    
  //! multiplication of a polynom by a scalar
  template<class T, class Allocator>
  template<class T0> inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator *=(const T0& alpha)
  {
    Mlt(alpha, *this);
    return *this;
  }
  

  //! multiplication of a polynom by a scalar
  template<class T, class Allocator>
  template<class T0> inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator /=(const T0& alpha)
  {
    Mlt(1.0/alpha, *this);
    return *this;
  }

  
  //! we add scalar to a polynom
  template<class T, class Allocator>
  template<class T0> inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator +=(const T0& alpha)
  {
    (*this)(0) += alpha;
    return *this;
  }
  
  
  //! we subtract scalar to a polynom
  template<class T, class Allocator>
  template<class T0> inline UnivariatePolynomial<T, Allocator>&
  UnivariatePolynomial<T, Allocator>::operator -=(const T0& alpha)
  {
    (*this)(0) -= alpha;
    return *this;
  }
   

  //! addition of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>
  operator +(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q)
  {
    UnivariatePolynomial<T, Allocator> R(P);
    Add(T(1), Q, R);
    return R;
  }
  
  
  //! difference of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>
  operator -(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q)
  {
    UnivariatePolynomial<T, Allocator> R(P);
    Add(T(-1), Q, R);
    return R;
  }
  

  //! opposite of a polynomial
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>
  operator -(const UnivariatePolynomial<T, Allocator>& P)
  {
    UnivariatePolynomial<T, Allocator> R(P);
    R *= -1.0;
    return R;
  }

  
  //! multiplication of two polynoms
  template<class T, class Allocator>
  inline UnivariatePolynomial<T, Allocator>
  operator *(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q)
  {
    UnivariatePolynomial<T, Allocator> R;
    MltAdd(T(1), P, Q, T(0), R);
    return R;
  }
  
    
    //! multiplication of a polynom by a scalar
  template<class T, class Allocator, class T0>
  inline UnivariatePolynomial<T, Allocator>
  operator *(const UnivariatePolynomial<T, Allocator>& P, const T0& alpha)
  {
    UnivariatePolynomial<T, Allocator> Q(P);
    Mlt(alpha, Q);
    return Q;
  }
  

  //! division of a polynom by a scalar
  template<class T, class Allocator, class T0>
  inline UnivariatePolynomial<T, Allocator>
  operator /(const UnivariatePolynomial<T, Allocator>& P, const T0& alpha)
  {
    UnivariatePolynomial<T, Allocator> Q(P);
    Mlt(T(1)/alpha, Q);
    return Q;
  }

  
  //! we add scalar to a polynom
  template<class T, class Allocator, class T0>
  inline UnivariatePolynomial<T, Allocator>
  operator +(const UnivariatePolynomial<T, Allocator>& P, const T0& alpha)
  {
    UnivariatePolynomial<T, Allocator> Q(P);
    Q(0) += alpha;
    return Q;
  }
  
  
  //! we subtract scalar to a polynom
  template<class T, class Allocator, class T0>
  inline UnivariatePolynomial<T, Allocator>
  operator -(const UnivariatePolynomial<T, Allocator>& P, const T0& alpha)
  {
    UnivariatePolynomial<T,Allocator> Q(P);
    Q(0) -= alpha;
    return Q;
  }

  
  //! Multiplication of a polynom by a scalar
  template<class T0, class T, class Allocator> inline UnivariatePolynomial<T, Allocator>
  operator *(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P)
  {
    UnivariatePolynomial<T,Allocator> Q = P;
    Mlt(alpha, Q); return Q;
  }
  

  //! Addition of a polynom and a scalar
  template<class T0, class T, class Allocator> inline UnivariatePolynomial<T, Allocator>
  operator +(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P)
  {
    UnivariatePolynomial<T,Allocator> Q = P;
    Q(0) += alpha; return Q;
  }
  
  
  //! Subtraction of a polynom by a scalar
  template<class T0, class T, class Allocator> inline UnivariatePolynomial<T, Allocator>
  operator -(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P)
  {
    UnivariatePolynomial<T, Allocator> Q = P;
    Q(0) -= alpha; Mlt(T(-1), Q); return Q;
  }
  
  
  //! derivation of 1-D polynom
  template<class T, class Allocator1, class Allocator2>
  inline void DerivatePolynomial(const UnivariatePolynomial<T, Allocator1>& P,
				UnivariatePolynomial<T, Allocator2>& Q)
  {
    int n = P.GetOrder(); 
    if (n == 0)
      Q.SetOrder(0);
    else
      Q.SetOrder(n-1);

    Q(0) = 0.0;
    for (int i = 1; i <= n; i++)
      Q(i-1) = P(i)*i;
      
  }
  
  
  //! integration of 1-D polynom
  template<class T, class Allocator1, class Allocator2>
  inline void IntegratePolynomial(const UnivariatePolynomial<T, Allocator1>& P,
				 UnivariatePolynomial<T, Allocator2>& Q)
  {
    int n = P.GetOrder(); Q.SetOrder(n+1);
    for (int i = 0; i <= n; i++)
      Q(i+1) = P(i)/(i+1);
    
  }
  
  
  //! R = beta*R
  template<class T0, class T, class Allocator>
  inline void Mlt(const T0& alpha, UnivariatePolynomial<T, Allocator>& R)
  {
    if (alpha == T0(0))
      R.Fill(T0(0));
    else
      for (int i = 0; i <= R.GetOrder(); i++)
	R(i) *= T(alpha);
  }
  

  template<class T0, class T, class Allocator1, class Allocator2> inline
  void Add(const T0& alpha, const UnivariatePolynomial<T, Allocator1>& P,
           UnivariatePolynomial<T, Allocator2>& Q)
  {
    int m = P.GetOrder();
    int n = Q.GetOrder();
    if (m > n)
      Q.ResizeOrder(m);
    
    for (int i = 0; i <= m; i++)
      Q(i) += alpha*P(i);
  }
  
}

#define MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_INLINE_CXX
#endif
