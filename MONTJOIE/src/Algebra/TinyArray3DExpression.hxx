#ifndef SELDON_FILE_TINY_ARRAY3D_EXPRESSION_HXX

namespace Seldon
{
  
  //! Expression between vectors
  template<class T, int m, int n, int p, class E>
  class TinyArray3DExpression
  {
  public:
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
    inline operator E&() { return static_cast<E&>(*this); }
    inline operator E const&() const { return static_cast<const E&>(*this); }
    
  };


  //! Difference between two expressions
  template<class T, int m, int n, int p, class E1, class E2>
  class TinyArray3DDifference : public TinyArray3DExpression<T, m, n, p, TinyArray3DDifference<T, m, n, p, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyArray3DDifference(const TinyArray3DExpression<T1, m, n, p, E1>&,
			  const TinyArray3DExpression<T2, m, n, p, E2>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Sum between two expressions
  template<class T, int m, int n, int p, class E1, class E2>
  class TinyArray3DSum : public TinyArray3DExpression<T, m, n, p, TinyArray3DSum<T, m, n, p, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyArray3DSum(const TinyArray3DExpression<T1, m, n, p, E1>&,
		   const TinyArray3DExpression<T2, m, n, p, E2>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Product between two expressions (elementwise)
  template<class T, int m, int n, int p, class E1, class E2>
  class TinyArray3DProduct : public TinyArray3DExpression<T, m, n, p, TinyArray3DProduct<T, m, n, p, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyArray3DProduct(const TinyArray3DExpression<T1, m, n, p, E1>&,
		       const TinyArray3DExpression<T2, m, n, p, E2>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Division between two expressions (elementwise)
  template<class T, int m, int n, int p, class E1, class E2>
  class TinyArray3DDivision : public TinyArray3DExpression<T, m, n, p, TinyArray3DDivision<T, m, n, p, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    template<class T1, class T2>
    TinyArray3DDivision(const TinyArray3DExpression<T1, m, n, p, E1>&,
			const TinyArray3DExpression<T2, m, n, p, E2>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Multiplication between a scalar and a 3-D array
  template<class T, int m, int n, int p, class T0, class E>
  class TinyArray3DScaled : public TinyArray3DExpression<T, m, n, p, TinyArray3DScaled<T, m, n, p, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyArray3DScaled(const T0&, const TinyArray3DExpression<T1, m, n, p, E>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Division between a scalar and a 3-D array
  template<class T, int m, int n, int p, class T0, class E>
  class TinyArray3DScalDiv : public TinyArray3DExpression<T, m, n, p, TinyArray3DScalDiv<T, m, n, p, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyArray3DScalDiv(const T0&, const TinyArray3DExpression<T1, m, n, p, E>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Addition between a scalar and a 3-D array
  template<class T, int m, int n, int p, class T0, class E>
  class TinyArray3DScalSum : public TinyArray3DExpression<T, m, n, p, TinyArray3DScalSum<T, m, n, p, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyArray3DScalSum(const T0&, const TinyArray3DExpression<T1, m, n, p, E>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Subtraction between a scalar and a 3-D array
  template<class T, int m, int n, int p, class T0, class E>
  class TinyArray3DScalDiff : public TinyArray3DExpression<T, m, n, p, TinyArray3DScalDiff<T, m, n, p, T0, E> >
  {
    const T0 alpha_;
    const E& u_;
    
  public:
    template<class T1>
    TinyArray3DScalDiff(const T0&, const TinyArray3DExpression<T1, m, n, p, E>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };


  //! Opposite of a 3-D array
  template<class T, int m, int n, int p, class E>
  class TinyArray3DOpposite : public TinyArray3DExpression<T, m, n, p, TinyArray3DOpposite<T, m, n, p, E> >
  {
    const E& u_;
    
  public:
    TinyArray3DOpposite(const TinyArray3DExpression<T, m, n, p, E>&);
    
    const T operator()(int, int, int) const;
    const T Val(int) const;
    
  };
  
  
  /*************
   * Operators *
   *************/  
  
  
  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DSum<T, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DSum<complex<T>, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DSum<complex<T>, m, n, p, E1, E2> 
  operator +(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDifference<T, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  inline const TinyArray3DDifference<complex<T>, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DDifference<complex<T>, m, n, p, E1, E2>
  operator -(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DProduct<T, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DProduct<complex<T>, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DProduct<complex<T>, m, n, p, E1, E2> 
  operator *(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DDivision<T, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DDivision<complex<T>, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<complex<T>, m, n, p, E1>& u,
	     const TinyArray3DExpression<T, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E1, class E2>
  const TinyArray3DDivision<complex<T>, m, n, p, E1, E2>
  operator /(const TinyArray3DExpression<T, m, n, p, E1>& u,
	     const TinyArray3DExpression<complex<T>, m, n, p, E2>& v);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<T, m, n, p, T, E>
  operator *(const T& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator *(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, T, E>
  operator *(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<T, m, n, p, T, E>
  operator *(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator *(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, T, E>
  operator *(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalDiv<T, m, n, p, T, E>
  operator /(const T& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalDiv<complex<T>, m, n, p, complex<T>, E>
  operator /(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalDiv<complex<T>, m, n, p, T, E>
  operator /(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u);
  
  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<T, m, n, p, T, E>
  operator /(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, complex<T>, E>
  operator /(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScaled<complex<T>, m, n, p, T, E> 
  operator /(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<T, m, n, p, T, E>
  operator +(const T& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator +(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator +(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<T, m, n, p, T, E>
  operator +(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator +(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator +(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<T, m, n, p, T, E>
  operator -(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, complex<T>, E>
  operator -(const TinyArray3DExpression<T, m, n, p, E>& u,
	     const complex<T>& alpha);

  template<class T, int m, int n, int p, class E>
  const TinyArray3DScalSum<complex<T>, m, n, p, T, E>
  operator -(const TinyArray3DExpression<complex<T>, m, n, p, E>& u,
	     const T& alpha);

  template<class T, int m, int n, int p, class E> 
  const TinyArray3DScalDiff<T, m, n, p, T, E>
  operator -(const T& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E> 
  const TinyArray3DScalDiff<complex<T>, m, n, p, complex<T>, E>
  operator -(const complex<T>& alpha,
	     const TinyArray3DExpression<T, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E> 
  const TinyArray3DScalDiff<complex<T>, m, n, p, T, E>
  operator -(const T& alpha,
	     const TinyArray3DExpression<complex<T>, m, n, p, E>& u);

  template<class T, int m, int n, int p, class E>
  inline const TinyArray3DOpposite<T, m, n, p, E>
  operator-(const TinyArray3DExpression<T, m, n, p, E>& u);
  
}

#define SELDON_FILE_TINY_ARRAY3D_EXPRESSION_HXX
#endif
