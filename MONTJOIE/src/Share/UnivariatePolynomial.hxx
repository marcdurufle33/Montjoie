#ifndef MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_HXX


namespace Montjoie
{
  //! polynom class of one variable x
  template<class T, class Allocator =
           typename SeldonDefaultAllocator<VectFull, T>::allocator >
  class UnivariatePolynomial : public Vector<T, Vect_Full, Allocator>
  {
    int n_; //!< order of the polynom
    // vector which contains the coefficients a_0, a_1, a_2 ...
    // P = a_0 + a_1 x + a_2 x^2 + a_3 x^3 + ... + a_n x^n 
        
    // basic methods
  public:
    
    // constructors
    UnivariatePolynomial();
    UnivariatePolynomial(int n);

    ~UnivariatePolynomial();
    
    // basic functions
    int GetOrder() const;
    void SetOrder(int order);
    void ResizeOrder(int order);
    
    // convenient functions
    template<class T0>
    T0 Evaluate(const T0& x) const;
    
    // operators
    UnivariatePolynomial<T, Allocator>&
    operator +=(const UnivariatePolynomial<T, Allocator>& P);
    
    UnivariatePolynomial<T, Allocator>&
    operator -=(const UnivariatePolynomial<T, Allocator>& P);
    
    UnivariatePolynomial<T, Allocator>&
    operator *=(const UnivariatePolynomial<T, Allocator>& P);
    
    template<class T0>
    UnivariatePolynomial<T, Allocator>& operator *=(const T0& alpha);
    
    template<class T0>
    UnivariatePolynomial<T, Allocator>& operator /=(const T0& alpha);
  
    template<class T0>
    UnivariatePolynomial<T, Allocator>& operator +=(const T0& alpha);
    
    template<class T0>
    UnivariatePolynomial<T, Allocator>& operator -=(const T0& alpha);

  };
  
  // other operators for polynoms
  template<class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator *(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q);
  
  template<class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator +(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q);
  
  template<class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator -(const UnivariatePolynomial<T, Allocator>& P,
	     const UnivariatePolynomial<T, Allocator>& Q);

  template<class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator -(const UnivariatePolynomial<T, Allocator>& P);
  
  template<class T0, class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator *(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P);
  
  template<class T0, class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator +(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P);
  
  template<class T0, class T, class Allocator> UnivariatePolynomial<T, Allocator>
  operator -(const T0& alpha, const UnivariatePolynomial<T, Allocator>& P);
    
  // other functions
  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const UnivariatePolynomial<T, Allocator>& P);

  template<class T, class Allocator>
  UnivariatePolynomial<T, Allocator>
  Pow(const UnivariatePolynomial<T, Allocator>& U, int r);
  
  template<class T, class Allocator1, class Allocator2>
  void DerivatePolynomial(const UnivariatePolynomial<T, Allocator1>& P,
			 UnivariatePolynomial<T, Allocator2>& Q);

  template<class T, class Allocator1, class Allocator2>
  void IntegratePolynomial(const UnivariatePolynomial<T, Allocator1>& P,
			  UnivariatePolynomial<T, Allocator2>& Q);
  
  //! Legendre polynomial
  template<class T>
  class LegendrePolynomial
  {
  protected :
    //! coefficients of the recurrence relation
    Vector<T> coef_a, coef_b;
  
  public :
    LegendrePolynomial(int n);
    
    void EvaluateDerivative(int n, const T& cos_teta, const T& sin_teta,
                            Vector<T>& Pn0, Vector<T>& Pn1, Vector<T>& dPn1_dteta) const;

    void EvaluatePn(int n, const T& cos_teta, Vector<T>& Pn0) const;
    
    void EvaluateDerive(int n, const T& x, Vector<T>& Pn, Vector<T>& dPn) const;
    
  };  

  
  //! associated Legendre polynomials (normalized)
  /*!
    P_l^m satisfies the relation of recurrence :
    (l-m) P_l^m = x (2l -1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    the class computes \bar{P}_l^m = \sqrt{\frac{(2l+1) (l-m)!}{4 pi (l+m)!}} P_l^m(cos \theta)
   */
  template<class T>
  class AssociatedLegendrePolynomial
  {
  protected :
    //! coefficients of the recurrence relation
    Vector<Vector<T> > coef_a, coef_b;
    //! coefficients for initial values
    Vector<T> coef_c, coef_d;
    
  public :
    AssociatedLegendrePolynomial();
    
    void Init(int lmax);
    
    void EvaluatePnm(int lmax, int mmax, const T& teta, Vector<Vector<T> >& Pnm);
    
  };

  template<class T, class Allocator1, class Allocator2>
  void ComputeLegendrePolynome(int begin_order, int n,
			       Vector<UnivariatePolynomial<T, Allocator1>,
			       Vect_Full, Allocator2>& Pn );

  template<class T0, class T, class Allocator1, class Allocator2>
  void Add(const T0& alpha, const UnivariatePolynomial<T, Allocator1>& P,
           UnivariatePolynomial<T, Allocator2>& Q);
  
  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, UnivariatePolynomial<T, Allocator>& R);

  template<class T0, class T1, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const UnivariatePolynomial<T, Allocator1>& P,
	      const UnivariatePolynomial<T, Allocator2>& Q,
	      const T1& beta, UnivariatePolynomial<T, Allocator3>& R);

  template<class T>
  void SolvePolynomialEquation(const UnivariatePolynomial<T>& P, Vector<T>& R, Vector<T>& Rimag);

  template<class T>
  void SolvePolynomialEquation(const UnivariatePolynomial<complex<T> >& P, Vector<complex<T> >& R);

}

#define MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_HXX
#endif
