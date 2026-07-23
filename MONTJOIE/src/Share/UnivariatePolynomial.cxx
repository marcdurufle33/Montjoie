#ifndef MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_CXX

#include "UnivariatePolynomialInline.cxx"

namespace Montjoie
{
  
  //! modification of the order
  /*!
    \warning previous data is removed
  */
  template<class T, class Allocator>
  void UnivariatePolynomial<T, Allocator>::SetOrder(int order)
  { 
    if (order != n_)
      {
	n_ = order;
	// we allocate the array coef and initialize to zero
	this->Reallocate(n_+1); this->Fill(T(0));
      }
  }
  
  
  //! modification of the order while keeping old coefficients
  template<class T, class Allocator>
  void UnivariatePolynomial<T, Allocator>::ResizeOrder(int order)
  {
    if (order != n_)
      {
	n_ = order;
	size_t old_m = this->m_;
	// we resize the array coef
	Vector<T, Vect_Full, Allocator>::Resize(n_+1);
	// zeros on additional coefficients
	for (size_t i = old_m; i < this->m_; i++)
	  this->data_[i] = T(0);
      }
  }
  
  
  //! displays the polynom 
  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const UnivariatePolynomial<T, Allocator>& P)
  {
    int n = P.GetOrder();
    if (n == 0)
      {
	out << P(0) << endl;
	return out;
      }
    
    if (P(0) == T(0))
      {
	if ( P(1) == T(1) )
	  out << "x ";
	else if ( P(1) != T(0) )
	  out << P(1) << " x ";
      }
    else
      {
	if ( P(1) == T(1) )
	  out << P(0) << " + x ";
	else if ( P(1) == T(0) )
	  out << P(0) << " ";
	else
	  out << P(0) << " + " << P(1) << " x ";
      }
    
    for (int i = 2; i <= n; i++)
      {
	if (P(i) == T(1))
	  out << "+ x^" << i << " ";
	else if (P(i) != T(0))
	  out << "+ " << P(i) << " x^" << i << " ";
      }

    return out;
  }
  
  
  //! R = beta*R + alpha*P*Q
  template<class T0, class T1, class T, class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const UnivariatePolynomial<T, Allocator1>& P,
	      const UnivariatePolynomial<T, Allocator2>& Q,
	      const T1& beta, UnivariatePolynomial<T, Allocator3>& R)
  {
    int m = P.GetOrder(), n = Q.GetOrder();
    int new_order = max(R.GetOrder(), m+n);
    Mlt(beta, R);
    R.ResizeOrder(new_order);
    T sum;
    for (int i = 0; i <= (m+n); i++)
      {
	sum = T(0);
	for (int k = max(0, i-n); k <= min(m,i); k++)
          sum += P(k)*Q(i-k);
          	
	R(i) += T(alpha)*sum;
      }
  }
  

  //! power of a polynom
  template<class T, class Allocator>
  UnivariatePolynomial<T, Allocator>
  Pow(const UnivariatePolynomial<T, Allocator>& U, int r)
  {
    UnivariatePolynomial<T, Allocator> P = U;
    UnivariatePolynomial<T, Allocator> Q;
    if (r == 0)
      {
        Q.SetOrder(0);
	Q(0) = 1.0;
      }
    else
      {
        Q = U;
        for (int k = 0; k < r-1; k++)
          {  
	    MltAdd(1, U, P, 0, Q);
	    P = Q;
          }
      }
    
    return Q;
  }

  
  //! computation of Legendre polynoms
  /*!
    \warning this explicit computation of coefficients for Legendre polynoms
    is highly instable, don't use for order greater than 50
   */
  template<class T, class Allocator1, class Allocator2>
  void ComputeLegendrePolynome(int begin_order, int n,
			       Vector<UnivariatePolynomial<T, Allocator1>,
			       Vect_Full, Allocator2>& Pn )
  {
    if (n <= 0)
      return ;
    
    if (n >= 50)
      {
	cout<<"unstable computation "<<endl;
	abort();
      }
    
    Pn.Reallocate(n);
    if (n == 1)
      {
	Pn(0).SetOrder(0);
	Pn(0)(0) = 1.0;
      }
    else
      {
	// the two first polynoms are 1 and x
	Pn(0).SetOrder(0);
	Pn(0)(0) = 1.0;
	Pn(1).SetOrder(1);
	Pn(1)(1) = 1.0;
	// for the others, the following recurrence relation is used
	// L_{k} = (2*k-1)/k x L_{k-1} - (k-1)/k L_{k-2}
	for (int k = 2; k < n; k++)
	  {
	    Pn(k).SetOrder(k);
	    for (int i = 0; i <= Pn(k-1).GetOrder(); i++)
	      Pn(k)(i+1) = Pn(k-1)(i)*(2*k-1)/k;
	    
	    for (int i = 0; i <= Pn(k-2).GetOrder(); i++)
	      Pn(k)(i) -= Pn(k-2)(i)*(k-1)/k;
	  }
      }
  }

  
  //! default constructor
  template<class T>
  LegendrePolynomial<T>::LegendrePolynomial(int n)
  {
    // coefficients of the recurrence are computed up to order n
    coef_a.Reallocate(n); coef_b.Reallocate(n);
    coef_a.Fill(0); coef_b.Fill(0);
    for (int k = 2; k < n; k++)
      {
        coef_a(k) = T(2*k-1)/k;
        coef_b(k) = -T(k-1)/k;
      }
  }

  
  //! computes P_i(cos_teta) for i =0 .. n-1
  template<class T>
  void LegendrePolynomial<T>::EvaluatePn(int n, const T& cos_teta, Vector<T>& Pn0) const
  {
    if (n <= 0)
      {
        Pn0.Clear();
        return;
      }

    Pn0.Reallocate(n);
    Pn0(0) = 1.0;
    if (n >= 2)
      Pn0(1) = cos_teta;
    
    for (int i = 2; i < n; i++)
      Pn0(i) = coef_a(i)*cos_teta*Pn0(i-1) + coef_b(i)*Pn0(i-2);    
  }
  

  //! computes P_i(x) and P_i'(x) for i =0 .. n-1
  template<class T>
  void LegendrePolynomial<T>::EvaluateDerive(int n, const T& x,
					     Vector<T>& Pn, Vector<T>& dPn) const
  {
    if (n <= 0)
      {
        Pn.Clear();
        dPn.Clear();
        return;
      }
    
    Pn.Reallocate(n);
    dPn.Reallocate(n);
    
    Pn(0) = 1.0;
    dPn(0) = 0.0;
    
    if (n >= 2)
      {
	Pn(1) = x;
	dPn(1) = 1.0;
      }
    
    for (int i = 2; i < n; i++)
      {
        dPn(i) = coef_a(i)*(Pn(i-1) + x*dPn(i-1)) + coef_b(i)*dPn(i-2);
        Pn(i) = coef_a(i)*x*Pn(i-1) + coef_b(i)*Pn(i-2);
      }
  }
  
  
  //! computes Pn(cos_teta), Pn_1 = sin_teta Pn'(cos_teta) and derivative of Pn_1
  template<class T>
  void LegendrePolynomial<T>::EvaluateDerivative(int n, const T& cos_teta,
                                                 const T& sin_teta, Vector<T>& Pn0,
                                                 Vector<T>& Pn1, Vector<T>& dPn1_dteta) const
  {
    if (n <= 0)
      {
        Pn0.Clear();
        Pn1.Clear();
        dPn1_dteta.Clear();
        return;
      }
    
    Pn0.Reallocate(n);
    Pn1.Reallocate(n);
    dPn1_dteta.Reallocate(n);
    Pn0(0) = 1.0;
    Pn1(0) = 0.0;
    dPn1_dteta(0) = 0.0;
    
    if (n >= 2)
      {
        Pn0(1) = cos_teta;
        Pn1(1) = sin_teta;
        dPn1_dteta(1) = cos_teta;        
      }
    
    T derPn, der2Pn, derPn_prev, der2Pn_prev, derPn_next, der2Pn_next;
    derPn_prev = 0.0; der2Pn_prev = 0.0;
    derPn = 1.0; der2Pn = 0.0;
    T sin2_teta = sin_teta*sin_teta;
    for (int i = 2; i < n; i++)
      {
        derPn_next = coef_a(i)*(Pn0(i-1) + cos_teta*derPn) + coef_b(i)*derPn_prev;
        der2Pn_next = coef_a(i)*(2.0*derPn + cos_teta*der2Pn) + coef_b(i)*der2Pn_prev;
        Pn0(i) = coef_a(i)*cos_teta*Pn0(i-1) + coef_b(i)*Pn0(i-2);
        Pn1(i) = sin_teta*derPn_next;
        dPn1_dteta(i) = cos_teta*derPn_next - sin2_teta*der2Pn_next;
        derPn_prev = derPn; der2Pn_prev = der2Pn;
        derPn = derPn_next; der2Pn = der2Pn_next;
      }
  }
  
  
  //! default constructor
  template<class T>
  AssociatedLegendrePolynomial<T>::AssociatedLegendrePolynomial()
  {
  }
  
  
  //! computes coefficients of the recurrence
  template<class T>
  void AssociatedLegendrePolynomial<T>::Init(int lmax)
  {
    if (lmax < 0)
      {
        coef_a.Clear(); coef_b.Clear();
        coef_c.Clear(); coef_d.Clear();
        return;
      }
    
    if (lmax > 1000)
      {
        cout << "The current implementation of spherical harmonics suffer from underflow for l > 1000" << endl;
        abort();
      }
    
    T one, u, v, w; SetComplexOne(one);

    coef_a.Reallocate(lmax+1);
    coef_b.Reallocate(lmax+1);
    for (int l = 0; l <= lmax; l++)
      {
        coef_a(l).Reallocate(l+1);
        coef_b(l).Reallocate(l+1);
      }
    
    for (int l = 0; l <= lmax; l++)
      for (int m = 0; m < l-1; m++)
        {
          u = T(2*l+1)/T(l-m);
          v = T(2*l-1)/T(l+m);
          coef_a(l)(m) = sqrt(u*v);
          
          v = T(l-m-1)/T(l+m);
          w = T(l+m-1)/T(2*l-3);
          coef_b(l)(m) = sqrt(u*v*w);
        }

    // P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
    // => \bar{P}_m^m = c_m (sin teta)^m with c_m = -c_{m-1} \sqrt{\frac{2m+1}{2m}}
    // P_{m+1}^m = (2m+1) x P_m^m
    // => \bar{P}_{m+1}^m = d_m cos teta \bar{P}_m^m with d_m = \sqrt{2m+3}
    coef_c.Reallocate(lmax+1);
    coef_d.Reallocate(lmax+1);
    Real_wp poids = sqrt(one/(4.0*pi_wp));
    coef_c(0) = poids; coef_d(0) = sqrt(T(3));
    for (int m = 1; m <= lmax; m++)
      {
        coef_c(m) = -coef_c(m-1)*sqrt(T(2*m+1)/T(2*m));
        coef_d(m) = sqrt(T(2*m+3));
      }
  }
  
  
  //! evaluates normalized associated Legendre polynomials up to order lmax
  template<class T>
  void AssociatedLegendrePolynomial<T>::EvaluatePnm(int lmax, int mmax,
                                                    const T& teta, Vector<Vector<T> >& P)
  {
    if ((lmax < 0) || (mmax > lmax))
      {
        P.Clear();
        return;
      }
    
    if (P.GetM() <= lmax)
      P.Reallocate(lmax+1);
    
    for (int l = 0; l <= lmax; l++)
      {
        if (P(l).GetM() <= l)
          P(l).Reallocate(min(mmax+1, l+1));
      }
    
    T cos_teta = cos(teta), sin_teta = sin(teta);
    T sin_teta_m = 1.0;
    for (int m = 0; m <= mmax; m++)
      {
        P(m)(m) = coef_c(m)*sin_teta_m;
        if (m < lmax)
          P(m+1)(m) = coef_d(m)*cos_teta*P(m)(m);
        
        for (int l = m+2; l <= lmax; l++)
          P(l)(m) = cos_teta*coef_a(l)(m)*P(l-1)(m) - coef_b(l)(m)*P(l-2)(m);
        
        sin_teta_m *= sin_teta;
      }
  }
  

  //! solves a real polynomial equation (zero-valued solutions are dropped)
  template<class T>
  void SolvePolynomialEquation(const UnivariatePolynomial<T>& P, Vector<T>& R, Vector<T>& Rimag)
  {
    T one(1), zero(0);
    
    // we detect the first non-null coefficient
    int first_coef = P.GetOrder();
    for (int i = 0; i < P.GetM(); i++)
      if (P(i) != zero)
	{
	  first_coef = i;
	  break;
	}
    
    // and the last non-null coefficient
    int last_coef = 0;
    for (int i = P.GetM()-1; i >= 0; i--)
      if (P(i) != zero)
	{
	  last_coef = i;
	  break;
	}
    
    // degree of the actual polynomial to solve
    int N = last_coef - first_coef;
    if (N <= 0)
      {
	R.Clear();
	Rimag.Clear();
	return;
      }
    
    // finding eigenvalues of companion matrix
    Matrix<T> A, V;
    A.Reallocate(N, N);
    A.Fill(zero);
    for (int i = 0; i < N-1; i++)
      A(i+1, i) = one;
    
    T coef = one/P(last_coef);
    for (int i = first_coef; i < last_coef; i++)
      A(i-first_coef, N-1) = -coef*P(i);
    
    GetEigenvaluesEigenvectors(A, R, Rimag, V);
  }


  //! solves a complex polynomial equation (zero-valued solutions are dropped)
  template<class T>
  void SolvePolynomialEquation(const UnivariatePolynomial<complex<T> >& P, Vector<complex<T> >& R)
  {
    complex<T> one(1, 0), zero(0, 0);
    
    // we detect the first non-null coefficient
    int first_coef = P.GetOrder();
    for (int i = 0; i < P.GetM(); i++)
      if (abs(P(i)) != zero)
	{
	  first_coef = i;
	  break;
	}
    
    // and the last non-null coefficient
    int last_coef = 0;
    for (int i = P.GetM()-1; i >= 0; i--)
      if (abs(P(i)) != zero)
	{
	  last_coef = i;
	  break;
	}
    
    // degree of the actual polynomial to solve
    int N = last_coef - first_coef;
    if (N <= 0)
      {
	R.Clear();
	return;
      }
    
    // finding eigenvalues of companion matrix
    Matrix<complex<T> > A;
    A.Reallocate(N, N);
    A.Fill(zero);
    for (int i = 0; i < N-1; i++)
      A(i+1, i) = one;
    
    complex<T> coef = one/P(last_coef);
    for (int i = first_coef; i < last_coef; i++)
      A(i-first_coef, N-1) = -coef*P(i);
    
    GetEigenvalues(A, R);
  }
  
}

#define MONTJOIE_FILE_UNIVARIATE_POLYNOMIAL_CXX
#endif
