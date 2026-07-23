#ifndef MONTJOIE_FILE_GAUSS_JACOBI_CXX

namespace Montjoie
{
  
  //! computes relation recurrence for orthogonal polynomials 
  /*!
    \param[out] ab coefficients defining recurrence relation satisfied by polynomials
    \param[in] n number of polynomials to compute
    \param[in] alpha exponent in weight (1-x)^alpha
    \param[in] beta exponent in weight (1+x)^beta
    we compute polynomials Pn orthogonal relatively to Jacobi weight (1-x)^alpha (1+x)^beta
    \int_{-1}^1  P_n(x) P_m(x) (1-x)^alpha (1+x)^beta dx = C_m \delta_{n, m}
    We have P_0 = 1, P_1 = x - ab(0, 0), and then :
    P_k(x) = (x - ab(k-1, 0)) P_{k-1}(x) - ab(k-1, 1) P_{k-2}(x)
    ab(0, 1) is the value of \int_0^1 \, (1-x)^alpha \, x^beta dx
    */
  template<class T>
  void GetJacobiPolynomial(Matrix<T>& ab, int n, const T& alpha, const T& beta)
  {  
    Vector<T> nab;

    if (n <= 0) return;

    if (n > 2)
      nab.Reallocate(n-2);
    
    ab.Reallocate(n, 2);
  
    ab(0,0) = (beta-alpha)/(alpha + beta + 2);
    ab(0,1) = tgamma(alpha+1) * tgamma(beta+1) / tgamma(alpha+beta+2);
    if (n > 1)
      {
	ab(1,0) = (beta*beta-alpha*alpha)/(( 2+alpha+beta)*((2+alpha+beta)+2));
	ab(1,1) = 4*(alpha+1)*(beta+1)/(square(alpha+beta+2)*(alpha+beta+3));
      }
    
    for (int i = 0; i < n-2; i++) 
      {
        nab(i) = 2*(i+2)+alpha+beta;
        ab(i+2,0) = (beta*beta-alpha*alpha)/(nab(i)*(nab(i)+2));
        ab(i+2,1) = 4*(i+2+alpha)*(i+2+beta)*(i+2)*(i+2+alpha+beta)
          /(square(nab(i))*(nab(i)+1)*(nab(i)-1));
      } 
  }
  
  
  //! evaluation of Jacobi polynomials defined by relation reccurrence ab
  /*!
    We have P_0 = 1, P_1 = x - ab(0, 0), and then :
    P_k(x) = (x - ab(k-1, 0)) P_{k-1}(x) - ab(k-1, 1) P_{k-2}(x)
   */
  template<class T>
  void EvaluateJacobiPolynomial(const Matrix<T>& ab, int n,
                                const T& x, Vector<T>& Pn)
  {  
    if (n < 0)
      return;
    
    Pn.Reallocate(n+1);
    // evaluating P_1
    Pn(0) = 1.0;
    if (n > 0)
      Pn(1) = x - ab(0, 0);
    
    for (int i = 1; i < n; i++)
      Pn(i+1) = (x - ab(i, 0))*Pn(i) - ab(i, 1)*Pn(i-1);
    
  }
  
  
  //! evaluation of a Jacobi polynomial defined by relation reccurrence ab
  /*!
    We have P_0 = 1, P_1 = x - ab(0, 0), and then :
    P_k(x) = (x - ab(k-1, 0)) P_{k-1}(x) - ab(k-1, 1) P_{k-2}(x)
   */
  template<class T>
  T EvaluateJacobiPolynomial(const Matrix<T>& ab, int n, const T& x)
  {  
    if (n < 0)
      return 0;
    
    if (n == 0)
      return 1.0;
    
    // evaluating P_1
    T P0 = 1.0, P1, Pn;
    P1 = x - ab(0, 0);
    
    for (int i = 1; i < n; i++)
      {
	Pn = (x - ab(i, 0))*P1 - ab(i, 1)*P0;
	P0 = P1; P1 = Pn;
      }
    
    return P1;    
  }
  
  
  //! evaluation of Jacobi polynomials (and derivatives) defined by relation reccurrence ab
  /*!
    We have P_0 = 1, P_1 = x - ab(0, 0), and then :
    P_k(x) = (x - ab(k-1, 0)) P_{k-1}(x) - ab(k-1, 1) P_{k-2}(x)
   */
  template<class T>
  void EvaluateJacobiPolynomial(const Matrix<T>& ab, int n, const T& x,
				Vector<T>& Pn, Vector<T>& dPn)
  {  
    if (n < 0)
      return ;
    
    Pn.Reallocate(n+1); dPn.Reallocate(n+1);
    Pn(0) = 1.0; dPn(0) = 0.0;
    // evaluating P_1
    if (n > 0)
      {
	Pn(1) = x - ab(0, 0);
	dPn(1) = 1.0;
      }
    
    for (int i = 1; i < n; i++)
      {
	Pn(i+1) = (x - ab(i, 0))*Pn(i) - ab(i, 1)*Pn(i-1);
	dPn(i+1) = Pn(i) + (x - ab(i, 0))*dPn(i) - ab(i, 1)*dPn(i-1);
      }
  }
  
  
  //! finds the roots of orthogonal polynomial, and associated weight quadrature
  /*!
    \param[out] x roots of the last orthogonal polynomial defined by relations ab
    \param[out] w quadrature weights
    \param[in] n polynomial number
    \param[in] alpha exponent in weight (1-x)^alpha
    \param[in] beta exponent in weight (1+x)^beta
    \param[in] ab array filled by GetJacobiPolynomial
    This method shifts points and weights on interval [0, 1] instead of [-1,1]
   */
  template<class T>
  void SolveGauss(Vector<T>& x, Vector<T>& w, int n, const Matrix<T>& ab)
  {  
    x.Reallocate(n); w.Reallocate(n);
    Vector<T> xreal(n), xim(n);
    Matrix<T> J(n,n), vp(n,n);
    J.Fill(0);
  
    for (int i = 0; i < n; i++) 
      J(i,i) = ab(i,0);
    
    for (int i = 1; i < n; i++) 
      {
        J(i,i-1) = Sqrt(ab(i,1));
        J(i-1,i) = J(i,i-1);
      }

    GetEigenvaluesEigenvectors(J, xreal, xim, vp);
  
    w.Fill(0);
    for (int i = 0; i < n; i++)
      {
        x(i) = (xreal(i) + 1)/2;
        w(i) = ab(0,1)*square(vp(0, i));
      }
    Sort(x,w);
  }

  
  // Gauss formulas : \int_0^1 f(x) dx = \sum_{j=1}^N wj f(xj)
  template<class T>
  void ComputeGaussLegendre(Vector<T>& x, Vector<T>& w, int r)
  {
    ComputeGaussFormulas(r+1, x, w, false);
  }


  // Gauss-Lobatto formulas : \int_0^1 f(x) dx = \sum_{j=2}^{N-1} wj f(xj) + w1 f(0) + wN f(1)
  template<class T>
  void ComputeGaussLobatto(Vector<T>& x, Vector<T>& w, int r)
  {
    ComputeGaussFormulas(r+1, x, w, true);
  }
  
  
  // Gauss-Jacobi : \int_0^1 (1−x)^α x^β f(x)dx = \sum_{j=1}^N wj f(xj)
  template<class T>
  void ComputeGaussJacobi(Vector<T>& x, Vector<T>& w, int r,
                          const T& alpha, const T& beta)
  {
    int n = r+1;
    Matrix<T> ab;
    GetJacobiPolynomial(ab, n, alpha, beta);
    SolveGauss(x, w, n, ab);
    Sort(x, w);
  }

  
  // Gauss-Radau-Jacobi : 
  // \int_0^1 (1−x)^α x^β f(x)dx = \sum_{j=1}^{N-1} wj f(xj) + wN f(1)  (ext = true)
  // \int_0^1 (1−x)^α x^β f(x)dx = w1 f(0) + \sum_{j=2}^N wj f(xj)      (ext = false)
  template<class T>
  void ComputeGaussRadauJacobi(Vector<T>& x, Vector<T>& w, int r,
			       const T& alpha, const T& beta, bool ext)
  {
    int n = r+1;
    Matrix<T> ab;
    GetJacobiPolynomial(ab, n, alpha, beta);
    
    if (!ext)
      ab(n-1,0) = -1+2*(n-1)*(n-1+alpha)/((2*(n-1)+alpha+beta)*(2*(n-1)+alpha+beta+1));
    else
      ab(n-1,0) = 1-2*(n-1)*(n-1+beta)/((2*(n-1)+alpha+beta)*(2*(n-1)+alpha+beta+1));
    
    SolveGauss(x, w, n, ab);
    Sort(x, w);
  }

  
  // Gauss-Lobatto-Jacobi : 
  // \int_0^1 (1−x)^α x^β f(x)dx = w1 f(0)  + \sum_{j=2}^{N-1} wj f(xj) + wN f(1)
  template<class T>
  void ComputeGaussLobattoJacobi(Vector<T>& x, Vector<T>& w, int r,
				 const T& alpha, const T& beta)
  {
    int n = r+1;
    Matrix<T> ab;
    GetJacobiPolynomial(ab, n, alpha, beta);

    ab(n-1,0) = (alpha-beta)/(2*(n-2)+alpha+beta+2);
    ab(n-1,1) = 4*((n-2)+alpha+1)*((n-2)+beta+1)*((n-2)+alpha+beta+1)
      /((2*(n-2)+alpha+beta+1)*square(T(2*(n-2)+alpha+beta+2)));
    
    SolveGauss(x, w, n, ab);
    Sort(x, w);
  }

  
  // Gauss-Chebyshev formulas
  template<class T>
  void ComputeGaussChebyshev(Vector<T>& x, Vector<T>& w, int r)
  {
    x.Reallocate(r+1);
    w.Reallocate(r+1);
    T one(1);
    for (int i = 0; i <= r; i++)
      {
	x(i) = (one + cos(T(2*i+1)/T(2*(r+1))*pi_wp))/2;
	w(i) = pi_wp/(2*(r+1));
      }
  }

  
  // Gauss formulas for weight t^a log(1/t)
  template<class T>
  void ComputeGaussLogarithmic(Vector<T>& x, Vector<T>& w, int r, const T& a)
  {
    Matrix<T> abm;
    GetJacobiPolynomial(abm, 2*(r+1), T(0), T(0));
    
    T one(1);
    for (int i = 0; i < 2*(r+1); i++)
      abm(i, 0) = (one+abm(i, 0))/2;
    
    abm(0, 1) = abm(0, 1)/2;
    for (int i = 1; i < 2*(r+1); i++)
      abm(i, 1) = abm(i, 1)/4;
    
    // function mm = mm_log(N, a)
    T c(1); Vector<T> mm(2*(r+1));
    for (int n = 0; n < 2*(r+1); n++)
      {
	if ( (a-floor(a) == T(0)) && (a < n))
	  {
	    T prod_p(one);
	    for (int p = n-toInteger(a); p <= n+toInteger(a)+1; p++)
	      prod_p *= p;
	    
	    mm(n) = one/prod_p;
	    if ((n-toInteger(a))%2 == 1)
	      mm(n) = -mm(n);
	    
	    mm(n) *= square(tgamma(a+1));
	  }
	else
	  {
	    if (n==0)
	      mm(0) = one/square(a+1);
	    else 
	      {
		T sum(0), prod(1);
		for (int k = 1; k <= n; k++)
		  {
		    sum += one/(a+1+k) - one/(a+1-k);
		    prod *= (a+1-k)/(a+1+k);
		  }
		
		mm(n) = (one/(a+1)+sum)*prod/(a+1);
	      }
	  }
	
	mm(n) = c*mm(n);
	c = 0.5*(n+1)*c/(2*n+1);
      }
    
    // function ab = chebyshev(N,mom,abj);
    Matrix<T> ab(r+1, 2);
    ab(0, 0) = abm(0, 0) + mm(1)/mm(0);
    ab(0, 1) = mm(0);
    if (r == 0)
      return;
    
    Matrix<T> sig(2*(r+1), 2*(r+1));
    sig.Fill(0);
    for (int i = 0; i < 2*(r+1); i++)
      sig(1, i) = mm(i);
    
    for (int n = 3; n <= (r+1)+1; n++)
      {
	for (int m = n-1; m <= 2*(r+1)-n+2; m++)
	  {
	    sig(n-1,m-1) = sig(n-2, m) - (ab(n-3, 0)-abm(m-1, 0))*sig(n-2, m-1)
	      - ab(n-3, 1)*sig(n-3, m-1) + abm(m-1, 1)*sig(n-2,m-2);
	  }

	ab(n-2, 0) = abm(n-2, 0) + sig(n-1, n-1)/sig(n-1, n-2)
	  - sig(n-2, n-2)/ sig(n-2,n-3);
	
	ab(n-2, 1) = sig(n-1, n-2)/sig(n-2,n-3);
      }
    
    // deducing points and weights
    SolveGauss(x, w, r+1, ab);
    Sort(x, w);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = 2*x(i)-one;
  }


  // Gauss formulas for weight t^a exp(-t)
  template<class T>
  void ComputeGaussLaguerre(Vector<T>& x, Vector<T>& w, int r, const T& a)
  {
    Matrix<T> ab(r+1, 2);

    T nu = a+1, mu = tgamma(a+1);
    ab(0, 0) = nu; ab(0, 1) = mu;
    for (int n = 1; n <= r; n++)
      {
	ab(n, 0) = 2*n+a+1;
	ab(n, 1) = n*(n+a);
      }
    
    // deducing points and weights
    SolveGauss(x, w, r+1, ab);
    Sort(x, w);

    T one(1);
    for (int i = 0; i < x.GetM(); i++)
      x(i) = 2*x(i)-one;    
  }

  
  // Gauss formulas for weight |t|^a exp(-t^2)
  template<class T>
  void ComputeGaussHermite(Vector<T>& x, Vector<T>& w, int r, const T& a)
  {
    Matrix<T> ab(r+1, 2);

    T one(1), zero(0);
    T m0 = tgamma((a+one)/2);
    ab(0, 0) = zero; ab(0, 1) = m0;
    for (int n = 1; n <= r; n++)
      {
	ab(n, 0) = zero;
	ab(n, 1) = T(n)/2;
	if (n%2 == 1)
	  ab(n, 1) += a/2;
      }
    
    // deducing points and weights
    SolveGauss(x, w, r+1, ab);
    Sort(x, w);

    for (int i = 0; i < x.GetM(); i++)
      x(i) = 2*x(i)-one;    
  }
  
}

#define MONTJOIE_FILE_GAUSS_JACOBI_CXX
#endif
