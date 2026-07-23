#ifndef MONTJOIE_FILE_GAUSS_JACOBI_HXX

namespace Montjoie
{
  template<class T>
  void GetJacobiPolynomial(Matrix<T>& ab, int n, const T& alpha, const T& beta);
  
  template<class T>
  void EvaluateJacobiPolynomial(const Matrix<T>& ab, int n,
                                const T& x, Vector<T>& Pn);

  template<class T>
  T EvaluateJacobiPolynomial(const Matrix<T>& ab, int n, const T& x);
  
  template<class T>
  void EvaluateJacobiPolynomial(const Matrix<T>& ab, int n, const T& x,
                                Vector<T>& Pn, Vector<T>& dPn);

  template<class T>
  void SolveGauss(Vector<T>& x, Vector<T>& w, int n, const Matrix<T>& ab);

  template<class T>
  void ComputeGaussLegendre(Vector<T>& x, Vector<T>& w, int r);
  
  template<class T>
  void ComputeGaussLobatto(Vector<T>& x, Vector<T>& w, int r);
  
  template<class T>
  void ComputeGaussJacobi(Vector<T>& x, Vector<T>& w, int r,
                          const T& alpha, const T& beta);
  
  template<class T>
  void ComputeGaussRadauJacobi(Vector<T>& x, Vector<T>& w, int r,
			       const T& alpha, const T& beta, bool ext);
  
  template<class T>
  void ComputeGaussLobattoJacobi(Vector<T>& x, Vector<T>& w, int r,
				 const T& alpha, const T& beta);

  template<class T>
  void ComputeGaussChebyshev(Vector<T>& x, Vector<T>& w, int r);

  template<class T>
  void ComputeGaussLogarithmic(Vector<T>& x, Vector<T>& w, int r, const T& alpha);
  
  template<class T>
  void ComputeGaussLaguerre(Vector<T>& x, Vector<T>& w, int r, const T& a);

  template<class T>
  void ComputeGaussHermite(Vector<T>& x, Vector<T>& w, int r, const T& a);

} // end namespace

#define MONTJOIE_FILE_GAUSS_JACOBI_HXX
#endif



