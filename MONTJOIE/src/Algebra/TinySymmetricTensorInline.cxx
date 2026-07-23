#ifndef SELDON_FILE_TINY_SYMMETRIC_TENSOR_INLINE_CXX

#include "TinySymmetricTensor.hxx"

namespace Seldon
{

  /***********************
   * TinySymmetricTensor *
   ***********************/

  
  //! default constructor : coefficients are set to 0
  template<class T, int m>
  inline TinySymmetricTensor<T, m>::TinySymmetricTensor() : C()
  {
    id.Fill();
  }
  
  
  //! all the coefficients are set to a given value x
  template<class T, int m> template<class T0>
  inline void TinySymmetricTensor<T, m>::Fill(const T0& x)
  {
    C.Fill(x);
  }
  
  
  //! returns C(i, j) (an elastic tensor can be viewed as a matrix)
  template<class T, int m>
  inline const T& TinySymmetricTensor<T, m>::operator()(int i, int j) const
  {
    return C(i, j);
  }
  
  
  //! returns C(i, j) (an elastic tensor can be viewed as a matrix)
  template<class T, int m>
  inline T& TinySymmetricTensor<T, m>::operator()(int i, int j)
  {
    return C(i, j);
  }
  
  
  //! returns C_{i, j, k, l}
  template<class T, int m>
  inline const T& TinySymmetricTensor<T, m>::operator()(int i, int j, int k, int l) const
  {
    return C(id(i, j), id(k, l));
  }
  
  
  //! returns C_{i, j, k, l}
  template<class T, int m>
  inline T& TinySymmetricTensor<T, m>::operator()(int i, int j, int k, int l)
  {
    return C(id(i, j), id(k, l));
  }
  
  
  //! returns the symmetric tensor under the form of a symmetric matrix
  /*template<class T, int m>
  inline const TinyMatrix<T, Symmetric, TinySymmetricTensor<T, m>::N,
			  TinySymmetricTensor<T, m>::N>&
  TinySymmetricTensor<T, m>::GetTensor() const
  {
    return C;
    }*/


  //! Computation of dV = C dU, assuming that C is orthotropic
  template<class T, int m> template<class T1>
  inline void TinySymmetricTensor<T, m>
  ::MltOrthotrope(const TinyVector<TinyVector<T1, m>, m>& dU,
                  TinyVector<TinyVector<T1, m>, m> & dV) const
  {
    TinySymTensorLoop<m-1>::MltOrtho(*this, dU, dV);      
  }


  //! Computation of dV = C dU
  template<class T, int m> template<class T1>
  inline void TinySymmetricTensor<T, m>
  ::Mlt(const TinyVector<TinyVector<T1, m>, m>& dU,
        TinyVector<TinyVector<T1, m>, m> & dV) const
  {
    FillZero(dV);
    TinySymTensorLoop<m-1>::MltAdd(*this, dU, dV);      
    TinySymTensorLoop<m-1>::Sym(dV);
  }


  //! Computation of dV = C dU, assuming that C is orthotropic
  template<class T, int m> template<class T1>
  inline void TinySymmetricTensor<T, m>
  ::MltOrthotrope(const TinyVector<T1, m*m>& dU, TinyVector<T1, m*m> & dV) const
  {
    TinySymTensorLoop<m-1>::MltOrtho(*this, dU, dV);      
  }


  //! Computation of dV = C dU
  template<class T, int m> template<class T1>
  inline void TinySymmetricTensor<T, m>
  ::Mlt(const TinyVector<T1, m*m>& dU, TinyVector<T1, m*m> & dV) const
  {
    FillZero(dV);
    TinySymTensorLoop<m-1>::MltAdd(*this, dU, dV);      
    TinySymTensorDoubleLoop<m, m-1>::Sym(dV);
  }


  //! product of a tensor by a scalar this = this * a where a is a scalar
  template<class T, int m> template<class T1>
  inline TinySymmetricTensor<T, m>& TinySymmetricTensor<T, m>::operator *=(const T1& a )
  {
    C *= a;
    return *this;
  }


  //! returns maximum modulus of coefficients stored in the tensor
  template<class T, int m>
  inline typename ClassComplexType<T>::Treal TinySymmetricTensor<T, m>
  ::GetNormInf() const
  {
    return MaxAbs(C);
  }
  

  /*************************
   * TinySymTensorQuadLoop *
   *************************/

  
  template<int p, int q, int r, int s> template<class T, int m, class T1>
  inline void TinySymTensorQuadLoop<p, q, r, s>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const T1& eps, TinyVector<TinyVector<T1, m>, m> & dV)
  {
    dV(r)(s) += C(TinyTensorNode<m, p, q>::p,
		  TinyTensorNode<m, r, s>::p) * eps;
    
    TinySymTensorQuadLoop<p, q, r, s-1>::MltNonDiag(C, eps, dV);
  }
  
  
  template<int p, int q, int r, int s> template<class T, int m, class T1>
  inline void TinySymTensorQuadLoop<p, q, r, s>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const T1& eps, TinyVector<T1, m*m> & dV)
  {
    dV(r*m + s) += C(TinyTensorNode<m, p, q>::p,
		     TinyTensorNode<m, r, s>::p) * eps;
    
    TinySymTensorQuadLoop<p, q, r, s-1>::MltNonDiag(C, eps, dV);
  }
  

  /***************************
   * TinySymTensorTripleLoop *
   ***************************/

  
  template<int p, int q, int r> template<class T, int m, class T1>
  inline void TinySymTensorTripleLoop<p, q, r>
  ::MltDiag(const TinySymmetricTensor<T, m>& C,
	    const TinyVector<TinyVector<T1, m>, m>& dU,
	    TinyVector<TinyVector<T1, m>, m> & dV)
  {
    dV(q)(r) += C(TinyTensorNode<m, p, p>::p,
		  TinyTensorNode<m, q, r>::p) * dU(p)(p);
    
    TinySymTensorTripleLoop<p, q, r-1>::MltDiag(C, dU, dV);
  }
  
  template<int p, int q, int r> template<class T, int m, class T1>
  inline void TinySymTensorTripleLoop<p, q, r>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const T1& eps, TinyVector<TinyVector<T1, m>, m> & dV)
  {      
    TinySymTensorQuadLoop<p, q, r, r>::MltNonDiag(C, eps, dV);
    TinySymTensorTripleLoop<p, q, r-1>::MltNonDiag(C, eps, dV);
  }
    
  
  template<int p, int q, int r> template<class T, int m, class T1>
  inline void TinySymTensorTripleLoop<p, q, r>
  ::MltDiag(const TinySymmetricTensor<T, m>& C,
	    const TinyVector<T1, m*m>& dU,
	    TinyVector<T1, m*m> & dV)
  {
    dV(q*m+r) += C(TinyTensorNode<m, p, p>::p,
		   TinyTensorNode<m, q, r>::p) * dU(p*m+p);
    
    TinySymTensorTripleLoop<p, q, r-1>::MltDiag(C, dU, dV);
  }

  template<int p, int q, int r> template<class T, int m, class T1>
  inline void TinySymTensorTripleLoop<p, q, r>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const T1& eps, TinyVector<T1, m*m> & dV)
  {      
    TinySymTensorQuadLoop<p, q, r, r>::MltNonDiag(C, eps, dV);
    TinySymTensorTripleLoop<p, q, r-1>::MltNonDiag(C, eps, dV);
  }
  
  template<int p, int q, int r> template<class T1>
  inline void TinySymTensorTripleLoop<p, q, r>::Sym(TinyVector<T1, p*p> & dV)
  {
    dV(r*p+q) = dV(q*p+r);
    TinySymTensorTripleLoop<p, q, r-1>::Sym(dV);
  }


  /***************************
   * TinySymTensorDoubleLoop *
   ***************************/

  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltOrtho(const TinySymmetricTensor<T, m>& C,
	     const TinyVector<TinyVector<T1, m>, m>& dU,
	     TinyVector<TinyVector<T1, m>, m> & dV)
  {
    dV(p)(q) = C(TinyTensorNode<m, p, q>::p,
		 TinyTensorNode<m, p, q>::p) * (dU(p)(q) + dU(q)(p));
    
    dV(q)(p) = dV(p)(q);
    
    TinySymTensorDoubleLoop<p, q-1>::MltOrtho(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
		 const TinyVector<TinyVector<T1, m>, m>& dU,
		 TinyVector<TinyVector<T1, m>, m> & dV)
  {
    dV(p)(p) += C(TinyTensorNode<m, p, p>::p,
		  TinyTensorNode<m, q, q>::p) * dU(q)(q);
    
    TinySymTensorDoubleLoop<p, q-1>::MltOrthoDiag(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltDiag(const TinySymmetricTensor<T, m>& C,
	    const TinyVector<TinyVector<T1, m>, m>& dU,
	    TinyVector<TinyVector<T1, m>, m> & dV)
  {
    TinySymTensorTripleLoop<p, q, q>::MltDiag(C, dU, dV);
    TinySymTensorDoubleLoop<p, q-1>::MltDiag(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const TinyVector<TinyVector<T1, m>, m>& dU,
	       TinyVector<TinyVector<T1, m>, m> & dV)
  {
    T1 eps = dU(p)(q) + dU(q)(p);
    TinySymTensorTripleLoop<p, q, m-1>::MltNonDiag(C, eps, dV);
    TinySymTensorDoubleLoop<p, q-1>::MltNonDiag(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltOrtho(const TinySymmetricTensor<T, m>& C,
	     const TinyVector<T1, m*m>& dU,
	     TinyVector<T1, m*m> & dV)
  {
    dV(p*m+q) = C(TinyTensorNode<m, p, q>::p,
		  TinyTensorNode<m, p, q>::p) * (dU(p*m+q) + dU(q*m+p));
    
    dV(q*m+p) = dV(p*m+q);
    
    TinySymTensorDoubleLoop<p, q-1>::MltOrtho(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
		 const TinyVector<T1, m*m>& dU,
		 TinyVector<T1, m*m> & dV)
  {
    dV(p*m+p) += C(TinyTensorNode<m, p, p>::p,
		   TinyTensorNode<m, q, q>::p) * dU(q*m+q);
    
    TinySymTensorDoubleLoop<p, q-1>::MltOrthoDiag(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltDiag(const TinySymmetricTensor<T, m>& C,
	    const TinyVector<T1, m*m>& dU,
	    TinyVector<T1, m*m> & dV)
  {
    TinySymTensorTripleLoop<p, q, q>::MltDiag(C, dU, dV);
    TinySymTensorDoubleLoop<p, q-1>::MltDiag(C, dU, dV);
  }
  
  template<int p, int q> template<class T, int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>
  ::MltNonDiag(const TinySymmetricTensor<T, m>& C,
	       const TinyVector<T1, m*m>& dU,
	       TinyVector<T1, m*m> & dV)
  {
    T1 eps = dU(p*m+q) + dU(q*m+p);
    TinySymTensorTripleLoop<p, q, m-1>::MltNonDiag(C, eps, dV);
    TinySymTensorDoubleLoop<p, q-1>::MltNonDiag(C, dU, dV);
  }
  
  template<int p, int q> template<int m, class T1>
  inline void TinySymTensorDoubleLoop<p, q>::Sym(TinyVector<TinyVector<T1, m>, m> & dV)
  {
    dV(q)(p) = dV(p)(q);
    TinySymTensorDoubleLoop<p, q-1>::Sym(dV);
  }
  
  template<int p, int q> template<class T1>
  inline void TinySymTensorDoubleLoop<p, q>::Sym(TinyVector<T1, p*p> & dV)
  {
    TinySymTensorTripleLoop<p, q, p-1>::Sym(dV);
    
    TinySymTensorDoubleLoop<p, q-1>::Sym(dV);
  }


  /*********************
   * TinySymTensorLoop *
   *********************/

  
  template<int p> template<class T, int m, class T1>
  inline void TinySymTensorLoop<p>::MltOrtho(const TinySymmetricTensor<T, m>& C,
					     const TinyVector<TinyVector<T1, m>, m>& dU,
					     TinyVector<TinyVector<T1, m>, m> & dV)
  {
    // term dV(i)(j) with i different from j
    TinySymTensorDoubleLoop<p, p-1>::MltOrtho(C, dU, dV);
    
    // term dV(i)(i)
    SetComplexZero(dV(p)(p));
    TinySymTensorDoubleLoop<p, m-1>::MltOrthoDiag(C, dU, dV);
    
    TinySymTensorLoop<p-1>::MltOrtho(C, dU, dV);
  }
  
  template<int p> template<class T, int m, class T1>
  inline void TinySymTensorLoop<p>::MltAdd(const TinySymmetricTensor<T, m>& C,
					   const TinyVector<TinyVector<T1, m>, m>& dU,
					   TinyVector<TinyVector<T1, m>, m> & dV)
  {
    // term dV(i)(j) with i different from j
    TinySymTensorDoubleLoop<p, p-1>::MltNonDiag(C, dU, dV);
    
    // term dV(i)(i)
    TinySymTensorDoubleLoop<p, m-1>::MltDiag(C, dU, dV);
    
    TinySymTensorLoop<p-1>::MltAdd(C, dU, dV);
  }

  template<int p> template<class T, int m, class T1>
  inline void TinySymTensorLoop<p>::MltOrtho(const TinySymmetricTensor<T, m>& C,
					     const TinyVector<T1, m*m>& dU,
					     TinyVector<T1, m*m> & dV)
  {
    // term dV(i)(j) with i different from j
    TinySymTensorDoubleLoop<p, p-1>::MltOrtho(C, dU, dV);
    
    // term dV(i)(i)
    SetComplexZero(dV(p*m+p));
    TinySymTensorDoubleLoop<p, m-1>::MltOrthoDiag(C, dU, dV);
    
    TinySymTensorLoop<p-1>::MltOrtho(C, dU, dV);
  }
  
  template<int p> template<class T, int m, class T1>
  inline void TinySymTensorLoop<p>::MltAdd(const TinySymmetricTensor<T, m>& C,
					   const TinyVector<T1, m*m>& dU,
					   TinyVector<T1, m*m> & dV)
  {
    // term dV(i)(j) with i different from j
    TinySymTensorDoubleLoop<p, p-1>::MltNonDiag(C, dU, dV);
    
    // term dV(i)(i)
    TinySymTensorDoubleLoop<p, m-1>::MltDiag(C, dU, dV);
    
    TinySymTensorLoop<p-1>::MltAdd(C, dU, dV);
  }

  template<int p> template<int m, class T1>
  inline void TinySymTensorLoop<p>::Sym(TinyVector<TinyVector<T1, m>, m> & dV)
  {
    TinySymTensorDoubleLoop<p, p-1>::Sym(dV);
    
    TinySymTensorLoop<p-1>::Sym(dV);
  }


  /*************
   * Functions *
   *************/

  
  template<class T1, class T, int m> 
  inline void Mlt(const T1& alpha, TinySymmetricTensor<T, m>& C)
  {
    C *= alpha;
  }

  template<class T, int p>
  inline typename ClassComplexType<T>::Treal MaxAbs(const TinySymmetricTensor<T, p>& C)
  {
    return C.GetNormInf();
  }
  
  
  template<class T, int m>
  inline ostream& operator <<(ostream& out, const TinySymmetricTensor<T, m>& A)
  {
    out << A.C;
    return out;
  }
  

  /**************************
   * TinySymmetricVecTensor *
   **************************/

  
  //! default constructor : coefficients are set to 0
  template<class T, int m>
  inline TinySymmetricVecTensor<T, m>::TinySymmetricVecTensor() : d()
  {
    id.Fill();
  }
  
  
  //! all the coefficients are set to a given value x
  template<class T, int m> template<class T0>
  inline void TinySymmetricVecTensor<T, m>::Fill(const T0& x)
  {
    d.Fill(x);
  }
  
  
  //! returns d_{k, i, j}
  template<class T, int m>
  inline const T& TinySymmetricVecTensor<T, m>::operator()(int k, int i, int j) const
  {
    return d(k, id(i, j));
  }
  
  
  //! returns d_{k, i, j}
  template<class T, int m>
  inline T& TinySymmetricVecTensor<T, m>::operator()(int k, int i, int j)
  {
    return d(k, id(i, j));
  }


  template<class T, int m>
  inline ostream& operator <<(ostream& out, const TinySymmetricVecTensor<T, m>& A)
  {
    out << A.d;
    return out;
  }

} // end namespace

#define SELDON_FILE_TINY_SYMMETRIC_TENSOR_INLINE_CXX
#endif
