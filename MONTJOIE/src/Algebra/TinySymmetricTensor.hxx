#ifndef SELDON_FILE_TINY_SYMMETRIC_TENSOR_HXX

namespace Seldon
{

  template<int p>
  class TinySymTensorLoop;

  template<int p, int q>
  class TinySymTensorDoubleLoop;
  
  //! base class for a symmetric tensor (for elastics)
  template<class T, int m>
  class TinySymmetricTensor
  {
    template<int p>
    friend class TinySymTensorLoop;

    template<int p, int q>
    friend class TinySymTensorDoubleLoop;
    
  protected :
    //! N = size of the associated matrix
    enum {N = m*(m+1)/2};    
    
    //! values of the symmetric tensor
    TinyMatrix<T, Symmetric, N, N> C;
    
    //! numbering for the associated matrix
    static TinyMatrix<int, Symmetric, m, m> id;
    
  public :
    TinySymmetricTensor();

    // inline methods
    template<class T0>
    void Fill(const T0& x);
    
    const T& operator()(int i, int j) const;
    T& operator()(int i, int j);
    
    const T& operator()(int i, int j, int k, int l) const;
    T& operator()(int i, int j, int k, int l);
    
    static inline int GetSize() { return TinyMatrix<T, Symmetric, N, N>::GetSize(); }
    inline const TinyMatrix<T, Symmetric, N, N>& GetTensor() const { return C; }

    template<class T1>
    void MltOrthotrope(const TinyVector<TinyVector<T1, m>, m>& dU,
                       TinyVector<TinyVector<T1, m>, m> & dV) const;

    template<class T1>
    void Mlt(const TinyVector<TinyVector<T1, m>, m>& dU,
             TinyVector<TinyVector<T1, m>, m> & dV) const;

    template<class T1>
    void MltOrthotrope(const TinyVector<T1, m*m>& dU,
                       TinyVector<T1, m*m> & dV) const;

    template<class T1>
    void Mlt(const TinyVector<T1, m*m>& dU,
             TinyVector<T1, m*m> & dV) const;

    template<class T1>
    TinySymmetricTensor<T, m>& operator *=(const T1& a );
    
    typename ClassComplexType<T>::Treal GetNormInf() const;

    // non-inline methods
    void FillIsotrope(const T& lambda, const T& mu);
    void AddIsotrope(const T& lambda, const T& mu);
    
    template<class T0>
    void ApplyRotation(const TinyMatrix<T0, General, m, m>& matQ);

    void GetInverse();
    
    template<class T1, int p>
    friend ostream& operator <<(ostream& out, const TinySymmetricTensor<T1, p>& A);
        
  };
  
  //! numbering used for the matrix associated with the symmetric tensor
  template<class T, int m>
  TinyMatrix<int, Symmetric, m, m> TinySymmetricTensor<T, m>::id;  
  
  //! numbering used for the matrix associated with the symmetric tensor
  template<int n, int i, int j>
  class TinyTensorNode
  {
  public :
    
    enum {p = (j > i)*( i * n - (i * (i + 1)) / 2 + j)
          + (j <= i)*( j * n - (j * (j + 1)) / 2 + i) };
    
  };
  
  //! class for unrolling quadruple loops of symmetric tensors
  template<int p, int q, int r, int s>
  class TinySymTensorQuadLoop
  {
  public :
    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const T1& eps, TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const T1& eps, TinyVector<T1, m*m> & dV);
    
  };


  //! class for terminating quadruple loops of symmetric tensors
  template<int p, int q, int r>
  class TinySymTensorQuadLoop<p, q, r, -1>
  {
  public :
    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const T1& eps,
                                  TinyVector<TinyVector<T1, m>, m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const T1& eps,
                                  TinyVector<T1, m*m> & dV)
    {
    }
  };


  //! class for unrolling triple loops of symmetric tensors  
  template<int p, int q, int r>
  class TinySymTensorTripleLoop
  {
  public :
    template<class T, int m, class T1>
    static void MltDiag(const TinySymmetricTensor<T, m>& C,
			const TinyVector<TinyVector<T1, m>, m>& dU,
			TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const T1& eps,
			   TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltDiag(const TinySymmetricTensor<T, m>& C,
			const TinyVector<T1, m*m>& dU,
			TinyVector<T1, m*m> & dV);

    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const T1& eps,
			   TinyVector<T1, m*m> & dV);

    template<class T1>
    static void Sym(TinyVector<T1, p*p> & dV);
    
  };

  
  //! class for terminating triple loops of symmetric tensors
  template<int p, int q>
  class TinySymTensorTripleLoop<p, q, -1>
  {
  public :
    template<class T, int m, class T1>
    static inline void MltDiag(const TinySymmetricTensor<T, m>& C,
                               const TinyVector<TinyVector<T1, m>, m>& dU,
                               TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const T1& eps,
                                  TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltDiag(const TinySymmetricTensor<T, m>& C,
                               const TinyVector<T1, m*m>& dU,
                               TinyVector<T1, m*m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const T1& eps,
                                  TinyVector<T1, m*m> & dV)
    {
    }

    template<class T1>
    static inline void Sym(TinyVector<T1, p*p> & dV)
    {}
    
  };


  //! class for unrolling double loops of symmetric tensors  
  template<int p, int q>
  class TinySymTensorDoubleLoop
  {
  public :
    template<class T, int m, class T1>
    static void MltOrtho(const TinySymmetricTensor<T, m>& C,
			 const TinyVector<TinyVector<T1, m>, m>& dU,
			 TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
			     const TinyVector<TinyVector<T1, m>, m>& dU,
			     TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltDiag(const TinySymmetricTensor<T, m>& C,
			const TinyVector<TinyVector<T1, m>, m>& dU,
			TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const TinyVector<TinyVector<T1, m>, m>& dU,
			   TinyVector<TinyVector<T1, m>, m> & dV);
    
    template<class T, int m, class T1>
    static void MltOrtho(const TinySymmetricTensor<T, m>& C,
			 const TinyVector<T1, m*m>& dU,
			 TinyVector<T1, m*m> & dV);
    
    template<class T, int m, class T1>
    static void MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
			     const TinyVector<T1, m*m>& dU,
			     TinyVector<T1, m*m> & dV);
    
    template<class T, int m, class T1>
    static void MltDiag(const TinySymmetricTensor<T, m>& C,
			const TinyVector<T1, m*m>& dU,
			TinyVector<T1, m*m> & dV);

    template<class T, int m, class T1>
    static void MltNonDiag(const TinySymmetricTensor<T, m>& C,
			   const TinyVector<T1, m*m>& dU,
			   TinyVector<T1, m*m> & dV);
    
    template<int m, class T1>
    static void Sym(TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T1>
    static void Sym(TinyVector<T1, p*p> & dV);
    
  };


  //! class for terminating double loops of symmetric tensors
  template<int p>
  class TinySymTensorDoubleLoop<p, -1>
  {
  public :
    template<class T, int m, class T1>
    static inline void MltOrtho(const TinySymmetricTensor<T, m>& C,
                                const TinyVector<TinyVector<T1, m>, m>& dU,
                                TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
                                    const TinyVector<TinyVector<T1, m>, m>& dU,
                                    TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltDiag(const TinySymmetricTensor<T, m>& C,
                               const TinyVector<TinyVector<T1, m>, m>& dU,
                               TinyVector<TinyVector<T1, m>, m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const TinyVector<TinyVector<T1, m>, m>& dU,
                                  TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltOrtho(const TinySymmetricTensor<T, m>& C,
                                const TinyVector<T1, m*m>& dU,
                                TinyVector<T1, m*m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltOrthoDiag(const TinySymmetricTensor<T, m>& C,
                                    const TinyVector<T1, m*m>& dU,
                                    TinyVector<T1, m*m> & dV)
    {
    }
    
    template<class T, int m, class T1>
    static inline void MltDiag(const TinySymmetricTensor<T, m>& C,
                               const TinyVector<T1, m*m>& dU,
                               TinyVector<T1, m*m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltNonDiag(const TinySymmetricTensor<T, m>& C,
                                  const TinyVector<T1, m*m>& dU,
                                  TinyVector<T1, m*m> & dV)
    {
    }
    
    template<int m, class T1>
    static inline void Sym(TinyVector<TinyVector<T1, m>, m> & dV) {}

    template<class T1>
    static inline void Sym(TinyVector<T1, p*p> & dV)
    {
    }
    
  };


  //! class for unrolling simple loops of symmetric tensors  
  template<int p>
  class TinySymTensorLoop
  {
  public :
    template<class T, int m, class T1>
    static void MltOrtho(const TinySymmetricTensor<T, m>& C,
			 const TinyVector<TinyVector<T1, m>, m>& dU,
			 TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltAdd(const TinySymmetricTensor<T, m>& C,
		       const TinyVector<TinyVector<T1, m>, m>& dU,
		       TinyVector<TinyVector<T1, m>, m> & dV);

    template<class T, int m, class T1>
    static void MltOrtho(const TinySymmetricTensor<T, m>& C,
			 const TinyVector<T1, m*m>& dU,
			 TinyVector<T1, m*m> & dV);

    template<class T, int m, class T1>
    static void MltAdd(const TinySymmetricTensor<T, m>& C,
		       const TinyVector<T1, m*m>& dU,
		       TinyVector<T1, m*m> & dV);

    template<int m, class T1>
    static void Sym(TinyVector<TinyVector<T1, m>, m> & dV);
    
  };

  //! class for terminating simple loops of symmetric tensors
  template<>
  class TinySymTensorLoop<-1>
  {
  public :
    template<class T, int m, class T1>
    static inline void MltOrtho(const TinySymmetricTensor<T, m>& C,
                                const TinyVector<TinyVector<T1, m>, m>& dU,
                                TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltAdd(const TinySymmetricTensor<T, m>& C,
                              const TinyVector<TinyVector<T1, m>, m>& dU,
                              TinyVector<TinyVector<T1, m>, m> & dV)
    {}

    template<class T, int m, class T1>
    static inline void MltOrtho(const TinySymmetricTensor<T, m>& C,
                                const TinyVector<T1, m*m>& dU,
                                TinyVector<T1, m*m> & dV)
    {
    }

    template<class T, int m, class T1>
    static inline void MltAdd(const TinySymmetricTensor<T, m>& C,
                              const TinyVector<T1, m*m>& dU,
                              TinyVector<T1, m*m> & dV)
    {}
    
    template<int m, class T1>
    static inline void Sym(TinyVector<TinyVector<T1, m>, m> & dV) {}

  };
  
  template<class T1, class T, int m> 
  void Mlt(const T1& alpha, TinySymmetricTensor<T, m>& C);

  template<class T, int p>
  typename ClassComplexType<T>::Treal MaxAbs(const TinySymmetricTensor<T, p>& C);
  
  // computation of eigenvalues of a symmetric tensor C
  template<class T, class T0, int p>
  void GetEigenvalues(const TinySymmetricTensor<T, p>& C,
                      TinyVector<T0, p>& k, TinyVector<T, p>& lambda);

  template<class T, int m>
  ostream& operator <<(ostream& out, const TinySymmetricTensor<T, m>& A);


  //! base class for a symmetric tensor to vector (piezo)
  template<class T, int m>
  class TinySymmetricVecTensor
  {
  protected :
    //! N = size of the associated matrix
    enum {N = m*(m+1)/2};    
    
    //! values of d_{k, i, j}
    TinyMatrix<T, General, m, N> d;
    
    //! numbering for the associated matrix
    static TinyMatrix<int, Symmetric, m, m> id;
    
  public :
    TinySymmetricVecTensor();

    // inline methods
    template<class T0>
    void Fill(const T0& x);
    
    const T& operator()(int k, int i, int j) const;
    T& operator()(int k, int i, int j);

    template<class T1, int p>
    friend ostream& operator <<(ostream& out, const TinySymmetricVecTensor<T1, p>& A);

  };

  //! numbering used for the matrix associated with the symmetric tensor
  template<class T, int m>
  TinyMatrix<int, Symmetric, m, m> TinySymmetricVecTensor<T, m>::id;  

  template<class T, int m>
  ostream& operator <<(ostream& out, const TinySymmetricVecTensor<T, m>& A);

}

#define SELDON_FILE_TINY_SYMMETRIC_TENSOR_HXX
#endif
