#ifndef MONTJOIE_FILE_MULTIVARIATE_POLYNOMIAL_HXX

namespace Montjoie
{
  
  //! Homogeneous polynome with multiple variables
  template<class T, class Allocator = 
           typename SeldonDefaultAllocator<VectFull, T>::allocator >
  class MultivariatePolynomial : public Vector<T, VectFull, Allocator>
  {
    // Attributes
  protected:
    //! number of variables of the polynom
    /*!
       \f$ (x_1, ... ,x_{\mbox{dimension}} ) \f$
    */
    unsigned dimension;
    //! degree of the polynom
    /*!
       we have \f$ \mbox{deg}(x_1) + \mbox{deg}(x_2) ...
       + \mbox{deg}(x_d) <= n_ \f$
     */
    int n_;
    // vector which contains the coefficients a_0, a_1, a_2 ...
    // P = a_0 + a_1 y + ... + a_n y^n + x * ( a_{n+1} + a_{n+2} y
    // + ... + a_{2*n-1} y^(n-1) ) + ... + a_{m_} x^n  
    // by example in dimension 2
    
    //! arrays used to find quickly monomial position
    IVect NumTenso, InvNumTenso;
    // NumTenso(i) returns N^{d-1} i_1 + N^{d-2} i_2 + ... + i_d
    // where i is associated with the monomial x_1^{i_1} \, x_2^{i_2} ... x_d^{i_d}
    
    // basic methods
  public:
    
    // constructors
    MultivariatePolynomial();
    MultivariatePolynomial(unsigned dimension_, int n);
    MultivariatePolynomial(const MultivariatePolynomial<T, Allocator>& P);

    ~MultivariatePolynomial();
    
    /*******************
     * Basic Functions *
     *******************/
    
    int GetM() const;
    size_t GetMemorySize() const;
    unsigned int GetDimension() const;
    int GetOrder() const;
    void SetOrder(unsigned d, int order);
    void ResizeOrder(unsigned d, int order);
    void Copy(const MultivariatePolynomial<T, Allocator>& Q);
    
  protected :
    
    /**********************
     * Internal Functions *
     **********************/
    
    void ConstructIndexArrays();
    void SubstituteLastVariable(const T& value);
    void SubstituteVariable(int num_variable, const T& value);
    
  public :
    
    
    /************************
     * Convenient Functions *
     ************************/
    
    int PowersToIndex(const IVect& powers) const;
    void IndexToPowers(int pos, IVect& powers) const;

    template<class TinyVector_>
    T Evaluate(const TinyVector_& x) const;
    
    void PrintRational(const T& epsilon) const;
    void PrintMonomial(int i, ostream& out) const;
        
    /*************
     * Operators *
     *************/
    
    T& operator() (int i);
    const T& operator() (int i) const;
    T& operator() (int i, int j);
    const T& operator() (int i, int j) const;
    T& operator() (int i, int j, int k);
    const T& operator() (int i, int j, int k) const;
    
    MultivariatePolynomial<T, Allocator>&
    operator =(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>
    operator +(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>
    operator -(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>
    operator *(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>&
    operator +=(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>&
    operator -=(const MultivariatePolynomial<T, Allocator>& P);
    MultivariatePolynomial<T, Allocator>&
    operator *=(const MultivariatePolynomial<T, Allocator>& P);

    template<class T0>
    MultivariatePolynomial<T, Allocator> operator *(const T0& alpha);
    template<class T0>
    MultivariatePolynomial<T, Allocator> operator +(const T0& alpha);
    template<class T0>
    MultivariatePolynomial<T, Allocator> operator -(const T0& alpha);
    template<class T0>
    MultivariatePolynomial<T, Allocator>& operator *=(const T0& alpha);
    template<class T0>
    MultivariatePolynomial<T, Allocator>& operator +=(const T0& alpha);
    template<class T0>
    MultivariatePolynomial<T, Allocator>& operator -=(const T0& alpha);

    bool operator ==(const MultivariatePolynomial<T, Allocator>& P);
    bool operator !=(const MultivariatePolynomial<T, Allocator>& P);
    
  };
  
  // other operators
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator *(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P);
  
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator +(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P);
  
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator -(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P);

  template<class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator -(const MultivariatePolynomial<T, Allocator>& P);

  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const MultivariatePolynomial<T, Allocator>& P);
  
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>
  Pow(const MultivariatePolynomial<T, Allocator>& U, int r);

  template<class T, class Allocator1, class Allocator2>
  void DerivatePolynomial(const MultivariatePolynomial<T, Allocator1>& P,
			 MultivariatePolynomial<T, Allocator2>& Q, int num_variable);
  
  template<class T, class Allocator> MultivariatePolynomial<T, Allocator>
  Derivate(const MultivariatePolynomial<T, Allocator>& P, int num_variable);

  template<class T, class Allocator1, class Allocator2>
  void GetCurlPolynomial(const TinyVector<MultivariatePolynomial<T, Allocator1>, 3>& P,
                         TinyVector<MultivariatePolynomial<T, Allocator2>, 3>& Q);

  template<class T, class Allocator1, class Allocator2>
  void GetDivPolynomial(const TinyVector<MultivariatePolynomial<T, Allocator1>, 3>& P,
                        MultivariatePolynomial<T, Allocator2>& Q);
  
  template<class T, class Allocator1, class Allocator2>
  void IntegratePolynomial(const MultivariatePolynomial<T, Allocator1>& P,
			  MultivariatePolynomial<T, Allocator2>& Q, int num_variable);

  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, MultivariatePolynomial<T, Allocator>& R);

  template<class T0, class T1, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const MultivariatePolynomial<T, Allocator1>& P,
	      const MultivariatePolynomial<T, Allocator2>& Q,
	      const T1& beta, MultivariatePolynomial<T, Allocator3>& R);

  template<class T, class Allocator>
  T IntegrateOnSimplex(const MultivariatePolynomial<T, Allocator>& P);

  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const T0& alpha, const MultivariatePolynomial<T1, Allocator1>& P,
	   MultivariatePolynomial<T2, Allocator2>& Q);

  template<class T, class Allocator>
  void ScalePolynomVariable(MultivariatePolynomial<T, Allocator>& P, const T& scale);

  template<class T, class Allocator>
  void CompressPolynom(MultivariatePolynomial<T, Allocator>& P, const T& threshold);

  template<class T, class Allocator>
  void OppositeX(MultivariatePolynomial<T, Allocator>& P);

  template<class T, class Allocator>
  void OppositeY(MultivariatePolynomial<T, Allocator>& P);

  template<class T, class Allocator>
  void PermuteXY(MultivariatePolynomial<T, Allocator>& P);
  
  template<class T1, class Allocator1, class Allocator2,class Allocator3, class Allocator4>
  void GenerateSymPol(const MultivariatePolynomial<T1, Allocator1>& P,
		      MultivariatePolynomial<T1, Allocator2>& Q, 
		      MultivariatePolynomial<T1, Allocator3>& R,
		      MultivariatePolynomial<T1, Allocator4>& S);

  template<class T1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GenerateSym4_Vertices(const MultivariatePolynomial<T1, Allocator1>& P,
			     MultivariatePolynomial<T1, Allocator2>& Q, 
			     MultivariatePolynomial<T1, Allocator3>& R,
			     MultivariatePolynomial<T1, Allocator4>& S);
  
  template<class T1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GenerateSym4_Edges(const MultivariatePolynomial<T1, Allocator1>& P,
			  MultivariatePolynomial<T1, Allocator2>& Q, 
			  MultivariatePolynomial<T1, Allocator3>& R,
			  MultivariatePolynomial<T1, Allocator4>& S);

  template<class T1, class Allocator1, class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5, class Allocator6,
	   class Allocator7, class Allocator8>
  void GenerateSym8(const MultivariatePolynomial<T1, Allocator1>& P,
		    MultivariatePolynomial<T1, Allocator2>& Q, 
		    MultivariatePolynomial<T1, Allocator3>& R,
		    MultivariatePolynomial<T1, Allocator4>& S, 
		    MultivariatePolynomial<T1, Allocator5>& T,
		    MultivariatePolynomial<T1, Allocator6>& U,
		    MultivariatePolynomial<T1, Allocator7>& V,
		    MultivariatePolynomial<T1, Allocator8>& W);
  
  template<class T>
  inline void FillZero(MultivariatePolynomial<T>& x){ x.Fill(0); }
  
}

#define MONTJOIE_FILE_MULTIVARIATE_POLYNOMIALHXX
#endif
