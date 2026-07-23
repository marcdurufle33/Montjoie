#ifndef MONTJOIE_FILE_MULTIVARIATE_POLYNOMIAL_CXX

namespace Montjoie
{
  
  //! default constructor, empty polynom
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>::MultivariatePolynomial()
    : Vector<T, Vect_Full, Allocator>()
  {
    dimension = 0; n_ = 0;
    NumTenso.Clear(); InvNumTenso.Clear();
  }
  
  
  //! polynome of order n
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>::MultivariatePolynomial(unsigned dimension_, int n)
    : Vector<T, Vect_Full, Allocator>()
  {
    // number of monomes is equal to C_{n+d}^n where d is the dimension
    dimension = dimension_; int nb_coef = 1; n_ = n;
    for (unsigned i = 1; i <= dimension ; i++)
      nb_coef = nb_coef*(n+i)/i;
    
    this->Reallocate(nb_coef); 
    // we set null coefficients
    T zero(0);
    this->Fill(zero);
    ConstructIndexArrays();
  }
  
  
  //! constructor by copy
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>::
  MultivariatePolynomial(const MultivariatePolynomial<T, Allocator>& P)
    : Vector<T, Vect_Full, Allocator>(P)
  {
    *this = P;
  }
  
  
  //! default destructor
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>::~MultivariatePolynomial() {}
  
  
  /*******************
   * Basic Functions *
   *******************/
  
  
  //! returns the number of coefficients
  template<class T, class Allocator>
  int MultivariatePolynomial<T, Allocator>::GetM() const
  {
    return this->m_;
  }
  
  
  template<class T, class Allocator>
  size_t MultivariatePolynomial<T, Allocator>::GetMemorySize() const
  {
    size_t taille = Vector<T, VectFull, Allocator>::GetMemorySize();
    taille += NumTenso.GetMemorySize() + InvNumTenso.GetMemorySize()+2*sizeof(int);
    return taille;
  }


  //! returns the dimension (number of variables)
  template<class T, class Allocator>
  unsigned int MultivariatePolynomial<T, Allocator>::GetDimension() const
  {
    return dimension;
  }
  
  
  //! returns the order of the polynom
  template<class T, class Allocator>
  int MultivariatePolynomial<T, Allocator>::GetOrder() const
  {
    return n_;
  }
  
  
  //! change the order of the polynom and dimension
  /*
    \warning the previous data is removed
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::SetOrder(unsigned d, int order)
  {
    if ((order!=n_)||(dimension!=d))
      {
	// number of variables
	dimension = d;
	// number of monomes is equal to C_{n+d}^n where d is the dimension
	int nb_coef = 1; n_ = order;
	for (unsigned i = 1; i <= dimension; i++)
	  nb_coef = nb_coef*(n_+i)/i;
	
	this->Reallocate(nb_coef);
	T zero(0); this->Fill(zero);
	
	ConstructIndexArrays();
      }
  }
  
  
  //! changes the order of the polynom and dimension
  /*!
    keeping old coefficients
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::ResizeOrder(unsigned d, int order)
  {
    if ((order!=n_)||(dimension!=d))
      {
	MultivariatePolynomial<T, Allocator> old_polynom(*this);
	// number of variables
	dimension = d;
	// number of monomes is equal to C_{n+d}^n where d is the dimension
	int nb_coef = 1; n_ = order;
	for (unsigned i = 1;  i <= dimension; i++)
	  nb_coef = (nb_coef*(n_+i))/i;
	
	this->Reallocate(nb_coef);
	T zero(0); this->Fill(zero);
	
	ConstructIndexArrays();
	
	IVect old_index(old_polynom.GetDimension()); old_index.Zero();
	IVect index(dimension); index.Zero();
	int D = min(old_polynom.GetDimension(), dimension);
	for (int i = 0; i < old_polynom.GetM(); i++)
	  {
	    old_polynom.IndexToPowers(i, old_index);
	    bool index_to_copy = true;
	    for (int k = dimension; k < old_index.GetM(); k++)
	      if (old_index(k) > 0)
		index_to_copy = false;
	    
	    for (int j = 0; j < D; j++)
	      if (old_index(j) > order)
		index_to_copy = false;
	    
	    if (index_to_copy)
	      {
		for (int j = 0; j < D; j++)
		  index(j) = old_index(j);
		
		int k = this->PowersToIndex(index);
		if ( k >= 0)
		  (*this)(k) = old_polynom(i);
	      }
	  }
	
      }
  }
  
  
  //! copy the polynome Q
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>
  ::Copy(const MultivariatePolynomial<T, Allocator>& Q)
  {
    Vector<T,Vect_Full,Allocator>::Copy(Q);
    n_ = Q.n_; NumTenso = Q.NumTenso; InvNumTenso = Q.InvNumTenso;
    dimension = Q.dimension;
  }
  
  
  /**********************
   * Internal Functions *
   **********************/
  
  
  //! internal method to reconstruct indexes
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::ConstructIndexArrays()
  {
    IVect sum_powers(dimension+1), powers(dimension+1); 
    powers.Zero(); sum_powers.Zero();
    
    NumTenso.Reallocate(this->m_);
    int nb_tenso = 1;
    for (unsigned j= 0; j < dimension; j++)
      nb_tenso *= n_+1;
    
    InvNumTenso.Reallocate(nb_tenso); InvNumTenso.Fill(-1);
    
    for (size_t i = 0; i < this->m_; i++)
      {
	int num_tenso = 0;
	int offset = 1;
	for (unsigned j = 0; j < dimension; j++)
	  {
	    num_tenso += offset*powers(j);
	    offset *= n_+1;
	  }
	
	NumTenso(i) = num_tenso;
	InvNumTenso(num_tenso) = i;
	
	// we update powers indices
	for (int j = (dimension-1); j >= 0; j--)
	  {
	    if (sum_powers(j)+powers(j) < n_)
	      {
		// we can still increase the index
		powers(j) ++;
		// we update sum 
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_powers(p) ++;
		break;
	      }
	    else
	      {
		// we couldn't increase the index of variable j
		// so we zero this index
		// and try to increase variable j-1 
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_powers(p) -= powers(j);
		powers(j) = 0;
	      }
	  } // end loop on indices
      }
    
    // DISP(NumTenso); DISP(InvNumTenso);
  }
  
  
  //! use of the Horner Algorithm
  /*!
    we substitute the last variable x_d by value
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::SubstituteLastVariable(const T& value)
  {
    dimension--;
    int nb_coef = 1;
    for (unsigned i = 1; i <= dimension; i++)
      nb_coef = nb_coef*(n_+i)/i;
    
    Vector<T, Vect_Full, Allocator> res(nb_coef);
    // nb_coef coefficients will be computed
    // indexes of the (dim-1) first variables
    IVect indexes(dimension+1),sum_indexes(dimension+1); 
    indexes.Zero(); sum_indexes.Zero();
    int degree_last_variable = 0;
    int k = 0; T coef; int index_data;
    // main loop
    for (int i = 0; i < nb_coef; i++)
      {
	degree_last_variable = n_ - sum_indexes(dimension);
	
	// now the following polynom is evaluated
	//  a_k + a_{k+1} x + ... a_{d} x^{d}
	// where d is degree_last_variable = n-\sum i_k 
	// with (i_k) the indices of the first (dim-1) variables
	index_data = k + degree_last_variable;
	coef = this->data_[index_data];
	int l = degree_last_variable;
	while (l > 0)
	  {
	    coef = coef*value + this->data_[--index_data];
	    l--;
	  }
	
	// we fill res
	res(i) = coef; k += degree_last_variable+1;
	
	// increment of indexes
	for (int j = (dimension-1); j >= 0; j--)
	  {
	    if (sum_indexes(j)+indexes(j) < n_)
	      {
		// we can still increase the index
		indexes(j) ++;
		// we update sum 
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_indexes(p) ++;
		break;
	      }
	    else
	      {
		// we couldn't increase the index of j-ieme variable
		// so we zero this index, and look for increasing (j-1) i�me variable 
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_indexes(p) -= indexes(j);
		indexes(j) = 0;
	      }
	  }
      }
    
    // we replace data_ by res.data_
    Vector<T, Vect_Full, Allocator>::Copy(res);
  }
  
  
  //! substitution of any variable
  /* 
     \warning the first variable is variable 0 (C-convention)
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::SubstituteVariable(int num_variable,
                                                                const T& value)
  {
    if (num_variable == int(dimension)-1)
      SubstituteLastVariable(value);
    
    // the final result
    MultivariatePolynomial res(dimension-1,n_);
    // indices
    IVect indexes(dimension+1),sum_indexes(dimension+1); 
    IVect new_indexes(dimension);
    // powers of value (1, x ... x^n_)
    Vector<T, Vect_Full, Allocator> powers_value(n_+1);
    powers_value(0) = 1.0;
    for (int i = 0; i < n_; i++)
      powers_value(i+1) = powers_value(i)*value;
    
    indexes.Zero(); sum_indexes.Zero();
    // main loop
    for (size_t i = 0; i < this->m_; i++)
      {
	// new_indexes are like the array indexes
	// where the num_variable value has been removed
	for (int j = 0; j < num_variable; j++)
	  new_indexes(j) = indexes(j);
	for (unsigned j = (num_variable+1); j <= dimension; j++)
	  new_indexes(j-1) = indexes(j);
	
	int pos = res.PowersToIndex(new_indexes); 
	res(pos) += powers_value(indexes(num_variable))*this->data_[i];
	
	// increment of indexes
	for (int j = (dimension-1); j >= 0; j--)
	  {
	    if (sum_indexes(j)+indexes(j) < n_)
	      {
		// we can still increase the index
		indexes(j) ++;
		// we update sum 
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_indexes(p)++;
		
		break;
	      }
	    else
	      {
		// we couldn't increase the index of variable j
		// we zero this index, and look for increasing variable j-1
		for (unsigned p = (j+1); p <= dimension; p++)
		  sum_indexes(p) -= indexes(j);
		indexes(j) = 0;
	      }
	  }
      }
    
    // we replace data_ by res.data_
    Copy(res);
  }
    
  
  /************************
   * Convenient Functions *
   ************************/
  
  
  //! cheap function, which uses index arrays
  /*!
    \param[in] powers power for each variable
    \return scalar subscript of the corresponding coefficient
  */
  template<class T, class Allocator>
  int MultivariatePolynomial<T, Allocator>::PowersToIndex(const IVect& powers) const
  {
    int num_tenso = 0;
    int N = n_+1, offset = 1;
    for (unsigned i = 0; i < dimension; i++)
      {
	num_tenso += powers(i)*offset;
	offset *= N;
      }
    
    if (num_tenso < InvNumTenso.GetM())
      return InvNumTenso(num_tenso);
    
    return -1;
  }
  
  
  //! cheap function, which uses index arrays
  /*!
    \param[in] pos scalar subscript of a coefficient
    \param[out] powers power for each variable
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::IndexToPowers(int pos, IVect& powers) const
  {
    pos = NumTenso(pos);
    powers.Reallocate(dimension);
    int N = n_+1;
    for (unsigned i = 0; i < dimension; i++)
      {
	powers(i) = pos%N;
	pos = pos/N;
      }
  }
  
  
  //! evaluation of the value of the polynome for x given
  /*! we compute \f$ P(x) = a_0 + a_1 y + ... + a_n y^n
    + x * ( a_{n+1} + a_{n+2} y + ... 
    + a_{2*n-1} y^(n-1) ) + ... + a_{m_} x^n \f$  
    with the Horner algorithm
    \f$ P(x) = a_0 + y*( a_1 + y*( a_2 + y*(a_3 ...+y*a_n) ) )
    + x * ( a_{n+1}+y*( a_{n+2} ...  ) ) ) \f$
  */
  template<class T, class Allocator> template<class TinyVector_>
  T MultivariatePolynomial<T, Allocator>::Evaluate(const TinyVector_& x) const
  {
    int d = dimension;
    MultivariatePolynomial<T, Allocator> Big_poly;
    Big_poly.Copy(*this);
    
    // we substitute each variable by its value
    while (d > 0)
      {
	Big_poly.SubstituteLastVariable(x(d-1));
        d--;
      }
    
    // finally, we have a 0-variable 0-degree polynom !
    return Big_poly(0);
  }
  
  
  //! displays polynome coefficients with rational coefficients
  /*
    \param[in] epsilon precision of the rational approximation
  */
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::PrintRational(const T& epsilon) const
  {
    // For each coefficient C, we look for continued fraction Q,
    // so that abs(Q-C) < epsilon
    // vector<int> decomp_coef;
    // int num, denom; int nb = 0;
    Vector<int> decomp_coef;
    int num, denom; int nb = 0;
    bool first_writing = true;
    for (size_t i = 0; i < this->m_; i++)
      {
	if (nb%10==1)
	  cout<<endl;
	
	if ((*this)(i) != T(0))
	  {
	    if (!first_writing)
	      cout<<" + ";
	    
	    DecomposeContinuedFraction((*this)(i), epsilon, decomp_coef);
	    GetNumeratorDenominator(decomp_coef, num, denom);
	    cout<<num<<"/"<<denom<<" ";
	    PrintMonomial(i, cout);
	    first_writing = false; nb++;
	  }
      }
  }
  
  
  //! displays the monomial associated with coef i
  template<class T, class Allocator>
  void MultivariatePolynomial<T, Allocator>::PrintMonomial(int i, ostream& out) const
  {
    IVect powers;
    IndexToPowers(i, powers);
    string xi, yj, zk;
    if (dimension <= 3)
      {
        if (powers.GetM() > 0)
          {
            if (powers(0)==0)
              xi = string("");
            else if (powers(0)==1)
              xi = string("x ");
            else
              xi = string("x^")+to_str(powers(0))+" ";
            
            if (dimension >= 2)
              {
                if (powers(1)==0)
                  yj = string("");
                else if (powers(1)==1)
                  yj = string("y ");
                else
                  yj = string("y^")+to_str(powers(1))+" ";
              }
            
            if (dimension >= 3)
              {
                if (powers(2)==0)
                  zk = string("");
                else if (powers(2)==1)
                  zk = string("z");
                else
                  zk = string("z^")+to_str(powers(2));
              }
            out<<xi<<yj<<zk;
          }
      }
    else
      {
	for (unsigned j = 0; j < dimension; j++)
	  {
	    if (powers(j) == 1)
	      out<<"x"<<j<<" ";
	    else if (powers(j) > 1)
	      out<<"x"<<j<<"^"<<powers(j)<<" ";
	  }
      }
  }
  
  
  /*************
   * Operators *
   *************/
  
  
  //! Element access and affectation.
  template<class T, class Allocator>
  T& MultivariatePolynomial<T, Allocator>::operator() (int i)
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "MultivariatePolynomial");
#endif
    
    return this->data_[i];
  }
  
  
  //! Element access 
  template<class T, class Allocator>
  const T& MultivariatePolynomial<T, Allocator>::operator() (int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->m_, "MultivariatePolynomial");
#endif
    
    return this->data_[i];
  }
  
  
  //! Element access and affectation for 2-D polynoms
  template<class T, class Allocator>
  T& MultivariatePolynomial<T, Allocator>::operator() (int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if ((i < 0)||(i >= long(this->m_))||(dimension != 2))
      throw WrongIndex("MultivariatePolynomial<Vect_Full>::operator()",
		       "Index should be in [0, "
		       + to_str(this->m_-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    int num_tenso = i + j*(n_+1);
    
    return this->data_[InvNumTenso(num_tenso)];
  }
  
  
  //! Element access for 2-D polynoms
  template<class T, class Allocator>
  const T& MultivariatePolynomial<T, Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if ((i < 0)||(i >= long(this->m_))||(dimension != 2))
      throw WrongIndex("MultivariatePolynomial<Vect_Full>::operator()",
		       "Index should be in [0, "
		       + to_str(this->m_-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    int num_tenso = i + j*(n_+1);
    
    return this->data_[InvNumTenso(num_tenso)];
  }
  
  
  //! Element access and affectation for 3-D polynoms
  template<class T, class Allocator>
  T& MultivariatePolynomial<T, Allocator>::operator() (int i, int j, int k)
  {
#ifdef SELDON_CHECK_BOUNDS
    if ((i < 0)||(i >= long(this->m_))||(dimension != 3))
      throw WrongIndex("MultivariatePolynomial<Vect_Full>::operator()",
		       "Index should be in [0, "
		       + to_str(this->m_-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    int num_tenso = i + (n_+1)*(j + k*(n_+1));
    
    return this->data_[InvNumTenso(num_tenso)];
  }
  
  
  //! Element access for 3-D polynoms
  template<class T, class Allocator>
  const T& MultivariatePolynomial<T, Allocator>::operator() (int i, int j, int k) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if ((i < 0)||(i >= long(this->m_))||(dimension != 3))
      throw WrongIndex("MultivariatePolynomial<Vect_Full>::operator()",
		       "Index should be in [0, "
		       + to_str(this->m_-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    int num_tenso = i + (n_+1)*(j + k*(n_+1));
    
    return this->data_[InvNumTenso(num_tenso)];
  }
  
  
  //! equality between two polynoms
  template<class T, class Allocator> MultivariatePolynomial<T, Allocator>&
  MultivariatePolynomial<T, Allocator>::
  operator =(const MultivariatePolynomial<T, Allocator>& P)
  {
    Vector<T, Vect_Full, Allocator>::Copy(static_cast<const Vector<T,
					  Vect_Full, Allocator>& >(P));
    dimension = P.dimension; 
    n_ = P.n_;
    NumTenso = P.NumTenso; InvNumTenso = P.InvNumTenso;
    return *this;
  }
  
  
  //! addition of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>
  MultivariatePolynomial<T, Allocator>::
  operator +(const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
    Add(T(1), P, Q);
    return Q;
  }
  
  
  //! difference of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>
  MultivariatePolynomial<T, Allocator>::
  operator -(const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
      Add(T(-1), P, Q);
      return Q;
  }
  
  
  //! multiplication of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator> MultivariatePolynomial<T, Allocator>::
  operator *(const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q;
    MltAdd(T(1), *this, P, T(0), Q);
    return Q;
  }
  
  
  //! addition of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>& MultivariatePolynomial<T, Allocator>::
  operator +=(const MultivariatePolynomial<T, Allocator>& P)
  {
    Add(T(1), P, *this);
    return *this;
  }
  

  //! difference of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>& MultivariatePolynomial<T, Allocator>::
  operator -=(const MultivariatePolynomial<T, Allocator>& P)
  {
    Add(T(-1), P, *this);
    return *this;
  }
  
  
  //! multiplication of two polynoms
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>& MultivariatePolynomial<T, Allocator>::
  operator *=(const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
    MltAdd(T(1), Q, P, T(0), *this);
    return *this;
  }
  
  
  //! multiplication of a polynom by a scalar
  template<class T, class Allocator> template<class T0>
  MultivariatePolynomial<T, Allocator> MultivariatePolynomial<T, Allocator>::
  operator *(const T0& alpha)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
    Mlt(alpha, Q);
    return Q;
  }
  
  
  //! we add scalar to a polynom
  template<class T, class Allocator>
  template<class T0> MultivariatePolynomial<T, Allocator>
  MultivariatePolynomial<T, Allocator>::operator +(const T0& alpha)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
    Q(0) += alpha;
    return Q;
  }
  
  
  //! we subtract scalar to a polynom
  template<class T, class Allocator>
  template<class T0> MultivariatePolynomial<T, Allocator>
  MultivariatePolynomial<T, Allocator>::operator -(const T0& alpha)
  {
    MultivariatePolynomial<T, Allocator> Q(*this);
    Q(0) -= alpha;
    return Q;
  }
  
  
  //! multiplication of a polynom by a scalar
  template<class T, class Allocator>
  template<class T0> MultivariatePolynomial<T, Allocator>&
  MultivariatePolynomial<T, Allocator>::operator *=(const T0& alpha)
  {
    Mlt(alpha, *this);
    return *this;
  }
  
  
  //! we add scalar to a polynom
  template<class T, class Allocator>
  template<class T0> MultivariatePolynomial<T, Allocator>&
  MultivariatePolynomial<T, Allocator>::operator +=(const T0& alpha)
  {
    (*this)(0) += alpha;
    return *this;
  }
  
  
  //! we subtract scalar to a polynom
  template<class T, class Allocator>
  template<class T0> MultivariatePolynomial<T, Allocator>&
  MultivariatePolynomial<T, Allocator>::operator -=(const T0& alpha)
  {
    (*this)(0) -= alpha;
    return *this;
  }
  
  
  //! comparison of two polynomials
  template<class T, class Allocator>
  bool MultivariatePolynomial<T, Allocator>::
  operator ==(const MultivariatePolynomial<T, Allocator>& P)
  {
    int m = P.GetOrder();
    bool eq = true;
    for (int i = 0; i <= m; i++)
      for (int j = 0; j <= m-i; j++)
        for (int k = 0; k <= m-i-j; k++)
          eq = (eq && (abs(P(i,j,k)-(*this)(i,j,k))<10E-9));

    return eq;
  }

  
  //! comparison of two polynomials
  template<class T, class Allocator>
  bool MultivariatePolynomial<T, Allocator>::
  operator !=(const MultivariatePolynomial<T, Allocator>& P)
  {
    int m = P.GetOrder();
    bool eq = true;
    for (int i = 0; i <= m; i++)
      for (int j = 0; j <= m-i; j++)
        for (int k = 0; k <= m-i-j; k++)
          eq = (eq && (abs(P(i,j,k)-(*this)(i,j,k))<10E-9));

    return !eq;
  }

  
  //! multiplication between a scalar and a polynom
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator *(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q = P;
    Mlt(alpha, Q); return Q;
  }
  
  
  //! addition of a scalar and a polynom
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator +(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q = P;
    Q(0) += alpha; return Q;
  }
  
  
  //! difference between a scalar and a polynom
  template<class T0, class T, class Allocator> MultivariatePolynomial<T, Allocator>
  operator -(const T0& alpha, const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q = P;
    Q(0) -= alpha; Mlt(T(-1), Q); return Q;
  }


  //! opposite of a polynome
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>
  operator -(const MultivariatePolynomial<T, Allocator>& P)
  {
    MultivariatePolynomial<T, Allocator> Q(P);
    Mlt(T(-1), Q);
    return Q;
  }

  
  //! displays a polynome
  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const MultivariatePolynomial<T, Allocator>& P)
  {
    bool first_writing = true;
    int nb = 0;
    bool IsZero = true;
    for (int i = 0; i < P.GetM(); i++)
      {	
	if (abs(P(i)-T(0))>10e-12)
	  {
	    IsZero = false;
	    if (!first_writing)
	      {
		if (P(i) >= 0)
		  out<<" + ";
		else
		  out <<" - ";
	      }
	    else
	      {
		if (P(i) < 0)
		  out<<"-";
	      }
	    
	    if ((abs(abs(P(i))-T(1))>10e-12)||(i == 0))
	      out <<abs(P(i))<<" ";
	    
	    P.PrintMonomial(i, out);
	    first_writing = false; nb++;
	    if (nb%10 == 0)
	      out<<endl;
	  }
      }

    if (IsZero)
	out<<" 0 ";
    
    if (first_writing)
      out << " 0 ";
    
    return out;
  }
  

  //! power of a polynomial
  template<class T, class Allocator>
  MultivariatePolynomial<T, Allocator>
  Pow(const MultivariatePolynomial<T, Allocator>& U, int r)
  {
    MultivariatePolynomial<T, Allocator> P = U;
    MultivariatePolynomial<T, Allocator> Q;
    if (r == 0)
      {
        Q.SetOrder(1,0);
	Q(0) = 1;
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
  
  
  //! derivation of a polynome
  template<class T, class Allocator1, class Allocator2>
  void DerivatePolynomial(const MultivariatePolynomial<T, Allocator1>& P,
			 MultivariatePolynomial<T, Allocator2>& Q, int num_variable)
  {
    int n = P.GetOrder(); int dimension = P.GetDimension();
    if (n==0)
      {
	Q.SetOrder(P.GetDimension(),0);
	Q(0) = 0.0;
	return;
      }
    else
      {
	Q.SetOrder(P.GetDimension(),n-1);
	IVect indexes(dimension+1),sum_indexes(dimension+1); 
	IVect new_indexes(dimension+1);
	
	indexes.Zero(); sum_indexes.Zero(); new_indexes.Zero();
	// main loop
	for (int i=0 ; i<P.GetM() ; i++)
	  {
	    // new_indexes are like the array indexes
	    // where the num_variable value has been decreased
	    for (int j=0 ; j<num_variable ; j++)
	      new_indexes(j) = indexes(j);
	    for (int j=(num_variable+1) ; j<=dimension ; j++)
	      new_indexes(j) = indexes(j);
	    new_indexes(num_variable) = indexes(num_variable) - 1;
	    
	    if (new_indexes(num_variable)>=0)
	      {
		int pos = Q.PowersToIndex(new_indexes); 
		Q(pos) = indexes(num_variable)*P(i);
	      }
	    // increment of indexes
	    for (int j=(dimension-1) ; j>=0 ; j--)
	      {
		if (sum_indexes(j)+indexes(j) < n)
		  {
		    // we can still increase the index
		    indexes(j) ++;
		    // we update sum 
		    for (int p=(j+1) ; p<=dimension ; p++)
		      sum_indexes(p) ++;
		    break;
		  }
		else
		  {
		    // we couldn't increase the index of j-ieme variable
		    // so we zero this index, and look for increasing j-th variable 
		    for (int p=(j+1) ; p<=dimension ; p++)
		      sum_indexes(p) -= indexes(j);
		    indexes(j) = 0;
		  }
	      }
	  } // end main loop
      }
  }
  
  
  //! returns derivate of a polynom
  template<class T, class Allocator> MultivariatePolynomial<T, Allocator>
  Derivate(const MultivariatePolynomial<T, Allocator>& P, int num_variable)
  {
    MultivariatePolynomial<T, Allocator> Q;
    DerivatePolynomial(P, Q, num_variable);
    return Q;
  }
  
  
  //! computes Q = curl P
  template<class T, class Allocator1, class Allocator2>
  void GetCurlPolynomial(const TinyVector<MultivariatePolynomial<T, Allocator1>, 3>& P,
                         TinyVector<MultivariatePolynomial<T, Allocator2>, 3>& Q)
  {
    MultivariatePolynomial<T, Allocator1> R1, R2;
    DerivatePolynomial(P(2), R1, 1);
    DerivatePolynomial(P(1), R2, 2);
    Q(0) = R1 - R2;
    
    DerivatePolynomial(P(0), R1, 2);
    DerivatePolynomial(P(2), R2, 0);
    Q(1) = R1 - R2;
    
    DerivatePolynomial(P(1), R1, 0);
    DerivatePolynomial(P(0), R2, 1);
    Q(2) = R1 - R2;
  }
  

  //! computes Q = div P
  template<class T, class Allocator1, class Allocator2>
  void GetDivPolynomial(const TinyVector<MultivariatePolynomial<T, Allocator1>, 3>& P,
                        MultivariatePolynomial<T, Allocator2>& Q)
  {
    MultivariatePolynomial<T, Allocator1> R1, R2, R3;
    DerivatePolynomial(P(0), R1, 0);
    DerivatePolynomial(P(1), R2, 1);
    DerivatePolynomial(P(2), R3, 2);
    Q = R1 + R2 + R3;
  }

  
  //! integration of a polynomial
  template<class T, class Allocator1, class Allocator2>
  void IntegratePolynomial(const MultivariatePolynomial<T, Allocator1>& P,
			  MultivariatePolynomial<T, Allocator2>& Q, int num_variable)
  {
    int n = P.GetOrder(); 
    Q.SetOrder(P.GetDimension(),n+1); 
    // not implemented
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  //! R = alpha*R
  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, MultivariatePolynomial<T, Allocator>& R)
  {
    if (alpha == T0(0))
      R.Fill(T0(0));
    else
      for (int i = 0; i < R.GetM(); i++)
	R(i) *= T(alpha);
  }
  
  
  //! R = beta*R + alpha*P*Q
  template<class T0, class T1, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const MultivariatePolynomial<T, Allocator1>& P,
	      const MultivariatePolynomial<T, Allocator2>& Q,
	      const T1& beta, MultivariatePolynomial<T, Allocator3>& R)
  {
    int m = P.GetOrder(), n = Q.GetOrder();
    int new_order = max(R.GetOrder(), m+n);
    unsigned new_dim = max(P.GetDimension(), Q.GetDimension());
    new_dim = max(R.GetDimension(), new_dim);
    R.ResizeOrder(new_dim, new_order);
    Mlt(beta, R);
    IVect index_P, index_Q, index_R(new_dim);
    for (int i = 0; i < P.GetM(); i++)
      {
	P.IndexToPowers(i, index_P);
	for (int j = 0; j < Q.GetM(); j++)
	  {
	    Q.IndexToPowers(j, index_Q);
	    // now we add indices
	    index_R.Zero();
	    for (int n = 0; n < index_P.GetM(); n++)
	      index_R(n) += index_P(n);
	    for (int n = 0; n < index_Q.GetM(); n++)
	      index_R(n) += index_Q(n);
	    
	    // we get position of the monomial
	    int pos = R.PowersToIndex(index_R);
	    
	    // we add contribution
	    R(pos) += T(alpha)*P(i)*Q(j);
	  }
      }
  }
  
  
  //! integral of polynom in unit simplex (triangle in 2-D, tet in 3-D)
  template<class T, class Allocator>
  T IntegrateOnSimplex(const MultivariatePolynomial<T, Allocator>& P)
  {
    T res(0), coef;
    int dimension = P.GetDimension();
    IVect powers(dimension);
    for (int i = 0; i < P.GetM(); i++)
      {
	P.IndexToPowers(i, powers);
	coef = T(1); int factor = powers(0)+1;
	for (int j = 1; j < dimension; j++)
	  {
	    for (int k = 1; k <= powers(j); k++)
	      {
		coef *= T(k)/T(factor);
		factor++;
	      }
	    coef *= T(1)/T(factor); factor++;
	  }
	
	coef *= T(1)/T(factor);
		
	res += P(i)*coef;
      }
    
    return res;
  }
  
  
  //! Q = Q + alpha*P
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const T0& alpha, const MultivariatePolynomial<T1, Allocator1>& P,
	   MultivariatePolynomial<T2, Allocator2>& Q)
  {
    int order = max(P.GetOrder(), Q.GetOrder());
    int dim = max(P.GetDimension(), Q.GetDimension());
    Q.ResizeOrder(dim, order);
    IVect index_P, index_Q(dim); index_Q.Zero();
    for (int i = 0; i < P.GetM(); i++)
      {
	P.IndexToPowers(i, index_P);
	for (int j = 0; j < index_P.GetM(); j++)
	  index_Q(j) = index_P(j);
	
	int pos = Q.PowersToIndex(index_Q);
	Q(pos) += T2(alpha*P(i));
      }
  }
  
  
  //! scaling of all variables x1, x2, etc
  template<class T, class Allocator>
  void ScalePolynomVariable(MultivariatePolynomial<T, Allocator>& P, const T& scale)
  {
    int dim = P.GetDimension();
    IVect index_P(dim);
    for (int i = 0; i < P.GetM(); i++)
      {
	P.IndexToPowers(i, index_P);
	int sum = 0;
	for (int k = 0; k < dim ; k++)
	  sum += index_P(k);
	
	if (sum > 0)
	  P(i) *= pow(scale, T(sum));
      }
  }
  
  
  //! reajusts order and dimension of polynom if null elements are found
  template<class T, class Allocator>
  void CompressPolynom(MultivariatePolynomial<T, Allocator>& P, const T& threshold)
  {
    int dim = P.GetDimension();
    int order = 0;
    IVect index_P(dim);
    for (int i = 0; i < P.GetM(); i++)
      {
	P.IndexToPowers(i, index_P);
	if (abs(P(i)) > threshold)
	  {
	    int sum = 0;
	    for (int k = 0; k < dim ; k++)
	      sum += index_P(k);
	    
	    order = max(sum, order);
	  }
      }
    
    P.ResizeOrder(dim, order);
    for (int i = 0; i < P.GetM(); i++)
      if (abs(P(i)) <= threshold)
	P(i) = T(0);
  }
  
  
  //! changes x into -x
  template<class T, class Allocator>
  void OppositeX(MultivariatePolynomial<T, Allocator>& P)
  {
    int m = P.GetOrder();
    for (int i = 1; i <= m; i+=2)
      for (int j = 0; j <= m-i; j++)
        for (int k = 0; k <= m-i-j; k++)
	  P(i,j,k) = -P(i,j,k); 
  }
  
  
  //! changes y into -y
  template<class T, class Allocator>
  void OppositeY(MultivariatePolynomial<T, Allocator>& P)
  {
      int m = P.GetOrder();
      for (int i = 0; i <= m; i++)
       for (int j = 1; j <= m-i; j+=2)
        for (int k = 0; k <= m-i-j; k++)
	  P(i,j,k) = -P(i,j,k);              
  }
  
  
  //! changes x into y and y into x
  template<class T, class Allocator>
  void PermuteXY(MultivariatePolynomial<T, Allocator>& P)
  {
    int m = P.GetOrder();
    T temp;
    for (int i = 0; i <= m; i++)
      for (int j = i+1; j <= m-i; j++)
        for (int k = 0; k <= m-i-j; k++)
	  {
	    temp = P(i,j,k);
	    P(i,j,k) = P(j,i,k);
	    P(j,i,k) = temp;
	  } 
  }
  
  //! generates four basis functions respecting the symmetry
  //! of the reference pyramid with one function
  template<class T1, class Allocator1, class Allocator2,class Allocator3, class Allocator4>
  void GenerateSymPol(const MultivariatePolynomial<T1, Allocator1>& P,
		     MultivariatePolynomial<T1, Allocator2>& Q, 
		     MultivariatePolynomial<T1, Allocator3>& R,
		     MultivariatePolynomial<T1, Allocator4>& S)
  {
    MultivariatePolynomial<T1> T = P;
    OppositeY(T); PermuteXY(T); Q = T;
    OppositeY(T); PermuteXY(T); R = T;
    OppositeY(T); PermuteXY(T); S = T;
  }
  
  //! generates four basis functions respecting the symmetry
  //! of the reference pyramid with one function for a vertex
  template<class T1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GenerateSym4_Vertices(const MultivariatePolynomial<T1, Allocator1>& P,
		    MultivariatePolynomial<T1, Allocator2>& Q, 
		    MultivariatePolynomial<T1, Allocator3>& R,
		    MultivariatePolynomial<T1, Allocator4>& S)
  {
    MultivariatePolynomial<T1> T = P;
    OppositeX(T); Q = T; OppositeY(T); R = T; OppositeX(T); S = T;
  }
  
  
  //! generates four basis functions respecting the symmetry
  //! of the reference pyramid with one function for a dof on an edge
  template<class T1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GenerateSym4_Edges(const MultivariatePolynomial<T1, Allocator1>& P,
		    MultivariatePolynomial<T1, Allocator2>& Q, 
		    MultivariatePolynomial<T1, Allocator3>& R,
		    MultivariatePolynomial<T1, Allocator4>& S)
  {
    MultivariatePolynomial<T1> T = P;
    OppositeY(T); R = T; PermuteXY(T); Q = T; OppositeX(T); S = T;
  }
  
  
  //! generates eight basis functions respecting the symmetry 
  //! of the reference pyramid with one function
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
		    MultivariatePolynomial<T1, Allocator8>& W)
  {
     MultivariatePolynomial<T1> A = P;
     OppositeX(A); Q = A; OppositeY(A); T = A; OppositeX(A); U = A;
     PermuteXY(A); R = A; 
     OppositeY(A); S = A; OppositeX(A); V = A; OppositeY(A); W = A;
  }

}

#define MONTJOIE_FILE_MULTIVARIATE_POLYNOMIAL_CXX
#endif
