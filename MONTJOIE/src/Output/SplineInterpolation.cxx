#ifndef MONTJOIE_FILE_SPLINE_INTERPOLATION_CXX

namespace Montjoie
{

  /***********************
   * SplineInterpolation *
   ***********************/
  
  
  //! default constructor
  template<class T>
  SplineInterpolation<T>::SplineInterpolation()
  {
    threshold = 1e-12;
    first = false;
  }
  
  
  
  //! sets the threshold used to include extremities
  template<class T>
  void SplineInterpolation<T>::SetThreshold(const Real_wp& eps)
  {
    threshold = eps;
  }
  

  //! returns second derivatives
  template<class T>
  const Vector<T>& SplineInterpolation<T>::GetSecondDerivatives() const
  {
    return ypp;
  }
    

  //! returns the memory used by the object in bytes
  template<class T>
  size_t SplineInterpolation<T>::GetMemorySize() const
  {
    size_t taille = x.GetMemorySize() + y.GetMemorySize() + ypp.GetMemorySize() + coef_h.GetMemorySize();
    return taille;
  }


  //! compute tridiagonal matrix needed to find values of y''
  template<class T>
  void SplineInterpolation<T>::ComputeTridiagonalMatrix(const VectReal_wp& h,
							TinyBandMatrix<Real_wp, 1>& A) const
  {
    int N = x.GetM();
    A.Reallocate(N-2, N-2);
    
    for (int i = 0; i < N-2; i++)
      A.Get(i, i) = (h(i)+h(i+1))/3;
    
    for (int i = 0; i < N-3; i++)
      {
	A.Get(i, i+1) = h(i+1)/6;
	A.Get(i+1, i) = h(i+1)/6;
      }
  }
  
  
  //! compute values of y'' from values of x and y
  /*!
    if first_ is true, y'' are set to 0 so that first-order polynomials
    are used instead of third-order polynomials
  */
  template<class T>
  void SplineInterpolation<T>
  ::Init(const Vector<Real_wp>& xdiv, const Vector<T>& ydiv, bool first_)
  {
    if (xdiv.GetM() < 2)
      {
	cout << "Not enough points to perform interpolation" << endl;
	abort();
      }
    
    int N = xdiv.GetM();
    x = xdiv;
    y = ydiv;
    
    // space steps h_i = x_{i+1} - x_i
    VectReal_wp h(N-1); coef_h.Reallocate(N-1);
    for (int i = 0; i < h.GetM(); i++)
      {
	h(i) = x(i+1) - x(i);
	coef_h(i) = square(h(i))/6.0;
      }

    if (xdiv.GetM() == 2)
      first_ = true;
    
    if (first_)
      {
	first = first_;
	return;
      }
    
    // right hand side of linear system
    Vector<T> rhs(N-2);
    for (int i = 0; i < N-2; i++)
      rhs(i) = 1.0/h(i+1)*(y(i+2) - y(i+1)) - 1.0/h(i)*(y(i+1) - y(i));

    // tridiagonal linear system
    TinyBandMatrix<Real_wp, 1> A;
    ComputeTridiagonalMatrix(h, A);

    // solving linear system
    A.Factorize();
    A.Solve(rhs);
    
    // then storing ypp
    ypp.Reallocate(N);
    SetComplexZero(ypp(0));
    SetComplexZero(ypp(N-1));
    for (int i = 0; i < N-2; i++)
      ypp(i+1) = rhs(i);

    threshold = 100.0*epsilon_machine*(x(N-1) - x(0));
  }
  
  
  //! finds the integers a, b such that xinterp belongs to [x_a, x_b]
  template<class T>
  void SplineInterpolation<T>::FindInterval(const Real_wp& xinterp, int& a, int& b) const
  {
    if (x.GetM() < 2)
      {
	cout << "Not enough points to perform interpolation" << endl;
	abort();
      }
    
    // searching interval with bisection procedure
    a = 0; b = x.GetM()-1;
    if ((xinterp < x(a)-threshold) || (xinterp > x(b)+threshold))
      {
	DISP(xinterp); DISP(x(a)); DISP(x(b));
	cout << "Point outside interval range" << endl;
	abort();
      }
    
    while (b > a+1)
      {
	int c = (a+b)/2;
	if (xinterp > x(c))
	  a = c;
	else
	  b = c;
      }      
  }
  
  
  //! computes interpolation on a given point xinterp
  template<class T>
  T SplineInterpolation<T>::Evaluate(const Real_wp& xinterp) const
  {
    // we find interval [x_a, x_b] where the point is
    int a, b;
    FindInterval(xinterp, a, b);
    
    // then using ypp for a cubic interpolation
    Real_wp zeta = (xinterp-x(a)) / (x(b) - x(a));
    T yinterp;
    if (first)
      yinterp = (1.0-zeta)*y(a) + zeta*y(b);
    else
      yinterp = (1.0-zeta)*y(a) + zeta*y(b)
	+ coef_h(a)*( (1.0-zeta)*(square(1.0-zeta)-1.0) * ypp(a)
		      + zeta*(zeta*zeta - 1.0) * ypp(b));
    
    return yinterp;
  }
  
  
  //! computes interpolation on a given point xinterp
  //! and provides first and second derivative
  template<class T>
  void SplineInterpolation<T>::EvaluateDerivative(const Real_wp& xinterp, T& yinterp,
						  T& dy_interp, T& d2y_interp) const
  {
    // we find interval [x_a, x_b] where the point is
    int a, b;
    FindInterval(xinterp, a, b);
    
    // then using ypp for a cubic interpolation
    Real_wp invH = 1.0/ (x(b) - x(a));
    Real_wp zeta = (xinterp-x(a))*invH;
    yinterp = (1.0-zeta)*y(a) + zeta*y(b);
    
    dy_interp = invH*(y(b) - y(a));
    d2y_interp = 0;
    
    if (!first)
      {
	yinterp += coef_h(a)*( (1.0-zeta)*(square(1.0-zeta)-1.0) * ypp(a)
			       + zeta*(zeta*zeta - 1.0) * ypp(b));
	
	dy_interp += invH*(coef_h(a)*((1.0-3.0*square(1.0-zeta))*ypp(a)
				      + (3.0*zeta*zeta-1.0)*ypp(b)));
	
	d2y_interp = (1.0-zeta)*ypp(a) + zeta*ypp(b);
      }
  }
  
  
  //! evalues basis functions of the spline
  /*!
    This function should not be used because spline are non-local
    so that vec_phi is always dense
  */
  template<class T>
  void SplineInterpolation<T>::EvaluateFctBasis(const Real_wp& point, VectReal_wp & vec_phi) const
  {
    int a, b;
    FindInterval(point, a, b);
    vec_phi.Fill(0.0);
    Real_wp invH = 1.0/ (x(b) - x(a));
    Real_wp zeta = (point-x(a))*invH;
    vec_phi(a) = (1-zeta);
    vec_phi(b) = zeta;
    
    if (!first)
      {
	int N = x.GetM();
	VectReal_wp h(N-1);
	for (int i = 0; i < h.GetM(); i++)
	  h(i) = x(i+1) - x(i);
	
	TinyBandMatrix<Real_wp, 1> A;
	ComputeTridiagonalMatrix(h, A);
	
	A.Factorize();
	
	VectReal_wp ones(N-2);
	
	if ((a >0) && (a < N-1))
	  {
	    // we retrieve the row  a-1 of A^{-1} B
	    Vector<Real_wp> rowA(N);
	    ones.Fill(0); ones(a-1) = 1.0;
	    A.Solve(ones);
	    rowA.Fill(0);
	    for (int i = 0; i < N-2; i++)
	      {
		rowA(i) += ones(i)/h(i);
		rowA(i+1) -= ones(i)*(1.0/h(i) + 1.0/h(i+1));
		rowA(i+2) += ones(i)/h(i+1);
	      }
	    
	    // then vec_phi is modified
	    Real_wp coef = coef_h(a)*(1.0-zeta)*(square(1.0-zeta)-1.0);
	    for (int j = 0; j < N; j++)
	      vec_phi(j) += coef*rowA(j);
	  }
	
	if ((b >0) && (b < N-1))
	  {
	    // we retrieve the row  b-1 of A^{-1} B
	    Vector<Real_wp> rowB(N);
	    ones.Fill(0); ones(b-1) = 1.0;
	    A.Solve(ones);
	    rowB.Fill(0);
	    for (int i = 0; i < N-2; i++)
	      {
		rowB(i) += ones(i)/h(i);
		rowB(i+1) -= ones(i)*(1.0/h(i) + 1.0/h(i+1));
		rowB(i+2) += ones(i)/h(i+1);
	      }
	    
	    // then vec_phi is modified
	    Real_wp coef = coef_h(a)*zeta*(zeta*zeta - 1.0);
	    for (int j = 0; j < N; j++)
	      vec_phi(j) += coef*rowB(j);
	  }
	
      }
  }
  

  //! retrieves cubic polynomials on each interval
  /*template<class T>
  void SplineInterpolation<T>::FillPolynomialBasis(Vector<UnivariatePolynomial<T> >& P)
  {    
    P.Reallocate(x.GetM()-1);
    UnivariatePolynomial<T> one;
    one.SetOrder(0); SetComplexOne(one(0));
    for (int i = 0; i < x.GetM()-1; i++)
      {
        UnivariatePolynomial<T> zeta;
        zeta.SetOrder(1);
        Real_wp invHi = 1.0/(x(i+1) - x(i));
        zeta(0) = -x(i) * invHi;
        zeta(1) = invHi;

        P(i) = (one-zeta)*y(i) + zeta*y(i+1) + coef_h(i)*((one-zeta)*(square(one-zeta)-one)*ypp(i)
                                                          + zeta*(zeta*zeta - one)*ypp(i+1));
      }
      }*/
  
  
  //! multiplies values of y and y'' by a same coefficient
  template<class T>
  void SplineInterpolation<T>::Mlt(const T& coef)
  {
    Seldon::Mlt(coef, y); 
    Seldon::Mlt(coef, ypp);
  }
  
  
  //! releases memory
  template<class T>
  void SplineInterpolation<T>::Clear()
  {
    x.Clear(); y.Clear(); ypp.Clear(); coef_h.Clear();
  }


  /*****************************
   * BSplineInterpolation_Base *
   *****************************/

  
  //! default constructor
  BSplineInterpolation_Base::BSplineInterpolation_Base()
  {
    order = 1;
    nb_basis_fct = 0;    
    threshold = 1e-12;
  }
  

  //! returns the order of approximation
  int BSplineInterpolation_Base::GetOrder() const
  {
    return order;
  }


  //! returns the number of basis functions
  int BSplineInterpolation_Base::GetNbDof() const
  {
    return nb_basis_fct;
  }

    
  //! returns the vector containing the final knots
  const VectReal_wp& BSplineInterpolation_Base::GetPoints() const
  {
    return t;
  }

  
  //! sets the threshold used to include extremities
  void BSplineInterpolation_Base::SetThreshold(const Real_wp& eps)
  {
    threshold = eps;
  }


  //! the point xinterp is found in the interval [t_a, t_b]
  void BSplineInterpolation_Base::FindInterval(const Real_wp& xinterp, int& a, int& b) const
  {
    if (t.GetM() < 2)
      {
	cout << "Not enough points to perform bisection algorithm" << endl;
	abort();
      }
    
    // searching interval with bisection procedure
    a = 0; b = t.GetM()-1;
    if ((xinterp < t(a)-threshold) || (xinterp > t(b)+threshold))
      {
	DISP(xinterp); DISP(t(a)); DISP(t(b));
	cout << "Point outside interval range" << endl;
	abort();
      }
    
    while (b > a+1)
      {
	int c = (a+b)/2;
	if (xinterp > t(c))
	  a = c;
	else
	  b = c;
      }      
  }
  
  
  //! evaluates basis functions at a given point
  int BSplineInterpolation_Base::EvaluateFctBasis(const Real_wp& point, VectReal_wp& vec_phi, bool init_vec) const
  {
    if (t.GetM() <= 1)
      {
	vec_phi.Clear();
	return 0;
      }

    // we find the interval where the point is
    int pos, b;
    FindInterval(point, pos, b);

    // skipping the first knots (that are added to have a clamped extremity)
    if (pos <= order)
      {
	pos = order;
	b = pos+1;
      }
    
    // naive algorithm : we start from (0, 0, ..., 1, 0, ..., 0)
    // and use the recurrence formula    
    int N = t.GetM() - 2*order;
    if (init_vec)
      {
	vec_phi.Reallocate(this->nb_basis_fct);
	vec_phi.Zero();
      }
    
    // only accessed values are initialized
    for (int i = pos-order; i < pos; i++)
      vec_phi(i) = 0;
    
    if (pos < N+order-2)
      vec_phi(pos+1) = 0;
    
    vec_phi(pos) = Real_wp(1);
        
    for (int k = 1; k <= order; k++)
      {
	// only indexes for which vec_phi is different from 0 are considered
	for (int i = pos-k; i <= pos; i++)
	  {
	    // B_{i, k} = (x - t_i) / (t_{i+k} - t_i) B_{i, k-1} + (t_{i+k+1} - x) / (t_{i+k+1} - t_{i+1}) B_{i+1, k-1}
	    // when a knot is multiple, we use the convention 0/0 = 0
	    if (t(i) != t(i+k))
	      vec_phi(i) *= (point - t(i)) / (t(i+k) - t(i));
	    else
	      vec_phi(i) = Real_wp(0);
	    
	    if (t(i+1) != t(i+k+1))
	      vec_phi(i) += (t(i+k+1) - point) / (t(i+k+1) - t(i+1)) * vec_phi(i+1);
	  }
      }
    
    return pos;
  }


  //! evaluates basis functions and derivatives at a given point
  int BSplineInterpolation_Base::EvaluateDeriveFctBasis(const Real_wp& point, VectReal_wp& vec_phi,
							VectReal_wp& der_phi, bool init_vec) const
  {
    if (t.GetM() <= 1)
      {
	der_phi.Clear();
	vec_phi.Clear();
	return 0;
      }

    // we find the interval where the point is
    int pos, b;
    FindInterval(point, pos, b);

    // skipping the first/last knots (that are added to have a clamped extremity)
    if (pos <= order)
      {
	pos = order;
	b = pos+1;
      }

    int N = t.GetM() - 2*order;
    if (pos >= N+order-2)
      {
	pos = N+order-2;
	b = pos+1;
      }
    
    // naive algorithm : we start from (0, 0, ..., 1, 0, ..., 0)
    // and use the recurrence formula    
    if (init_vec)
      {
	vec_phi.Reallocate(this->nb_basis_fct);
	der_phi.Reallocate(this->nb_basis_fct);
	vec_phi.Zero();
	der_phi.Zero();
      }
    
    // only accessed values are initialized
    for (int i = pos-order; i < pos; i++)
      vec_phi(i) = 0;
    
    if (pos < N+order-2)
      vec_phi(pos+1) = 0;
    
    vec_phi(pos) = Real_wp(1);
        
    for (int k = 1; k <= order; k++)
      {
	if (k == order)
	  {
	    // derivatives of basis functions can be computed here
	    for (int i = pos-k; i <= pos; i++)
	      {
		if (t(i) != t(i+k))
		  der_phi(i) = Real_wp(k) / (t(i+k) - t(i)) * vec_phi(i);
		else
		  der_phi(i) = Real_wp(0);

		if (t(i+1) != t(i+k+1))
		  der_phi(i) -= Real_wp(k) / (t(i+k+1) - t(i+1)) * vec_phi(i+1);
	      }
	  }

	// only indexes for which vec_phi is different from 0 are considered
	for (int i = pos-k; i <= pos; i++)
	  {
	    // B_{i, k} = (x - t_i) / (t_{i+k} - t_i) B_{i, k-1} + (t_{i+k+1} - x) / (t_{i+k+1} - t_{i+1}) B_{i+1, k-1}
	    // when a knot is multiple, we use the convention 0/0 = 0
	    if (t(i) != t(i+k))
	      vec_phi(i) *= (point - t(i)) / (t(i+k) - t(i));
	    else
	      vec_phi(i) = Real_wp(0);
	    
	    if (t(i+1) != t(i+k+1))
	      vec_phi(i) += (t(i+k+1) - point) / (t(i+k+1) - t(i+1)) * vec_phi(i+1);
	  }

      }

    return pos;
  }
  
  
  //! Constructs QR factorisation of VDM matrix needed to compute the weights
  void BSplineInterpolation_Base::Init(const Vector<Real_wp>& xdiv, int r,
				       const Vector<Real_wp>& knots, Matrix<Real_wp>& VDM, VectReal_wp& tau)
  {
    order = r;
    
    // we assume here that knots are already sorted
    // extremities are clamped
    int N = knots.GetM();
    t.Reallocate(knots.GetM() + 2*order);
    for (int i = 0; i < order; i++)
      {
	t(i) = knots(0);
	t(N+order+i) = knots(N-1);
      }
    
    for (int i = 0; i < N; i++)
      t(i+order) = knots(i);

    nb_basis_fct = N+order-1;      
    
    // basis functions are evaluated for given points in xdiv
    value_phi.Reallocate(nb_basis_fct);
    value_dphi.Reallocate(nb_basis_fct);
    value_dphi.Zero();
    VDM.Reallocate(xdiv.GetM(), nb_basis_fct);
    for (int i = 0; i < xdiv.GetM(); i++)
      {
	Real_wp points = xdiv(i);
	EvaluateFctBasis(points, value_phi);
	for (int j = 0; j < nb_basis_fct; j++)
	  VDM(i, j) = value_phi(j);
      }

    // Least-squares problem is solved with QR factorisation
    GetQR(VDM, tau);    
  }


  /************************
   * BSplineInterpolation *
   ************************/
  

  //! Returns the weights associated with the basis functions
  template<class T>
  const Vector<T>& BSplineInterpolation<T>::GetWeights() const
  {
    return weights;
  }
    

  //! Returns the size in bytes used by the object
  template<class T>
  size_t BSplineInterpolation<T>::GetMemorySize() const
  {
    return weights.GetMemorySize() + this->value_phi.GetMemorySize() +
      this->value_dphi.GetMemorySize() + this->t.GetMemorySize();
  }
  
  
  //! Initializes the spline with given points (xdiv, ydiv)
  /*
    This method performs a least-squares fitting of a B-spline at the points (xdiv, ydiv)
    The order of approximation is given by r, and the knots are the desired knots given by the user
   */
  template<class T>
  void BSplineInterpolation<T>::Init(const Vector<Real_wp>& xdiv, const Vector<T>& ydiv, int r,
				     const Vector<Real_wp>& knots)
  {
    Matrix<Real_wp> VDM;
    VectReal_wp tau;
    BSplineInterpolation_Base::Init(xdiv, r, knots, VDM, tau);

    weights = ydiv;
    SolveQR(VDM, tau, weights);
  }

  
  //! Initializes the B-spline directly with knots and weights
  template<class T>
  void BSplineInterpolation<T>::SetData(const VectReal_wp& knots, const Vector<T>& weights_, int r)
  {
    this->t = knots;   
    this->order = r;
    this->nb_basis_fct = this->t.GetM() - this->order - 1;
    this->value_phi.Reallocate(this->nb_basis_fct);
    this->value_dphi.Reallocate(this->nb_basis_fct);
    this->value_phi.Zero();
    this->value_dphi.Zero();

    // last values of weights_ are ignored
    this->weights.Reallocate(nb_basis_fct);
    for (int i = 0; i < nb_basis_fct; i++)
      this->weights(i) = weights_(i);
  }

  
  //! Evaluates the spline at a given point
  template<class T>
  T BSplineInterpolation<T>::Evaluate(const Real_wp& point) const
  {
    // solution is interpolated on a single point
    
    // basis functions are evaluated for this point
    int pos = this->EvaluateFctBasis(point, this->value_phi, false);
    
    // yinterp is a combination with weights
    T yinterp;
    SetComplexZero(yinterp);

    // the combination is performed with only non-null values
    for (int j = pos-this->order; j <= pos; j++)
      yinterp += this->value_phi(j)*weights(j);
    
    return yinterp;
  }
  
  
  //! Evaluates the spline and its derivative at a given point
  template<class T>
  void BSplineInterpolation<T>::EvaluateDerivative(const Real_wp& point, T& y, T& dy) const
  {
    // basis functions and derivatives are evaluated for this point
    int pos = this->EvaluateDeriveFctBasis(point, this->value_phi, this->value_dphi, false);
    
    // y and dy are a combination with weights
    SetComplexZero(y);
    SetComplexZero(dy);
    // the combination is performed with only non-null values
    for (int j = pos-this->order; j <= pos; j++)
      {
	y += this->value_phi(j)*weights(j);
	dy += this->value_dphi(j)*weights(j);
      }
  }
    
  
  //! weights are multiplied by coef
  template<class T>
  void BSplineInterpolation<T>::Mlt(const T& coef)
  {
    weights *= coef;
  }
  

  //! releases memory used by the spline
  template<class T>
  void BSplineInterpolation<T>::Clear()
  {
    weights.Clear();
    this->t.Clear();
    this->value_phi.Clear();
    this->value_dphi.Clear();
    order = 1;
    nb_basis_fct = 0;
  }
  
}

#define MONTJOIE_FILE_SPLINE_INTERPOLATION_CXX
#endif
