#ifndef MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_CXX

namespace Montjoie
{

  /****************************
   * VariableParameter1D_Base *
   ****************************/


  //! finds the interval where the point xinterp is
  template<class T>
  int VariableParameter1D_Base<T>::FindInterval(const Real_wp& xinterp) const
  {
    // searching interval with bisection procedure
    int a = 0, b = points.GetM()-1;
    if ((xinterp < points(a)) || (xinterp > points(b)))
      {
	DISP(xinterp); DISP(points(a)); DISP(points(b));
	cout << "Point outside interval range" << endl;
	abort();
      }
    
    while (b > a+1)
      {
	int c = (a+b)/2;
	if (xinterp > points(c))
	  a = c;
	else
	  b = c;
      }      
    
    return a;
  }


  //! Specifies the approximation to use for the 1-D parameter
  template<class T>
  VariableParameter1D_Base<T>* VariableParameter1D_Base<T>
  ::GetNewApproximation(const Real_wp& a, const Real_wp& b, const Vector<string>& param)
  {
    if (param.GetM() <= 1)
      {
	cout << "Enter the discretization method for c and number of points" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }

    VariableParameter1D_Base<T>* var;
    
    // regular subdivision
    int nb_intervals = to_num<int>(param(1));
    VectReal_wp points_reg;
    Linspace(a, b, nb_intervals+1, points_reg);
    
    // method to use
    if (param(0) == "Polynomial")
      {
	PolynomialInterpolationIndex<T>* pol_approx;
	pol_approx = new PolynomialInterpolationIndex<T>();
	pol_approx->SetApproximation(param, points_reg);
	var = pol_approx;
      }
    else if (param(0) == "PiecewiseLinear")
      {
	SplineInterpolationIndex<T>* spline_approx;
	spline_approx = new SplineInterpolationIndex<T>();
	spline_approx->SetApproximation(points_reg, 1);
	var = spline_approx;
      }
    else if (param(0) == "CubicSpline")
      {
	SplineInterpolationIndex<T>* spline_approx;
	spline_approx = new SplineInterpolationIndex<T>();
	spline_approx->SetApproximation(points_reg, 3);
	var = spline_approx;
      }   
    else if (param(0) == "BSpline")
      {
	BSplineInterpolationIndex<T>* spline_approx;
	spline_approx = new BSplineInterpolationIndex<T>();
	spline_approx->SetApproximation(param, points_reg);
	var = spline_approx;
      }
    else
      {
	cout << "Unknown approximation " << param(0) << endl;
	abort();
      }

    return var;
  }


  //! evaluates f, df and d^2 f
  template<class T>
  void VariableParameter1D_Base<T>::EvaluateSecondDerivative(const Real_wp& x, T&, T&, T&) const
  {
    cout << "Second derivative not implemented for this index" << endl;
    abort();
  }
  

  //! Projects a function f to the approximation space
  template<class T>
  void VariableParameter1D_Base<T>
  ::ProjectFunction(T (&f)(Real_wp), Vector<T>& projF)
  {
    // default discretization is assumed to be interpolatory
    projF.Reallocate(nb_points_grid_c);
    for (int i = 0; i < all_points.GetM(); i++)
      projF(i) = f(all_points(i));
  }
  

  /****************************
   * SplineInterpolationIndex *
   ****************************/
  
  
  //! Inits the index with values written on a file
  template<class T>
  void SplineInterpolationIndex<T>::Init(const string& data_file, T& cte)
  {
    // values are read on a file
    Matrix<T> data;
    data.ReadText(data_file);

    // first column : abscissae x_i
    // second column : values of the index at at points x_i
    int N = data.GetM();
    Vector<Real_wp> x(N);
    Vector<T> y(N);
    for (int i = 0; i < N; i++)
      {
        x(i) = realpart(data(i, 0));
        y(i) = data(i, 1);
      }

    Init(x, y, cte, false);
  }


  //! Inits the index with points x_i and associated values y_i
  template<class T>
  void SplineInterpolationIndex<T>::Init(const Vector<Real_wp>& x,
					 const Vector<T>& y, T& cte,
					 bool first_order)
  {
    this->points = x;
    this->all_points = x;
    this->nb_points_grid_c = x.GetM();
    
    if (first_order)
      this->order = 1;
    else
      this->order = 3;
    
    int N = y.GetM();
    spline.Init(x, y, first_order);
    // the constant is set to the last value (constant at infinity)
    cte = y(N-1);
  }


  template<class T>
  void SplineInterpolationIndex<T>::SetApproximation(const VectReal_wp& points_reg, int r)
  {
    int nb_intervals = points_reg.GetM()-1;
    this->order = r;
    this->nb_points_grid_c = nb_intervals+1;
    this->points = points_reg;
    this->all_points = points_reg;

    Vector<T> Ones(this->nb_points_grid_c); Ones.Fill(1);
    bool first_order = (this->order == 1);
    spline.Init(this->points, Ones, first_order);
  }
  
  
  /*****************************
   * BSplineInterpolationIndex *
   *****************************/


  //! Inits the index with knots and associated weights
  template<class T>
  void BSplineInterpolationIndex<T>::Init(const string& data_file, int r, T& cte)
  {
    Matrix<T> data;
    data.ReadText(data_file);
    int N = data.GetM();
    
    VectReal_wp x(N);
    Vector<T> y(N);
    for (int i = 0; i < N; i++)
      {
	x(i) = realpart(data(i, 0));
	y(i) = data(i, 1);
      }

    cte = y(N-1);
    this->order = r;
    knots = x;
    this->nb_points_grid_c = x.GetM();
    this->points.Clear();
    this->all_points.Clear();
    bspline.SetData(x, y, this->order);
  }


  //! Projects a function f to the approximation space
  template<class T>
  void BSplineInterpolationIndex<T>
  ::ProjectFunction(T (&f)(Real_wp), Vector<T>& projF)
  {
    projF.Reallocate(this->all_points.GetM());
    for (int i = 0; i < this->all_points.GetM(); i++)
      projF(i) = f(this->all_points(i));
    
    bspline.Init(this->all_points, projF, this->order, this->points);
    projF = bspline.GetWeights();    
  }


  template<class T>
  void BSplineInterpolationIndex<T>::SetApproximation(const Vector<string>& param, const VectReal_wp& points_reg)
  {
    if (param.GetM() <= 2)
      {
	cout << "Enter the polynomial order, example : ApproximationC = BSpline nb_points order" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }
    
    this->order = to_num<int>(param(2));
    this->points = points_reg;
    
    knots.Reallocate(this->points.GetM() + 2*this->order);
    for (int i = 0; i < this->order; i++)
      {
	knots(i) = this->points(0);
	knots(this->points.GetM() + this->order + i) = this->points(this->points.GetM()-1);
      }
    
    for (int i = 0; i < this->points.GetM(); i++)
      knots(this->order+i) = this->points(i);
    
    this->nb_points_grid_c = knots.GetM() - this->order - 1;
    Vector<T> Ones(this->nb_points_grid_c); Ones.Fill(1);
    bspline.SetData(knots, Ones, this->order);

    Real_wp a = points_reg(0);
    Real_wp b = points_reg(points_reg.GetM()-1);
    int nb_intervals = points_reg.GetM();
    
    // minimization with regular points
    Linspace(a, b, nb_intervals + 3*this->order+2, this->all_points);
  }


  /********************************
   * PolynomialInterpolationIndex *
   ********************************/
  

  //! Default constructor
  template<class T>
  PolynomialInterpolationIndex<T>::PolynomialInterpolationIndex()
  {
    offset_per_elt = 0;
    discontinuous = false;
  }
  

  //! returns the size of the object in bytes
  template<class T>
  inline size_t PolynomialInterpolationIndex<T>::GetMemorySize() const
  {
    return yn.GetMemorySize() + lob_basis.GetMemorySize();
  }

  
  //! computes basis functions at a given point x
  template<class T>
  void PolynomialInterpolationIndex<T>::ComputeValuesPhiRef(const Real_wp& x, VectReal_wp& phi) const
  {
    VectReal_wp phi_loc(this->order+1);
    int pos = this->FindInterval(x);
    Real_wp L = (x - this->points(pos)) / (this->points(pos+1) - this->points(pos));
    lob_basis.ComputeValuesPhiRef(L, phi_loc);
    
    int offset = pos*offset_per_elt;
    phi.Reallocate(this->nb_points_grid_c);
    phi.Zero();
    for (int i = 0; i <= this->order; i++)
      phi(offset + i) = phi_loc(i);
  }
  

  //! Sets values associated with interpolation points
  template<class T>
  void PolynomialInterpolationIndex<T>::Init(const Vector<T>& cn)
  {
    yn = cn;
  }
  

  //! Multiplies values of the index by a coefficient
  template<class T>
  inline void PolynomialInterpolationIndex<T>::Mlt(const T& coef)
  {
    yn *= coef;
  }
  
  
  //! returns the value of the index at a given point x
  template<class T>
  T PolynomialInterpolationIndex<T>::Evaluate(const Real_wp& x) const
  {
    int pos = this->FindInterval(x);
    VectReal_wp phi(this->order+1);
    Real_wp L = (x - this->points(pos)) / (this->points(pos+1) - this->points(pos));
    lob_basis.ComputeValuesPhiRef(L, phi);
    int offset = pos*this->offset_per_elt;
    T rho; SetComplexZero(rho);
    for (int j = 0; j < phi.GetM(); j++)
      rho += phi(j) * yn(offset + j);
    
    return rho;      
  }
  

  //! evaluate the index and its derivative at a given point
  template<class T>
  void PolynomialInterpolationIndex<T>::EvaluateDerivative(const Real_wp& x, T&, T&) const
  {
    cout << "not implemented" << endl;
    abort();
  }
  

  //! returns a duplicata of the current index
  template<class T>
  VariableParameter1D_Base<T>* PolynomialInterpolationIndex<T>::GetDuplicate() const
  {
    return new PolynomialInterpolationIndex<T>(*this);
  }


  //! sets parameters of the polynomial approximation
  template<class T>
  void PolynomialInterpolationIndex<T>::SetApproximation(const Vector<string>& param, const VectReal_wp& points_reg)
  {
    if (param.GetM() <= 2)
      {
	cout << "Enter the polynomial order, example : ApproximationC = Polynomial nb_points order" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }
    
    this->order = to_num<int>(param(2));
    this->discontinuous = false;
    if (param.GetM() > 3)
      {
	if (param(3) == "Discontinuous")
	  this->discontinuous = true;
      }

    this->points = points_reg;
    int nb_intervals = points_reg.GetM()-1;
    if (this->discontinuous)
      {
	this->offset_per_elt = this->order+1;
	// Gauss points are used as interpolation on each interval
	lob_basis.ConstructQuadrature(this->order, lob_basis.QUADRATURE_GAUSS);
	
	this->nb_points_grid_c = (this->order+1)*nb_intervals;
      }
    else
      {
	this->offset_per_elt = this->order;
	// Gauss-Lobatto points are used as interpolation on each interval
	lob_basis.ConstructQuadrature(this->order, lob_basis.QUADRATURE_LOBATTO);
	this->nb_points_grid_c = this->order*nb_intervals + 1;
      }
    
    this->all_points.Reallocate(this->nb_points_grid_c);
    for (int i = 0; i < nb_intervals; i++)
      {
	Real_wp dx = this->points(i+1) - this->points(i);
	for (int j = 0; j < this->offset_per_elt; j++)
	  this->all_points(i*offset_per_elt + j) = this->points(i) + dx*lob_basis.Points(j);
      }
	
    // the two extremities
    if (!this->discontinuous)
      {
	this->all_points(0) = points_reg(0);
	this->all_points(this->nb_points_grid_c-1) = points_reg(nb_intervals);
      }
  }
  
  
  /************************
   * PhysicalVaryingMedia *
   ************************/
  

  //! default constructor
  template<class T>
  PhysicalVaryingMedia<Dimension1, T>::PhysicalVaryingMedia()
  {
    type = CONSTANT;
    SetComplexZero(cte_coef);
    compute_grad = false;
    index = NULL;    
  }


  //! Destructor
  template<class T>
  PhysicalVaryingMedia<Dimension1, T>::~PhysicalVaryingMedia()
  {
    if (index != NULL)
      delete index;
  }


  //! Copy constructor
  template<class T>
  PhysicalVaryingMedia<Dimension1, T>
  ::PhysicalVaryingMedia(const PhysicalVaryingMedia<Dimension1, T>& var)
  {
    index = NULL;
    *this = var;
  }


  //! Assignment operator
  template<class T>
  PhysicalVaryingMedia<Dimension1, T>& PhysicalVaryingMedia<Dimension1, T>
  ::operator=(const PhysicalVaryingMedia<Dimension1, T>& var)
  {
    type = var.type;
    eval_coef = var.eval_coef;
    grad_coef = var.grad_coef;
    cte_coef = var.cte_coef;
    offset_eval = var.offset_eval;
    compute_grad = var.compute_grad;

    // previous index is deleted if present
    if (index != NULL)
      delete index;

    // arrays are duplicated
    if (var.index != NULL)
      index = var.index->GetDuplicate();
    
    return *this;
  }  
  

  //! returns the size of the object in bytes
  template<class T>
  inline size_t PhysicalVaryingMedia<Dimension1, T>::GetMemorySize() const
  {
    size_t taille = eval_coef.GetMemorySize() + grad_coef.GetMemorySize() + offset_eval.GetMemorySize();
    if (index != NULL)
      taille += index->GetMemorySize();
    
    return taille;
  }

    
  //! reads the index from parameters of the data file
  template<class T>
  void PhysicalVaryingMedia<Dimension1, T>
  ::SetInputData(const VectString& param, int& nb)
  {
    // previous index is removed if present
    if (index != NULL)
      {
	delete index;
	index = NULL;
      }
    
    if (param(nb) == "SPLINE")
      {
        nb++;
        type = VARIABLE;
        string nom = param(nb++);
	SplineInterpolationIndex<T>* index_spline;
	index_spline = new SplineInterpolationIndex<T>();
        index_spline->Init(nom, cte_coef);
	index = index_spline;
      }
    else if (param(nb) == "BSPLINE")
      {
        nb++;
        type = VARIABLE;
	int order = to_num<int>(param(nb++));
        string nom = param(nb++);
	BSplineInterpolationIndex<T>* index_bspline;
	index_bspline = new BSplineInterpolationIndex<T>();
	index_bspline->Init(nom, order, cte_coef);
	index = index_bspline;
      }
    else if (param(nb) == "FILE")
      {
        nb++;
        type = USER;
        string nom = param(nb++);
        eval_coef.ReadText(nom);
      }
    else
      {
        // numerical value
        type = CONSTANT;
        cte_coef = to_num<T>(param(nb++));
      }
  }
    

#ifdef MONTJOIE_WITH_ONE_DIM
  //! allocates arrays needed to evaluate the coefficient on quadrature points
  template<class T>
  void PhysicalVaryingMedia<Dimension1, T>
  ::InitArray(int ref, const VarProblem_1D& var)
  {
    // counting the number of elements of reference ref
    int nb_elt = 0;
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      if (var.mesh.Element(i).GetReference() == ref)
        nb_elt++;

    // constructing offset_eval
    offset_eval.Reallocate(nb_elt+1);
    nb_elt = 0;
    int offset = 0;
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      if (var.mesh.Element(i).GetReference() == ref)
        {
          offset_eval(nb_elt) = offset;
          // internal quadrature points
          offset += var.GetReferenceElement(i).GetNbPointsQuadratureInside();
          // the two extremities
          offset += 2;
          nb_elt++;
        }

    offset_eval(nb_elt) = offset;
    if (type == VARIABLE)
      {
        eval_coef.Reallocate(offset);
        eval_coef.Zero();
	if (compute_grad)
	  {
	    grad_coef.Reallocate(offset);
	    grad_coef.Zero();	    
	  }
      }
  }
#endif
  
  
  //! computes and stores coefficient at point (i, j) if needed
  template<class T>
  void PhysicalVaryingMedia<Dimension1, T>
  ::ComputeCoefficient(const VarPhysicalProblem& var, int i, int j, const Real_wp& x)
  {
    if (compute_grad)
      {
	T f, df;
	SetComplexZero(f);
	SetComplexZero(df);
	if (type == VARIABLE)
	  {
	    index->EvaluateDerivative(x, f, df);    
	    
	    eval_coef(GetOffsetPoint(var, i) + j) = f;
	    grad_coef(GetOffsetPoint(var, i) + j) = df;
	  }
      }
    else
      {
	if (type == VARIABLE)
	  eval_coef(GetOffsetPoint(var, i) + j) = index->Evaluate(x);    
      }
  }
  
  
  //! Evaluates the coefficient at point x
  template<class T>
  const T PhysicalVaryingMedia<Dimension1, T>::EvaluateCoefficient(const Real_wp& x)
  {
    if (type == VARIABLE)
      return index->Evaluate(x);
    
    return cte_coef;
  }
  

  //! sets coefficients on quadrature points
  template<class T>
  void PhysicalVaryingMedia<Dimension1, T>
  ::SetUserIndex(const Vector<T>& user_coef, const T& cte)
  {
    eval_coef = user_coef;
    cte_coef = cte;
    type = USER;
  }
  
  
}

#define MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_CXX
#endif
