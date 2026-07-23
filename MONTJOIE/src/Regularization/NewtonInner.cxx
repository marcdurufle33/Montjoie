#ifndef ITREG_FILE_NEWTON_INNER_CXX

namespace itreg
{
  
  //! Default constructor
  template<class T>
  NewtonInnerP<T>::NewtonInnerP()
  {
    type = eInexact;
    rho = 0.8;
    maxInnerSteps = 50;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void NewtonInnerP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "rho")
      rho = to_num<T>(value);
    else if (keyword == "maxInnerSteps")
      maxInnerSteps = to_num<int>(value);
    else
      RegMethodP<T>::SetParameter(keyword, value);
  }
  
  
  //! Default constructor
  template<class T>
  NewtonLwP<T>::NewtonLwP()
  {
    omega = 0.9;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void NewtonLwP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "scaling:omega")
      omega = to_num<T>(value);
    else
      NewtonInnerP<T>::SetParameter(keyword, value);
  }
  

  //! Default constructor
  template<class T>
  NewtonNuP<T>::NewtonNuP()
  {
    omega = 0.9;
    nu = 0.5;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void NewtonNuP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "scaling:omega")
      omega = to_num<T>(value);
    else if (keyword == "nu")
      nu = to_num<T>(value);
    else
      NewtonInnerP<T>::SetParameter(keyword, value);
  }
  

  /***************
   * NewtonInner *
   ***************/
  
  
  //! Constructor with parameters and forward operator
  template<class T, class VectorRhs, class VectorSol>
  NewtonInner<T, VectorRhs, VectorSol>
  ::NewtonInner(const NewtonInnerP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    innerStep = 0;
    justInitialized = false;
    mustInit = true;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void NewtonInner<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& residualNorm)
  {    
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, residualNorm);
    sum_inner_steps += innerStep;
    if (param.print_level >= 1)
      cout << "\t NrSteps=" << innerStep << endl;
  }
  

  template<class T, class VectorRhs, class VectorSol>
  void NewtonInner<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    T one; SetComplexOne(one);
    VectorRhs yn(ydelta), residual(ydelta);
    VectorSol x0(xn), update(xn);
    
    this->F.Init(xn);
    
    // initial residual is computed, residual = ydelta - F(xn)
    this->F.Evaluate(yn);
    residual = ydelta - yn;
    
    this->residualNormInit = this->F.GetNorm2_Rhs(residual);
    T residualNorm = this->residualNormInit;

    // first residual is displayed
    int step = 0;
    mustInit = true;
    justInitialized = false;
    sum_inner_steps = 0;
    this->OutputStep(step, xn, residualNorm);
    
    // main loop
    while (!this->Stop(step, delta, xn, residualNorm))
      {
	if (param.type == NewtonInnerP<T>::eInexact)
	  update.Zero();
	else
	  update = x0 - xn; 
	
	// inner iteration to solve DF* DF update = DF* residual
	InnerIteration(residual, residualNorm, update);
	
	// Gauss-Newton iterate
	// x^n+1 = x^n + (DF* DF)^{-1} DF* residual
	xn += update;
	
	step++;
	
	mustInit = this->MustInit();
	if (mustInit)
	  {	    
	    this->F.Init(xn);
	    this->F.Evaluate(yn);
	    justInitialized = true;
	  }
	else
	  {
	    this->F.Evaluate2(xn, yn);
	    justInitialized = false;
	  }

	// residual is updated = ydelta - F(xn)
	residual = ydelta - yn;
	residualNorm = this->F.GetNorm2_Rhs(residual);
	
	// display informations about current iteration
	this->OutputStep(step, xn, residualNorm);
	
	// int test_input; cout << "Waiting" << endl; cin >> test_input;
      }
  }


  /*******************
   * NewtonLandweber *
   *******************/
  

  template<class T, class VectorRhs, class VectorSol>  
  NewtonLandweber<T, VectorRhs, VectorSol>
  ::NewtonLandweber(const NewtonLwP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonInner<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonLandweber<T, VectorRhs, VectorSol>
  ::InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    T scaling = param.omega / this->ComputeNorm(z, x);
    VectorRhs zn(z); VectorSol aux(x);
    
    this->innerStep = 0;
    zn = -z;
    do
      {
	this->innerStep++;
	this->F.Adjoint(zn, aux);
	Add(-scaling, aux, x);
	this->F.Derivative(x, zn);
	zn -= z;
      }
    while ((this->F.GetNorm2_Rhs(zn) >= param.rho*residualNorm) && (this->innerStep <= param.maxInnerSteps));
  }
  
  
  /******************
   * NewtonNuMethod *
   ******************/
  

  template<class T, class VectorRhs, class VectorSol>  
  NewtonNuMethod<T, VectorRhs, VectorSol>
  ::NewtonNuMethod(const NewtonNuP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonInner<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
    
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonNuMethod<T, VectorRhs, VectorSol>
  ::InnerIteration(const VectorRhs& y, const T& residualNorm, VectorSol& x)
  {
    VectorSol start(x);
    VectorRhs z(y), zp(y), zpp(y), resi(y);
    T one(1), w, mu, nu = param.nu;
    
    this->innerStep = 1;
    T scaling = param.omega / this->ComputeNorm(y, x);
    
    zp.Zero();
    z *= (2*nu + 2) / (4*nu + 1);
    this->F.Adjoint(z, x);
    x = scaling*x + start;
    this->F.Derivative(x, resi);
    resi -= y;
    
    while ((this->F.GetNorm2_Rhs(resi) >= param.rho*residualNorm) && (this->innerStep <= param.maxInnerSteps))
      {
	this->innerStep++;
	zpp = zp;
	zp = z;
	mu = (this->innerStep-1)*(2*this->innerStep-3)*(2*this->innerStep+2*nu-1);
	mu /= (this->innerStep+2*nu-1)*(2*this->innerStep+4*nu-1)*(2*this->innerStep+2*nu-3);
	w = 4*(2*this->innerStep+2*nu-1)*(this->innerStep+nu-1);
	w /= (this->innerStep+2*nu-1)*(2*this->innerStep+4*nu-1);
	
	this->F.Derivative(x, resi);
	resi -= y;
	z = -w*resi + (one+mu)*zp - mu*zpp;
	this->F.Adjoint(z, x);
	x = scaling*x + start;
      }	   
  }


  /************
   * NewtonCG *
   ************/
  
  
  template<class T, class VectorRhs, class VectorSol>  
  NewtonCG<T, VectorRhs, VectorSol>
  ::NewtonCG(const NewtonInnerP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonInner<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCG<T, VectorRhs, VectorSol>
  ::InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    VectorRhs q(z), r(z);
    VectorSol s(x), d(x);
    
    this->F.Derivative(x, q);
    
    r = z - q;
    this->F.Adjoint(r, s);
    d = s;
    
    T normsOld = this->F.GetNorm2_Sol(s);
    this->innerStep = 1;
    do
      {
	this->F.Derivative(d, q);
	T alpha = square(normsOld / this->F.GetNorm2_Rhs(q));
	x += alpha*d;
	r -= alpha*q;
	
	this->F.Adjoint(r, s);
	T norms = this->F.GetNorm2_Sol(s);
	T beta = square(norms / normsOld);
	normsOld = norms;
	d = beta*d + s;
	this->innerStep++;
      }
    while ((this->F.GetNorm2_Rhs(r) >= this->rho*residualNorm) && (this->innerStep <= this->maxInnerSteps));
    
  }


  /*******************
   * NewtonCGReOrtho *
   *******************/
  

  template<class T, class VectorRhs, class VectorSol>  
  NewtonCGReOrtho<T, VectorRhs, VectorSol>
  ::NewtonCGReOrtho(const NewtonInnerP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonInner<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  

  template<class T, class VectorRhs, class VectorSol>  
  bool NewtonCGReOrtho<T, VectorRhs, VectorSol>::MustInit() const
  {
    if ((this->innerStep > 5) && !this->justInitialized)
      return true;
    
    return false;
  }
  
    
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCGReOrtho<T, VectorRhs, VectorSol>
  ::InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    VectorRhs q(z), d(x), s(x), r(z);

    SvdOrtho<T> H(param.maxInnerSteps, true);
    this->innerStep = 0;
    this->F.Derivative(x, q);
    r = z - q;
    this->F.Adjoint(r, s);
    d = s;

    T alpha, beta, normq;
    T normsOld = H.ReOrtho(s), norms = normsOld;
    do
      {
	this->F.Derivative(d, q);
	normq = this->F.GetNorm2_Rhs(q);
	alpha = square(normsOld / normq);
	x += alpha*d;
	r -= alpha*q;
	this->F.Adjoint(r, s);
	if (this->innerStep <= param.maxInnerSteps-2)
	  norms = H.ReOrtho(s);
	else
	  norms = this->F.GetNorm2_Sol(s);
	
	beta = square(norms / normsOld);
	normsOld = norms;
	d = beta*d + s;
	this->innerStep++;
      }
    while ( this->F.GetNorm2_Rhs(r) >= param.rho*residualNorm && this->innerStep <= param.maxInnerSteps);
  }

  
  /*******************
   * NewtonCGPrecond *
   *******************/
  
  
  template<class T, class VectorRhs, class VectorSol>  
  NewtonCGPrecond<T, VectorRhs, VectorSol>
  ::NewtonCGPrecond(const NewtonInnerP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonInner<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    singValues.Reallocate(param.maxInnerSteps);
    singVectors.Reallocate(param.maxInnerSteps);
    singValues.Zero();
  }
 

  template<class T, class VectorRhs, class VectorSol>  
  bool NewtonCGPrecond<T, VectorRhs, VectorSol>::MustInit() const
  {
    if ((this->innerStep > 5) && !this->justInitialized)
      return true;
    
    return false;
  }

  
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCGPrecond<T, VectorRhs, VectorSol>
  ::InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    if (this->mustInit)
      IterateReInit(z, residualNorm, x);
    else
      IterateWithPrecond(z, residualNorm, x);
  }


  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCGPrecond<T, VectorRhs, VectorSol>
  ::SolvePreconditioning(const VectorSol& y, VectorSol& res)
  {
    T aux;
    res = y;
    for (int i = 0; i < nrSingValues; i++)
      {
	aux = DotProd(y, singVectors(i));
	aux *= T(1)/singValues(i) - T(1);
	Add(aux, singVectors(i), res);
      }
  }

  
  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCGPrecond<T, VectorRhs, VectorSol>
  ::IterateWithPrecond(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    SvdOrtho<T> H(param.maxInnerSteps);
    T alpha, beta, normq, norms, normsOld;
    VectorSol d(x), s(x), aux(x);
    VectorRhs r(z), q(z);
   
    this->innerStep = 0;
    this->F.Derivative(x, q);
    r = z - q;
    this->F.Adjoint(r, aux);
    this->SolvePreconditioning(aux, s);
    d = s;
    normsOld = H.ReOrtho(s);
    
    do
      {
	this->SolvePreconditioning(d, aux);
	this->F.Derivative(aux, q);
	normq = this->F.GetNorm2_Rhs(q);
	alpha = square(normsOld / normq);
	Add(alpha, d, x);
	Add(-alpha, q, r);
	this->F.Adjoint(r, aux);
	
	this->SolvePreconditioning(aux, s);
	if (this->innerStep <= param.maxInnerSteps-2)
	  norms = H.ReOrtho(s);
	else
	  norms = this->F.GetNorm2_Sol(s);
	
	beta = square(norms / normsOld);
	normsOld = norms;
	
	d = beta*d + s;
	H.D(this->innerStep) = T(1)/sqrt(alpha);
	H.UD(this->innerStep) = -sqrt(beta / alpha);
	this->innerStep++;
      }
    while ( this->F.GetNorm2_Rhs(r) >= param.rho*residualNorm && this->innerStep < param.maxInnerSteps);
    
    this->SolvePreconditioning(x, aux);
    x = aux;
  }


  template<class T, class VectorRhs, class VectorSol>  
  void NewtonCGPrecond<T, VectorRhs, VectorSol>
  ::IterateReInit(const VectorRhs& z, const T& residualNorm, VectorSol& x)
  {
    T eps=1e-14;
    int N = param.maxInnerSteps;
    SvdOrtho<T> H(N);

    T alpha, beta, normq, norms, normsOld;
    VectorSol d(x), s(x);
    VectorRhs r(z), q(z);

    this->innerStep = 0;
    this->F.Derivative(x, q);
    r = z - q;
    this->F.Adjoint(r, s);
    d = s;
    normsOld = H.ReOrtho(s);
    
    do
      {
	this->F.Derivative(d, q);
	normq = this->F.GetNorm2_Rhs(q);
	alpha = square(normsOld / normq);
	Add(alpha, d, x);
	Add(-alpha, q, r);
	this->F.Adjoint(r, s);
	
	if (this->innerStep <= N-2)
	  norms = H.ReOrtho(s);
	else
	  norms = this->F.GetNorm2_Sol(s);
	
	beta = square(norms / normsOld);
	normsOld = norms;
	
	d = beta*d + s;
	H.D(this->innerStep) = T(1)/sqrt(alpha);
	H.UD(this->innerStep) = -sqrt(beta / alpha);
	this->innerStep++;
      }
    while ( (this->F.GetNorm2_Rhs(r) >= param.rho*residualNorm) && (this->innerStep < param.maxInnerSteps));
    
    nrSingValues = this->innerStep;
    // iterate further without updating x in order to get a better SVD
    if (this->innerStep < param.maxInnerSteps)
      {
	do
	  {
	    this->F.Derivative(d, q);
	    normq = this->F.GetNorm2_Rhs(q);
	    alpha = square(normsOld / normq);
	    Add(-alpha, q, r);
	    this->F.Adjoint(r, s);
	    
	    if (this->innerStep <= N-2)
	      norms = H.ReOrtho(s);
	    else
	      norms = this->F.GetNorm2_Sol(s);
	    
	    beta = square(norms / normsOld);
	    normsOld = norms;
	    
	    d = beta*d + s;
	    H.D(this->innerStep) = T(1)/sqrt(alpha);
	    H.UD(this->innerStep) = -sqrt(beta / alpha);
	    this->innerStep++;
	  }
	while ( normsOld > eps && this->innerStep < param.maxInnerSteps);
      }
    
    N = this->innerStep;
    Matrix<T> V;
    H.ComputeSVD(N, V);
    
    for (int j = 0; j < nrSingValues; j++)
      singValues(j) = H.D(j);
    
    for (int i = 0; i < nrSingValues; i++)
      {
	singVectors(i).Reallocate(V.GetM());
	for (int j = 0; j < V.GetM(); j++)
	  singVectors(i)(j) = V(j, i);
      }
  }
  

  /************
   * SvdOrtho *
   ************/
  

  template<class T>
  SvdOrtho<T>::SvdOrtho(int size, bool ortho)
  {
    H.Reallocate(size);
    U.Reallocate(size);
    D.Reallocate(size);
    UD.Reallocate(size);
    D.Zero(); UD.Zero();
    onlyReOrtho = ortho;
    k = 0;
  }
  
  
  template<class T>
  T SvdOrtho<T>::ReOrtho(Vector<T>& r)
  {
    Vector<T> w(r);
    
    // Apply previous Householder transformations to w
    for (int l = 0; l <= k-1; l++)
      ApplyHouseholderTransformation(H(l), l, w);
    
    // Householder transformation is computed for w(k:end)
    H(k).Reallocate(w.GetM() - k);
    for (int i = k; i < w.GetM(); i++)
      H(k)(i-k) = w(i);
    
    T v1 = GetHouseholderNormale(H(k));
    
    // Householder transformations are applied to the canonical vector
    w.Zero();
    SetComplexOne(w(k));
    for (int l = k; l >= 0; l--)
      ApplyHouseholderTransformation(H(l), l, w);
    
    r = v1*w;
    if (!onlyReOrtho)
      U(k) = w;
    
    k++;
    return abs(v1);
  }
  

  template<class T>
  void SvdOrtho<T>::ComputeSVD(int n, Matrix<T>& V)
  {
#ifdef SELDON_WITH_LAPACK
    int m = U(0).GetM();
    Vector<T> E(n), work(4*n-4);
    V.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	V(i, j) = U(j)(i);
    
    E.Zero();
    for (int i = 0; i < n-1; i++)
      E(i) = UD(i);
    
    int nru = 0, ncc = 0, ldu = 1, ldc = 1;
    char upper('U'); int info;
    dbdsqr_(&upper, &n, &m, &nru, &ncc, D.GetData(), E.GetData(),
	    V.GetData(), &n, NULL, &ldu, NULL, &ldc, work.GetData(), &info);
    
    if (info != 0)
      {
	cout << "Error during ComputeSVD" << endl;
	abort();
      }
#else
    cout << "Recompile Montjoie with Lapack" << endl;
    abort();
#endif
  }
  
}

#define ITREG_FILE_NEWTON_INNER_CXX
#endif
