#ifndef ITREG_FILE_IRGNM_CG_CXX

namespace itreg
{

  template<class T>
  NewtonCG_with_G_P<T>::NewtonCG_with_G_P()
  {
    weight_G = 1.0;
    type_G = 1;
    rho = 0.8;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void NewtonCG_with_G_P<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "rho")
      rho = to_num<T>(value);
    else if (keyword == "weight_G")
      weight_G = to_num<T>(value);
    else if (keyword == "type_G")
      type_G = to_num<int>(value);
    else
      IRGNMP<T>::SetParameter(keyword, value);
  }
  

  /************
   * IRGNM_CG *
   ************/
  

  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG<T, VectorRhs, VectorSol>
  ::IRGNM_CG(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)      
  {
    initialNewtonCGPhase = (param.alpha0 < 0 ? true : false);
    weight_G = T(1);
    G = NULL;
    nrNoiseCopies = 10;
    change_NewtonCG_IRGNMCG = 0;
  }

  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    VectorRhs yn(ydelta), rhs1(ydelta);   
    VectorSol x0(xn), rhs2(xn), h(xn);
    h.Zero();
    
    int step = 0;
    sqrtRegParG = weight_G*sqrt(this->regPar(step));
    
    innerSteps = 0;
    sum_inner_steps = 0;
    bool mustInit = true;
    frozenAt = 0;
    if (G != NULL)
      G->Init(xn);
    
    // residual is evaluated
    this->F.Init(xn);
    this->F.Evaluate(yn);
    step_type = 0;
    if (G != NULL)
      {
	G->Init(xn);
	G->Evaluate(rhsG);
	rhsG *= -sqrtRegParG;
      }
    
    rhs1 = ydelta - yn;
    this->residualNormInit = this->F.GetNorm2_Rhs(rhs1);
    this->residualNorm = this->residualNormInit;
    
    rhs2 = x0 - xn;
    sqrtRegPar = 0;
    this->OutputStep(step, xn, this->residualNorm);    

    // main loop
    while (!this->Stop(step, delta, xn, this->residualNorm))
      {
	step++;
	sqrtRegPar = sqrt(this->regPar(step));
	rhs2 *= sqrtRegPar;
	
	// descent vector is computed
	step_type = 0;
	ComputeUpdate(step, x0, xn, rhs1, rhs2, h, mustInit);
	
	// xn is updated
	if (param.stoppingRule != RegMethodP<T>::eDiscrepancy)
	  this->EstimateRisk(this->rhoLepskij(step), this->expectedResidual(step), delta, this->nrNoiseCopies);
	
	xn += h;

	if (initialNewtonCGPhase && (innerSteps > 5))
	  {
	    initialNewtonCGPhase = false;
	    param.alpha0 = FindAlpha0(step);
	    this->change_NewtonCG_IRGNMCG = step;
	  }
	
	// and residual also
	mustInit = this->MustInit(step);
	if (mustInit)
	  {
	    this->F.Init(xn);
	    this->F.Evaluate(yn);
	    frozenAt = step;

	    if (G != NULL)
	      {
		G->Init(xn);
		G->Evaluate(rhsG);
		rhsG *= -sqrtRegParG;
	      }	 
	    
	    sqrtRegParG = weight_G*sqrt(this->regPar(step));
	  }
	else
	  {
	    this->F.Evaluate2(xn, yn);
	    if (G != NULL)
	      {
		G->Evaluate2(xn, rhsG);
		rhsG *= -sqrtRegParG;
	      }	    
	  }

	rhs1 = ydelta - yn;
	this->residualNorm = this->F.GetNorm2_Rhs(rhs1);
	rhs2 = x0 - xn;

	// display informations about current iteration
	this->OutputStep(step, xn, residualNorm);	
      }

    if (param.stoppingRule != RegMethodP<T>::eDiscrepancy)
      this->FindBestIterate(delta);
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>
  ::EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries)
  {
    cout << "Not implemented " << endl;
    abort();
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    int Nadd = this->rhsG.GetM();
    Vector<T> rhs(rhs1.GetM() + rhs2.GetM() + Nadd);
    rhs.Zero();
    for (int i = 0; i < rhs1.GetM(); i++)
      rhs(i) = rhs1(i);

    int offset = rhs1.GetM();
    for (int i = 0; i < rhs2.GetM(); i++)
      rhs(offset+i) = rhs2(i);
    
    h.Zero();
    innerSteps = Cgne(rhs, h, T(1e-1), 100);    
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG<T, VectorRhs, VectorSol>::FindAlpha0(int step)
  {
    return -param.alpha0; 
  }
  

  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG<T, VectorRhs, VectorSol>::MustInit(int step) const
  {
    return true;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    return (initialNewtonCGPhase ? 0 : param.alpha0 * pow(T(1)/param.rAlpha, k));
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>::OutputStep(int step, const VectorSol& xn, const T& res)
  {
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, res);
    if (param.print_level >= 1)
      cout << "(" << res / this->residualNormInit << ")\t sqrta= " << this->sqrtRegPar
	   <<"\tstp=" << this->innerSteps << endl;    
  }
  

  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>
  ::ExtractVector(const Vector<T>& x, int Nsol, Vector<T>& x1, Vector<T>& x2, Vector<T>& x3)
  {
    int Nrhs = x.GetM() - Nsol;
    if (G != NULL)
      Nrhs -= this->rhsG.GetM();
    
    T* x_ptr = x.GetData();
    x1.SetData(Nrhs, x_ptr); x_ptr += Nrhs;
    x2.SetData(Nsol, x_ptr); x_ptr += Nsol;
    if (G != NULL)
      x3.SetData(this->rhsG.GetM(), x_ptr);
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>::MltVector(const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> y1, y2, y3;
    this->ExtractVector(y, x.GetM(), y1, y2, y3);
    
    this->F.Derivative(x, y1);
    y2 = this->sqrtRegPar*x;
    if (G != NULL)
      {
	G->Derivative(x, y3);
	y3 *= this->sqrtRegParG;
      }

    y1.Nullify(); y2.Nullify(); y3.Nullify();
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG<T, VectorRhs, VectorSol>::MltVectorTranspose(const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> x1, x2, x3;
    this->ExtractVector(x, y.GetM(), x1, x2, x3);

    this->F.Adjoint(x1, y);
    y += this->sqrtRegPar*x2;
    if (G != NULL)
      {
	Vector<T> aux;
	G->Adjoint(x3, aux);
	y += this->sqrtRegParG*aux;
      }

    x1.Nullify(); x2.Nullify(); x3.Nullify();
  }


  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG<T, VectorRhs, VectorSol>
  ::GetScalarProduct_SolCgne(const Vector<T>& x, const Vector<T>& y)
  {
    return DotProd(x, y);
  }


  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG<T, VectorRhs, VectorSol>::GetNorm2_SolCgne(const VectorSol& r)
  {
    return Norm2(r);
  }


  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG<T, VectorRhs, VectorSol>::GetNorm2_RhsCgne(const Vector<T>& r)
  {
    return Norm2(r);
  }


  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG<T, VectorRhs, VectorSol>
  ::StopCgne(const T& alpha, const T& norm, const VectorSol& x0, const T& eps)
  {
    if (initialNewtonCGPhase)
      return (norm <= 0.8*this->residualNorm);
    
    return (alpha <= square(this->sqrtRegPar)*this->GetNorm2_SolCgne(x0)*eps);
  }


  template<class T, class VectorRhs, class VectorSol>
  int IRGNM_CG<T, VectorRhs, VectorSol>
  ::Cgne(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter)
  {
    int m = rhs.GetM(), n = x0.GetM();
    Vector<T> d(m), r(n), q(m), p(n);
    
    this->MltVector(x0, d);
    d = rhs - d;
    this->MltVectorTranspose(d, r);
    
    T alpha = this->GetNorm2_SolCgne(r);
    if (this->StopCgne(alpha, this->GetNorm2_RhsCgne(d), x0, eps))
      return 1;
    
    p = r;
    int k = 0;
    T normRes, help, beta;
    
    do
      {
	k++;
	this->MltVector(p, q);
	normRes = this->GetNorm2_SolCgne(r);
	alpha = square(normRes / this->GetNorm2_RhsCgne(q));
	x0 += alpha*p;
	d -= alpha*q;
	
	this->MltVectorTranspose(d, r);
	help = this->GetNorm2_SolCgne(r);
	beta = square(help / normRes);
	p = beta*p + r;
      }
    while (!this->StopCgne(help, this->GetNorm2_RhsCgne(d), x0, eps) && (k < nb_max_iter-1));
    
    return k+1;
  }
  
  
  /*******************
   * IRGNM_CG_with_G *
   *******************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_with_G<T, VectorRhs, VectorSol>
  ::IRGNM_CG_with_G(IRGNM_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG<T, VectorRhs, VectorSol>(param_, op), param(param_)      
  {
    this->G = &op;
    this->weight_G = param.weight_G;
  }


  /********************
   * IRGNM_CG_Reortho *
   ********************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_Reortho<T, VectorRhs, VectorSol>
  ::IRGNM_CG_Reortho(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG<T, VectorRhs, VectorSol>(param_, op), param(param_)      
  {
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Reortho<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    Vector<T> rhs(rhs1.GetM() + rhs2.GetM());
    for (int i = 0; i < rhs1.GetM(); i++)
      rhs(i) = rhs1(i);

    for (int i = rhs1.GetM(); i < rhs.GetM(); i++)
      rhs(i) = rhs2(i - rhs1.GetM());
    
    h.Zero();
    bool init = this->initialNewtonCGPhase;
    this->initialNewtonCGPhase = false;
    this->innerSteps = CgneFullReortho(rhs, h, T(1)/3, 100);    
    this->initialNewtonCGPhase = init;
  }


  template<class T, class VectorRhs, class VectorSol>
  int IRGNM_CG_Reortho<T, VectorRhs, VectorSol>
  ::CgneFullReortho(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter)
  {
    int m = rhs.GetM(), n = x0.GetM();
    Vector<T> d(m), r(n), q(m), p(n);
    this->MltVector(x0, d);
    d = rhs - d;    
    this->MltVectorTranspose(d, r);
    
    T alpha = this->GetNorm2_SolCgne(r);
    if (this->StopCgne(alpha, this->GetNorm2_RhsCgne(d), x0, eps))
      return 1;
    
    Vector<Vector<T> > H(nb_max_iter+1), Q(nb_max_iter+1);
    Q(0) = r;
    GetReorthogonalization(H, 0, Q(0));
    ApplyHouseholderTransformation(H(0), 0, r);
    
    Q(0) *= r(0);
    p = Q(0);
    r = Q(0);
    
    int k = 0;
    T normRes, help, beta;
    
    do
      {
	k++;
	this->MltVector(p, q);
	normRes = this->GetNorm2_SolCgne(Q(k-1));
	alpha = square(normRes / this->GetNorm2_RhsCgne(q));
	x0 += alpha*p;
	d -= alpha*q;
	
	this->MltVectorTranspose(d, r);
	Q(k) = r;
	if (k != nb_max_iter)
	  {
	    for (int i = 0; i < k; i++)
	      {
		ApplyHouseholderTransformation(H(i), i, Q(k));
		ApplyHouseholderTransformation(H(i), i, r);
	      }

	    GetReorthogonalization(H, k, Q(k));
	    ApplyHouseholderTransformation(H(k), k, r);
	    Q(k) *= r(k);
	    r = Q(k);
	  }
	else
	  break;

	help = this->GetNorm2_SolCgne(r);
	beta = square(help / normRes);
	p = beta*p + r;
      }
    while (!this->StopCgne(help, this->GetNorm2_RhsCgne(d), x0, eps) && (k < nb_max_iter-1));
    
    return k+1;
  }  


  /********************
   * Newton_CG_with_G *
   ********************/

  
  template<class T, class VectorRhs, class VectorSol>
  NewtonCG_with_G<T, VectorRhs, VectorSol>
  ::NewtonCG_with_G(NewtonCG_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    this->G = &op;
    this->weight_G = param.weight_G;
  }

  
  template<class T, class VectorRhs, class VectorSol>
  bool NewtonCG_with_G<T, VectorRhs, VectorSol>::MustInit(int step) const
  {
    return true;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void NewtonCG_with_G<T, VectorRhs, VectorSol>::MltVector(const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> y1, y2, y3;
    this->ExtractVector(y, x.GetM(), y1, y2, y3);
    
    this->F.Derivative(x, y1);
    y2.Zero();
    if (this->G != NULL)
      {
	this->G->Derivative(x, y3);
	y3 *= this->sqrtRegParG;
      }

    y1.Nullify(); y2.Nullify(); y3.Nullify();
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void NewtonCG_with_G<T, VectorRhs, VectorSol>::MltVectorTranspose(const Vector<T>& x, Vector<T>& y)
  {
    Vector<T> x1, x2, x3;
    this->ExtractVector(x, y.GetM(), x1, x2, x3);

    this->F.Adjoint(x1, y);
    if (this->G != NULL)
      {
	Vector<T> aux;
	this->G->Adjoint(x3, aux);
	y += this->sqrtRegParG*aux;
      }
    
    x1.Nullify(); x2.Nullify(); x3.Nullify();
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void NewtonCG_with_G<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    this->sqrtRegPar = sqrt(this->regPar(step));
    this->innerSteps = 1;

    int Nadd = this->rhsG.GetM();
    Vector<T> rhs(rhs1.GetM() + rhs2.GetM() + Nadd);
    rhs.Zero();
    for (int i = 0; i < rhs1.GetM(); i++)
      rhs(i) = rhs1(i);

    int offset = rhs1.GetM();
    for (int i = 0; i < rhs2.GetM(); i++)
      rhs(offset+i) = rhs2(i);

    offset += rhs2.GetM();
    for (int i = 0; i < this->rhsG.GetM(); i++)
      rhs(offset+i) = this->rhsG(i);

    int m = rhs.GetM(), n = xn.GetM();
    Vector<T> d(m), q(m);
    Vector<T> p(n), s(n), aux(n);
    
    SvdOrtho<T> H(param.maxInnerSteps, true);

    h.Zero();
    d = rhs;
    MltVectorTranspose(d, s);
    p = s;
    T normsOld = H.ReOrtho(s), norms(0);
    // T stoptest = square(this->sqrtRegPar) * this->F.GetNorm2_Sol(h);
    do
      {
	MltVector(p, q);
	T normq = Norm2(q);
	T alpha = square(normsOld / normq);
	h += alpha*p;
	d -= alpha*q;
	MltVectorTranspose(d, s);
	
	if (this->innerSteps <= param.maxInnerSteps-1)
	  norms = H.ReOrtho(s);
	else
	  norms = Norm2(s);
	
	T beta = square(norms/normsOld);
	normsOld = norms;
	p = beta*p + s;
	// stoptest = square(this->sqrtRegPar) * this->F.GetNorm2_Sol(h);
	this->innerSteps++;
      }
    while (Norm2(d) >= this->residualNorm*param.rho && this->innerSteps < param.maxInnerSteps);
    
  }
  

  /********************
   * IRGNM_CG_Precond *
   ********************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::IRGNM_CG_Precond(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    nrSingValues = 0;
    singValues.Reallocate(param.maxInnerSteps);
    singValues.Zero();
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    int Nadd = this->rhsG.GetM();
    Vector<T> rhs(rhs1.GetM() + rhs2.GetM() + Nadd);
    rhs.Zero();
    for (int i = 0; i < rhs1.GetM(); i++)
      rhs(i) = rhs1(i);
    
    int offset = rhs1.GetM();
    for (int i = 0; i < rhs2.GetM(); i++)
      rhs(offset+i) = rhs2(i);
    
    h.Zero();
    if (mustInit)
      {
	int nb_max_iter = 8;
	if (!this->initialNewtonCGPhase)
	  nb_max_iter = param.maxInnerSteps;
	
	this->innerSteps = CgneFullReorthoWithEigApp(rhs, h, T(1e-9), nb_max_iter,
						     this->singValues, this->singVectors);
	
	int dim = singValues.GetM();
	nrSingValues = 0;
	for (int i = 0; i < dim; i++)
	  {
	    T value = -singValues(i);
	    if ((value > 1.1*square(this->sqrtRegPar)) && 
		((this->sqrtRegPar > 0) || (nrSingValues < dim-1)) )
	      {
		singValues(nrSingValues) = sqrt(value - square(this->sqrtRegPar));
		singVectors(nrSingValues) = singVectors(i);
		nrSingValues++;
	      }
	  }
	
	if (param.print_level >= 1)
	  cout << "discarded " << dim - nrSingValues << endl;
      }
    else
      this->innerSteps = Pcgne(rhs, h, T(1)/3, param.maxInnerSteps);
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  T IRGNM_CG_Precond<T, VectorRhs, VectorSol>::FindAlpha0(int step)
  {
    T res = 0.01*square(singValues(nrSingValues)) * pow(param.rAlpha, step);
    return res;
  }
  

  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MustInit(int step) const
  {
    if ( (sqrt(T(step))-sqrt(T(this->frozenAt))>=0.999999999 && this->innerSteps>8)
       || nrSingValues == 0 
       || this->initialNewtonCGPhase 
       || (!this->initialNewtonCGPhase && step==this->change_NewtonCG_IRGNMCG))
      return true;
    
    return false;    
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVector(const Vector<T>& x, Vector<T>& y)
  {
    IRGNM_CG<T, VectorRhs, VectorSol>::MltVector(x, y);
  }
  

  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVectorTranspose(const Vector<T>& x, Vector<T>& y)
  {
    IRGNM_CG<T, VectorRhs, VectorSol>::MltVectorTranspose(x, y);
  }

  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::SolvePreconditioning(const Vector<T>& x, Vector<T>& y)
  {
    T aux, alpha = square(this->sqrtRegPar);
    y = T(1)/alpha * x;
    for (int i = 0; i < nrSingValues; i++)
      {
	aux = DotProd(x, singVectors(i));
	aux *= -T(1)/alpha + T(1) / (square(singValues(i)) + alpha);
	y += aux*singVectors(i);
      }    
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::SolveSqrtPrec(const Vector<T>& x, Vector<T>& y)
  {
    T aux, alpha = this->sqrtRegPar;
    y = T(1)/alpha * x;
    for (int i = 0; i < nrSingValues; i++)
      {
	aux = DotProd(x, singVectors(i));
	aux *= -T(1)/alpha + T(1) / sqrt(square(singValues(i)) + alpha*alpha);
	y += aux*singVectors(i);
      }    
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::MltSqrtPrec(const Vector<T>& x, Vector<T>& y)
  {
    T aux, alpha = this->sqrtRegPar;
    y = alpha * x;
    for (int i = 0; i < nrSingValues; i++)
      {
	aux = DotProd(x, singVectors(i));
	aux *= -alpha + sqrt(square(singValues(i)) + alpha*alpha);
	y += aux*singVectors(i);
      }    
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Precond<T, VectorRhs, VectorSol>::OutputStep(int step, const VectorSol& xn, const T& res)
  {
    IRGNM_CG<T, VectorRhs, VectorSol>::OutputStep(step, xn, res);
  }
  

  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::Stop(int step, const T& delta, const VectorSol& xn, const T& res)
  {
    return IRGNM_CG<T, VectorRhs, VectorSol>::Stop(step, delta, xn, res);    
  }



  template<class T, class VectorRhs, class VectorSol>
  int IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::Pcgne(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter)
  {
    int m = rhs.GetM(), n = x0.GetM();
    Vector<T> d(m), r(n), q(m), p(n), z(n);
    
    this->MltVector(x0, d);
    d = rhs - d;
    this->MltVectorTranspose(d, r);
    
    T alpha = this->GetNorm2_SolCgne(r);
    if (this->StopCgne(alpha, this->GetNorm2_RhsCgne(d), x0, eps))
      return 1;
    
    this->SolvePreconditioning(r, z);    
    p = z;
    
    int k = 0;
    T normRes, help, beta;
    
    do
      {
	k++;
	this->MltVector(p, q);
	normRes = this->GetScalarProduct_SolCgne(r, z);
	alpha = normRes / square(this->GetNorm2_RhsCgne(q));
	x0 += alpha*p;
	d -= alpha*q;
	
	this->MltVectorTranspose(d, r);
	this->SolvePreconditioning(r, z);
	help = this->GetScalarProduct_SolCgne(r, z);
	beta = help / normRes;
	p = beta*p + z;
      }
    while (!this->StopCgne(sqrt(help), this->GetNorm2_RhsCgne(d), x0, eps) && (k < nb_max_iter-1));
    
    return k+1;
  }


  template<class T, class VectorRhs, class VectorSol>
  int IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  ::CgneFullReorthoWithEigApp(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter,
			      Vector<T>& eigen_values, Vector<Vector<T> >& eigen_vectors)
  {
    int m = rhs.GetM(), n = x0.GetM();
    Vector<T> d(m), r(n), q(m), p(n);
    this->MltVector(x0, d);
    d = rhs - d;    
    this->MltVectorTranspose(d, r);
    
    T alpha = this->GetNorm2_SolCgne(r);
    if (this->StopCgne(alpha, this->GetNorm2_RhsCgne(d), x0, eps))
      return 1;
    
    Vector<Vector<T> > H(nb_max_iter+1), Q(nb_max_iter+1);
    Q(0) = r;
    GetReorthogonalization(H, 0, Q(0));
    ApplyHouseholderTransformation(H(0), 0, r);
    
    Q(0) *= r(0);
    p = Q(0);
    r = Q(0);
    
    int k = 0;
    T normRes, help, beta, betaOverAlpha(0);
    
    Matrix<T, Symmetric, RowSymPacked> tridiag(nb_max_iter, nb_max_iter);
    tridiag.Zero();
    
    Vector<Vector<T> > residual(nb_max_iter);
    help = 1.0/this->GetNorm2_SolCgne(r);
    residual(0) = help*r;
    
    do
      {
	k++;
	this->MltVector(p, q);
	normRes = this->GetNorm2_SolCgne(Q(k-1));
	alpha = square(normRes / this->GetNorm2_RhsCgne(q));
	x0 += alpha*p;
	d -= alpha*q;
	
	this->MltVectorTranspose(d, r);
	Q(k) = r;
	if (k != nb_max_iter)
	  {
	    for (int i = 0; i < k; i++)
	      {
		ApplyHouseholderTransformation(H(i), i, Q(k));
		ApplyHouseholderTransformation(H(i), i, r);
	      }

	    GetReorthogonalization(H, k, Q(k));
	    ApplyHouseholderTransformation(H(k), k, r);
	    Q(k) *= r(k);
	    r = Q(k);
	  }

	help = this->GetNorm2_SolCgne(r);
	if (k < nb_max_iter)
	  residual(k) = 1.0/help*r;
	
	beta = square(help / normRes);
	p = beta*p + r;

	if (k == 1)
	  {
	    tridiag.Set(0, 0, T(1)/alpha);
	    tridiag.Set(0, 1, -sqrt(beta) / alpha);
	  }
	else
	  {
	    tridiag.Set(k-1, k-1, T(1)/alpha + betaOverAlpha);
	    if (k < nb_max_iter)
	      tridiag.Set(k-1, k, -sqrt(beta) / alpha);
	  }
	
	betaOverAlpha = beta/alpha;
      }
    while (!this->StopCgne(help, this->GetNorm2_RhsCgne(d), x0, eps) && (k < nb_max_iter-1));
    
    cout << "k = " << k << endl;
    tridiag.Resize(k, k);
    Matrix<T, General, RowMajor> eigen_vect;
    GetEigenvaluesEigenvectors(tridiag, eigen_values, eigen_vect);
    
    eigen_vectors.Reallocate(k);
    for (int j = 0; j < k; j++)
      eigen_vectors(j).Reallocate(n);
    
    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
	{
	  SetComplexZero(help);
	  for (int l = 0; l < k; l++)
	    help += residual(l)(i)*eigen_vect(l, j);
	  
	  eigen_vectors(j)(i) = help;
	}
    
    eigen_vect.Clear();
    residual.Clear();

    // eigenvalues are reverted 
    eigen_values = -eigen_values;
    Sort(eigen_values, eigen_vectors);
    
    return k+1;
  }


  /***************************
   * IRGNM_CG_Precond_with_G *
   ***************************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_Precond_with_G<T, VectorRhs, VectorSol>
  ::IRGNM_CG_Precond_with_G(IRGNM_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG_Precond<T, VectorRhs, VectorSol>(param_, op), param(param_)      
  {
    this->G = &op;
    this->weight_G = param.weight_G;
  }


  /********************
   * IRGNM_CG_Stefan1 *
   ********************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_Stefan1<T, VectorRhs, VectorSol>
  ::IRGNM_CG_Stefan1(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG_Precond<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG_Stefan1<T, VectorRhs, VectorSol>::MustInit(int step) const
  {
    if ( (sqrt(T(step))-sqrt(T(this->frozenAt))>=0.999999999 && this->update_Prec) 
       || this->nrSingValues == 0 
       || this->initialNewtonCGPhase 
       || (!this->initialNewtonCGPhase && step==this->change_NewtonCG_IRGNMCG))
      return true;
    
    return false;
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan1<T, VectorRhs, VectorSol>
  ::MltVector(const Vector<T>& x, Vector<T>& y)
  {
    if (this->step_type != 2)
      return IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVector(x, y);
    
    Vector<T> x2(x);
    this->SolveSqrtPrec(x, x2);
    IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVector(x2, y);
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan1<T, VectorRhs, VectorSol>
  ::MltVectorTranspose(const Vector<T>& x, Vector<T>& y)
  {
    if (this->step_type != 2)
      return IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVectorTranspose(x, y);

    Vector<T> y2(y);
    IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVectorTranspose(x, y2);
    this->SolveSqrtPrec(y2, y);
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan1<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    Vector<T> rhs(rhs1.GetM() + rhs2.GetM());
    for (int i = 0; i < rhs1.GetM(); i++)
      rhs(i) = rhs1(i);

    for (int i = rhs1.GetM(); i < rhs.GetM(); i++)
      rhs(i) = rhs2(i - rhs1.GetM());

    h.Zero();
    if (mustInit)
      {
	this->step_type = 0;
	int nb_max_iter = 8;
	if (!this->initialNewtonCGPhase)
	  nb_max_iter = param.maxInnerSteps;
	
	this->innerSteps = this->CgneFullReorthoWithEigApp(rhs, h, T(1e-9), nb_max_iter,
							   this->singValues, this->singVectors);
	
	int dim = this->singValues.GetM();
	this->nrSingValues = 0;
	for (int i = 0; i < dim; i++)
	  {
	    T value = -this->singValues(i);
	    if ((value > 1.1*square(this->sqrtRegPar)) && 
		((this->sqrtRegPar > 0) || (this->nrSingValues < dim-1)) )
	      {
		this->singValues(this->nrSingValues) = sqrt(value - square(this->sqrtRegPar));
		this->singVectors(this->nrSingValues) = this->singVectors(i);
		this->nrSingValues++;
	      }
	  }
	
	this->singValues.Resize(this->nrSingValues);
	this->singVectors.Resize(this->nrSingValues);
	
	if (param.print_level >= 1)
	  cout << "discarded " << dim - this->nrSingValues << endl;

	update_counter = 0;
	update_Prec = false;
      }
    else
      {
	if (!update_Prec)
	  {
	    this->step_type = 1;
	    this->innerSteps = this->Pcgne(rhs, h, T(1)/3, param.maxInnerSteps);
	    update_counter++;
	    if ((this->innerSteps > 5 && update_counter >= 4)
		|| (this->innerSteps>10 && update_counter >=3 ))
	      update_Prec = true;	    	    
	  }
	else
	  {
	    Vector<T> new_eigen_values; Vector<Vector<T> > new_eigen_vectors;
	    VectorSol cgne_aux_X1(h);
	    cgne_aux_X1.Zero();
	    this->step_type = 2;
	    this->innerSteps = this->CgneFullReorthoWithEigApp(rhs, cgne_aux_X1, T(1e-2), param.maxInnerSteps,
							       new_eigen_values, new_eigen_vectors);
	    this->MltSqrtPrec(cgne_aux_X1, h);

	    int dim = new_eigen_values.GetM();
	    for (int i = 0; i < dim; i++)
	      {
		T value = -new_eigen_values(i);
		if (value > 1.1)
		  {
		    this->singValues.PushBack(this->sqrtRegPar*sqrt(value - 1.0));
		    this->singVectors.PushBack(new_eigen_vectors(i));
		    this->nrSingValues++;
		  }
	      }

	    // int test_input; cout << "Waiting" << endl; cin >> test_input;
	    update_counter = 0;
	    update_Prec = false;
	  }
      }
  }


  /********************
   * IRGNM_CG_Stefan2 *
   ********************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_CG_Stefan2<T, VectorRhs, VectorSol>
  ::IRGNM_CG_Stefan2(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM_CG_Precond<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    SingVectors = NULL;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  bool IRGNM_CG_Stefan2<T, VectorRhs, VectorSol>::MustInit(int step) const
  {
    if ( sqrt(T(step)) < this->innerSteps-1 && step > this->frozenAt+1 && step > 30)
      return true;

    return false;
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan2<T, VectorRhs, VectorSol>
  ::MltVector(const Vector<T>& x, Vector<T>& y)
  {
    if (this->step_type != 2)
      return IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVector(x, y);
    
    Vector<T> x2(x);
    this->SolveSqrtPrec(x, x2);
    IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVector(x2, y);
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan2<T, VectorRhs, VectorSol>
  ::MltVectorTranspose(const Vector<T>& x, Vector<T>& y)
  {
    if (this->step_type != 2)
      return IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVectorTranspose(x, y);

    Vector<T> y2(y);
    IRGNM_CG_Precond<T, VectorRhs, VectorSol>::MltVectorTranspose(x, y2);
    this->SolveSqrtPrec(y2, y);
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_CG_Stefan2<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, VectorSol& h, bool mustInit)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
}

#define ITREG_FILE_IRGNM_CG_CXX
#endif
