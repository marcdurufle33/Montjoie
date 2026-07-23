#ifndef ITREG_FILE_NEWTON_TIKH_CXX

namespace itreg
{
  
  //! Default constructor
  template<class T>
  IRGNMP<T>::IRGNMP()
  {
    maxInnerSteps = 100;
    alpha0 = 10.0;
    rAlpha = 1.5;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void IRGNMP<T>::SetParameter(const string& keyword, const string& value)
  {    
    if (keyword == "alpha0")
      alpha0 = to_num<T>(value);
    else if (keyword == "alphaR")
      rAlpha = to_num<T>(value);
    else
      RegMethodP<T>::SetParameter(keyword, value);
  }


  //! Default constructor
  template<class T>
  IRGNM_with_G_P<T>::IRGNM_with_G_P()
  {
    weight_G = 1.0;
    type_G = 1;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void IRGNM_with_G_P<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "weight_G")
      weight_G = to_num<T>(value);
    else if (keyword == "type_G")
      type_G = to_num<int>(value);
    else
      IRGNMP<T>::SetParameter(keyword, value);
  }
  

  //! Default constructor
  template<class T>
  LM2P<T>::LM2P()
  {
    this->alpha0 = 0.001;
    this->rAlpha = 1.2;
    rho = 0.8;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void LM2P<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "rho")
      rho = to_num<T>(value);
    else
      IRGNMP<T>::SetParameter(keyword, value);
  }
  
  
  //! Default constructor
  template<class T>
  IterTikhP<T>::IterTikhP()
  {
    nrInnerSteps = 2;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void IterTikhP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "nrInnerSteps")
      nrInnerSteps = to_num<int>(value);
    else
      IRGNMP<T>::SetParameter(keyword, value);
  }
  

  //! Default constructor
  template<class T>
  NLTikhP<T>::NLTikhP()
  {
    alpha = 0.001;
  }
  

  //! modifies a parameter with a keyword and associated value
  template<class T>
  void NLTikhP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "alpha")
      alpha = to_num<T>(value);
    else
      IRGNMP<T>::SetParameter(keyword, value);
  }
  
  
  /**************
   * NewtonTikh *
   **************/
  
  
  //! Constructor with parameters and forward operator
  template<class T, class VectorRhs, class VectorSol>
  NewtonTikh<T, VectorRhs, VectorSol>
  ::NewtonTikh(const IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    nrNoiseCopies = 25;
  }
  
  
  //! Solves F(X) = Y with regularized Gauss-Newton (with Tikhonov regularization)
  template<class T, class VectorRhs, class VectorSol>
  void NewtonTikh<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    VectorSol x0(xn), h(xn), rhs2(xn);
    VectorRhs yn(ydelta), rhs1(ydelta);
    Matrix<T> A;
    rhs2.Zero();
    
    this->F.Init(xn);

    // residual is computed
    this->F.Evaluate(yn);
    rhs1 -= yn;
    
    int step = 0;
    bool mustInit = true;
    
    this->residualNormInit = this->F.GetNorm2_Rhs(rhs1);
    T residualNorm = this->residualNormInit;
    this->OutputStep(step, xn, residualNorm);
    
    // main loop
    while (!this->Stop(step, delta, xn, residualNorm))
      {
	if (mustInit)
	  this->F.FullDerivative(xn, yn, JacobianMatrix);
	
	if (param.stoppingRule == RegMethodP<T>::eLepskij_ran)
	  EstimateRisk(this->rhoLepskij(step), this->expectedResidual(step), delta, nrNoiseCopies);
	
	// descent vector is computed
	ComputeUpdate(step, x0, xn, rhs1, rhs2, A, h);
	
	step++;
	xn += h;
	
	mustInit = MustInit(step);
	
	// residual is updated
	if (mustInit)
	  {
	    this->F.Init(xn);	    
	    this->F.Evaluate(yn);
	  }
	else
	  this->F.Evaluate2(xn, yn);
	
	rhs1 = ydelta - yn;
	residualNorm = this->F.GetNorm2_Rhs(rhs1);

	// displays informations about current iteration
	this->OutputStep(step, xn, residualNorm);

	// int test_input; cout << "Waiting..." << endl; cin >> test_input;
      }    
    
    if (param.stoppingRule == RegMethodP<T>::eLepskij_ran)
      {
	if (mustInit)
	  this->F.FullDerivative(xn, yn, JacobianMatrix);
	
	EstimateRisk(this->rhoLepskij(step), this->expectedResidual(step), delta, nrNoiseCopies);
      }
    
    if (param.stoppingRule != RegMethodP<T>::eDiscrepancy)
      this->FindBestIterate(delta);
  }

  
  //! Constructs the linear system to be solved
  template<class T, class VectorRhs, class VectorSol>
  void NewtonTikh<T, VectorRhs, VectorSol>
  ::SetupEqs(const T& factor, const VectorRhs& rhs1, const VectorSol& rhs2,
	     Matrix<T>& A, Vector<T>& rhs)
  {
    int m = this->JacobianMatrix.GetM();
    int n = this->JacobianMatrix.GetN();
    A.Reallocate(m+n, n);
    rhs.Reallocate(m+n);
    A.Zero();
    
    // Jacobian part
    for (int i = 0; i < m; i++)
      {
	rhs(i) = rhs1(i);
	for (int j = 0; j < n; j++)
	  A(i, j) = this->JacobianMatrix(i, j);	  
      }

    // Gram-matrix part
    // Gram matrix is assumed to be the identity matrix
    for (int i = 0; i < n; i++)
      {
	rhs(i+m) = factor*rhs2(i);
	A(i+m, i) = factor;
      }
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  T NewtonTikh<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    T one(1);
    return param.alpha0*pow(one/param.rAlpha, k);
  }


  template<class T, class VectorRhs, class VectorSol>
  void NewtonTikh<T, VectorRhs, VectorSol>
  ::EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries)
  {
    cout << "Not implemented " << endl;
    abort();
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void NewtonTikh<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& residualNorm)
  {
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, residualNorm);
    if (param.print_level >= 1)
      cout << "(" << residualNorm/this->residualNormInit << ")"
	   << "\t sqrtalpha = " << sqrt(regPar(step)) << endl;    
  }
  

  /*********
   * IRGNM *
   *********/
  
  
  //! Constructor with parameters and a forward operator
  template<class T, class VectorRhs, class VectorSol>
  IRGNM<T, VectorRhs, VectorSol>
  ::IRGNM(const IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonTikh<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  //! Inner iteration
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, Matrix<T>& A, VectorSol& h)
  {
    T factor = sqrt(this->regPar(step));
    rhs2 = x0 - xn;

    Vector<T> tau, rhs;    
    this->SetupEqs(factor, rhs1, rhs2, A, rhs);
    
    GetQR(A, tau);
    SolveQR(A, tau, rhs);

    for (int i = 0; i < h.GetM(); i++)
      h(i) = rhs(i);
  }


  /****************
   * IRGNM_with_G *
   ****************/


  template<class T, class VectorRhs, class VectorSol>
  IRGNM_with_G<T, VectorRhs, VectorSol>::
  IRGNM_with_G(const IRGNM_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    G = &op;
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_with_G<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    VectorSol x0(xn), h(xn), rhs2(xn);
    VectorRhs yn(ydelta), rhs1(ydelta);
    Matrix<T> A;
    rhsG.Reallocate(yn.GetM());
    
    this->F.Init(xn);
    G->Init(xn);
    
    bool mustInit = true;
    this->F.Evaluate(yn);
    G->Evaluate(rhsG);
    rhsG = -rhsG;    
    rhs1 = ydelta - yn;
    T normYdelta = this->F.GetNorm2_Rhs(ydelta);
    
    int step = 0;
    this->residualNormInit = this->F.GetNorm2_Rhs(rhs1) / normYdelta;
    T residualNorm = this->residualNormInit;
    this->OutputStep(step, xn, residualNorm);
    
    while (!this->Stop(step, delta, xn, residualNorm))
      {
	if (mustInit)
	  {
	    this->F.FullDerivative(xn, yn, this->JacobianMatrix);
	    G->FullDerivative(xn, rhsG, DG);
	  }
	
	this->ComputeUpdate(step, x0, xn, rhs1, rhs2, A, h);
	xn += h;
	
	step++;
	mustInit = this->MustInit(step);
	if (mustInit)
	  {
	    G->Init(xn);
	    this->F.Init(xn);
	    this->F.Evaluate(yn);
	    G->Evaluate(rhsG);
	    rhsG = -rhsG;
	    rhs2.Zero();
	  }
	else
	  {
	    this->F.Evaluate2(xn, yn);
	    G->Evaluate2(xn, rhsG);
	    rhsG = -rhsG;
	  }

	rhs1 = ydelta - yn;
	residualNorm = this->F.GetNorm2_Rhs(rhs1) / normYdelta;

	// displays informations about current iteration
	this->OutputStep(step, xn, residualNorm);
      }

    if (param.stoppingRule != RegMethodP<T>::eDiscrepancy)
      this->FindBestIterate(delta);
  }


  template<class T, class VectorRhs, class VectorSol>
  void IRGNM_with_G<T, VectorRhs, VectorSol>
  ::SetupEqs(const T& factorGram, const VectorRhs& rhs1, const VectorSol& rhs2,
	     Matrix<T>& A, Vector<T>& rhs)
  {
    T factorG = factorGram * param.weight_G;
    A.Reallocate(this->JacobianMatrix.GetM() + DG.GetM() + rhs2.GetM(), this->JacobianMatrix.GetN());
    rhs.Reallocate(A.GetM());
    A.Zero();
    for (int i = 0; i < this->JacobianMatrix.GetM(); i++)
      {
	rhs(i) = rhs1(i);
	for (int j = 0; j < this->JacobianMatrix.GetN(); j++)
	  A(i, j) = this->JacobianMatrix(i, j);
      }
    
    int offsetG = this->JacobianMatrix.GetM();
    for (int i = 0; i < this->DG.GetM(); i++)
      {
	rhs(offsetG+i) = factorG*rhsG(i);
	for (int j = 0; j < this->JacobianMatrix.GetN(); j++)
	  A(offsetG+i, j) = factorG * this->DG(i, j);
      }
    
    offsetG += this->DG.GetM();
    for (int i = 0; i < rhs2.GetM(); i++)
      {
	rhs(offsetG+i) = factorGram*rhs2(i);
	A(offsetG+i, i) = factorGram;
      }
  }


  /****************
   * LevMarquardt *
   ****************/

  
  template<class T, class VectorRhs, class VectorSol>
  LevMarquardt<T, VectorRhs, VectorSol>
  ::LevMarquardt(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonTikh<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }


  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardt<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, Matrix<T>& A, VectorSol& h)
  {
    T factor = sqrt(this->regPar(step));
    rhs2.Zero();
    
    Vector<T> tau, rhs;    
    this->SetupEqs(factor, rhs1, rhs2, A, rhs);
    
    GetQR(A, tau);
    SolveQR(A, tau, rhs);
    for (int i = 0; i < h.GetM(); i++)
      h(i) = rhs(i);

    if (step == 0)
      {
	VectorRhs rhs1aux(rhs1);
	Mlt(this->JacobianMatrix, h, rhs1aux);

	rhs1aux -= rhs1;
	T quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
	if (param.print_level >= 1)
	  cout << "quot = " << quot << endl;
	
	while (quot > 0.8)
	  {
	    param.alpha0 /= param.rAlpha;
	    factor = sqrt(this->regPar(step));
	    this->SetupEqs(factor, rhs1, rhs2, A, rhs);
	    
	    GetQR(A, tau);
	    SolveQR(A, tau, rhs);
	    for (int i = 0; i < h.GetM(); i++)
	      h(i) = rhs(i);
	    
	    Mlt(this->JacobianMatrix, h, rhs1aux);
	    rhs1aux -= rhs1;

	    quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
	    if (param.print_level >= 1)
	      cout << "quot = " << quot << endl;
	  }
      }
  }
  

  /*******
   * LM2 *
   *******/  

  
  template<class T, class VectorRhs, class VectorSol>
  LM2<T, VectorRhs, VectorSol>
  ::LM2(LM2P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonTikh<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {    
  }
    

  template<class T, class VectorRhs, class VectorSol>
  void LM2<T, VectorRhs, VectorSol>
  ::ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, Matrix<T>& A, VectorSol& h)
  {
    VectorRhs rhs1aux(rhs1);
    rhs2.Zero();
    T& alphaU = param.alpha0;
    T alpha = alphaU;
    T alphaL = 0;
    
    Vector<T> tau, rhs;    
    this->SetupEqs(sqrt(alpha), rhs1, rhs2, A, rhs);
    
    GetQR(A, tau);
    SolveQR(A, tau, rhs);
    for (int i = 0; i < h.GetM(); i++)
      h(i) = rhs(i);

    Mlt(this->JacobianMatrix, h, rhs1aux);
    rhs1aux -= rhs1;

    T quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
    while (quot < param.rho)
      {
	alphaL = alpha;
	alphaU *= 2;
	
	this->SetupEqs(sqrt(alphaU), rhs1, rhs2, A, rhs);

	GetQR(A, tau);
	SolveQR(A, tau, rhs);
	for (int i = 0; i < h.GetM(); i++)
	  h(i) = rhs(i);
	
	Mlt(this->JacobianMatrix, h, rhs1aux);
	rhs1aux -= rhs1;

	quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
      }

    while (abs(quot-param.rho) > 0.001)
      {
	alpha = (alphaU + alphaL)/2;
	
	this->SetupEqs(sqrt(alpha), rhs1, rhs2, A, rhs);

	GetQR(A, tau);
	SolveQR(A, tau, rhs);
	for (int i = 0; i < h.GetM(); i++)
	  h(i) = rhs(i);
	
	Mlt(this->JacobianMatrix, h, rhs1aux);
	rhs1aux -= rhs1;

	quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
	if (quot < param.rho)
	  alphaL = alpha;
	else
	  alphaU = alpha;
      }    
  }


  template<class T, class VectorRhs, class VectorSol>
  T LM2<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    return param.alpha0;
  }
  

  /************
   * IterTikh *
   ************/


  template<class T, class VectorRhs, class VectorSol>
  IterTikh<T, VectorRhs, VectorSol>
  ::IterTikh(IterTikhP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonTikh<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  

  template<class T, class VectorRhs, class VectorSol>
  T IterTikh<T, VectorRhs, VectorSol>::LepskijNoiseFct(int iterIndex, const T& delta)
  {
    T res = 0;
    if (param.stoppingRule == RegMethodP<T>::eLepskij)
      res = param.nrInnerSteps*delta /sqrt(this->regPar(iterIndex));// deterministic estimate
    else if (param.stoppingRule == RegMethodP<T>::eLepskij_ran)
      res = this->rhoLepskij(iterIndex); // statistical estimate
    
    return res;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void IterTikh<T, VectorRhs, VectorSol>
  ::EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries)
  {
    cout << "Not implemented " << endl;
    abort();    
  }
    

  template<class T, class VectorRhs, class VectorSol>
  void IterTikh<T, VectorRhs, VectorSol>::
  ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		  VectorSol& rhs2, Matrix<T>& A, VectorSol& h)
  {
    T factor = sqrt(this->regPar(step));
    rhs2 = x0 - xn;
    
    Vector<T> tau, rhs;    
    this->SetupEqs(factor, rhs1, rhs2, A, rhs);

    GetQR(A, tau);
    SolveQR(A, tau, rhs);
    for (int i = 0; i < h.GetM(); i++)
      h(i) = rhs(i);
    
    int m = rhs1.GetM(), n = xn.GetM();
    for (int innerStep = 2; innerStep <= param.nrInnerSteps; innerStep++)
      {
	rhs.Reallocate(m+n);
	for (int i = 0; i < m; i++)
	  rhs(i) = rhs1(i);

	for (int i = 0; i < n; i++)
	  {
	    rhs(m+i) = factor*h(i);
	    rhs2(i) = h(i);
	  }
	
	SolveQR(A, tau, rhs);
	for (int i = 0; i < h.GetM(); i++)
	  h(i) = rhs(i);
      }

    if (step == 0)
      {
	VectorRhs rhs1aux(rhs1);
	Mlt(this->JacobianMatrix, h, rhs1aux);
	rhs1aux -= rhs1;
	T quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);
	if (param.print_level >= 1)
	  cout << "quot = " << quot << endl;
	
	while (quot > 0.8)
	  {
	    param.alpha0 /= param.rAlpha;
	    rhs2 = x0 - xn;
	    factor = sqrt(this->regPar(step));	    
	    this->SetupEqs(factor, rhs1, rhs2, A, rhs);
	    
	    GetQR(A, tau);
	    SolveQR(A, tau, rhs);
	    for (int i = 0; i < h.GetM(); i++)
	      h(i) = rhs(i);
	    
	    for (int innerStep = 2; innerStep <= param.nrInnerSteps; innerStep++)
	      {
		rhs.Reallocate(m+n);
		for (int i = 0; i < m; i++)
		  rhs(i) = rhs1(i);
		
		for (int i = 0; i < n; i++)
		  {
		    rhs(m+i) = factor*h(i);
		    rhs2(i) = h(i);
		  }
				
		SolveQR(A, tau, rhs);
		for (int i = 0; i < h.GetM(); i++)
		  h(i) = rhs(i);
	      }
	    
	    Mlt(this->JacobianMatrix, h, rhs1aux);
	    rhs1aux -= rhs1;
	    quot = this->F.GetNorm2_Rhs(rhs1aux) / this->F.GetNorm2_Rhs(rhs1);

	    if (param.print_level >= 1)
	      cout << "quot = " << quot << endl;	    
	  }
      }
  }


  /**********
   * NLTikh *
   **********/

  
  template<class T, class VectorRhs, class VectorSol>
  NLTikh<T, VectorRhs, VectorSol>
  ::NLTikh(const NLTikhP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : NewtonTikh<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    residualNormOld = 1.0;
    alphaNLTikh = 0;
    //param.alpha * param.deltaMax;
  }
    

  template<class T, class VectorRhs, class VectorSol>
  T NLTikh<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    T aux = 10*exp(log(T(2)/3)*k);
    return (aux > this->alphaNLTikh ? aux : this->alphaNLTikh);
  }
  

  template<class T, class VectorRhs, class VectorSol>
  bool NLTikh<T, VectorRhs, VectorSol>
  ::Stop(int step, const T& delta, const VectorSol& xn, const T& residualNorm)
  {
    if (step == 0)
      residualNormOld = residualNorm;
    else
      {
	if ((residualNormOld - residualNorm) / residualNormOld < 1e-6)
	  return true;
	
	residualNormOld = residualNorm;
	return NewtonTikh<T, VectorRhs, VectorSol>::Stop(step, delta, xn, residualNorm);
      }
    
    return false;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void NLTikh<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& normResidual)
  {
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, normResidual);
    if (param.print_level >= 1)
      cout << "\t alpha = " << this->regPar(step) << endl;    
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void NLTikh<T, VectorRhs, VectorSol>::
  ComputeUpdate(int step, const VectorSol& x0, const VectorSol& xn, const VectorRhs& rhs1,
		VectorSol& rhs2, Matrix<T>& A, VectorSol& h)
  {
    rhs2 = x0 - xn;        
    T factor = sqrt(this->regPar(step));	    

    Vector<T> rhs, tau;
    this->SetupEqs(factor, rhs1, rhs2, A, rhs);
    
    GetQR(A, tau);
    SolveQR(A, tau, rhs);
    for (int i = 0; i < h.GetM(); i++)
      h(i) = rhs(i);
  }
    

  /***************
   * IRGNMFrozen *
   ***************/

  
  template<class T, class VectorRhs, class VectorSol>
  IRGNMFrozen<T, VectorRhs, VectorSol>
  ::IRGNMFrozen(const IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : IRGNM<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  

  /**********************
   * LevMarquardtFrozen *
   **********************/


  template<class T, class VectorRhs, class VectorSol>
  LevMarquardtFrozen<T, VectorRhs, VectorSol>
  ::LevMarquardtFrozen(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : LevMarquardt<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }


  /***********************
   * LevMarquardtGenuine *
   ***********************/


  template<class T, class VectorRhs, class VectorSol>
  LevMarquardtGenuine<T, VectorRhs, VectorSol>
  ::LevMarquardtGenuine(const LM2P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
    ydelta_ptr = NULL;
  }


  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardtGenuine<T, VectorRhs, VectorSol>
  ::EvaluateF(const VectorSol& x, VectorRhs& feval)
  {
    this->F.Init(x);
    this->F.Evaluate(feval);
    feval -= *ydelta_ptr;    
  }

  

  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardtGenuine<T, VectorRhs, VectorSol>
#ifdef SELDON_WITH_MPI
  ::EvaluateJacobian(const VectorSol& x, VectorRhs& feval, DistributedMatrix<T, General, RowMajor>& fjac)
#else
  ::EvaluateJacobian(const VectorSol& x, VectorRhs& feval, Matrix<T>& fjac)
#endif
  {
    this->F.EvaluateJacobian(x, feval, fjac);
    feval -= *ydelta_ptr;    
  }

  
  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardtGenuine<T, VectorRhs, VectorSol>
  ::GiveIterate(int step, VectorSol& x, const T& normResidual)
  {
    this->OutputStep(step, x, normResidual);
    if (param.print_level >= 1)
      cout << endl;
  }


#ifdef SELDON_WITH_MPI
  template<class T, class VectorRhs, class VectorSol>
  T LevMarquardtGenuine<T, VectorRhs, VectorSol>::GetNorm2_Rhs(const VectorRhs& f)
  {
    T norm = DotProd(f, f), res = norm;
    Vector<int64_t> xtmp;
    MpiAllreduce(this->comm, &norm, xtmp, &res, 1, MPI_SUM);    
    return sqrt(res);
  }

  
  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardtGenuine<T, VectorRhs, VectorSol>
  ::InitJacobian(DistributedMatrix<T, General, RowMajor>& A)
  {
    A.SetCommunicator(this->comm);
  }
#endif
    
  
  template<class T, class VectorRhs, class VectorSol>
  void LevMarquardtGenuine<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {  
    ydelta_ptr = const_cast<VectorRhs*>(&ydelta);
    this->m_ = ydelta.GetM();
    this->n_ = xn.GetM();
    Montjoie::SolveLeastSquaresLvm(*this, xn, param.delta, param.delta, param.maxSteps, 0, param.alpha0);
  }
  
}

#define ITREG_FILE_NEWTON_TIKH_CXX
#endif

