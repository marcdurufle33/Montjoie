#ifndef ITREG_FILE_LANDWEBER_CXX

namespace itreg
{
  
  template<class T>
  LandweberP<T>::LandweberP()
  {
    omega = 0.8;
  }

  
  //! modifies a parameter with a keyword and associated value
  template<class T>
  void LandweberP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "scaling:omega")
      omega = to_num<T>(value);
    else
      RegMethodP<T>::SetParameter(keyword, value);
  }
  

  template<class T>
  ModLandweberP<T>::ModLandweberP()
  {
    omega = 0.001;
    l0 = 100;
    psi = 0.9;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void ModLandweberP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "scaling:omega")
      omega = to_num<T>(value);
    else if (keyword == "l0")
      l0 = to_num<T>(value);
    else if (keyword == "psi")
      psi = to_num<T>(value);
    else
      RegMethodP<T>::SetParameter(keyword, value);
  }
   

  template<class T, class VectorRhs, class VectorSol>
  ScaledLandweber<T, VectorRhs, VectorSol>
  ::ScaledLandweber(LandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void ScaledLandweber<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    VectorRhs yn(ydelta), residual(ydelta);
    VectorSol update(xn);

    // first residual y^delta - F(x^n) is computed
    this->F.Init(xn);
    Scaling(residual, update);
    this->F.Evaluate(yn);

    residual = ydelta - yn;
    this->residualNormInit = this->F.GetNorm2_Rhs(residual);
    T residualNorm = this->residualNormInit;
    
    // first residual is displayed
    int step = 0;
    this->OutputStep(step, xn, residualNorm);
    
    // main loop
    while (!this->Stop(step, delta, xn, residualNorm))
      {
	step++;
	
	this->F.Adjoint(residual, update);
	
	xn += param.omega*update;
	
	// residual is updated
	this->F.Init(xn);
	this->F.Evaluate(yn);
	residual = ydelta - yn;
	residualNorm = this->F.GetNorm2_Rhs(residual);
	
	// the residual is displayed
	this->OutputStep(step, xn, residualNorm);
      }
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void ScaledLandweber<T, VectorRhs, VectorSol>::Scaling(VectorSol& x, VectorRhs& y)
  {
    x.FillRand();
    T coef = T(0.5)/RAND_MAX;
    x = coef*x - T(0.5);
    param.omega /= this->ComputeNorm(y, x, true);
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void ScaledLandweber<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& residualNorm)
  {
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, residualNorm);
    if (param.print_level >= 1)
      cout << endl;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  Landweber<T, VectorRhs, VectorSol>
  ::Landweber(LandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : ScaledLandweber<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }


  /****************
   * ModLandweber *
   ****************/


  template<class T, class VectorRhs, class VectorSol>
  ModLandweber<T, VectorRhs, VectorSol>
  ::ModLandweber(ModLandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void ModLandweber<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
    VectorRhs yn(ydelta), residual(ydelta);
    VectorSol padj(xn), pdiff(xn), x0(xn);

    // first residual y^delta - F(x^n) is computed
    this->F.Init(xn);
    this->F.Evaluate(yn);

    residual = ydelta - yn;
    this->residualNormInit = this->F.GetNorm2_Rhs(residual);
    T residualNorm = this->residualNormInit;
    
    // first residual is displayed
    int step = 0;
    this->OutputStep(step, xn, residualNorm);
    
    // main loop
    while (!this->Stop(step, delta, xn, residualNorm))
      {
	step++;
	
	pdiff = x0 - xn;
	this->F.Adjoint(residual, padj);
	
	xn += param.omega*padj + this->regPar(step)*pdiff;
	
	// residual is updated
	this->F.Init(xn);
	this->F.Evaluate(yn);
	residual = ydelta - yn;
	residualNorm = this->F.GetNorm2_Rhs(residual);
	
	// the residual is displayed
	this->OutputStep(step, xn, residualNorm);
      }

  }
  

  template<class T, class VectorRhs, class VectorSol>
  T ModLandweber<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    return pow(k+param.l0, -param.psi);
  }
    
  
  template<class T, class VectorRhs, class VectorSol>
  void ModLandweber<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& residualNorm)
  {
    RegularizationMethod<T, VectorRhs, VectorSol>::OutputStep(step, xn, residualNorm);
    if (param.print_level >= 1)
      cout << "\t alpha=" << this->regPar(step) << endl;
  }

}

#define ITREG_FILE_LANDWEBER_CXX
#endif
