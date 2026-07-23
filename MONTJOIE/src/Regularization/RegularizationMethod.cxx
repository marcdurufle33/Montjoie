#ifndef ITREG_FILE_REGULARIZATION_METHOD_CXX

namespace itreg
{
  
  //! Default constructor
  template<class T>
  RegMethodP<T>::RegMethodP()
  {
    stoppingRule = eDiscrepancy;
    tau = 2.0;
    /* deltaMax = 0;
    deltaStep = 0.5;
    deltaMin = 1e-8;
    repetitions = 1;
    tauLepskij = 2.1;
    radiusDomainF = 100.0; */
    maxSteps = 50;    
    print_level = 0;
    delta = 0;
  }


  //! modifies a parameter with a keyword and associated value
  template<class T>
  void RegMethodP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "StoppingRule")
      {
	if (value == "discrepancy")
	  stoppingRule = eDiscrepancy;
	else
	  {
	    cout << "Other rules not implemented" << endl;
	    abort();
	  }
      }
    else if (keyword == "delta")
      delta = to_num<T>(value);
    else if (keyword == "tau")
      tau = to_num<T>(value);
    else if (keyword == "max_nr_steps")
      maxSteps = to_num<int>(value);    
    else if (keyword == "history")
      history_file = value;
    else
      {
	cout << "Unknown keyword = " << keyword << endl;
	abort();
      }
  }
  

  /************************
   * RegularizationMethod *
   ************************/
  

  //! constructor with parameters and given operator
  template<class T, class VectorRhs, class VectorSol>
  RegularizationMethod<T, VectorRhs, VectorSol>
  ::RegularizationMethod(const RegMethodP<T>& param_,
			 ForwardOperator<T, VectorRhs, VectorSol>& op)
    : param(param_), F(op)
  {
    SetComplexOne(residualNormInit);
    sigma_noise = T(0);
    nrStepsLepskij = 0;
    timer.Reset(Montjoie::VirtualTimer::ALL);
    timer.Start(Montjoie::VirtualTimer::ALL);

#ifdef SELDON_WITH_MPI
    comm = MPI_COMM_SELF;
#endif
    
  }
    

  //! Destructor
  template<class T, class VectorRhs, class VectorSol>
  RegularizationMethod<T, VectorRhs, VectorSol>::~RegularizationMethod()
  {
  }

#ifdef SELDON_WITH_MPI
  template<class T, class VectorRhs, class VectorSol>
  const MPI_Comm& RegularizationMethod<T, VectorRhs, VectorSol>::GetCommunicator() const
  {
    return comm;
  }
  
  template<class T, class VectorRhs, class VectorSol>
  void RegularizationMethod<T, VectorRhs, VectorSol>::SetCommunicator(const MPI_Comm& comm_input)
  {
    comm = comm_input;
  }
#endif


  template<class T, class VectorRhs, class VectorSol>
  void RegularizationMethod<T, VectorRhs, VectorSol>
  ::GetNewSolver(const Vector<string>& keywords, 
		 ForwardOperator<T, VectorRhs, VectorSol>& op,
		 RegMethodP<T>*& param_,
		 RegularizationMethod<T, VectorRhs, VectorSol>*& solver_)
  {
    if (keywords.GetM()%2 == 0)
      {
	cout << "No regularization method selected or incorrect parameters" << endl;
	cout << "Given parameters are " << endl << keywords << endl;
	abort();
      }
        
    if (keywords(0) == "IRGNM")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM<T>(*p, op);
      }
    else if (keywords(0) == "LevMarquardt")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eLevMarquardt;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new LevMarquardt<T>(*p, op);
      }
    else if (keywords(0) == "LevMarquardt2")
      {
	LM2P<T>* p = new LM2P<T>();
	p->mt = RegMethodP<T>::eLevMarquardt2;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new LM2<T>(*p, op);	
      }
    else if (keywords(0) == "IterTikh")
      {
	IterTikhP<T>* p = new IterTikhP<T>();
	p->mt = RegMethodP<T>::eIterTikh;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new IterTikh<T>(*p, op);
      }
    else if (keywords(0) == "NonlinearTikhonovReg")
      {
	NLTikhP<T>* p = new NLTikhP<T>();
	p->mt = RegMethodP<T>::eNLTikh;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NLTikh<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM_CG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new IRGNM_CG<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_Reortho")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM_CG_Reortho;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM_CG_Reortho<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_Precond")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM_CG_Precond;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new IRGNM_CG_Precond<T>(*p, op);
      }
    else if (keywords(0) == "NewtonLandweber")
      {
	NewtonLwP<T>* p = new NewtonLwP<T>();
	p->mt = RegMethodP<T>::eNewtonLandweber;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonLandweber<T>(*p, op);
      }
    else if (keywords(0) == "NewtonTschebychev")
      {
	NewtonNuP<T>* p = new NewtonNuP<T>();
	p->mt = RegMethodP<T>::eNewtonTschebychev;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonNuMethod<T>(*p, op);
      }
    else if (keywords(0) == "NewtonNuMethod")
      {
	NewtonNuP<T>* p = new NewtonNuP<T>();
	p->mt = RegMethodP<T>::eNewtonNuMethod;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonNuMethod<T>(*p, op);
      }
    else if (keywords(0) == "NewtonCG")
      {
	NewtonInnerP<T>* p = new NewtonInnerP<T>();
	p->mt = RegMethodP<T>::eNewtonCG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonCG<T>(*p, op);
      }
    else if (keywords(0) == "NewtonCGReOrtho")
      {
	NewtonInnerP<T>* p = new NewtonInnerP<T>();
	p->mt = RegMethodP<T>::eNewtonCGReOrtho;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonCGReOrtho<T>(*p, op);
      }
    else if (keywords(0) == "NewtonCGPrecond")
      {
	NewtonInnerP<T>* p = new NewtonInnerP<T>();
	p->mt = RegMethodP<T>::eNewtonCGPrecond;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new NewtonCGPrecond<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_Stefan1")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM_CG_Stefan1;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM_CG_Stefan1<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_Stefan2")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNM_CG_Stefan2;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new IRGNM_CG_Stefan2<T>(*p, op);
      }
    else if (keywords(0) == "BakushLandweber")
      {
	NewtonLwP<T>* p = new NewtonLwP<T>();
	p->mt = RegMethodP<T>::eBakushLandweber;	
	p->type = NewtonInnerP<T>::eBakush;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonLandweber<T>(*p, op);
      }
    else if (keywords(0) == "BakushTschebychev")
      {
	NewtonNuP<T>* p = new NewtonNuP<T>();
	p->mt = RegMethodP<T>::eBakushTschebychev;
	p->type = NewtonInnerP<T>::eBakush;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonNuMethod<T>(*p, op);
      }
    else if (keywords(0) == "BakushNuMethod")
      {
	NewtonNuP<T>* p = new NewtonNuP<T>();
	p->mt = RegMethodP<T>::eBakushNuMethod;
	p->type = NewtonInnerP<T>::eBakush;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new NewtonNuMethod<T>(*p, op);
      }
    else if (keywords(0) == "BakushCG")
      {
	NewtonInnerP<T>* p = new NewtonInnerP<T>();
	p->mt = RegMethodP<T>::eBakushCG;
	p->type = NewtonInnerP<T>::eBakush;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new NewtonCG<T>(*p, op);
      }
    else if (keywords(0) == "Landweber")
      {
	LandweberP<T>* p = new LandweberP<T>();
	p->mt = RegMethodP<T>::eLandweber;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new Landweber<T>(*p, op);
      }
    else if (keywords(0) == "ScaledLandweber")
      {
	LandweberP<T>* p = new LandweberP<T>();
	p->mt = RegMethodP<T>::eScaledLandweber;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new ScaledLandweber<T>(*p, op);
      }
    else if (keywords(0) == "ModLandweber")
      {
	ModLandweberP<T>* p = new ModLandweberP<T>();
	p->mt = RegMethodP<T>::eModLandweber;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new ModLandweber<T>(*p, op);
      }
    else if (keywords(0) == "IRGNMFrozen")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eIRGNMFrozen;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new IRGNMFrozen<T>(*p, op);
      }
    else if (keywords(0) == "LevMarquardtFrozen")
      {
	IRGNMP<T>* p = new IRGNMP<T>();
	p->mt = RegMethodP<T>::eLevMarquardtFrozen;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new LevMarquardtFrozen<T>(*p, op);
      }
    else if (keywords(0) == "LevMarquardtGenuine")
      {
	LM2P<T>* p = new LM2P<T>();
	p->mt = RegMethodP<T>::eLevMarquardtGenuine;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));

	param_ = p;
	solver_ = new LevMarquardtGenuine<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_with_G")
      {
	IRGNM_with_G_P<T>* p = new IRGNM_with_G_P<T>();
	p->mt = RegMethodP<T>::eIRGNMwithG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM_with_G<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_with_G")
      {
	IRGNM_with_G_P<T>* p = new IRGNM_with_G_P<T>();
	p->mt = RegMethodP<T>::eIRGNMCGwithG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM_CG_with_G<T>(*p, op);
      }
    else if (keywords(0) == "IRGNM_CG_Precond_with_G")
      {
	IRGNM_with_G_P<T>* p = new IRGNM_with_G_P<T>();
	p->mt = RegMethodP<T>::eIRGNMCGPrecondwithG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new IRGNM_CG_Precond_with_G<T>(*p, op);
      }
    else if (keywords(0) == "NewtonCG_with_G")
      {
	NewtonCG_with_G_P<T>* p = new NewtonCG_with_G_P<T>();
	p->mt = RegMethodP<T>::eNewtonCGwithG;
	for (int i = 1; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new NewtonCG_with_G<T>(*p, op);
      }
    else if (keywords(0) == "GSL")
      {
	GslRegMethodP<T>* p = new GslRegMethodP<T>();
	p->mt = RegMethodP<T>::eGSL;
	for (int i = 3; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new MinimizationGslReg<T>(*p, op);	
	Montjoie::VirtualMinimizedFunction<T>* sol_ = 
	  dynamic_cast<Montjoie::VirtualMinimizedFunction<T>* >(solver_);
	
	if (keywords(2) == "CG_PR")
	  sol_->SetGslAlgorithm(sol_->CG);
	else if (keywords(2) == "BFGS")
	  sol_->SetGslAlgorithm(sol_->BFGS2);
	else if (keywords(2) == "DESCENT")
	  sol_->SetGslAlgorithm(sol_->STEEPEST_DESCENT);
	else
	  sol_->SetGslAlgorithm(sol_->CG_FR);
      }
    else if (keywords(0) == "ALGLIB")
      {
	GslRegMethodP<T>* p = new GslRegMethodP<T>();
	p->mt = RegMethodP<T>::eALGLIB;
	for (int i = 3; i < keywords.GetM(); i += 2)
	  p->SetParameter(keywords(i), keywords(i+1));
	
	param_ = p;
	solver_ = new MinimizationAlglibReg<T>(*p, op);	
	Montjoie::VirtualMinimizedFunction<T>* sol_ = 
	  dynamic_cast<Montjoie::VirtualMinimizedFunction<T>* >(solver_);
	
	if (keywords(2) == "CG")
	  sol_->SetGslAlgorithm(sol_->CG);
	else if (keywords(2) == "BFGS")
	  sol_->SetGslAlgorithm(sol_->BFGS2);
      }
    else
      {
	cout << "Unknown regularization method : " << keywords(0) << endl;
	abort();
      }
  }

  
  template<class T, class VectorRhs, class VectorSol>
  int RegularizationMethod<T, VectorRhs, VectorSol>::FindBestIterate(const T& delta)
  {
    cout << "Not implemented" << endl;
    abort();
    return 0;
  }
    
  
  //! returns true if the iterative method should stop
  template<class T, class VectorRhs, class VectorSol>
  bool RegularizationMethod<T, VectorRhs, VectorSol>
  ::Stop(int step, const T& delta, const VectorSol& xn, const T& residualNorm)
  {
    if (param.stoppingRule != RegMethodP<T>::eDiscrepancy)
      {
	//previous_iterates(step) = new VectorSol(xn);
	//nrStepsLepskij = step;	
	cout << "Not implemented " << endl;
	abort();
      }

    if (residualNorm > 6*residualNormInit)
      {
	cout << "The method seems to diverge" << endl;	
	return true;
      }
    
    if (step >= param.maxSteps)
      {
	if (param.print_level >= 1)
	  cout << "Iterations stopped at maximal iteration number" << endl;
	
	return true;
      }
    
    if ((param.stoppingRule == RegMethodP<T>::eDiscrepancy)
	&& (residualNorm <= param.tau*delta))
      return true;

    if ((param.stoppingRule != RegMethodP<T>::eDiscrepancy)
	&& LepskijDone(step-1, delta))
      {
	if (param.print_level >= 1)
	  cout << "Maximum number of iterations for Lepskij performed " << endl;
	
	nrStepsLepskij = step-1;
	return true;
      }
    
    return false;
  }

  
  //! displays information about convergence of the method
  template<class T, class VectorRhs, class VectorSol>
  void RegularizationMethod<T, VectorRhs, VectorSol>
  ::OutputStep(int step, const VectorSol& xn, const T& residualNorm)
  {
    F.GiveIterate(step, xn);
    
    if (param.print_level >= 1)
      {
        T err_direct(0);
	if (this->F.ExactSolutionKnown())
	  {
            err_direct = this->F.DistanceToExactSolution(xn);
            cout << "Error = " << err_direct << '\t';
          }
	
	cout << "Residue at iteration " << step << " = " << residualNorm;

        if (param.history_file.size() > 1)
          {
            timer.Stop(Montjoie::VirtualTimer::ALL);
            double elap_time = timer.GetSeconds(Montjoie::VirtualTimer::ALL);
            
            ofstream file_out;
            file_out.precision(15);
            if (step <= 1)
              file_out.open(param.history_file.data());
            else
              file_out.open(param.history_file.data(), ios::app);
            
            file_out << elap_time << " " << err_direct << " " << residualNorm << '\n';
            file_out.close();

            timer.Start(Montjoie::VirtualTimer::ALL);
          }
      }
  }


  //! Returns an evaluation of spectral radius of DF DF* 
  template<class T, class VectorRhs, class VectorSol>
  T RegularizationMethod<T, VectorRhs, VectorSol>
  ::ComputeNorm(const VectorRhs& z, const VectorSol& x, bool transpose)
  {
    VectorRhs auxY(z); VectorSol auxX(x);
    T normaux = 0, normauxOld = 0;
    int i = 0;
    if (transpose)
      {
	normaux = this->F.GetNorm2_Sol(auxX);
	auxX *= T(1) / normaux;
      }
    else
      {
	normaux = this->F.GetNorm2_Rhs(auxY);
	auxY *= T(1) / normaux;
      }
    
    // basic power iterative method by starting with z
    do
      {
	i++;
	if (transpose)
	  {
	    this->F.Derivative(auxX, auxY);
	    this->F.Adjoint(auxY, auxX);
	  }
	else
	  {
	    this->F.Adjoint(auxY, auxX);
	    this->F.Derivative(auxX, auxY);
	  }
	
	normauxOld = normaux;
	if (transpose)
	  {
	    normaux = this->F.GetNorm2_Sol(auxX);
	    auxX *= T(1) / normaux; 
	  }
	else
	  {
	    normaux = this->F.GetNorm2_Rhs(auxY);
	    auxY *= T(1) / normaux; 
	  }
      }
    while ((normauxOld / normaux < T(0.99)) || (i < 3));
    
    if (param.print_level >= 1)
      cout << "Norm determined in " << i << " its: " << normaux << endl;
    
    return normaux;
  }
  

  template<class T, class VectorRhs, class VectorSol>
  bool RegularizationMethod<T, VectorRhs, VectorSol>::LepskijDone(int its, const T& delta)
  {
    cout << "Not implemented" << endl;
    abort();
    //if ((its >= 0) && (LepskijNoiseFct(its, delta) >= param.radiusDomainF/3))
    //return true;
    
    return false;
  }
  

  template<class T, class VectorRhs, class VectorSol>
  T RegularizationMethod<T, VectorRhs, VectorSol>::regPar(int k) const
  {
    return T(1);
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  T RegularizationMethod<T, VectorRhs, VectorSol>::LepskijNoiseFct(int iterIndex, const T& delta)
  {
    T res = 0;
    if (param.stoppingRule == RegMethodP<T>::eLepskij)
      res = 0.5*delta / sqrt(regPar(iterIndex)); // deterministic estimate
    else if (param.stoppingRule == RegMethodP<T>::eLepskij_ran)
      res = this->rhoLepskij(iterIndex); // statistical estimate;
    
    return res;
  }
     
}

#define ITREG_FILE_REGULARIZATION_METHOD_CXX
#endif

