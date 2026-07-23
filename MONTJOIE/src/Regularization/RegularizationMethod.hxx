#ifndef ITREG_FILE_REGULARIZATION_METHOD_HXX

namespace itreg
{
  
  //! Parameter class for class RegularizationMethod
  template<class T>
  class RegMethodP
  {
  public:
    RegMethodP();
    
    enum StoppingRuleT {eDiscrepancy, eLepskij, eLepskij_ran, eLcurve, eGCV};

    enum MethodTyp {eIRGNM, eLevMarquardt, eLevMarquardt2, eIterTikh, eNLTikh,
		    eIRGNM_CG, eIRGNM_CG_Reortho,eIRGNM_CG_Precond, 
		    eNewtonLandweber, eNewtonTschebychev, eNewtonNuMethod,
		    eNewtonCG, eNewtonCGReOrtho, eNewtonCGPrecond, 
		    eIRGNM_CG_Stefan1, eIRGNM_CG_Stefan2, 
		    eBakushLandweber, eBakushTschebychev,eBakushNuMethod,eBakushCG,
		    eLandweber, eScaledLandweber, eModLandweber,
		    eIRGNMFrozen,eLevMarquardtFrozen,eStateSpaceReg,
		    eIRGNMwithG, eIRGNMCGwithG, eIRGNMCGPrecondwithG,
		    eNewtonCGwithG, eLevMarquardtGenuine, eGSL, eALGLIB};
    
    int print_level;
    StoppingRuleT stoppingRule;
    MethodTyp mt;
    
    // T deltaMax, deltaStep, deltaMin;
    // int repetitions;
    // T tauLepskij, radiusDomainF;
    T tau, delta;
    int maxSteps;
    string history_file;
    
    virtual ~RegMethodP(){}
    virtual void SetParameter(const string& keyword, const string& value);
    
  };


  //! Base class for defining an iterative regularization method
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class RegularizationMethod
  {
  private:
    const RegMethodP<T>& param;

  protected:
    ForwardOperator<T, VectorRhs, VectorSol>& F;
    Vector<T> rhoLepskij, expectedResidual;
    T residualNormInit;
    T sigma_noise; // square root of norm of covariance operator of noise
    int nrStepsLepskij;
    Montjoie::MontjoieTimer timer;
#ifdef SELDON_WITH_MPI
    MPI_Comm comm;
#endif
    
  public:
    RegularizationMethod(const RegMethodP<T>& param_,
			 ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    virtual ~RegularizationMethod();
    
#ifdef SELDON_WITH_MPI
    const MPI_Comm& GetCommunicator() const;
    void SetCommunicator(const MPI_Comm&);
#endif

    static void GetNewSolver(const Vector<string>& keywords, 
			     ForwardOperator<T, VectorRhs, VectorSol>& op,
			     RegMethodP<T>*& param_,
			     RegularizationMethod<T, VectorRhs, VectorSol>*& solver_);
    
    virtual void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn) = 0;
    
    T ComputeNorm(const VectorRhs& z, const VectorSol& x, bool transpose = false);
    
    virtual bool LepskijDone(int, const T& delta);
    virtual T regPar(int k) const;
    virtual T LepskijNoiseFct(int, const T& delta);
    virtual int FindBestIterate(const T& delta);

  protected:
    virtual bool Stop(int step, const T& delta, const VectorSol& xn, const T&);
    virtual void OutputStep(int step, const VectorSol&, const T&);
        
  };

}

#define ITREG_FILE_REGULARIZATION_METHOD_HXX
#endif

