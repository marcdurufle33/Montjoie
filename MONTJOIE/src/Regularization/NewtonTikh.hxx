#ifndef ITREG_FILE_NEWTON_TIKH_HXX

namespace itreg
{
  
  //! Parameter class for IRGNM-like methods (Iteratively Regularized Gauss Newton)
  template<class T>
  class IRGNMP : public RegMethodP<T>
  {
  public:
    IRGNMP();
    
    int maxInnerSteps;
    T alpha0, rAlpha;
    // the regularization parameter in the n-th step is alpha0*rAlpha^(-n)

    void SetParameter(const string& keyword, const string& value);
  };


  //! Parameter class for class IRGNM_with_G (IRGNM with additional regularization term)
  template<class T>
  class IRGNM_with_G_P : public IRGNMP<T>
  {
  public:
    IRGNM_with_G_P();
    
    T weight_G;
    int type_G;

    void SetParameter(const string& keyword, const string& value);
  };


  //! Parameter class for class LM2 (Levenberg Marquardt algorithm 2)
  template<class T>
  class LM2P : public IRGNMP<T>
  {
  public:
    LM2P();
    
    T rho;

    void SetParameter(const string& keyword, const string& value);
  };

  
  //! Parameter class for class IterTikh (Newton method with iterated Tikhonov regularization)
  template<class T>
  class IterTikhP : public IRGNMP<T>
  {
  public:
    IterTikhP();
    
    int nrInnerSteps;

    void SetParameter(const string& keyword, const string& value);
  };


  //! Parameter class for class NLTikh (nonlinear Tikhonov regularization)
  template<class T>
  class NLTikhP : public IRGNMP<T>
  {
  public:
    NLTikhP();
    
    T alpha;

    void SetParameter(const string& keyword, const string& value);
  };

  
  //! Base class for regularization methods with Tikhonov regularization
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonTikh : public RegularizationMethod<T>
  {
  private:
    const IRGNMP<T>& param;
    
  protected:
    Matrix<T> JacobianMatrix;
    int nrNoiseCopies;
    
  public:
    NewtonTikh(const IRGNMP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    virtual void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);

    virtual void SetupEqs(const T& factor, const VectorRhs& rhs1, const VectorSol& rhs2,
			  Matrix<T>& A, Vector<T>& rhs);
    
    virtual inline bool MustInit(int step) const { return true; }
    virtual T regPar(int k) const;
    virtual void EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries);
    
  protected:
    virtual void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
			       VectorSol&, Matrix<T>&, VectorSol&) = 0;
    
    void OutputStep(int step, const VectorSol&, const T&);
    
  };

  
  //! Iteratively Regularized Gauss-Newton
  /*!
    see Bakushinskii (1992) for details
  */
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM : public NewtonTikh<T, VectorRhs, VectorSol>
  {
  private:
    const IRGNMP<T>& param;
    
  public:
    IRGNM(const IRGNMP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, Matrix<T>&, VectorSol&);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_with_G : public IRGNM<T, VectorRhs, VectorSol>
  {
  private:
    const IRGNM_with_G_P<T>& param;
    
  protected:
    ForwardOperator<T, VectorRhs, VectorSol>* G;
    Vector<T> rhsG; Matrix<T> DG;
    
  public:
    IRGNM_with_G(const IRGNM_with_G_P<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    virtual void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);

    void SetupEqs(const T& factor, const VectorRhs& rhs1, const VectorSol& rhs2,
		  Matrix<T>& A, Vector<T>& rhs);
        
  };


  
  //! Levenberg-Marquardt algorithm with a difference choice of regularization
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class LevMarquardt : public NewtonTikh<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  public:
    LevMarquardt(IRGNMP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    virtual void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
			       VectorSol&, Matrix<T>&, VectorSol&);
    
  };

  
  //! Levenberg-Marquardt algorithm with a difference choice of regularization
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class LM2 : public NewtonTikh<T, VectorRhs, VectorSol>
  {
  private:
    LM2P<T>& param;
    
  public:
    LM2(LM2P<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, Matrix<T>&, VectorSol&);

    T regPar(int k) const;
    
  };


  //! IRGNM with iterated Tikhonov regularization
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IterTikh : public NewtonTikh<T, VectorRhs, VectorSol>
  {
  private:
    IterTikhP<T>& param;
    
  public:
    IterTikh(IterTikhP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);

    T LepskijNoiseFct(int, const T& delta);
    void EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries);
    
  protected:
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, Matrix<T>&, VectorSol&);
    
  };


  //! Nonlinear Tikhonov regularization
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NLTikh : public NewtonTikh<T, VectorRhs, VectorSol>
  {
  private:
    const NLTikhP<T>& param;
    T residualNormOld, alphaNLTikh;
    
  public:
    NLTikh(const NLTikhP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);

    T regPar(int k) const;
    
  protected:
    bool Stop(int step, const T& delta, const VectorSol& xn, const T&);
    void OutputStep(int step, const VectorSol&, const T&);
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, Matrix<T>&, VectorSol&);
    
  };
  
  
  //! Similar class to IRGNM, but the derivative is always evaluated at the initial guess
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNMFrozen : public IRGNM<T, VectorRhs, VectorSol>
  {
  private:
    const IRGNMP<T>& param;
    
  public:
    IRGNMFrozen(const IRGNMP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);

    inline bool MustInit(int step) const { return (step%10 == 0 ? true : false);  }  
    
  };


  //! Similar class to LevMarquardt, but the derivative is always evaluated at the initial guess
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class LevMarquardtFrozen : public LevMarquardt<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  public:
    LevMarquardtFrozen(IRGNMP<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
  
    inline bool MustInit(int step) const { return (step%10 == 0 ? true : false);  }  
    
  };
  

  //! interface with SolveLeastSquaresLvm
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class LevMarquardtGenuine : public RegularizationMethod<T, VectorRhs, VectorSol>,
#ifndef SELDON_WITH_MPI
			      public Montjoie::VirtualLeastSquaresFunction<T, VectorSol, VectorRhs,
									   Matrix<T> >
#else
			      public Montjoie::VirtualLeastSquaresFunction<T, VectorSol, VectorRhs,
									   DistributedMatrix<T, General, RowMajor>,
                                                                           Matrix<T, Symmetric, RowSymPacked> >
#endif
  {
  private:
    const LM2P<T>& param;
    VectorRhs* ydelta_ptr;
    
  public:
    LevMarquardtGenuine(const LM2P<T>&, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void FindInitGuess(VectorSol& x) {};
    void EvaluateF(const VectorSol& x, VectorRhs& feval);
#ifdef SELDON_WITH_MPI
    void EvaluateJacobian(const VectorSol& x, VectorRhs& feval, DistributedMatrix<T, General, RowMajor>& fjac);    
#else
    void EvaluateJacobian(const VectorSol& x, VectorRhs& feval, Matrix<T>& fjac);    
#endif
    void GiveIterate(int step, VectorSol& x, const T&);
#ifdef SELDON_WITH_MPI
    T GetNorm2_Rhs(const VectorRhs& f);
    void InitJacobian(DistributedMatrix<T, General, RowMajor>& A);
#endif
    
    virtual void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    
  };
  
}

#define ITREG_FILE_NEWTON_TIKH_HXX
#endif

