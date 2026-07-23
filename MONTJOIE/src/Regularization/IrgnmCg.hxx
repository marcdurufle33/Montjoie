#ifndef ITREG_FILE_IRGNM_CG_HXX

namespace itreg
{

  template<class T>
  class NewtonCG_with_G_P : public IRGNMP<T>
  {
  public:
    NewtonCG_with_G_P();
    
    T rho;
    T weight_G;
    int type_G;

    void SetParameter(const string& keyword, const string& value);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG : public RegularizationMethod<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  protected:
    int sum_inner_steps, innerSteps;
    ForwardOperator<T, VectorRhs, VectorSol>* G;
    bool initialNewtonCGPhase;
    T weight_G;
    int frozenAt, step_type;
    
    Vector<T> rhsG;
    T residualNorm, sqrtRegPar, sqrtRegParG;
    int nrNoiseCopies, change_NewtonCG_IRGNMCG;

  public:
    
    IRGNM_CG(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    virtual void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);    
    
    virtual void EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries);
    
  protected:
    virtual void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
			       VectorSol&, VectorSol&, bool);
    
    virtual T FindAlpha0(int);
    virtual bool MustInit(int step) const;
    virtual T regPar(int k) const;
    
    virtual void OutputStep(int step, const VectorSol&, const T&);
    
    void ExtractVector(const Vector<T>& x, int Nsol, Vector<T>& x1, Vector<T>& x2, Vector<T>& x3);
    
    virtual void MltVector(const Vector<T>& x, Vector<T>& y);
    virtual void MltVectorTranspose(const Vector<T>& x, Vector<T>& y);
    
    virtual T GetScalarProduct_SolCgne(const Vector<T>& x, const Vector<T>& y);
    virtual T GetNorm2_RhsCgne(const Vector<T>& r);
    virtual T GetNorm2_SolCgne(const VectorSol& r);
    
    virtual bool StopCgne(const T& alpha, const T& norm, const VectorSol& x0, const T& eps);
    
    virtual int Cgne(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter);
    
  };



  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_with_G : public IRGNM_CG<T, VectorRhs, VectorSol>
  {
  private:
    IRGNM_with_G_P<T>& param;
    
  public:
    
    IRGNM_CG_with_G(IRGNM_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  };
  

  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_Reortho : public IRGNM_CG<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  public:
    
    IRGNM_CG_Reortho(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    int CgneFullReortho(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter);
    
  protected:
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, VectorSol&, bool);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonCG_with_G : public IRGNM_CG<T, VectorRhs, VectorSol>
  {
  private:
    NewtonCG_with_G_P<T>& param;
    
  public:
    
    NewtonCG_with_G(NewtonCG_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    bool MustInit(int step) const;

    void MltVector(const Vector<T>& x, Vector<T>& y);
    void MltVectorTranspose(const Vector<T>& x, Vector<T>& y);
    
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, VectorSol&, bool);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_Precond : public IRGNM_CG<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  protected:
    int nrSingValues;
    Vector<T> singValues;
    Vector<Vector<T> > singVectors;
    
  public:    
    IRGNM_CG_Precond(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void EstimateRisk(T& rhoLepskij, T& expectedResidual, const T& delta, int nrTries);
    
  protected:
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, VectorSol&, bool);

    T FindAlpha0(int);
    bool MustInit(int) const;
    
    void MltVector(const Vector<T>& x, Vector<T>& y);
    void MltVectorTranspose(const Vector<T>& x, Vector<T>& y);
    virtual void SolvePreconditioning(const Vector<T>& x, Vector<T>& y);
    virtual void SolveSqrtPrec(const Vector<T>& x, Vector<T>& y);
    virtual void MltSqrtPrec(const Vector<T>& x, Vector<T>& y);
    
    void OutputStep(int step, const VectorSol&, const T&);
    bool Stop(int step, const T& delta, const VectorSol& xn, const T&);
    
    int CgneFullReorthoWithEigApp(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter,
				  Vector<T>&, Vector<Vector<T> >&);
    
    int Pcgne(const Vector<T>& rhs, VectorSol& x0, const T& eps, int nb_max_iter);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_Precond_with_G : public IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  {
  private:
    IRGNM_with_G_P<T>& param;
    
  public:
    
    IRGNM_CG_Precond_with_G(IRGNM_with_G_P<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  };
  

  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_Stefan1 : public IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  protected:
    bool update_Prec;
    int update_counter;
    
  public:    
    IRGNM_CG_Stefan1(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    bool MustInit(int) const;

    void MltVector(const Vector<T>& x, Vector<T>& y);
    void MltVectorTranspose(const Vector<T>& x, Vector<T>& y);
    
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, VectorSol&, bool);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class IRGNM_CG_Stefan2 : public IRGNM_CG_Precond<T, VectorRhs, VectorSol>
  {
  private:
    IRGNMP<T>& param;
    
  protected:
    SvdOrtho<T>* SingVectors;
    
  public:    
    IRGNM_CG_Stefan2(IRGNMP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    bool MustInit(int) const;

    void MltVector(const Vector<T>& x, Vector<T>& y);
    void MltVectorTranspose(const Vector<T>& x, Vector<T>& y);
    
    void ComputeUpdate(int, const VectorSol&, const VectorSol&, const VectorRhs&,
		       VectorSol&, VectorSol&, bool);
    
  };
  
}

#define ITREG_FILE_IRGNM_CG_HXX
#endif
