#ifndef ITREG_FILE_NEWTON_INNER_HXX

namespace itreg
{
  
  template<class T>
  class NewtonInnerP : public RegMethodP<T>
  {
  public:
    NewtonInnerP();
   
    enum NewtonType {eInexact, eBakush};    
    NewtonType type;

    T rho;
    int maxInnerSteps;
    
    void SetParameter(const string& keyword, const string& value);
    
  };

  template<class T>
  class NewtonLwP : public NewtonInnerP<T>
  {
  public:
    NewtonLwP();
    
    T omega;

    void SetParameter(const string& keyword, const string& value);
  };
    
  template<class T>
  class NewtonNuP : public NewtonInnerP<T>
  {
  public:
    NewtonNuP();
    
    T omega, nu;

    void SetParameter(const string& keyword, const string& value);
  };

  
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonInner : public RegularizationMethod<T, VectorRhs, VectorSol>
  {
  private:
    const NewtonInnerP<T>& param;
    
  protected:
    int sum_inner_steps, innerStep;
    bool justInitialized, mustInit;
    
  public:
    
    NewtonInner(const NewtonInnerP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

    void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    
  protected:
    void OutputStep(int step, const VectorSol& xn, const T& residualNorm);
    
    virtual bool MustInit() const {return true;}

    virtual void InnerIteration(const VectorRhs& residual, const T& residualNorm,
				VectorSol& update) = 0;
    
  };

  
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonLandweber : public NewtonInner<T, VectorRhs, VectorSol>
  {
  private:
    const NewtonLwP<T>& param;
    
  public:
    NewtonLandweber(const NewtonLwP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);

  protected:
    void InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonNuMethod : public NewtonInner<T, VectorRhs, VectorSol>,
			 protected NewtonNuP<T>
  {
  private:
    const NewtonNuP<T>& param;
    
  public:
    NewtonNuMethod(const NewtonNuP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    void InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonCG : protected NewtonInnerP<T>,
		   public NewtonInner<T, VectorRhs, VectorSol>
		   
  {
  private:
    const NewtonInnerP<T>& param;
    
  public:
    NewtonCG(const NewtonInnerP<T>& param, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    void InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonCGReOrtho : protected NewtonInnerP<T>,
			  public NewtonInner<T, VectorRhs, VectorSol>			  
  {
  private:
    const NewtonInnerP<T>& param;
    
  public:
    NewtonCGReOrtho(const NewtonInnerP<T>& param, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:    
    bool MustInit() const;
    void InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x);
    
  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class NewtonCGPrecond : protected NewtonInnerP<T>,
			  public NewtonInner<T, VectorRhs, VectorSol>
			  
  {
  private:
    const NewtonInnerP<T>& param;
    
  public:
    NewtonCGPrecond(const NewtonInnerP<T>& param, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
  protected:
    int nrSingValues;
    Vector<T> singValues;
    Vector<Vector<T> > singVectors;
    
    bool MustInit() const;
    void InnerIteration(const VectorRhs& z, const T& residualNorm, VectorSol& x);
    
    void SolvePreconditioning(const VectorSol& y, VectorSol& res);
    void IterateWithPrecond(const VectorRhs& z, const T&, VectorSol& x);
    void IterateReInit(const VectorRhs& z, const T&, VectorSol& x);
    
  };


  template<class T>
  class SvdOrtho
  {
  protected:
    bool onlyReOrtho;
    Vector<Vector<T> > H, U;
    int k;
    
  public:
    Vector<T> D, UD;

    SvdOrtho(int size, bool ortho = false);

    T ReOrtho(Vector<T>& r);    
    void ComputeSVD(int n, Matrix<T>& V);
    
  };
  
}

#define ITREG_FILE_NEWTON_INNER_HXX
#endif

