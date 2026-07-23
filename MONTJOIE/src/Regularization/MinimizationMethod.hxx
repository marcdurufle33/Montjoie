#ifndef ITREG_FILE_MINIMIZATION_METHOD_HXX

namespace itreg
{

  template<class T>
  class GslRegMethodP : public RegMethodP<T>
  {
  public:
    T step_size, tol;
    
    GslRegMethodP();
    void SetParameter(const string& keyword, const string& value);
    
  };

  
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class MinimizationGslReg : public RegularizationMethod<T, VectorRhs, VectorSol>,
			     public Montjoie::VirtualMinimizedFunction<T>
  {
  protected:
    const GslRegMethodP<T>& param;
    Vector<T> yn, y_measured, residual;
    T residualNorm;

  public:
    MinimizationGslReg(const GslRegMethodP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    
    virtual void FindInitGuess(Vector<T>& param);
    virtual void EvaluateFunction(const Vector<T>& x, T& feval);
    virtual void EvaluateFunctionGradient(const Vector<T>& x,
					  T& feval, Vector<T>& fjac);
    
    virtual void GiveIntermediateSolution(int step, Vector<T>& x);

  };


  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class MinimizationAlglibReg : public MinimizationGslReg<T, VectorRhs, VectorSol>
  {
  public:
    MinimizationAlglibReg(const GslRegMethodP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    
  };

}

#define ITREG_FILE_MINIMIZATION_METHOD_HXX
#endif

