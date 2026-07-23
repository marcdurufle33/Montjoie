#ifndef ITREG_FILE_MINIMIZATION_METHOD_CXX

namespace itreg
{

  //! Default constructor
  template<class T>
  GslRegMethodP<T>::GslRegMethodP()
  {
    step_size = 0.01;
    tol = 1e-4;
  }
  

  //! modifies a parameter with a keyword and associated value
  template<class T>
  void GslRegMethodP<T>::SetParameter(const string& keyword, const string& value)
  {
    if (keyword == "step_size")
      step_size = to_num<T>(value);
    else if (keyword == "tol")
      tol = to_num<T>(value);
    else
      RegMethodP<T>::SetParameter(keyword, value);
  }
  

  //! Constructor with parameters and forward operator
  template<class T, class VectorRhs, class VectorSol>
  MinimizationGslReg<T, VectorRhs, VectorSol>
  ::MinimizationGslReg(const GslRegMethodP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : RegularizationMethod<T, VectorRhs, VectorSol>(param_, op), param(param_)
  {
  }
  

  template<class T, class VectorRhs, class VectorSol>
  void MinimizationGslReg<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
#ifdef MONTJOIE_WITH_GSL
    this->n = xn.GetM();
    this->step_size = double(param.step_size);
    this->step_line_minimization = double(param.tol);
    y_measured = ydelta;

    MinimizeParametersGsl(*this, xn, double(delta), param.maxSteps);
#else
    cout << "Recompile Montjoie with GSL" << endl;
    abort();
#endif
  }

  
  template<class T, class VectorRhs, class VectorSol>
  void MinimizationGslReg<T, VectorRhs, VectorSol>::FindInitGuess(Vector<T>& x)
  {
    // nothing, xn contains already the initial guess
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void MinimizationGslReg<T, VectorRhs, VectorSol>
  ::EvaluateFunction(const Vector<T>& x, T& feval)
  {
    this->F.Init(x);
    this->F.Evaluate(yn);
    residual = y_measured - yn;
    
    feval = this->F.GetNorm2_Rhs(residual);
    residualNorm = feval;
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void MinimizationGslReg<T, VectorRhs, VectorSol>
  ::EvaluateFunctionGradient(const Vector<T>& x,
			     T& feval, Vector<T>& fjac)
  {
    this->F.Init(x);
    this->F.Evaluate(yn);
    residual = y_measured - yn;

    feval = this->F.GetNorm2_Rhs(residual);
    residualNorm = feval;    
    
    this->F.Adjoint(residual, fjac);
    fjac *= -T(1)/feval;    
  }
  
  
  template<class T, class VectorRhs, class VectorSol>
  void MinimizationGslReg<T, VectorRhs, VectorSol>
  ::GiveIntermediateSolution(int step, Vector<T>& x)
  {
    this->OutputStep(step, x, residualNorm);
    if (this->param.print_level >= 1)
      cout << endl;
  }


  /*************************
   * MinimizationAlglibReg *
   *************************/

  
  //! Constructor with parameters and forward operator
  template<class T, class VectorRhs, class VectorSol>
  MinimizationAlglibReg<T, VectorRhs, VectorSol>
  ::MinimizationAlglibReg(const GslRegMethodP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op)
    : MinimizationGslReg<T, VectorRhs, VectorSol>(param_, op)
  {
  }
  

  template<class T, class VectorRhs, class VectorSol>
  void MinimizationAlglibReg<T, VectorRhs, VectorSol>
  ::Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn)
  {
#ifdef MONTJOIE_WITH_ALGLIB
    this->n = xn.GetM();
    this->step_size = double(this->param.step_size);
    this->step_line_minimization = double(this->param.tol);
    this->y_measured = ydelta;
    
    MinimizeParametersAlglib(*this, xn, double(delta), this->param.maxSteps);
#else
    cout << "Recompile Montjoie with ALGLIB" << endl;
    abort();
#endif
  }

}

#define ITREG_FILE_MINIMIZATION_METHOD_CXX
#endif

