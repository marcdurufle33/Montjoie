#ifndef ITREG_FILE_FORWARD_OPERATOR_HXX

namespace itreg
{
  
  //! Base class for defining a Forward operator
  /*!
    Inverse problems are solved by constructing the forward operator F(x)
    where x are the parameters and F(x) returns the observables
    Any class that define a new forward operator should derive from this class
    Regularization methods are taking a forward operator as argument
  */
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class ForwardOperator
  {
  public:
    virtual ~ForwardOperator();
    
    // inits the computation of F(x) or F*(x), etc    
    virtual void Init(const VectorSol& x) = 0;
    
    // evaluates the operator F(x), assuming that Init has been called with x
    virtual void Evaluate(VectorRhs& y) = 0;
    
    // creates synthetic data y 
    // from the exact solution x, this function
    // should be different from Evaluate to avoid inverse crime
    virtual void EvaluateDifferently(VectorRhs& y) = 0;
    
    // evaluates the derivative of F for a direction h 
    // it should compute res = DF(h) 
    virtual void Derivative(const VectorSol& h, VectorRhs& res) = 0;

    // evaluates the adjoint operator DF*(h)
    virtual void Adjoint(const VectorRhs& h, VectorSol& res) = 0;

    // evaluates the jacobian matrix (full DF)
    virtual void FullDerivative(const VectorSol&, const VectorRhs&, Matrix<T>&);

    // evaluates at the same time F and DF
    virtual void EvaluateJacobian(const VectorSol& x, VectorRhs& f, Matrix<T>& df);
    
    // additional Evaluate function for frozen Newton methods
    // evaluates F(x) at a different point from the point given in Init
    virtual void Evaluate2(const VectorSol& x, VectorRhs& y);
    
    // exact solution known ?
    virtual bool ExactSolutionKnown() const;
    
    // distance to the exact solution
    virtual T DistanceToExactSolution(const VectorSol& x);

    // L^2 norm of vectors
    virtual T GetNorm2_Rhs(const VectorRhs& x);
    virtual T GetNorm2_Sol(const VectorSol& x);
    
    // gives the iterate x^n (in order to draw it for instance)
    virtual void GiveIterate(int n, const VectorSol& x);
    
    void CheckDerivatives(const Vector<T>&, const T& h);
    
  };

}

#define ITREG_FILE_FORWARD_OPERATOR_HXX
#endif

