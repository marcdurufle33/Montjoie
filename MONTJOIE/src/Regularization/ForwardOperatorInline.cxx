#ifndef ITREG_FILE_FORWARD_OPERATOR_INLINE_CXX

namespace itreg
{

  //! Destructor
  template<class T, class VectorRhs, class VectorSol>
  inline ForwardOperator<T, VectorRhs, VectorSol>::~ForwardOperator()
  {
  }


  //! Evalues y = F(x)
  template<class T, class VectorRhs, class VectorSol>
  inline void ForwardOperator<T, VectorRhs, VectorSol>::Evaluate2(const VectorSol& x, VectorRhs& y)
  {
    // default choice : we call Evaluate
    Init(x);
    Evaluate(y);
  }
  

  //! Returns true if the exact solution is known
  template<class T, class VectorRhs, class VectorSol>
  inline bool ForwardOperator<T, VectorRhs, VectorSol>::ExactSolutionKnown() const
  {
    return false;
  }
    
  
  //! Returns the distance to the exact solution
  template<class T, class VectorRhs, class VectorSol>
  inline T ForwardOperator<T, VectorRhs, VectorSol>::DistanceToExactSolution(const VectorSol& x)
  {
    return T(1);
  }

  
  //! returns the L^2 norm of a solution vector
  template<class T, class VectorRhs, class VectorSol>
  inline T ForwardOperator<T, VectorRhs, VectorSol>::GetNorm2_Sol(const VectorSol& x)
  {
    return Norm2(x);
  }


  //! returns the L^2 norm of a right hand side vector
  template<class T, class VectorRhs, class VectorSol>
  inline T ForwardOperator<T, VectorRhs, VectorSol>::GetNorm2_Rhs(const VectorRhs& x)
  {
    return Norm2(x);
  }
  
    
  //! gives the iterate x^n (in order to draw it for instance)
  template<class T, class VectorRhs, class VectorSol>
  inline void ForwardOperator<T, VectorRhs, VectorSol>::GiveIterate(int n, const VectorSol& x)
  {    
  }

}

#define ITREG_FILE_FORWARD_OPERATOR_INLINE_CXX
#endif

