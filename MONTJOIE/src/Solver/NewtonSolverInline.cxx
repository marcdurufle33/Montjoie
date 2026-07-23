#ifndef MONTJOIE_FILE_NEWTON_SOLVER_INLINE_CXX

namespace Montjoie
{

  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::SetStoppingCriterion(Treal eps1, Treal eps2)
  {
    if (eps2 < 0)
      eps2 = eps1;
      
    tol = eps1;
    tol_convergence = eps2; 
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::SetMaxNumberOfIterations(int N)
  {
    nbitermax = N;
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::EnableConvergenceTest(bool test)
  {
    convergence_required = test; 
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::SetAbsoluteSchemeNorm(Treal norm)
  {
    threshold_norme_scheme = norm;
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::ForceReevaluationJacobian(bool force)
  {
    force_reevaluation_jacobian = force;
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::EnableBacktracking(bool yes)
  {
    backtracking = yes; 
  }
  
  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::SetPrintLevel(int lvl)
  {
    print_level = lvl;
  }    

  
  template<class T, class VectorSol, class VectorScale>
  inline void NewtonSolver<T, VectorSol, VectorScale>::EnableDecreasingResidue(bool yes)
  {
    decreasing_residual_required = yes;
  }

  
  template<class T, class VectorSol, class VectorScale>
  void NewtonSolver<T, VectorSol, VectorScale>::SetScalingFactors(const VectorScale& scale)
  {
    scale_newton_test = true;
    scale_newton = scale;
  }
  
}

#define MONTJOIE_FILE_NEWTON_SOLVER_INLINE_CXX
#endif
