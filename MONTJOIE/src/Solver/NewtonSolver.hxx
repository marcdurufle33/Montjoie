#ifndef MONTJOIE_FILE_NEWTON_SOLVER_HXX

namespace Montjoie
{

  //! base class to derive when solving a non-linear system with Newton algorithm
  template<class T, class VectorSol = Vector<T>,
	   class VectorScale = Vector<typename ClassComplexType<T>::Treal> >
  class NonLinearEquations_Newton
  {
  public :
    // computes y = F(x) where F is the non-linear system
    virtual void ComputeScheme(const VectorSol& x, VectorSol& y) = 0;
    
    // computes and factorizes jacobian matrix
    virtual void ComputeAndFactoriseDiff(const VectorSol&, const VectorScale&) = 0;
    
    // solves linear system DF x = b where DF has been previously computed and factorized
    virtual void SolveDifferential(const VectorSol&, VectorSol&) = 0;

    // evaluates norm2 of a vector
    virtual Real_wp GetNorm2Vector(const VectorSol&) = 0;    
    
  };
  
  
  //! Newton solver
  template<class T, class VectorSol = Vector<T>,
	   class VectorScale = Vector<typename ClassComplexType<T>::Treal> >
  class NewtonSolver
  {
  protected:
    typedef typename ClassComplexType<T>::Treal Treal;
    
    Treal tol, tol_convergence; // stopping criteria
    int nbitermax; // le nb max diter pour newton
    VectorSol x_memory, x_prec; // solutions at previous iterates
    VectorSol scheme, save_scheme;// scheme = F(x)        
    VectorScale scale_newton; // scaling vector
    Treal erreur, threshold_norme_scheme;
    bool scale_newton_test, convergence_required;
    bool force_reevaluation_jacobian;
    int print_level;
    bool backtracking;
    ofstream file_out_residual_newton;
    bool decreasing_residual_required;

  public:
    NewtonSolver();
    
    // Inline methods
    void SetStoppingCriterion(Treal eps1, Treal eps2 = Treal(-1));
    void SetMaxNumberOfIterations(int N);
    void EnableConvergenceTest(bool test = true);
    void SetAbsoluteSchemeNorm(Treal);
    void ForceReevaluationJacobian(bool force = true);
    void EnableBacktracking(bool yes = true);
    void SetPrintLevel(int lvl);
    void EnableDecreasingResidue(bool yes = true);
    void SetScalingFactors(const VectorScale& scale);
    
    // Other functions
    size_t GetMemorySize() const;

    void SetInputData(const string & description_field, const Vector<string> &parameters);
    
    void Init(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var, VectorSol& x_init);
    int Solve(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var, VectorSol& x);    
    void Newton(NonLinearEquations_Newton<T, VectorSol, VectorScale>& var, VectorSol & x, int &nb,
                const Treal& norme_ref, bool reevaluate_jacobian = false,
		bool scheme_already_evaluated = false);
    
  };
  
}

#define MONTJOIE_FILE_NEWTON_SOLVER_HXX
#endif
