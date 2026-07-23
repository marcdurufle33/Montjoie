#ifndef MONTJOIE_FILE_NONLINEAR_EQUATIONS_HXX

#ifdef MONTJOIE_WITH_ALGLIB
#include "optimization.h"
#endif

namespace Montjoie

{
  // Minpack functions
  
  //! base class to derive when solving a non-linear system with Minpack algorithm
  template<class T, class VectorSol = Vector<T>, class MatrixJac = Matrix<T> >
  class NonLinearEquations_Base
  {
  public :
    // computes y = F(x) where F is the non-linear system
    virtual void EvaluateFunction(const VectorSol& x, VectorSol& y) = 0;
    
    // computes Jac = DF(x)
    virtual void EvaluateJacobian(const VectorSol&, MatrixJac& Jac) = 0;
    
  };

  // solving non-linear equations knowing derivatives
  template<class T, class Vector, class Matrix>
  int SolveMinpack(NonLinearEquations_Base<T, Vector, Matrix>& eqns,
		   Vector& x_sol, Vector& fvec, Matrix& fjac, Vector& scale_eqn,
		   const IVect& Control, const VectReal_wp& RControl);
  
  template<class Vector, class Matrix>
  void Compute_QR_Factorisation(Matrix& A, bool pivot, IVect& ipvt,
				Vector& rdiag, Vector& acnorm, Vector& wa);

  template<class Vector, class Matrix>
  void QrSolve(Matrix& r, IVect& ipvt, Vector& diag, Vector& qtb,
	       Vector& x, Vector& sdiag, Vector& wa);

  template<class Vector, class Matrix>
  void qform(Matrix& q, Vector& wa);

  template<class real, class Vector, class Vector2>
  void dogleg(Vector2& r, Vector& diag, Vector& qtb,
	      real& delta, Vector& x, Vector& wa1, Vector& wa2);
  
  template<class Vector, class Matrix>
  void r1mpyq(Matrix& a, Vector& v, Vector& w);

  template<class Vector>
  void r1mpyq(Vector& a, Vector& v, Vector& w);
  
  template<class Vector2, class Vector>
  void r1updt(Vector2& s, Vector& u, Vector& v, Vector& w, bool& sing);

  //! base class to specify a function to minimize with gsl or mkl
  template<class T>
  class VirtualMinimizedFunction
  {
  protected :
    int n;
    int type_algo;
    double step_size;
    double step_line_minimization;
    int step_number;
    
  public :
    // available algorithms in Gsl
    // BFGS : Broyden-Fletcher-Goldfarb-Shanno algorithm 
    // BFGS2 : improved BFGS (algorithms 2.6.2 and 2.6.4 of Fletcher's Practical Methods of Optimization
    // CG : non-linear conjugate gradient with Polak-Ribiere formula
    // CG_FR : non-linear conjugate gradient with Fletcher-Reeves formula
    // STEEPEST_DESCENT : steepest descent as described in Gsl documentation
    // SIMPLEX : nmsimplex in Gsl (algorithm of simplex of Nelder and Mead)
    // SIMPLEX2 : nmsimplex2 in Gsl (algorithm of simplex of Nelder and Mead)
    // SIMPLEX2_RAND : nmsimplex2_rand in Gsl (algorithm of simplex of Nelder and Mead)
    enum{BFGS, CG, STEEPEST_DESCENT, SIMPLEX, SIMPLEX2, SIMPLEX2_RAND,
         CG_FR, BFGS2};
    
    inline VirtualMinimizedFunction() { n = 1; type_algo = BFGS; step_size = 0.01; step_line_minimization = 1e-4;
      step_number = 0; }
    
    inline int GetGslAlgorithm() const { return type_algo; }
    inline void SetGslAlgorithm(int type) { type_algo = type; }

    //! returns the size of the first trial step
    inline double GetStepSize() const { return step_size; }

    //! returns the tolerance for the line minimization procedure
    inline double GetLineMinimizationParameter() const { return step_line_minimization; }
    
    inline int GetM() const { return n; }
    inline void SetStepNumber(int m) { step_number = m; }
    
    virtual void FindInitGuess(Vector<T>& param) = 0;
    virtual void EvaluateFunction(const Vector<T>& x, T& feval) = 0;
    virtual void EvaluateFunctionGradient(const Vector<T>& x,
					  T& feval, Vector<T>& fjac) = 0;
    
    virtual inline void GiveIntermediateSolution(int step, Vector<T>& x) {}
    
#ifdef MONTJOIE_WITH_GSL
    static double my_f(const gsl_vector *v, void *params);
    static void my_df(const gsl_vector *v, void *params, gsl_vector *df);
    static void my_fdf(const gsl_vector *v, void *params, double* f, gsl_vector *df);
#endif

#ifdef MONTJOIE_WITH_ALGLIB
    static void alglib_fdf(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr);
    static void alglib_rep(const alglib::real_1d_array &x, double func, void* ptr); 
#endif
    
  };


#ifdef MONTJOIE_WITH_GSL
  template<class T>
  T MinimizeParametersGsl(VirtualMinimizedFunction<T>& fct, Vector<T>& xsol,
			  double epsilon = 1e-12, unsigned nb_max_iter = 5000);
#endif


#ifdef MONTJOIE_WITH_ALGLIB
  template<class T>
  T MinimizeParametersAlglib(VirtualMinimizedFunction<T>& fct, Vector<T>& xsol,
			     double epsilon = 1e-12, unsigned nb_max_iter = 5000);
#endif
  
} // end namespace

#define MONTJOIE_FILE_NONLINEAR_EQUATIONS_HXX
#endif
