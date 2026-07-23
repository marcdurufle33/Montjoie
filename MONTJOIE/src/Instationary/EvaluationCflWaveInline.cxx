#ifndef MONTJOIE_FILE_EVALUATION_CFL_WAVE_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class TypeEquation>
  inline EigenProblemTimeMontjoie
  ::EigenProblemTimeMontjoie(HyperbolicProblem<TypeEquation>& var)
    : var_problem(var.var_harmonic), var_boundary(var.var_harmonic),
      var_time(var), var_comm(var.var_harmonic)
  {
    InitDefaultValues();
  }


  //! evaluation of CFL 
  template<class TypeEquation>
  inline Real_wp EvaluateCFL(HyperbolicProblem<TypeEquation>& var_time,
			     const UnivariatePolynomial<Real_wp>& P, bool init_computation)
  {
    typedef typename TypeEquation::TypeEquationStationary TypeEquationW;
    EllipticProblem<TypeEquationW>& var = var_time.var_harmonic;
    EigenProblemTimeMontjoie pb_eig(var_time);
    
    return EvaluateCFL_Unsteady(var, var_time, pb_eig, P, init_computation);
  }
  
}

#endif
