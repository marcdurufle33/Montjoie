#ifndef MONTJOIE_FILE_TIME_ACOUSTIC_INLINE_CXX

namespace Montjoie
{
  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline TimeAcoustic_Base<Dimension>::TimeAcoustic_Base(HyperbolicProblem<TypeEquation>& var)
    : 
#ifdef MONTJOIE_WITH_TIME_REVERSAL 
    time_reversal(var),
#endif
    var_time(var), var_problem(var.var_harmonic), var_helm(var.var_harmonic),
    var_output(var.var_harmonic), var_boundary(var.var_harmonic), var_source(var.var_harmonic)
  {
  }


  //! constructor
  template<class Dimension> template<class TypeEquation>
  inline DiscontinuousMassMatrixAcoustic<Dimension>
  ::DiscontinuousMassMatrixAcoustic(HyperbolicProblem<TypeEquation>& var)
    : DiscontinuousUnsteadyMassMatrix<Dimension>(var),
      var_problem(var.var_harmonic), var_time(var), var_boundary(var.var_harmonic),
      var_helm(var.var_harmonic)
  {
  }

  
  //! constructor
  template<class Dimension> template<class TypeEquation>
  inline ContinuousMassMatrixAcoustic<Dimension>
  ::ContinuousMassMatrixAcoustic(HyperbolicProblem<TypeEquation>& var)
    : ContinuousUnsteadyMassMatrix<Dimension>(var),
    var_problem(var.var_harmonic), var_time(var), var_boundary(var.var_harmonic)
  {
  }
  
}

#define MONTJOIE_FILE_TIME_ACOUSTIC_INLINE_CXX
#endif

  
  
