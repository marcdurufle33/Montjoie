#ifndef MONTJOIE_FILE_TIME_REVERSAL_TIME_INLINE_CXX

namespace Montjoie
{

  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline TimeReversalSource<Dimension>
  ::TimeReversalSource(const EllipticProblem<TypeEquation>& var,
		       TimeAcousticReversal<Dimension>& time_rev)
    : VirtualSourceFEM<Real_wp, Dimension>(var), time_reversal(time_rev)
  {
    order = 4;
    t0_buffer = 0;
    t1_buffer = 0;
    n0_buffer = 0;
    n1_buffer = 0;
    t_current = 0;
  }


  template<class Dimension>
  inline bool TimeAcousticReversal<Dimension>::IsPresent() const
  {
    return presence_time_reversal;
  }
  
  
  template<class Dimension>
  inline int TimeAcousticReversal<Dimension>::GetSimulationType() const
  {
    return type_simulation;
  }
  
  
  template<class Dimension>
  inline Real_wp TimeAcousticReversal<Dimension>::GetInitialTime() const
  {
    return tbegin;
  }
  
  
  template<class Dimension>
  inline Real_wp TimeAcousticReversal<Dimension>::GetFinalTime() const
  {
    return tend; 
  }


  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline TimeAcousticReversal<Dimension>::TimeAcousticReversal(HyperbolicProblem<TypeEquation>& var)
    : src_reverse(var.var_harmonic, *this), var_problem(var.var_harmonic),
      var_source(var.var_harmonic), var_time(var)
  {
    presence_time_reversal = false;
    type_simulation = DIRECT;
    tbegin = 0;
    tend = 0;
    deltat = 0;
  }


  
}

#define MONTJOIE_FILE_TIME_REVERSAL_TIME_INLINE_CXX
#endif


