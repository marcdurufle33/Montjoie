#ifndef MONTJOIE_FILE_TIME_AERO_ACOUSTIC_CXX

namespace Montjoie
{

  template<class Dimension>
  void HyperbolicProblem<TimeLinearizedEulerEquation<Dimension> >
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    rho = Real_wp(1);
    sigma = this->var_harmonic.eval_sigma(i)(j);
    
    if (this->var_harmonic.type_model == this->var_harmonic.CONSERVATIVE)
      {
        rho *= this->var_harmonic.eval_rho0(i)(j);
        sigma *= this->var_harmonic.eval_rho0(i)(j);
      }
  }
  
  
  template<class Dimension>
  void HyperbolicProblem<TimeAeroAcousticEquation<Dimension> >
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    rho = Real_wp(1);
    sigma = this->var_harmonic.eval_sigma(i)(j);
  }
  
}

#define MONTJOIE_FILE_TIME_AERO_ACOUSTIC_CXX
#endif
