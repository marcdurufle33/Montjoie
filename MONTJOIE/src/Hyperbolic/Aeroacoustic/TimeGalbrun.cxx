#ifndef MONTJOIE_FILE_TIME_GALBRUN_CXX

namespace Montjoie
{

  template<class Dimension>
  void HyperbolicProblem<TimeGalbrunEquation<Dimension> >
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    rho = this->var_harmonic.eval_rho0(i)(j);
    sigma = this->var_harmonic.eval_rho0(i)(j)*this->var_harmonic.eval_sigma(i)(j);
  }
  
}

#define MONTJOIE_FILE_TIME_GALBRUN_CXX
#endif
