#ifndef MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HELIO_HelmholtzAxisymEquation>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension2, 1>(var), var_helm(var)
  {
    this->read_param_condition = false;
  }

}

#define MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_INLINE_CXX
#endif
