#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline VarTransmission_Maxwell3D::VarTransmission_Maxwell3D(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Base<Dimension3>(var), fimpedance(var),
      var_problem(var), var_maxwell(var), var_boundary(var)      
  {
    InitDefaultValues();
  }


  template<class TypeEquation>
  inline VarTransmission<HarmonicMaxwellEquation_3D>::VarTransmission(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Maxwell3D(var)
  {
  }
  
  
  template<class TypeEquation>
  inline ImpedanceTransmissionMaxwell3D
  ::ImpedanceTransmissionMaxwell3D(const EllipticProblem<TypeEquation>& vars)
    : ImpedanceFunction_Base<Complex_wp, Dimension3>(vars),
      var_maxwell(vars), var_transmission(vars.var_transmission)
  {
    coef_val = 0.0; coef_grad = 0.0; coef_u = 0.0;
  }

  
  inline bool ImpedanceTransmissionMaxwell3D::PresenceGradient() const
  {
    return true;
  }
  
}

#define MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_INLINE_CXX
#endif
