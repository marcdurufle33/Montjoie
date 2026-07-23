#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_INLINE_CXX

namespace Montjoie
{

  //! Default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarTransmission_Elas_Base<Dimension>::VarTransmission_Elas_Base(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Base<Dimension>(var), var_problem(var), var_boundary(var)
  {    
    coef_ch = 1.0;
    coef_sh = 1.0;
    curvilinear_case = false;
  }
  

  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarTransmission_Elas<Complexe, Dimension>::VarTransmission_Elas(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Elas_Base<Dimension>(var), var_problem(var), var_boundary(var)
  {
    rho_layer = 1.0;
    lambda_layer = 1.0;
    mu_layer = 1.0;
    sigma_layer = 0.0;
  }
  

  template<class Dimension> template<class TypeEquation>
  inline VarTransmission<HarmonicElasticEquation<Dimension> > ::VarTransmission(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Elas<Complex_wp, Dimension>(var)
  {
  }


  template<class Dimension> template<class TypeEquation>
  inline VarTransmission<ElasticEquation<Dimension> > ::VarTransmission(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Elas<Real_wp, Dimension>(var)
  {
  }

  
  template<class Complexe, class Dimension>
  inline bool ImpedanceTransmission_Elas<Complexe, Dimension>::PresenceGradient() const
  {
    return true;
  }
  
}

#define MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_INLINE_CXX
#endif

