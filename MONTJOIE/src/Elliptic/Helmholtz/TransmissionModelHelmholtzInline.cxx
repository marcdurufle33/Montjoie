#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{
  
  //! Default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarTransmission_Helm<Complexe, Dimension>
  ::VarTransmission_Helm(EllipticProblem<TypeEquation>& var)
    : VarTransmission_Base<Dimension>(var), var_problem(var), var_boundary(var),
      var_helm(var)
  {
    InitDefaultValues();
  }

  
  template<class Complexe, class Dimension>
  inline void VarTransmission_Helm<Complexe, Dimension>
  ::AddTransmissionTerms(const Real_wp& alpha,
                         const GlobalGenericMatrix<Real_wp>& nat_mat,
                         VirtualMatrix<Real_wp>& mat_sp, int, int)
  {
  }

  
  inline VarTransmission<HelmholtzEquation<Dimension2> > 
  ::VarTransmission(EllipticProblem<HelmholtzEquation<Dimension2> >& var)
    : VarTransmission_Helm<Complex_wp, Dimension2>(var)
  {
  }
  
}

#define MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_INLINE_CXX
#endif
