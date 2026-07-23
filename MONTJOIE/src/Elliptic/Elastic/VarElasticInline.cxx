#ifndef MONTJOIE_FILE_VAR_ELASTIC_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline VarElastic_Base::VarElastic_Base(EllipticProblem<TypeEquation>& var)
    : var_problem_base(var)
  {
    display_stress = true;
    omega_bar = 1.0;
    u_bar = 1.0;
  }  


  //! default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarElastic_Dim<Complexe, Dimension>::VarElastic_Dim(EllipticProblem<TypeEquation>& var)
    : VarElastic_Base(var), var_problem(var)
  {
  }
  

  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension>(var), var_problem(var), var_elas(var)
  {
  }


  //! constructor with a given problem
  template<class T, class Dimension> template<class TypeEquation>
  inline ImpedanceABC<T, ElasticEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, Dimension>(var), var_problem(var), var_elas(var)
  {
  }

  
  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<ElasticEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, ElasticEquation<Dimension> >(var_)
  {
  }
      

  template<class Dimension>
  inline FemMatrixFreeClass<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<HarmonicElasticEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HarmonicElasticEquation<Dimension> >(var_)
  {
  }
    
}

#define MONTJOIE_FILE_VAR_ELASTIC_INLINE_CXX
#endif
