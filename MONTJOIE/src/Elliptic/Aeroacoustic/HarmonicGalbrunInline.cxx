#ifndef MONTJOIE_FILE_HARMONIC_GALBRUN_INLINE_CXX

namespace Montjoie
{

  template<class Dimension> template<class TypeEquation>
  inline VarGalbrunIndex_Base<Dimension>
  ::VarGalbrunIndex_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var)
  {
    adjustment_neumann = false;
    compute_div_flow = false;
    compute_grad_flow = false;
    compute_grad_rho = false;
    compute_grad_sigma = false;
    compute_grad_c0 = false;
    compute_hessian_flow = false;
    compute_hessian_p0 = false;
    compute_gamma = false;
    compute_grad_gamma = false;
    store_grad_rho0_c0 = false;
  }
  

  //! Default constructor for Galbrun equation
  template<class TypeEquation>
  inline VarGalbrun_Base::VarGalbrun_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var)
  {
    apply_convective_derivate_source = false;
    drop_unstable_terms = DROP_NONE;
    coef_convective_term = 1.0;
  };
  
  
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarGalbrun_Dim<Complexe, Dimension>
  ::VarGalbrun_Dim(EllipticProblem<TypeEquation>& var)
    : VarGalbrun_Base(var), VarGalbrunIndex_Base<Dimension>(var)
  {
  }


  template<class Dimension>
  inline FemMatrixFreeClass<Real_wp, GalbrunStationaryEquationDG<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<GalbrunStationaryEquationDG<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Real_wp, GalbrunStationaryEquationDG<Dimension> >(var_)
  {
  }

}

#define MONTJOIE_FILE_HARMONIC_GALBRUN_INLINE_CXX
#endif

