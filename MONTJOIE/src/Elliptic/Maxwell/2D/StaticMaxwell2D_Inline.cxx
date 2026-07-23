#ifndef MONTJOIE_FILE_STATIC_MAXWELL_2D_INLINE_CXX

namespace Montjoie
{
  
  
  /****************************
   * StaticMaxwellEquation_2D *
   ****************************/
  
  
  //! penalization coefficient
  template<class TypeEquation, class T0>
  inline void StaticMaxwellEquation_2D::
  GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
		      int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef)
  {
    coef = var.delta_penalization*nat_mat.GetCoefDamping()
      /var.coefficient_impedance_absorbing(ref);
  }
  
  
  /*******************************
   * StaticMaxwellEquation_2D_DG *
   *******************************/

  
  //! not used
  template<class TypeEquation>
  inline void StaticMaxwellEquation_2D_DG
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int i, const ElementReference_Dim<Dimension>&)
  {
  }
  

  /****************
   * ImpedanceABC *
   ****************/
  

  template<class T> template<class TypeEquation>
  inline ImpedanceABC<T, StaticMaxwellEquation_2D>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, Dimension2>(var), var_problem(var), var_maxwell(var)
  {
  }
  
   
}

#define MONTJOIE_FILE_STATIC_MAXWELL_2D_INLINE_CXX
#endif

