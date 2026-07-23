#ifndef MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline ModalSourceMaxwell3D_Base
  ::ModalSourceMaxwell3D_Base(const EllipticProblem<TypeEquation>& var,
			      const ModalSourceBoundary_Base& mode)
    : var_problem(var), modal_source(mode)
  {
  }


  template<class T> template<class TypeEquation>
  inline ModalSourceMaxwell3D<T>::ModalSourceMaxwell3D(const EllipticProblem<TypeEquation>& var)
    : ModalSourceBoundary_Dim<T, Dimension3>(var), ModalSourceMaxwell3D_Base(var, *this),
      var_problem(var), var_maxwell(var)
  {
    SetComplexOne(voltage);
  }
  
}

#define MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_INLINE_CXX
#endif

  
