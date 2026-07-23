#ifndef MONTJOIE_FILE_MODAL_SOURCE_LAPLACE_INLINE_CXX

namespace Montjoie
{

  template<class Dimension> template<class TypeEquation>
  inline ModalSourceLaplace_Dim<Dimension>
  ::ModalSourceLaplace_Dim(const EllipticProblem<TypeEquation>& var,
			   const ModalSourceBoundary_Base& mode)
    : var_problem(var), modal_source(mode)
  {
  }


  //! constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline ModalSourceLaplace<Complexe, Dimension>
  ::ModalSourceLaplace(const EllipticProblem<TypeEquation>& var)
    : ModalSourceBoundary_Dim<Complexe, Dimension>(var),
      ModalSourceLaplace_Dim<Dimension>(var, *this),
      var_problem(var), var_laplace(var)
  {
  }
    
}

#define MONTJOIE_FILE_MODAL_SOURCE_LAPLACE_INLINE_CXX
#endif
