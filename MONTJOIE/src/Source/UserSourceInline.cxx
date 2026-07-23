#ifndef MONTJOIE_FILE_USER_SOURCE_INLINE_CXX

namespace Montjoie
{

  //! constructor with given problem
  template<class T, class Dimension> template<class TypeEquation>
  inline UserDefinedSource<T, Dimension>::UserDefinedSource(const EllipticProblem<TypeEquation>& var, const Vector<VectString>& param)
    : VirtualSourceFEM<T, Dimension>(var), source_space_param(param(0))
  {
  }  
  
  
  //! default constructor
  template<class T, class Dimension> template<class TypeEquation>
  inline InitialUserFunction<T, Dimension>
  ::InitialUserFunction(const EllipticProblem<TypeEquation>& var, const Vector<string>& param)
    : VirtualProjectorFEM<T, Dimension>(var), param_initial_condition(param)
  {
  }
  
}

#define MONTJOIE_FILE_USER_SOURCE_INLINE_CXX
#endif
