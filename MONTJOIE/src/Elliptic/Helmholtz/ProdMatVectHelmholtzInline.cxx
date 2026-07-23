#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{

  template<class Dimension>
  inline FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<HelmholtzEquation<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquation<Dimension> >(var_)
  {
  }


  template<class Dimension>
  inline FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension> >
  ::FemMatrixFreeClass(const EllipticProblem<HelmholtzEquationDG<Dimension> >& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquationDG<Dimension> >(var_)
  {
  }
  
}

#define MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_INLINE_CXX
#endif
