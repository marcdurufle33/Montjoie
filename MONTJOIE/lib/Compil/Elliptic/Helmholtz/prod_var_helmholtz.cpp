#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#include "Elliptic/Helmholtz/MontjoieLaplaceHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieLaplaceInline.hxx"

#include "Elliptic/Helmholtz/CellulairePeriodique.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/ProdMatVectHelmholtz.cxx"
#include "Computation/ProdMatVectScalarH1.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension3> >;

  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension3> >;

  SELDON_EXTERN template class FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension3> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension3> >;

  SELDON_EXTERN template class FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension3> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension2> >;
  SELDON_EXTERN template class FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension3> >;

}
