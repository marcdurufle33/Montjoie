#include "MontjoieFlag.hxx"

#include "Elliptic/Helmholtz/MontjoieHelmholtzHeader.hxx"
#include "Elliptic/Helmholtz/MontjoieHelmholtzInline.hxx"

#include "Elliptic/Helmholtz/CellulairePeriodique.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/CellulairePeriodique.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarHelmholtz_CellulairePeriodique<HelmholtzEquation_CellulairePeriodique<Dimension2> >;
  SELDON_EXTERN template class VarHelmholtz_CellulairePeriodique<HelmholtzEquation_CellulairePeriodique<Dimension3> >;

  SELDON_EXTERN template class EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension2> >;
  SELDON_EXTERN template class EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension3> >;

  SELDON_EXTERN template class MuSource<Dimension2>;
  SELDON_EXTERN template class MuSource<Dimension3>;

  SELDON_EXTERN template MuSource<Dimension2>::MuSource(const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension2> >&);
  SELDON_EXTERN template MuSource<Dimension3>::MuSource(const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension3> >&);

  SELDON_EXTERN template class SourceAbsorbante<Dimension2>;
  SELDON_EXTERN template class SourceAbsorbante<Dimension3>;

  SELDON_EXTERN template SourceAbsorbante<Dimension2>::SourceAbsorbante(const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension2> >&);
  SELDON_EXTERN template SourceAbsorbante<Dimension3>::SourceAbsorbante(const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension3> >&);

}
