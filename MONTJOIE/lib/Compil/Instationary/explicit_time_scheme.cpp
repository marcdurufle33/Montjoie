#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieUnsteadyHeader.hxx"
#include "Instationary/MontjoieUnsteadyInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/ExplicitTimeSchemes.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void RunTimeScheme(const Real_wp&, const Real_wp&, const Real_wp&,
					    VectReal_wp&, VirtualOdeSystem<Real_wp>&, VirtualTimeScheme<Real_wp>&, int);

  SELDON_EXTERN template void RunTimeScheme(const Real_wp&, const Real_wp&, const Real_wp&,
					    VectComplex_wp&, VirtualOdeSystem<Complex_wp>&, VirtualTimeScheme<Complex_wp>&, int);

  SELDON_EXTERN template class RungeKutta_Iterator<Real_wp>;
  SELDON_EXTERN template class RungeKutta_Iterator<Complex_wp>;

  SELDON_EXTERN template class LinearRungeKutta_Iterator<Real_wp>;
  SELDON_EXTERN template class LinearRungeKutta_Iterator<Complex_wp>;

  SELDON_EXTERN template class Talezer_Iterator<Real_wp>;
  SELDON_EXTERN template class Talezer_Iterator<Complex_wp>;
  
  SELDON_EXTERN template class AdamsBashforth_Moulton_Iterator<Real_wp>;
  SELDON_EXTERN template class AdamsBashforth_Moulton_Iterator<Complex_wp>;

  SELDON_EXTERN template class TaylorSeries_Iterator<Real_wp>;
  SELDON_EXTERN template class TaylorSeries_Iterator<Complex_wp>;

  SELDON_EXTERN template class LowStorageRK_Iterator<Real_wp>;
  SELDON_EXTERN template class LowStorageRK_Iterator<Complex_wp>;

  SELDON_EXTERN template class MultiStepButcher_Iterator<Real_wp>;
  SELDON_EXTERN template class MultiStepButcher_Iterator<Complex_wp>;

  SELDON_EXTERN template class SplitScheme_Iterator<Real_wp>;
  SELDON_EXTERN template class SplitScheme_Iterator<Complex_wp>;

  SELDON_EXTERN template class Nystrom_Iterator<Real_wp>;
  SELDON_EXTERN template class Nystrom_Iterator<Complex_wp>;

  SELDON_EXTERN template class RungeKuttaNystrom_Iterator<Real_wp>;
  SELDON_EXTERN template class RungeKuttaNystrom_Iterator<Complex_wp>;

  SELDON_EXTERN template class SymmetricMultistep_Iterator<Real_wp>;
  SELDON_EXTERN template class SymmetricMultistep_Iterator<Complex_wp>;

}
