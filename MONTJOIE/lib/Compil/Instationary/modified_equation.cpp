#include "MontjoieFlag.hxx"

#include "Instationary/MontjoieTimeHeader.hxx"
#include "Instationary/MontjoieTimeInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Instationary/ModifiedEquation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void RunSecondOrderScheme(const Real_wp&, const Real_wp&, const Real_wp&,
						   VectReal_wp&, VectReal_wp&, VirtualOdeSystem<Real_wp>&, VirtualTimeScheme<Real_wp>&);

  SELDON_EXTERN template void RunSecondOrderScheme(const Real_wp&, const Real_wp&, const Real_wp&,
					           VectComplex_wp&, VectComplex_wp&, VirtualOdeSystem<Complex_wp>&, VirtualTimeScheme<Complex_wp>&);

  SELDON_EXTERN template void RunFirstOrderScheme(const Real_wp&, const Real_wp&, const Real_wp&,
					          VectReal_wp&, VectReal_wp&, VirtualOdeSystem<Real_wp>&, VirtualTimeScheme<Real_wp>&);

  SELDON_EXTERN template void RunFirstOrderScheme(const Real_wp&, const Real_wp&, const Real_wp&,
					          VectComplex_wp&, VectComplex_wp&, VirtualOdeSystem<Complex_wp>&, VirtualTimeScheme<Complex_wp>&);

  SELDON_EXTERN template class ModifiedEquationIterator<Real_wp>;
  SELDON_EXTERN template class ModifiedEquationIterator<Complex_wp>;

  SELDON_EXTERN template class ModifiedEquationSystemIterator<Real_wp>;
  SELDON_EXTERN template class ModifiedEquationSystemIterator<Complex_wp>;

  SELDON_EXTERN template class OptimalModifiedEquationSecond_Iterator<Real_wp>;
  SELDON_EXTERN template class OptimalModifiedEquationSecond_Iterator<Complex_wp>;

  SELDON_EXTERN template class OptimalModifiedEquation_Iterator<Real_wp>;
  SELDON_EXTERN template class OptimalModifiedEquation_Iterator<Complex_wp>;
  
}
