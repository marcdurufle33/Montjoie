#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwell3D_Header.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D_Inline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/3D/ComputationModeMaxwell.cxx"
#include "Computation/ProdMatVectHcurl2D.cxx"
#include "Computation/ElementaryMatrixHcurl2D.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarHarmonic<ModeMaxwellEquation<Real_wp> >;
  SELDON_EXTERN template class VarModeMaxwell<ModeMaxwellEquation<Real_wp> >;
  SELDON_EXTERN template class VarHarmonic<ModeMaxwellEquation<Complex_wp> >;
  SELDON_EXTERN template class VarModeMaxwell<ModeMaxwellEquation<Complex_wp> >;

  SELDON_EXTERN template class ModalSourceMaxwell3D<Real_wp>;
  SELDON_EXTERN template class ModalSourceMaxwell3D<Complex_wp>;
  
}
