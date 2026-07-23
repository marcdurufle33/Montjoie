#include "MontjoieFlag.hxx"

#include "Elliptic/Maxwell/MontjoieMaxwellAxiHeader.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwellAxiInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Maxwell/Axi/DefineSourceMaxwellAxi.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VolumetricSource_MaxwellAxi<Real_wp>;
  SELDON_EXTERN template class VolumetricSource_MaxwellAxi<Complex_wp>;
}
