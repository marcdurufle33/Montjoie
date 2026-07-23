#include "MontjoieFlag.hxx"

#include "Elliptic/Elastic/MontjoieElasticHeader.hxx"
#include "Elliptic/Elastic/MontjoieElasticInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Elastic/TransmissionModelElastic.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class VarTransmission_Elas_Base<Dimension2>;
  SELDON_EXTERN template class VarTransmission_Elas_Base<Dimension3>;

  SELDON_EXTERN template class VarTransmission_Elas<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarTransmission_Elas<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarTransmission_Elas<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarTransmission_Elas<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarTransmission<HarmonicElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarTransmission<HarmonicElasticEquation<Dimension3> >;
  SELDON_EXTERN template class VarTransmission<ElasticEquation<Dimension2> >;
  SELDON_EXTERN template class VarTransmission<ElasticEquation<Dimension3> >;

  SELDON_EXTERN template class ImpedanceTransmission_Elas_Base<Dimension2>;
  SELDON_EXTERN template class ImpedanceTransmission_Elas_Base<Dimension3>;

  SELDON_EXTERN template class ImpedanceTransmission_Elas<Real_wp, Dimension2>;
  SELDON_EXTERN template class ImpedanceTransmission_Elas<Real_wp, Dimension3>;
  SELDON_EXTERN template class ImpedanceTransmission_Elas<Complex_wp, Dimension2>;
  SELDON_EXTERN template class ImpedanceTransmission_Elas<Complex_wp, Dimension3>;

}
