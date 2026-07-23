#include "MontjoieFlag.hxx"

#include "Quadrature/MontjoieQuadratureHeader.hxx"
#include "Quadrature/MontjoieQuadratureInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/GaussLobatto.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class Globatto<Real_wp>;
  
  SELDON_EXTERN template void Globatto<Real_wp>::AffectPoints(const Vector<Real_wp>&);

  SELDON_EXTERN template void Globatto<Real_wp>::AffectWeights(const Vector<Real_wp>&);
  

}

