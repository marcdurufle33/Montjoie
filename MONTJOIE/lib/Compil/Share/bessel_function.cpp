#include "MontjoieFlag.hxx"

#include "Share/MontjoieCommonHeader.hxx"
#include "Share/MontjoieCommonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Share/BesselFunctionsInterface.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template int ComputeOrder(const Real_wp&, const Real_wp&);
  
}
