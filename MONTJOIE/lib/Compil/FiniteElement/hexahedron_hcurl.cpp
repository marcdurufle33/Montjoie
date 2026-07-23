#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/Hexahedron/HexahedronHcurlFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlHpFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlOptimalHpFirstFamily.cxx"
#include "FiniteElement/Hexahedron/HexahedronHcurlLobatto.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void HexahedronHcurlLobatto::ApplyRhGen(const VectReal_wp&, VectReal_wp&) const;  
  SELDON_EXTERN template void HexahedronHcurlLobatto::ApplyRhGen(const VectComplex_wp&, VectComplex_wp&) const;  

  SELDON_EXTERN template void HexahedronHcurlLobatto::ApplyRhTransposeGen(const VectReal_wp&, VectReal_wp&) const;  
  SELDON_EXTERN template void HexahedronHcurlLobatto::ApplyRhTransposeGen(const VectComplex_wp&, VectComplex_wp&) const;  

}
