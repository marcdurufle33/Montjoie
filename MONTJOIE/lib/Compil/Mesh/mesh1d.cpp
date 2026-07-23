#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/Mesh1D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class Mesh<Dimension1>;  
}
