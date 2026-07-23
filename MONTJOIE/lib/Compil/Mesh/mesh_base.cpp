#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MeshBase.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class Mesh_Base<Dimension2>;
  SELDON_EXTERN template class Mesh_Base<Dimension3>;

}

