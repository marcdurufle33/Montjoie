#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MeshElement.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template bool operator==(const Face<Dimension2>&, const Face<Dimension2>&);
  SELDON_EXTERN template bool operator==(const Face<Dimension3>&, const Face<Dimension3>&);

  SELDON_EXTERN template bool operator<(const Face<Dimension2>&, const Face<Dimension2>&);
  SELDON_EXTERN template bool operator<(const Face<Dimension3>&, const Face<Dimension3>&);

  SELDON_EXTERN template bool operator>(const Face<Dimension2>&, const Face<Dimension2>&);
  SELDON_EXTERN template bool operator>(const Face<Dimension3>&, const Face<Dimension3>&);

}

