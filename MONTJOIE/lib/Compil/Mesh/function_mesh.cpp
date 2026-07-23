#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/FunctionsMesh.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void WriteElementMesh(const Mesh<Dimension2>&, const ElementGeomReference<Dimension2>&, const SetPoints<Dimension2>&, const string&, int);
  SELDON_EXTERN template void WriteElementMesh(const Mesh<Dimension3>&, const ElementGeomReference<Dimension3>&, const SetPoints<Dimension3>&, const string&, int);

}
