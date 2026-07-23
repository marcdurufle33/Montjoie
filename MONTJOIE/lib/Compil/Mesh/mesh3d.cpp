#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/Mesh3D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class Mesh<Dimension3>;

  SELDON_EXTERN template ostream& operator <<(ostream& out, const Mesh<Dimension3, Real_wp>& mesh);
  
}
