#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/Mesh2D.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class Mesh<Dimension2>;

  SELDON_EXTERN template ostream& operator <<(ostream& out, const Mesh<Dimension2, Real_wp>& mesh);
  
  SELDON_EXTERN template
  void ExtrudePMLLayer(Mesh<Dimension2, Real_wp>& mesh, const TinyVector<Real_wp, 2>& axis,
		       const Real_wp& scal_max, const Vector<int>& ref_cond, int ref_target,
		       const Real_wp& thickness, int nb_layers, bool);
  
}
