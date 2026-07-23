#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Mesh/MeshBoundaries.cxx"
#include "Mesh/Mesh2DBoundaries.cxx"
#include "Mesh/Mesh3DBoundaries.cxx"
#include "vector/Vector.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void Vector<Montjoie::PmlRegionParameter<Montjoie::Dimension2, Real_wp> >::Resize(size_t);
  SELDON_EXTERN template void Vector<Montjoie::PmlRegionParameter<Montjoie::Dimension3, Real_wp> >::Resize(size_t);


}

namespace Montjoie
{
  SELDON_EXTERN template class PmlRegionParameter_Base<Dimension2, Real_wp>;
  SELDON_EXTERN template class PmlRegionParameter_Base<Dimension3, Real_wp>;

  SELDON_EXTERN template class PmlRegionParameter<Dimension2, Real_wp>;
  SELDON_EXTERN template class PmlRegionParameter<Dimension3, Real_wp>;

  SELDON_EXTERN template class MeshBoundaries_Base<Dimension2, Real_wp>;
  SELDON_EXTERN template class MeshBoundaries_Base<Dimension3, Real_wp>;

  SELDON_EXTERN template void MeshBoundaries_Base<Dimension2, Real_wp>::SetBoundaryCondition(const Mesh<Dimension2>&);
  SELDON_EXTERN template void MeshBoundaries_Base<Dimension2, Real_wp>::SetBoundaryCondition(const Mesh<Dimension3>&);
  SELDON_EXTERN template void MeshBoundaries_Base<Dimension3, Real_wp>::SetBoundaryCondition(const Mesh<Dimension2>&);
  SELDON_EXTERN template void MeshBoundaries_Base<Dimension3, Real_wp>::SetBoundaryCondition(const Mesh<Dimension3>&);

  SELDON_EXTERN template class MeshBoundaries<Dimension2, Real_wp>;
  SELDON_EXTERN template class MeshBoundaries<Dimension3, Real_wp>;

  SELDON_EXTERN template void MeshBoundaries<Dimension2, Real_wp>::
  GetBoundaryMesh(int ref, Mesh<Dimension2>& mesh_surf, IVect& Index_EdgeSurf_to_EdgeRef,
		  IVect& Liste_Vertices, IVect&, IVect&,
		  const IVect& ref_cond, int ref_domain) const;

  SELDON_EXTERN template void MeshBoundaries<Dimension3, Real_wp>::
  GetBoundaryMesh(int ref, Mesh<Dimension3>& mesh_surf, IVect& Index_EdgeSurf_to_EdgeRef,
		  IVect& Liste_Vertices, IVect&, IVect&,
		  const IVect& ref_cond, int ref_domain) const;
  
}
