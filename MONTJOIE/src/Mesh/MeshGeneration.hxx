#ifndef MONTJOIE_FILE_MESH_GENERATION_CXX

namespace Montjoie
{

  void CreateTetrahedralMesh(Mesh<Dimension3>& mesh_tri, Mesh<Dimension3>& mesh_tet);
  
  void CreatePyramidalLayer(const Mesh<Dimension3>& mesh_surf, Mesh<Dimension3>& mesh_tri,
			    Mesh<Dimension3>& mesh_pyramid, const Real_wp& height_pyramid,
                            int new_ref);

  void CreateSquareMesh(int nb_points_x, int nb_points_y, const R3& ptA, const R3& ptB,
                        const R3& ptC, const R3& ptD, Mesh<Dimension3>& mesh);
  
  void CreateTriangularMesh(int n1, int n2, int n3, int n4,
			    const R3& ptA, const R3& ptB, const R3& ptC,
                            const R3& ptD, Mesh<Dimension3>& mesh);
  
}

#define MONTJOIE_FILE_MESH_GENERATION_CXX
#endif

