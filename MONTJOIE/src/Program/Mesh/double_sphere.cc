#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 5)
    {
      cout<<"Not enough arguments"<<endl;
      return -1;
    }
  
  // origin is the center
  R3 center;
  Real_wp a = to_num<Real_wp>(argv[1]);
  Real_wp b = to_num<Real_wp>(argv[2]);
  Real_wp h = to_num<Real_wp>(argv[3]);
  string nom_fichier(argv[4]);
  
  // on genere le maillage d'un demi-disque
  Mesh<Dimension2> mesh_s1, mesh_s2; Mesh<Dimension3> mesh1, mesh2;
  int nb_div_teta = int(ceil(pi_wp*a/h));
  Real_wp step_teta = pi_wp/nb_div_teta, teta;
  mesh_s1.ReallocateVertices(nb_div_teta+1);
  for (int i = 0; i <= nb_div_teta; i++)
    {
      teta = -0.5*pi_wp + i*step_teta;
      mesh_s1.Vertex(i)(0) = a*cos(teta);
      mesh_s1.Vertex(i)(1) = a*sin(teta);
    }
      
  mesh_s1.ReallocateBoundariesRef(nb_div_teta);
  for (int i = 0; i < nb_div_teta; i++)
    mesh_s1.BoundaryRef(i).Init(i, i+1, 1);
  
  // et on fait tourner, pour avoir la sphere complete
  IVect Index_FaceSurf_to_EdgeRef, Index_VertexSurf_to_Vertex;
  VectReal_wp AngleVertex;
  IVect ref_cond(mesh_s1.GetNbReferences()+1); ref_cond.Fill();
  mesh_s1.GenerateSurfaceOfRevolution(1, mesh1, Index_FaceSurf_to_EdgeRef,
				      Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
  
  // on decale la sphere
  for (int i = 0; i < mesh1.GetNbVertices(); i++)
    Add(1.0, center, mesh1.Vertex(i));

  nb_div_teta = int(ceil(pi_wp*b/h));
  step_teta = pi_wp/nb_div_teta;
  mesh_s2.ReallocateVertices(nb_div_teta+1);
  for (int i = 0; i <= nb_div_teta; i++)
    {
      teta = -0.5*pi_wp + i*step_teta;
      mesh_s2.Vertex(i)(0) = b*cos(teta);
      mesh_s2.Vertex(i)(1) = b*sin(teta);
    }
      
  mesh_s2.ReallocateBoundariesRef(nb_div_teta);
  for (int i = 0; i < nb_div_teta; i++)
    mesh_s2.BoundaryRef(i).Init(i, i+1, 3);
  
  // et on fait tourner, pour avoir la sphere complete
  ref_cond.Reallocate(mesh_s2.GetNbReferences()+1);
  ref_cond.Fill();
  mesh_s2.GenerateSurfaceOfRevolution(3, mesh2, Index_FaceSurf_to_EdgeRef,
				      Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
      
  // on decale la sphere
  for (int i = 0; i < mesh2.GetNbVertices(); i++)
    Add(1.0, center, mesh2.Vertex(i));

  
  mesh1.AppendMesh(mesh2, false);
  mesh1.Write(nom_fichier);
  
  return 0;
}
