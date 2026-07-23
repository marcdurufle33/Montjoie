#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout<<"Cette commande demande au moints trois arguments"<<endl;
      cout<<"extract_boundary toto.vf surf.mesh 1"<<endl;
      cout<<"the file surf.mesh is created"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  string extension_input = GetExtension(file_input);
  string extension_output = GetExtension(file_output);
  
  Mesh<Dimension3> mesh;
  SurfacicMesh<Dimension3> mesh_surf;
  
  // on lit le maillage
  mesh.Read(file_input);
  int ref_max = 0;
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    ref_max = max(ref_max, mesh.BoundaryRef(i).GetReference());
  
  // on redirige toutes les ref demandees a une condition de 1
  IVect ref_cond(ref_max+1); ref_cond.Fill(0);
  for (int i = 3; i < argc; i++)
    {
      int ref = to_num<int>(string(argv[i]));
      if (ref <= ref_max)
	ref_cond(ref) = 1;
    }
  
  // on extrait le maillage surfacique
  mesh.GetBoundaryMesh(1, mesh_surf, ref_cond);
  
  // on ecrit maillage
  mesh_surf.Write(file_output);
  
  return 0;
}  

