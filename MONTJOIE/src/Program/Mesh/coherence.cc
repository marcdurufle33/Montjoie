#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc != 3)
    {
      cout<<"Cette commande demande deux arguments"<<endl;
      cout<<"coherence toto.msh toto.mesh"<<endl;
      cout<<"the file toto.mesh is created"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  string extension_input = GetExtension(file_input);
  string extension_output = GetExtension(file_output);
  
  int dimension = GetDimensionMesh(file_input, extension_input);
  if (dimension == 2)
    {
      Mesh<Dimension2> mesh;
      
      // on lit le maillage
      mesh.Read(file_input);
      //DISP(mesh.GetNbVertices());
      
      // on force la coherence
      mesh.ForceCoherenceMesh();
      
      // on ecrit maillage
      mesh.Write(file_output);
    }
  else
    {
      Mesh<Dimension3> mesh;
      
      // on lit le maillage
      mesh.Read(file_input);
      
      // on force la coherence
      mesh.ForceCoherenceMesh();
      
      // on ecrit maillage
      mesh.Write(file_output);
    }
  
  return FinalizeMontjoie();
}  

