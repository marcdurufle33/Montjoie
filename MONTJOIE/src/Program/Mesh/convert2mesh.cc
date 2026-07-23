#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout<<"Cette commande demande au moins deux arguments"<<endl;
      cout<<"convert_mesh toto.vf toto.mesh"<<endl;
      cout<<"the file toto.mesh is created"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  string file_output(argv[2]);
  string extension_input = GetExtension(file_input);
  string extension_output = GetExtension(file_output);
  int r = 1;
  if (argc > 3)
    r = atoi(argv[3]);

  bool ascii = true, double_prec = true;
  if (argc > 4)
    {
      if (string(argv[4]) == "binary")
	ascii = false;
    }

  if (argc > 5)
    {
      if (string(argv[5]) == "float")
	double_prec = false;
    }
  
  int dimension = GetDimensionMesh(file_input, extension_input);
  DISP(dimension);
  if (dimension == 2)
    {
      Mesh<Dimension2> mesh;
      mesh.SetGeometryOrder(r);
      
      // on lit le maillage
      mesh.Read(file_input);
      
      // on ecrit maillage
      mesh.Write(file_output, double_prec, ascii);
    }
  else
    {
      Mesh<Dimension3> mesh;
      mesh.SetGeometryOrder(r);
      
      mesh.print_level = 3;
      // on lit le maillage
      mesh.Read(file_input);
      
      // on ecrit maillage
      mesh.Write(file_output, double_prec, ascii);
    }

  return FinalizeMontjoie();
}  

