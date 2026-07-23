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
      cout<<"texture_vtk toto.vtk coef"<<endl;
      cout<<"the file toto.vtk is scaled"<<endl;
      return -1;
    }
  
  string file_input(argv[1]);
  Real_wp coef = atof(argv[2]);

  bool ascii = false, double_prec = false;
  string ext = "vtk";
  
  Mesh<Dimension3> mesh;
  ifstream file_in(file_input.data());
  
  if (!file_in.is_open())
    {
      cout<<"Mesh file not found"<<endl;
      cout<<file_input<<endl;
      abort();
    }

  // on lit le maillage
  mesh.Read(file_in, ext);
  
  string chaine; int n;
  file_in >> chaine >> n;
  string name; VectReal_wp x_sol(n);
  ReadVtk(x_sol, name, file_in, ascii);
  
  file_in.close();  
  
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    mesh.Vertex(i)(2) *= coef;
  
  ofstream file_out(file_input.data());
  
  mesh.Write(file_out, ext, double_prec, ascii);
  
  file_out << "POINT_DATA " << n << endl;  
  //WriteVtk(x_sol, name, file_out, double_prec, ascii);
  
  file_out << "TEXTURE_COORDINATES coor 2 float" << endl;
  Vector<float> texture(2*n);
  Real_wp theta = -54.0*pi_wp/180.0;
  for (int i = 0; i < n; i++)
    {
      Real_wp x2 = 0.9*mesh.Vertex(i)(0);
      Real_wp y2 = 0.1+0.5*mesh.Vertex(i)(1);
      Real_wp x = cos(theta) * x2 - sin(theta)*y2;
      Real_wp y = sin(theta) * x2 + cos(theta)*y2;
      // texture(2*i) = x - floor(x);
      // texture(2*i+1) = y - floor(y);
      texture(2*i) = x;
      texture(2*i+1) = -0.45+0.55*y;
    }
  
  WriteBinaryDoubleOrFloat(texture, file_out, double_prec, false, true);
  file_out << '\n';
  
  file_out.close();
  
  return FinalizeMontjoie();
}  

