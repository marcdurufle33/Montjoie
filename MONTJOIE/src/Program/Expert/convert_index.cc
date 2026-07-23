#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

template<class Dimension>
void ConvertMeditData(int order, Mesh<Dimension>& mesh, MeshNumbering<Dimension>& mesh_num,
                      const string& name_mesh, const string& name_source, const string& name_dest)
{
  mesh_num.SetOrder(order);
  mesh.SetGeometryOrder(order);
  mesh.Read(name_mesh);
  
  VectReal_wp nu;
  ReadMedit(nu, name_source);
  
  mesh_num.number_map.SetOrder(mesh, order);
  mesh_num.NumberMesh();
  
  Vector<VectReal_wp> nu_elt(mesh.GetNbElt());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int nb_nodes = mesh_num.Element(i).GetNbDof();
      nu_elt(i).Reallocate(nb_nodes);
      for (int j = 0; j < nb_nodes; j++)
        nu_elt(i)(j) = nu(mesh_num.Element(i).GetNumberDof(j));            
    }

  WriteMeshData(nu_elt, name_dest);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 5)
    {
      cout << "Usage : convert_index mesh order source destination" << endl;
      abort();
    }
  
  string name_mesh(argv[1]), name_source(argv[3]), name_dest(argv[4]);
  int order = atoi(argv[2]);
  
  int dim = GetDimensionMesh(name_mesh);
  if (dim == 2)
    {
      Mesh<Dimension2> mesh;
      MeshNumbering<Dimension2> mesh_num(mesh);
      
      ConvertMeditData(order, mesh, mesh_num,
		       name_mesh, name_source, name_dest);
    }
  else
    {
      Mesh<Dimension3> mesh;
      MeshNumbering<Dimension3> mesh_num(mesh);
      
      ConvertMeditData(order, mesh, mesh_num,
		       name_mesh, name_source, name_dest);
    }
  
  return FinalizeMontjoie();
}
