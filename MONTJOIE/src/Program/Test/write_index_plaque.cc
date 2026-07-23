#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"
#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

typedef Real_wp Complexe;

void ComputeIndex(const Real_wp& x, const Real_wp& y, int ref, int n, Complexe& delta)
{
  switch (ref)
    {
    case 1: delta = 0.002*exp(-(x*x+y*y)); break;
    case 2 : delta = 0.005; break;
    case 5 : delta = 0.01*exp(-0.3*x*x); break;
    case 6 : delta = 0.007*exp(-0.5*y*y); break;
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout << "Usage : write_index nodal_points.dat maillage.mesh destination" << endl;
      abort();
    }
  
  string name_source(argv[1]);
  Vector<string> name_dest(argc-3);
  for (int i = 3; i < argc; i++)
    name_dest(i-3) = string(argv[i]);
  
  string name_mesh(argv[2]);
  Mesh<Dimension2> mesh;
  mesh.Read(name_mesh);
  
  ifstream file_in(name_source.data());
  int nb_elt, nb_nodes;
  file_in >> nb_elt >> nb_nodes;
  // on lit tous les points
  Vector<R2> points(nb_nodes);
  for (int i = 0; i < nb_nodes; i++)
    file_in >> points(i)(0) >> points(i)(1);
  
  if (mesh.GetNbElt() != nb_elt)
    {
      cout << "Nombre d'elements different entre le maillage et nodal_points.dat" << endl;
      abort();
    }
  
  // on lit les numeros des ddls
  Vector<IVect> Nodle(nb_elt);
  for (int i = 0; i < nb_elt; i++)
    {
      int nb_nodes_elt, num_node;
      file_in >> nb_nodes_elt;
      Nodle(i).Reallocate(nb_nodes_elt);
      for (int j = 0; j < nb_nodes_elt; j++)
        {
          file_in >> num_node;
          Nodle(i)(j) = num_node;
        }
    }
  
  file_in.close();
  
  // boucle sur tous les fichiers de sortie
  for (int j = 0; j < name_dest.GetM(); j++)
    {
      // on remplit l'indice par element
      Vector<Vector<Complexe> > nu_elt;
      nu_elt.Reallocate(nb_elt);
      
      for (int i = 0; i < nb_elt; i++)
        {
          int ref = mesh.Element(i).GetReference();
          nu_elt(i).Reallocate(Nodle(i).GetM());
          for (int j = 0; j < Nodle(i).GetM(); j++)
            {
              int num_node = Nodle(i)(j);
              ComputeIndex(points(num_node)(0), points(num_node)(1), ref, j, nu_elt(i)(j));
            }
        }
      
      WriteMeshData(nu_elt, name_dest(j));
    }
  
  return FinalizeMontjoie();
}
