#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

typedef Real_wp Complexe;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout << "Usage : regularize_data nodal_points.dat source destination" << endl;
      abort();
    }
  
  string name_nodal(argv[1]);
  string name_source(argv[2]);
  string name_dest(argv[3]);
  
  ifstream file_in(name_nodal.data());
  int nb_elt, nb_nodes;
  file_in >> nb_elt >> nb_nodes;
  // on lit tous les points
  for (int i = 0; i < nb_nodes; i++)
    {
      R2 point;
      file_in >> point(0) >> point(1);
    }
  
  // on lit tous les noeuds
  Vector<IVect> Nodle(nb_elt);
  for (int i = 0; i < nb_elt; i++)
    {
      int nb_nodes_elt = 0;
      file_in >> nb_nodes_elt;
      Nodle(i).Reallocate(nb_nodes_elt);
      for (int j = 0; j < nb_nodes_elt; j++)
        file_in >> Nodle(i)(j);
    }
  
  file_in.close();
  
  // on lit la source
  Vector<Vector<Complexe> > nu_elt;
  ReadMeshData(nu_elt, name_source);
  if (nu_elt.GetM() != nb_elt)
    {
      cout << "Source not compatible with given nodal points" << endl;
      abort();
    }
  
  // on forme la moyenne sur les noeuds
  Vector<Complexe> nu_node(nb_nodes);
  nu_node.Fill(0);
  Vector<int> nb_eval(nb_nodes);
  nb_eval.Zero();
  for (int i = 0; i < nb_elt; i++)
    {
      if (nu_elt(i).GetM() != Nodle(i).GetM())
        {
          cout << "Source not compatible with given nodal points" << endl;
          abort();
        }
      
      for (int j = 0; j < Nodle(i).GetM(); j++)
        {
          nu_node(Nodle(i)(j)) += nu_elt(i)(j);
          nb_eval(Nodle(i)(j))++;
        }
    }
  
  for (int i = 0; i < nb_nodes; i++)
    nu_node(i) /= Real_wp(nb_eval(i));
  
  // ensuite on redistribue les valeurs sur les elements
  for (int i = 0; i < nb_elt; i++)
    for (int j = 0; j < Nodle(i).GetM(); j++)
      nu_elt(i)(j) = nu_node(Nodle(i)(j));
  
  WriteMeshData(nu_elt, name_dest);
  
  return FinalizeMontjoie();
}
