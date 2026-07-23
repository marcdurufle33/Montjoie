#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

typedef Real_wp Complexe;

void ComputeIndex(const Real_wp& x, const Real_wp& y, int n, Complexe& coef)
{
  /* if (n == 0)
    coef = sqrt(1.0-x*x/20-y*y/32);
  else if (n == 1)
    coef = exp(-x*x/16-y*y/12);
  else
    {
      //Real_wp r = sqrt(x*x+y*y);
      //Real_wp f = r*(4.0-r)/10;
      Real_wp rho = sqrt(1.0-x*x/20-y*y/32);
      rho *= 10.0;
      if (n == 2)
	coef = x/rho;
      else
	coef = y/rho;
	} */
  
  if (n == 0)
    coef = 1.5 + 0.2*cos(pi_wp*x/4)*sin(pi_wp*y/2);
  else if (n == 1)
    coef = 0.8 + 0.16*sin(3.0*pi_wp*x/4)*cos(pi_wp*y/4);
  else
    {
      Real_wp rho = 1.5 + 0.2*cos(pi_wp*x/4)*sin(pi_wp*y/2);
      //if (n == 2)
      //coef = y/rho;
      //else
      //coef = -x/rho;
      
      if (n == 2)
        coef = (0.3+0.1*cos(pi_wp*y/4))/rho;
      else
        coef = (0.2+0.08*sin(pi_wp*x/4))/rho;
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Usage : write_index source destination" << endl;
      abort();
    }
  
  string name_source(argv[1]);
  Vector<string> name_dest(argc-2);
  for (int i = 2; i < argc; i++)
    name_dest(i-2) = string(argv[i]);
  
  ifstream file_in(name_source.data());
  int nb_elt, nb_nodes;
  file_in >> nb_elt >> nb_nodes;
  // on lit tous les points
  Vector<R2> points(nb_nodes);
  for (int i = 0; i < nb_nodes; i++)
    file_in >> points(i)(0) >> points(i)(1);
  
  // on evalue l'indice sur tous les points (cas d'un indice continu)
  Vector<Vector<Complexe> > nu_node(name_dest.GetM());
  for (int j = 0; j < name_dest.GetM(); j++)
    {
      nu_node(j).Reallocate(nb_nodes);  
      for (int i = 0; i < nb_nodes; i++)
	ComputeIndex(points(i)(0), points(i)(1), j, nu_node(j)(i));
    }
  
  // ensuite on remplit l'indice par element
  Vector<Vector<Vector<Complexe> > > nu_elt(name_dest.GetM());
  for (int j = 0; j < name_dest.GetM(); j++)
    nu_elt(j).Reallocate(nb_elt);
  
  int nb_nodes_elt, num_node;
  for (int i = 0; i < nb_elt; i++)
    {
      file_in >> nb_nodes_elt;
      for (int j = 0; j < name_dest.GetM(); j++)
	nu_elt(j)(i).Reallocate(nb_nodes_elt);
      
      for (int j = 0; j < nb_nodes_elt; j++)
	{
          file_in >> num_node;
          for (int k = 0; k < name_dest.GetM(); k++)
	    nu_elt(k)(i)(j) = nu_node(k)(num_node);
	}
    }
  
  file_in.close();
  
  for (int j = 0; j < name_dest.GetM(); j++)
    WriteMeshData(nu_elt(j), name_dest(j));
  
  return FinalizeMontjoie();
}
