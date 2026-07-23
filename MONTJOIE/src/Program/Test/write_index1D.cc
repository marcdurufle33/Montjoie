#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

typedef Complex_wp Complexe;

void ComputeIndex(const Real_wp& x, Vector<Complexe>& coef)
{
  // on evalue tous les indices et on remplit le tableau coef
  Real_wp alpha = 1000.0;
  Real_wp c = 8.0e3;
  Real_wp gamma = 1.25e4;
  Real_wp rho = exp(-alpha*(x-1.0));
  // coef(0) = rho_tilde = 1/(rho c^2)
  coef(0) = 1.0/(rho*c*c);
  
  // coef(1) = mu_tilde = 1 / rho
  coef(1) = 1.0/rho;
  
  // coef(2) = 2.0 gamma / (rho c^2)
  coef(2) = 2.0*gamma / (rho*c*c);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Entrez le nom du fichier nodal_points.dat et les fichiers de sortie" << endl;
      cout << "Usage : ./write.x nodal_points.dat rho.don mu.don sigma.don" << endl;
      abort();
    }
  
  string file_nodal(argv[1]);
  int nb_index = argc-2;
  Vector<string> name_dest(nb_index);
  for (int j = 0; j < nb_index; j++)
    name_dest(j) = string(argv[2+j]);
  
  Vector<Complexe> coef(nb_index);
  Vector<Vector<Complexe> > all_coef(nb_index);

  VectReal_wp points;
  points.ReadText(file_nodal);

  int nb_points = points.GetM();
  for (int j = 0; j < nb_index; j++)
    all_coef(j).Reallocate(nb_points);
  
  for (int i = 0; i < points.GetM(); i++)
    {
      ComputeIndex(points(i), coef);
      for (int j = 0; j < nb_index; j++)
        all_coef(j)(i) = coef(j);
    }
  
  for (int j = 0; j < nb_index; j++)
    all_coef(j).WriteText(name_dest(j));
  
  return FinalizeMontjoie();
}
