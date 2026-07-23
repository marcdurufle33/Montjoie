//pour la corde
#define MONTJOIE_WITH_ONE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Corde/MontjoieString.hxx"

using namespace Montjoie;

template<int nb_base, int nb_quad, class TypeEquation>
void RunAll(Shank<nb_base, nb_quad, TypeEquation>& var, const string& input_file)
{
  typedef Montjoie::Real_wp Real_wp;

  ReadInputFile(input_file, var);
  
  DISP(var.L);
  var.ConstructAll(0);
  
  int N = var.GetM();
  DISP(N);
  Matrix<Real_wp> val(N, 3);
  Vector<Real_wp> scheme(N);

  int nb_iterations = toInteger(var.t_end/var.Deltat);
  DISP(nb_iterations);
  
  var.Init(nb_iterations, val);

  

  // boucle en temps
  DISP(var.DOSSIER);
  WriteOnTheGoWithTinyBuffer<Real_wp, 4> output_energy;
  output_energy.Init(var.DOSSIER+"energy.dat", 10, true);
  TinyVector<Real_wp, 4> tmp;
  for (int nt = 0; nt < nb_iterations; nt++)
    {
      Real_wp t = var.Deltat*nt;
      var.t_courant = t;

      var.WriteSnapshots(nt, t, val);      
      var.Advance(nt, t, val);
      //var.ImposeDirichlet(val);
      
      int n = toInteger(round(t/(var.t_display)));
      if ((abs(t/(var.t_display) - n) <var.Deltat/2))
        {
          Real_wp e_cin, e_pot, e_int;
          Real_wp energy = var.GetEnergy(val, e_cin, e_pot, e_int);
          cout << "At time : " << t << ", ";
          cout << "theta = " << val(var.teta_dof_number, 0) << endl;
          cout << "energy = " << energy << endl;
          tmp(0) = e_cin; tmp(1) = e_pot; tmp(2) = e_int; tmp(3) = energy;
          output_energy.AddTinyVect(tmp);
        }
    }
  
  output_energy.CloseBuffer();
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  srand(0);
  
  if (argc < 2)
    {
      cout << "entrez le nom du fichier de donnees" << endl;
      abort();
    }
  
  Hammer hammer;
  Shank<5, 5, WaveEquationStiffString> var(hammer);  
  string input_file(argv[1]);
  
  RunAll(var, input_file);
  
  return FinalizeMontjoie();
}
