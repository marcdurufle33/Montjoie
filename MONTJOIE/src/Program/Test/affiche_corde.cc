#define MONTJOIE_WITH_ONE_DIM

#include "Corde/MontjoieString.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout << "Entrez le plan de cordes a utiliser" << endl;
      abort();
    }
  
  string plan_cordes(argv[1]);
  ParameterMultistring param_corde_data;

  param_corde_data.InitPlanCordes(plan_cordes);
  
  Real_wp L, A, rho, T0, E, I, G, k_prime;
  R3 amo, amoB2;
  Real_wp masse_marteau, exposant, Ks, Rs;
  Real_wp x0,y0;
  cout.setf(ios::scientific);
  for (int n = 3; n < 88; n++)
    {
      string note = param_corde_data.GetNoteString(n);
      param_corde_data.GetParamNote(note, L, A, rho, T0, E, I, G, k_prime, amo, amoB2);
      param_corde_data.GetParamHammer(note, masse_marteau, exposant, Ks, Rs); 
      param_corde_data.GetParamAttache(note,x0,y0),

      cout << "Note " << note << endl;
      cout << "L = " << L << endl;
      cout << "A = " << A << endl;
      cout << "rho = " << rho << endl;
      cout << "T0 = " << T0 << endl;
      cout << "E = " << E << endl;
      cout << "I = " << I << endl;
      cout << "G = " << G << endl;
      cout << "k_prime = " << k_prime << endl;
      cout << "amo constant = " << amo << endl;
      cout << "amo quadratique = " << amoB2 << endl;
      cout << "masse_marteau = " << masse_marteau << endl;
      cout << "exposant = " << exposant << endl;
      cout << "K^H = " << Ks << endl;
      cout << "R^H = " << Rs << endl;
      cout << "x0 = " << x0 << endl;
      cout << "y0 = " << y0 << endl;
      cout << endl;
    }
  
  return FinalizeMontjoie();
}
