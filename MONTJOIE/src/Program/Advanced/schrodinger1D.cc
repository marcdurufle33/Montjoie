#define MONTJOIE_WITH_ONE_DIM

#include "Harmonic/MontjoieHarmonic.hxx"
#include "Instationary/MontjoieTime.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.hxx"
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.hxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.cxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.cxx"
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.cxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  typedef Montjoie::Real_wp Real_wp;

  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout << "Fournissez un fichier de donnees " << endl;
      abort();
    }

  cout.setf(ios::scientific);
  string input_file(argv[1]);
  
  // on recupere le parametre d'adimensionalisation
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  // simulation principale
  KerrProblem var;
  
  ReadInputFile(input_file, var);
  
  // construction du numero de dossier
  // si le dossier n'est pas connu
  int num = -1;
  if (var.DOSSIER.size() == 0)
    {
      var.DOSSIER = "[STIFFOUT]/"; 
      EcritDossier(var.DOSSIER, input_file, num);
    }
  
  var.ConstructAll();
  
  int nb_iterations = toInteger(var.zmax/var.dz);
  string name_file = GetBaseString(var.output_grid_z.GetTotalFieldFile());

  // condition initiale
  VectComplex_wp E0;
  var.GetInitialCondition(E0); 
  var.time_scheme.SetInitialCondition(0.0, var.dz, E0, var);
  
  for (int nt = 0; nt <= nb_iterations; nt++)
    {
      glob_chrono.Start(VirtualTimer::ALL);
      
      Real_wp z = var.dz*nt;
      
      // on ecrit le snapshot
      if (var.output_grid_z.SnapshotToStore(z))
        {
          int num = var.output_grid_z.GetSnapshotNumber();
	  string name = name_file + NumberToString(num) + ".dat";
	  string name_time = name_file + NumberToString(num) + "_time.dat";
	  var.WriteOutputFile(name, name_time, z, var.time_scheme.GetIterate());
	  var.output_grid_z.IncrementSnapshot();
	}      
      
      // on donne une solution intermediaire a var
      var.GiveIterate(nt, z, var.time_scheme.GetIterate());
      
      // on avance le schema
      var.time_scheme.Advance(z, nt, var);
      
      // on modifie l'itere si necessaire
      var.ModifyIterate(nt, z+var.dz, var.time_scheme.GetIterate());
      
      if ( ((var.print_level >= 2) && (nt% 100 == 0))
           || (var.print_level >= 6))
	  {
	    cout << "At z = " << z*var.z0_adim << endl;
	    cout << "|| u || = " << var.GetEnergy(var.time_scheme.GetIterate()) << endl;
            
            if (var.print_level >= 4)
              {
                //cout.setf(ios::fixed);
                //cout << "Overal time spend in EvaluateFunction : " << glob_chrono.GetSeconds(VirtualTimer::JACOBIAN) << endl;
                cout << "Time spent for outputs : " << glob_chrono.GetSeconds(VirtualTimer::OUTPUT) << endl;
                cout << "Time spent to compute the phase : " << glob_chrono.GetSeconds(VirtualTimer::PML) << endl;
                cout << "Time spent in the FFTs : " << glob_chrono.GetSeconds(VirtualTimer::FLUX) << endl;
                cout << "Time spent in non-linear part : " << glob_chrono.GetSeconds(VirtualTimer::PROD) << endl;
                cout << "Time spent in computation of envelope : " << glob_chrono.GetSeconds(VirtualTimer::EXTRAPOL) << endl;
                cout << "Time spent in final expression of scheme : " << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
                cout << "Time spent in the time scheme : " << glob_chrono.GetSeconds(VirtualTimer::SCHEME) << endl;
                cout << "Time spent in tests : " << glob_chrono.GetSeconds(VirtualTimer::MASS) << endl;
                cout << "Time spent in the global simulation : " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
                cout << "Time spent in EvaluateFunction : " << glob_chrono.GetSeconds(VirtualTimer::JACOBIAN) << endl;
                //cout.unsetf(ios::fixed);
                //cout.setf(ios::scientific);
              }
	  }
      
      glob_chrono.Stop(VirtualTimer::ALL);
      //DISP(z);
      //var.time_scheme.GetIterate().Write("xref.dat");
      //int testinput; cout << "Waiting..." << endl; cin >> testinput;
    }

  return FinalizeMontjoie();
}
