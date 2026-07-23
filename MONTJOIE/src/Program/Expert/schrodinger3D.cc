#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Elliptic/Maxwell/MontjoieSchrodinger3D.hxx"

using namespace Montjoie;

//! main function
void RunAll(SchrodingerProblem3D& var, const string& input_file)
{
  cout.setf(ios::scientific);
  
  All_TimeScheme<Complex_wp> scheme;
  
  // on lit le fichier de donnees
  ReadInputFile(input_file, var);

  // construction du numero de dossier
  // si le dossier n'est pas connu
  int num = -1;
  if (var.DOSSIER.size() == 0)
    {
      var.DOSSIER = "[STIFFOUT]/"; 
      EcritDossier(var.DOSSIER, input_file, num);
    }
  
  // initialisation du schema en espace
  var.ConstructAll();
  
  // initialisation du schema en temps (z)  
  int Nz = toInteger(ceil(var.GetZmax() / var.GetDz()));
  string name_file = GetBaseString(var.output_grid_z.GetTotalFieldFile());
  if (Nz < 1)
    return;
  
  Vector<Complex_wp> E0;
  var.GetInitialCondition(E0);
  var.time_scheme.SetInitialCondition(0.0, var.GetDz(), E0, var);
  
  // boucle principale
  for (int n = 0; n <= Nz; n++)
    {
      glob_chrono.Start(VirtualTimer::ALL);
      
      Real_wp z = n*var.GetDz();
      
      // snapshot is written on the file
      if (var.output_grid_z.SnapshotToStore(z))
        {
          int num = var.output_grid_z.GetSnapshotNumber();
	  string name = name_file + NumberToString(num) + ".dat";
	  string name_time = name_file + NumberToString(num) + "_time.dat";
	  var.WriteOutputFile(name, name_time, z, var.time_scheme.GetIterate());
	  var.output_grid_z.IncrementSnapshot();
	}
      
      // the solution at z = z^n is given to the object var
      var.GiveIterate(n, z, var.time_scheme.GetIterate());
      
      // scheme is advanced to replace u^n by u^{n+1}
      var.time_scheme.Advance(z, n, var);
      
      // on modifie l'itere si necessaire
      var.ModifyIterate(n, z+var.dz, var.time_scheme.GetIterate());

      if ( ((var.print_level >= 2) && (n% 100 == 0))
           || (var.print_level >= 6))
	{
          cout << "Norme pour z = " << z*var.z0_adim << " : "
               << Norm2(var.time_scheme.GetIterate()) << endl;
          
          if (var.print_level >= 4)
	    {
	      cout.setf(ios::fixed);
	      cout << "Time spent for outputs : " << glob_chrono.GetSeconds(VirtualTimer::OUTPUT) << endl;
	      cout << "Time spent to compute the phase : " << glob_chrono.GetSeconds(VirtualTimer::PML) << endl;
	      cout << "Time spent in the FFTs : " << glob_chrono.GetSeconds(VirtualTimer::FLUX) << endl;
	      cout << "Time spent in non-linear part : " << glob_chrono.GetSeconds(VirtualTimer::PROD) << endl;
	      cout << "Time spent in final expression of scheme : " << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
	      cout << "Time spent in the time scheme : " << glob_chrono.GetSeconds(VirtualTimer::SCHEME) << endl;
	      cout << "Time spent in the global simulation : " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
	      cout.unsetf(ios::fixed);
              cout.setf(ios::scientific);
	    }
	}
      
      glob_chrono.Stop(VirtualTimer::ALL);
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      SchrodingerProblem3D var;
      string input_file(argv[1]);
      
      string type_elt, type_equation;
      getElement_Equation(input_file, type_elt, type_equation);
      
      RunAll(var, input_file);
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"schrodinger3D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
