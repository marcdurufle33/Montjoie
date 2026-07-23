//pour la corde
#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Corde/MontjoieString.hxx"

using namespace Montjoie;


int main(int argc, char** argv)
{
  
  InitMontjoie(argc,argv);
  
  if (argc == 1)
  {
	cout << "Cette commande demande au moins un argument " << endl;
	cout << "Par exemple, ./multistring.x example/data_file.ini " << endl;
	abort();
  }

  {  
  string name_data_file(argv[1]);

  // 5 fonctions de base par element => ordre 4 en espace
  enum{n=5};
  
  // (u,v,phi)
#ifdef PIANO_NONLINEAR_STIFF_STRING
   MultiString<n, n, WaveEquationStiffNLString> multi;
#endif

  //(u) --> parfait
#ifdef PIANO_LINEAR_STRING
   MultiString<n, n, WaveEquationString> multi;
#endif
    
  // (u,phi) --> parfait
#ifdef PIANO_STIFF_STRING
   MultiString<n, n, WaveEquationStiffString> multi;
#endif

  // (u,v) --> pas parfait mais en forcant le calcul de la jacob on a 1e-10
#ifdef PIANO_NONLINEAR_STRING
  MultiString<n, n, WaveEquationNLString> multi;
#endif

#ifdef PIANO_LINEAR_STIFF_STRING
  MultiString<n, n, WaveEquationStiffLINString> multi;
#endif

#ifdef PIANO_LINEAR_STRING_TL
  MultiString<n,n,WaveEquationStringLINTL> multi;
#endif

#ifdef PIANO_LINEAR_STRING_2T
  MultiString<n,n,WaveEquationStringLIN2T> multi;
#endif

#ifdef PIANO_NONLINEAR_STRING_2T
  MultiString<n, n, WaveEquationStringNL2T> multi;
#endif

#ifdef PIANO_NONLINEAR_STIFF_STRING_2T
  MultiString<n,n,WaveEquationStringStiffNL2T> multi;
#endif

#ifdef PIANO_LINEAR_STIFF_STRING_2T
  MultiString<n,n,WaveEquationStringStiffLIN2T> multi;
#endif
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
  int rank_proc(0), nb_proc(1);
#endif

  // (un truc lineaire dans U)
  //  MultiString<5,5,WaveEquationNLPipo> multi;

  ReadInputFile(name_data_file, multi);
  multi.name_data_file = name_data_file;
  
  // on essaie de lire le numero de l'iteration
  IVect num_iter;
  string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
  string prefix = multi.path_reprise;

  if (multi.load_reprise)
    {
      num_iter.Read(prefix+"IterationPianoSave"+suffix);
      if (num_iter.GetM() != 3)
        {
          // on pense alors que l'ecriture critique a ete interrompu,
          // on regarde l'autre fichier dans ce cas
          num_iter.Read(prefix+"IterationPianoSaveNew"+suffix);
                
          if (num_iter.GetM() != 3)
            {
              cout << "Echec de la reprise" << endl;
              cout << "Veuillez corriger le probleme" << endl;
              abort();
            }
        }
      
      multi.load_iter_reprise = num_iter(0);
    }

  int num_dossier = -1;
  //bool remove_file = true;
  Matrix<Real_wp> val;
  if (multi.load_reprise)
    {
      num_dossier = num_iter(1);
      //remove_file = false;

      suffix = to_str(num_iter(2)) + "_P" + suffix;
      val.Read(prefix+"ValPianoSave"+suffix);
      VectReal_wp tmp;
      tmp.Read(prefix+"ForceMarteauPianoSave"+suffix);
      multi.beginning_time_force = tmp(0);
      multi.last_time_force = tmp(1);
      multi.final_time_shank = tmp(2);
      multi.hammer.energy_free_hammer = tmp(3);
      multi.shank.free_kinetic_energy = tmp(4);
      multi.shank.free_potential_energy = tmp(5);
      multi.shank.free_internal_energy = tmp(6);
      multi.shank.free_shank_energy = tmp(7);
    }
  
  multi.ConstructAll(num_dossier);
  string file_name = multi.DOSSIER + multi.file_storage_param;
  ofstream file_out(file_name.data());
  multi.Display(file_out, true);
  file_out.close();
  
  if (multi.load_reprise)
    {
      int nb_dof_note = val.GetM();
      
      bool reprise_failed = false;
      if (multi.val.GetM() != nb_dof_note)
        reprise_failed = true;
      
      if (reprise_failed)
        {
          cout << "Echec de la reprise" << endl;
          abort();
        }

      multi.val = val;
    }
  
  multi.RunAll();
  } // multi is destructed
  
  return FinalizeMontjoie();
}
