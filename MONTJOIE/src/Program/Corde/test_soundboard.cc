//pour la corde
#define MONTJOIE_WITH_ONE_DIM
//#define SELDON_WITH_MUMPS
// #define SELDON_FILE_ITERATIVE_GMRES_CXX

// pour le soundboard
#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Corde/MontjoieSoundboard.hxx"

using namespace Montjoie;



int main(int argc, char** argv)
{
  //  ./Test_Soundboard.x fichier.ini Deltat NbPasDeTemps DOSSIER
  // ./Test_Soundboard.x example/reissner/plaquetest.ini 1e-4 100  TEST
  //80  0.3 0.3 0.1 0.1
  
  InitMontjoie(argc,argv);
  typedef Montjoie::Real_wp Real_wp;

  string name_data_file(argv[1]);
  Real_wp Deltat = to_num<Real_wp>(argv[2]);
  int Nt = to_num<int>(argv[3]);
  string DOSSIER(argv[4]);

  string type_element, type_equation;
  getElement_Equation(name_data_file, type_element, type_equation);
   
  DISP(Deltat);
  DISP(Nt);
  DISP(DOSSIER);
  
  // on cree l'objet
  Soundboard soundboard;
  
  Vector<string> lines_data_file;
#ifdef SELDON_WITH_MPI
  ReadLinesFile(name_data_file, lines_data_file, MPI_COMM_WORLD);
#else
  ReadLinesFile(name_data_file, lines_data_file);
#endif
  
  SourceChi fsrc;
  soundboard.init_file=name_data_file;
  soundboard.Deltat = Deltat;
  soundboard.RunSoundboardDiagonalization(string("./"), fsrc, lines_data_file, true);
  
  soundboard.ConstruitVectSecondMembreModal();
  
  // a voir si le second membre est bien ce qu'on croit (MF,U) doit donner F_j cad la composante de F dans la base des U_j...;
  soundboard.ConstruitMatricesResolExacte();
  soundboard.ConstruitLambda0etPoint();
  
  // prise en compte de la donnee ini	
  //Real_wp rhodelta = rho*soundboard.VarSoundboard.ref_delta(1).GetConstant();
  Real_wp rhodelta = 0.0;
  soundboard.Initialise(rhodelta);
  
  // pour garder en memoire la force (t)
  VectReal_wp VecteurTemps(Nt);
  VecteurTemps.Fill(0);
  
  VectReal_wp ener(Nt+1);
  ener.Fill(0);
  
  // pour la donnee ini
  ener(0) = soundboard.GetKineticEnergy();
  ener(0)+= soundboard.GetPotentialEnergy();
  DISP(ener(0));
  
  soundboard.VarSoundboard.InitOutput(Real_wp(0));
  
  for (int nt=0; nt<Nt; nt++)
    {
      Real_wp fonction_temps = 0.0;
      //-(2.0*square(pi_wp)*square(pulsation*nt*Deltat-1.0)-1.0)*exp(-square(pi_wp)*square(pulsation*nt*Deltat-1.0));
      VecteurTemps(nt) = fonction_temps;
      // sin(nt*Deltat*pulsation*M_PI/180);
      
      Real_wp frequency = soundboard.VarSoundboard.GetFrequency();
      //Real_wp pulse = sin(nt*Deltat*2.0*pi_wp*frequency);
      //Real_wp pulse = RickerTime(nt*Deltat, frequency);
      Real_wp t = nt*Deltat;
      Real_wp gaussian_term = exp(-square(pi_wp*(frequency*t - 1.0)));
      Real_wp pulse = gaussian_term*(-1.0 + 2.0*square(pi_wp*(frequency*t - 1.0)));
      //DISP(t); DISP(pulse);
      
      soundboard.ActualiseLambda(fonction_temps, pulse);
      soundboard.WriteOutputFEM(nt, nt*Deltat, soundboard.VarSoundboard);
      
      ener(nt+1) = soundboard.GetKineticEnergy();
      ener(nt+1)+= soundboard.GetPotentialEnergy();
      //DISP(ener(nt+1));
    }
  
  VecteurTemps.WriteText(DOSSIER + "/SourceTemps.dat");
  ener.WriteText(DOSSIER + "/Energie.dat");
  
  return 0;
}
