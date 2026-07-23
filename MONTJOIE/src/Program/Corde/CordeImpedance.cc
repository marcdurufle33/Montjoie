//pour la corde
#define MONTJOIE_WITH_ONE_DIM
//#define SELDON_WITH_MUMPS

// pour le soundboard
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Corde/MontjoieSoundboard.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  if (argc==1)
    {
      cout << "file name missing"<<endl;
      abort();
    }
  
  string data_file = argv[1];
  typedef Montjoie::Real_wp Real_wp;
  
  CouplagePiano Piano(argc,argv);
  
  Piano.ConstructAllImpedance(data_file);
  Piano.Display();
  
  Real_wp t0 = Piano.note.t_begin;
  Real_wp dt = Piano.note.Deltat;
  Real_wp tf = Piano.note.t_end;  
  
  // then main loop in time
  int nb_max_iter;
  if (dt != Real_wp(0))
    nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
  else
    {
      cout << "Enter a time step different from 0 " << endl;
      abort();
    }
  
  
  for (int nt = 0; nt < nb_max_iter; nt++)
    {
      Real_wp temps = t0+nt*dt;
      Piano.AdvanceStringPlate(nt,temps);
      
      Piano.WriteEnergy();
    }
        
  Piano.Display();
  
  FinalizeMontjoie();
  return(0);  
}   
