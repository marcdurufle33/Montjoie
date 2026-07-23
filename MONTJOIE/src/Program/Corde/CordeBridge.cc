//pour la corde
#define MONTJOIE_WITH_ONE_DIM

// pour le soundboard
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Corde/MontjoieBridge.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc==1)
    {
      cout << "file name missing"<<endl;
      abort();
    }
 
  typedef Montjoie::Real_wp Real_wp; 
  string data_file = argv[1];
  
  CouplageBridge string_bridge(data_file);
  
  string_bridge.ConstructAll(data_file);
  string_bridge.Display();
  //~ cout << "Verify all parameters. Press any keys to start calculation." << endl;
  //~ getchar();
  
  Real_wp t0 = string_bridge.note.t_begin;
  Real_wp dt = string_bridge.note.Deltat;
  Real_wp tf = string_bridge.note.t_end;  
  
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
      string_bridge.AdvanceStringBridge(nt, temps);
      
      string_bridge.WriteEnergy();
    }
        
  string_bridge.Display();
  
  return FinalizeMontjoie();
}
