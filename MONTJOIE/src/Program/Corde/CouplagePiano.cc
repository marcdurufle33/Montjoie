//pour la corde
#define MONTJOIE_WITH_ONE_DIM

// pour le soundboard
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
// #define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_TRANSMISSION

#include "Corde/MontjoiePiano.hxx"

using namespace Montjoie;



int main(int argc, char** argv)
{
  InitMontjoie(argc,argv);
  
  if (argc != 2)
    {
      cout << "file name missing"<<endl;
      abort();
    }

  typedef Montjoie::Real_wp Real_wp;  
  string data_file = argv[1];

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif  

  
  if (rank_proc == 0)
    {
      int size_tot; MPI_Comm_size(MPI_COMM_WORLD, &size_tot);
      cout << "Simulation launched on " << size_tot << " procs. " << endl;
    }
  string type_element, type_equation;
  getElement_Equation(data_file, type_element, type_equation);
  
  CouplagePiano Piano(argc,argv);
  
  Piano.ConstructAll(data_file);
  	
  Real_wp t0 = Piano.var_fluide.var_vibro.GetInitialTime();
  Real_wp dt = Piano.var_fluide.var_vibro.GetTimeStep();
  Real_wp tf = Piano.var_fluide.var_vibro.GetFinalTime();  

  // then main loop in time
  int nb_max_iter;
  if (dt != Real_wp(0))
    nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
  else
    {
      cout << "Enter a time step different from 0 " << endl;
      abort();
    }
  
  if (rank_proc == 0)
    {  
      
      cout << "--Celerity: " << Piano.var_fluide.var_vibro.GetCelerity() << "m.s-1" << endl << endl;
      cout << "--Frequency: " <<  Piano.var_fluide.var_vibro.GetFrequency() << "Hz" <<endl << endl;
      cout << "--Time Step dt: "<<  Piano.var_fluide.var_vibro.GetTimeStep() << endl << endl;
      cout << "Frequence supportee par le maillage :" << ( Piano.var_fluide.var_vibro.GetCelerity()*
                                                           Piano.var_fluide.var_vibro.GetOrder())
        / (10* Piano.var_fluide.var_vibro.GetMeshSize()) << endl << endl;
      cout << "---------------------------------------" << endl;
    }
  
  int nt0 = 0;
  bool save_snapshot = true;
  if (Piano.var_fluide.var_vibro.var_volume.load_reprise)
    {
      nt0 = Piano.var_fluide.var_vibro.var_volume.load_iter_reprise;
      save_snapshot = false;
    }
	
  for (int nt = nt0; nt < nb_max_iter; nt++)
    {
      Real_wp temps = t0+nt*dt;
      
      // on sauvegarde les instantanes si besoin est
      if (nt > nt0+10)
        Piano.SaveDatas(nt);
      
      Piano.Advance(nt, temps, save_snapshot);
      
      if (rank_proc == 0)
        {
          Piano.chrono.Start(VirtualTimer::ALL);
          Piano.WriteEnergy();          
          Piano.chrono.Stop(VirtualTimer::ALL);
        }      
      
      save_snapshot = true;
    }
  
  return FinalizeMontjoie();  
}
