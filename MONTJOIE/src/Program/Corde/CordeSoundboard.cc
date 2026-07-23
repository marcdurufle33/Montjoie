//pour la corde
#define MONTJOIE_WITH_ONE_DIM

// pour le soundboard
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Corde/MontjoieSoundboard.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc,argv);
  if (argc==1) {
    cout << "file name missing"<<endl;
    abort();
  }

  string data_file = argv[1];
  typedef Montjoie::Real_wp Real_wp;  

  string type_element, type_equation;
  getElement_Equation(data_file, type_element, type_equation);
  

  CouplagePiano Piano(argc,argv);
  glob_chrono.Start(VirtualTimer::EXTRAPOL);
  Piano.ConstructAllStringSoundboard(data_file);
  

  Piano.Display();
  glob_chrono.Stop(VirtualTimer::EXTRAPOL);

  Real_wp t0 = Piano.note.t_begin;
  Real_wp dt = Piano.note.Deltat;
  Real_wp tf = Piano.note.t_end;  
  Real_wp t_display = Piano.note.t_display;
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  Vector<int64_t> t_tmp;
  MpiBcast(MPI_COMM_WORLD, &t0, t_tmp, 1, 0);
  MpiBcast(MPI_COMM_WORLD, &dt, t_tmp, 1, 0);
  MpiBcast(MPI_COMM_WORLD, &tf, t_tmp, 1, 0);
  MpiBcast(MPI_COMM_WORLD, &t_display, t_tmp, 1, 0);
#else
  int rank_proc(0);
#endif
  

  // then main loop in time
  int nb_max_iter;
  if (dt != Real_wp(0))
    nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
  else
    {
      cout << "Enter a time step different from 0 " << endl;
      abort();
    }
  
  for (int nt=0; nt<nb_max_iter; nt++)
    {
      glob_chrono.Start(VirtualTimer::PML);
      Real_wp ener_tot;

      glob_chrono.Start(VirtualTimer::PROD);
      Real_wp temps = t0+nt*dt;

      if(rank_proc == 0)
        Piano.AdvanceStringPlate(nt,temps);
      glob_chrono.Stop(VirtualTimer::PROD);
      
      glob_chrono.Start(VirtualTimer::STIFFNESS);
      if(rank_proc == 0)
        ener_tot = Piano.WriteEnergy();
      
      int n = toInteger(round(temps/(t_display)));
      if(rank_proc == 0)
        {
          if (abs(temps/t_display - n) < dt/2)		
            {
              DISP(glob_chrono.GetSeconds(VirtualTimer::EXTRAPOL));
              DISP(glob_chrono.GetSeconds(VirtualTimer::PML));
              DISP(glob_chrono.GetSeconds(VirtualTimer::OUTPUT));
              DISP(glob_chrono.GetSeconds(VirtualTimer::STIFFNESS));
              DISP(glob_chrono.GetSeconds(VirtualTimer::PROD));
              DISP(chrono_string.GetSeconds(Piano.CHRONO_NOTE));
              DISP(chrono_string.GetSeconds(Piano.CHRONO_COUPLAGE));
              DISP(chrono_string.GetSeconds(Piano.CHRONO_PLATE));
              DISP(chrono_string.GetSeconds(Piano.CHRONO_ENERGY));
            }
        }
      glob_chrono.Stop(VirtualTimer::STIFFNESS);
      
      glob_chrono.Stop(VirtualTimer::PML);
    }
  
  Piano.Display();
  
  FinalizeMontjoie();
  return(0);
  
}
