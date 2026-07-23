#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_TIME_REVERSAL
#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc>1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

#ifdef MONTJOIE_WITH_REAL_TIMING      
      RealTimer chrono;
      chrono.Reset(0);
      chrono.Start(0);
#endif
      
      if ((!type_equation.compare("ACOUSTIC_DG")) || (!type_equation.compare("ACOUSTIC_HDG"))
	   || (!type_equation.compare("HELMHOLTZ_HDG")) )
	{
	  HyperbolicProblem<AcousticEquationDG<Dimension2> > Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" Acoustic Solver with LDG formulation"<<endl; 
	  
          Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}
      else
	{
          HyperbolicProblem<AcousticEquation<Dimension2> > Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" Acoustic Solver with continuous elements"<<endl; 
	  
          Vars.RunAll(file_name_data, type_element, type_equation);
	  
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}

#ifdef MONTJOIE_WITH_REAL_TIMING      
      chrono.Stop(0);
      int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);

      if (rank_proc == 0)
        cout << "Simulation lasted " << chrono.GetSeconds(0) << " s" << endl;
#endif
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"acous2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
