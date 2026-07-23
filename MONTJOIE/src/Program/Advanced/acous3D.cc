#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define SELDON_WITH_PRECONDITIONING
#define MONTJOIE_WITH_TIME_REVERSAL

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);

#ifdef SELDON_WITH_MPI
      int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
      int rank_proc(0);
#endif
      
      if (rank_proc == 0)
	cout << " Data file " << file_name_data << endl;
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      if (rank_proc == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;

      if ((!type_equation.compare("ACOUSTIC_DG")) || (!type_equation.compare("ACOUSTIC_HDG"))
	  || (!type_equation.compare("HELMHOLTZ_HDG")) )
	{
	  HyperbolicProblem<AcousticEquationDG<Dimension3> > Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
	    cout<<" Acoustic Solver with Gauss-Lobatto hexahedra"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<rank_proc<<" we destroy the variables "<<endl;
	}
      else
	{
	  HyperbolicProblem<AcousticEquation<Dimension3> > Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
	    cout<<" Acoustic Solver with Discontinuous Galerkin and Gauss-Lobatto hexahedral"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"acous3D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
