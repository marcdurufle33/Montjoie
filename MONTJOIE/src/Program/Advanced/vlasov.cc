// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_NODAL_DG

#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Maxwell/MontjoieVlasovMaxwell.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

#ifdef SELDON_WITH_MPI
      int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
      int rank_proc(0);
#endif
      
      if (rank_proc == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;

      int dim_N = 3;
      if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
	dim_N = 2;
      
      if (dim_N == 2)
	{
	  if (type_equation == "TIME_MAXWELL_2D")
	    {
	      HyperbolicProblem<VlasovMaxwellEquation2D> Vars;
	      Vars.RunAll(file_name_data, type_element, type_equation);
	      if (rank_proc == 0)
		cout<<" we destroy the variables "<<endl;
	    }
	  else
	    {
	      if (rank_proc == 0)
		cout<<"Maxwell Solver with DG Lobatto "<<endl; 
	      
	      HyperbolicProblem<VlasovMaxwellEquation_2D_DG> Vars;
	      Vars.RunAll(file_name_data, type_element, type_equation);
	      if (rank_proc == 0)
		cout<<" we destroy the variables "<<endl;
	    }
	}
      else
	{
	  if (rank_proc == 0)
	    cout<<"Maxwell Solver with hexahedras of Nedelec's second family "<<endl; 
	  
	  HyperbolicProblem<VlasovMaxwellEquation3D> Vars;
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (rank_proc == 0)
            cout<<" we destroy the variables "<<endl;
	}
    }
  
  return FinalizeMontjoie();
}
