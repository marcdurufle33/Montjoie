// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define SELDON_WITH_PRECONDITIONING

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_HP_HCURL
#define MONTJOIE_WITH_NODAL_HCURL

#include "Elliptic/Maxwell/MontjoieMaxwell2D.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"

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

      int dim_N = 3;
      if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
	dim_N = 2;
      
      if (dim_N == 2)
	{
	  EllipticProblem<StaticMaxwellEquation_2D> Vars;
          if (Vars.GetRankProcMode() == 0)
	    cout<<"Maxwell Solver with triangles of Nedelec's first family "<<endl;
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
	    cout<<" we destroy the variables "<<endl;
	}
      else
	{	  
	  if ((type_equation == "STATIC_MAXWELL3D_DG")
              || (type_equation == "TIME_MAXWELL3D_DG"))
            {
              EllipticProblem<StaticMaxwellEquation_3D_DG> Vars;
              if (Vars.GetRankProcMode() == 0)
                cout<<"Maxwell Solver with tetrahedra of Nedelec's first family "<<endl; 

              Vars.RunAll(file_name_data, type_element, type_equation);

              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
            }
          else
            {
              EllipticProblem<StaticMaxwellEquation_3D> Vars;
              if (Vars.GetRankProcMode() == 0)
                cout<<"Maxwell Solver with tetrahedra of Nedelec's first family "<<endl;  
	      
	      Vars.RunAll(file_name_data, type_element, type_equation);

              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
            }
	}
    }
  
  return FinalizeMontjoie();
}
