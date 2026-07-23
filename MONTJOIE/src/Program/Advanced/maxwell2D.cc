#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#include "Elliptic/Maxwell/MontjoieMaxwell2D.hxx"

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
      
      EllipticProblem<HarmonicMaxwellEquation_2D> Vars;      
      if (Vars.GetRankProcMode() == 0)
	cout<<"Maxwell Solver with edge elements "<<endl; 
      
      Vars.RunAll(file_name_data, type_element, type_equation);
      if (Vars.GetRankProcMode() == 0)
	cout<<" we destroy the variables "<<endl;
    }

  return FinalizeMontjoie();
}
