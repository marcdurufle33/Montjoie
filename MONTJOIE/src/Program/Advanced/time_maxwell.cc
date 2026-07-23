// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_HCURL_AXI

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_NODAL_HCURL

#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Maxwell/MontjoieTimeMaxwellAxi.hxx"

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
      
      if (type_equation == "TIME_MAXWELL_AXI_DG")
	{
          HyperbolicProblem<TimeMaxwellEquation_Axi_DG> Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<"Time Maxwell Solver with DG Gauss"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;	      
	}
      else
	{
          HyperbolicProblem<TimeMaxwellEquation_HcurlAxi> Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<"Time Maxwell Solver with H-curl axi"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;	      
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"time.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
