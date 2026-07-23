#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Hyperbolic/Elastic/MontjoieTimeElastic.hxx"

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
      
      if (!type_equation.compare("HARMONIC_REISSNER_MINDLIN"))
	{
	  EllipticProblem<HarmonicReissnerMindlinEquation> Vars;
	  if (Vars.GetRankProcMode() == 0)
	    cout<<"Reissner Mindlin Solver with Gauss-Lobatto quadrilaterals in time-harmonic domain"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.GetRankProcMode() == 0)
	    cout<<" we destroy the variables "<<endl;
	}
      else if (!type_equation.compare("STATIC_REISSNER_MINDLIN"))
	{
	  EllipticProblem<ReissnerMindlinEquation> Vars;
	  if (Vars.GetRankProcMode() == 0)
	    cout<<"Reissner Mindlin Solver with Gauss-Lobatto quadrilaterals in static domain"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.GetRankProcMode() == 0)
	    cout<<" we destroy the variables "<<endl;
	}
      else
	{
	  HyperbolicProblem<TimeReissnerMindlinEquation> Vars;
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
	    cout<<"Elastic Solver with Gauss-Lobatto quadrilaterals"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  if (Vars.var_harmonic.GetRankProcMode() == 0)
	    cout<<" we destroy the variables "<<endl;
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
