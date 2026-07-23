#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"
#include "Hyperbolic/WaterWaves/CauchyPoisson.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc>1)
    {
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
          HyperbolicProblem<CauchyPoissonEquation<Dimension2> > Vars;
          if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<"Solver with Gauss-Lobatto quadrangles Qk"<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
	  
          if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}
      else
	{
          HyperbolicProblem<CauchyPoissonEquation<Dimension3> > Vars;
          if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<"Solver with classical triangles Pk "<<endl; 
	  
	  Vars.RunAll(file_name_data, type_element, type_equation);
          
          if (Vars.var_harmonic.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"cauchy_poisson.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
