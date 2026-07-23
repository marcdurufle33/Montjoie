#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define SELDON_WITH_PRECONDITIONING

#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

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
   
      if ((type_equation == "LAPLACE_DG") || (type_equation == "LAPLACE_HDG"))
	{
	  if (dim_N == 2)
	    {
              EllipticProblem<LaplaceEquationDG<Dimension2> > Vars;	  
              if (Vars.GetRankProcMode() == 0)
                cout<<"Laplace solver with LDG formulation "<<endl; 
                	    
	      Vars.RunAll(file_name_data, type_element, type_equation);

              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
	    }
	  else
	    {
	      EllipticProblem<LaplaceEquationDG<Dimension3> > Vars;	  
              if (Vars.GetRankProcMode() == 0)
                cout<<"Laplace solver with LDG formulation "<<endl; 
              
	      Vars.RunAll(file_name_data, type_element, type_equation);

              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
	    }          
	}
      else
	{	  
	  if (dim_N == 2)
	    {
	      EllipticProblem<LaplaceEquation<Dimension2> > Vars;
              if (Vars.GetRankProcMode() == 0)
                cout<<"Poisson solver with continuous elements"<<endl; 
	      
              Vars.RunAll(file_name_data, type_element, type_equation);

              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
	    }
	  else
	    {
	      EllipticProblem<LaplaceEquation<Dimension3> > Vars;
              if (Vars.GetRankProcMode() == 0)
                cout<<"Poisson solver with continuous elements"<<endl; 
              
              Vars.RunAll(file_name_data, type_element, type_equation);
              
              if (Vars.GetRankProcMode() == 0)
                cout<<" we destroy the variables "<<endl;
	    }
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"laplace.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
