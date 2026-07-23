// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#define SELDON_WITH_PRECONDITIONING

#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"
// #include "FiniteElement/Pyramid/PyramidHcurlOther.hxx"
// #include "FiniteElement/Pyramid/PyramidHcurlOther.cxx"


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
      
      EllipticProblem<StaticMaxwellEquation_3D> Vars;
      if (Vars.GetRankProcMode() == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (Vars.GetRankProcMode() == 0)
	cout<<"Maxwell Solver with tetrahedra of Nedelec's first family "<<endl; 
      
      Vars.RunAll(file_name_data, type_element, type_equation);

      if (Vars.GetRankProcMode() == 0)
	cout<<" we destroy the variables "<<endl;
      
    }
  
  return FinalizeMontjoie();
}
