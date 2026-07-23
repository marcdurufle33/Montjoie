#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL
#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG

// #define MONTJOIE_WITH_WIRES
#define MONTJOIE_WITH_TWO_DIM

#include "Hyperbolic/Maxwell/MontjoieTimeMaxwell3D.hxx"

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
      
      //HyperbolicProblem<TimeMaxwellEquation_3D_DG> Vars;
      HyperbolicProblem<TimeMaxwellEquation_3D> Vars;

      if (Vars.var_harmonic.GetRankProcMode() == 0)
        cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;

      if (Vars.var_harmonic.GetRankProcMode() == 0)
        cout<<" Maxwell Solver with DG ortho"<<endl; 
      
      Vars.RunAll(file_name_data, type_element, type_equation);
      if (Vars.var_harmonic.GetRankProcMode() == 0)
        cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file as an argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
