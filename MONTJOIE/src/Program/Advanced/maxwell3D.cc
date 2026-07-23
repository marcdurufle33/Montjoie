// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
//#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#define SELDON_WITH_PRECONDITIONING
#define MONTJOIE_WITH_TRANSMISSION

//#define MONTJOIE_WITH_WIRES

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

      if (type_equation == "HARMONIC_MAXWELL_HDG")
        {
          EllipticProblem<HarmonicMaxwellEquationHdg_3D> Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else if (type_equation == "HARMONIC_MAXWELL_DG")
        {
          EllipticProblem<HarmonicMaxwellEquation_3D_DG> Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else
        {
          EllipticProblem<HarmonicMaxwellEquation_3D> Vars;     
          Vars.RunAll(file_name_data, type_element, type_equation);

          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"maxwell3D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
