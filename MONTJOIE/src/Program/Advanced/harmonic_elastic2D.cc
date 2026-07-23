#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_TRANSMISSION
#define SELDON_WITH_PRECONDITIONING

#include "Elliptic/Elastic/MontjoieElastic.hxx"

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
      
      if (type_equation == "HARMONIC_ELASTIC_DG")
        {
          EllipticProblem<HarmonicElasticEquationDG<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Harmonic Elastic Solver with discontinuous finite elements"<<endl; 
          
          Vars.RunAll(file_name_data, type_element, type_equation);
          
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
      else
        {
          EllipticProblem<HarmonicElasticEquation<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Harmonic Elastic Solver with continuous finite elements"<<endl; 
      
          Vars.RunAll(file_name_data, type_element, type_equation);

          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"harmonic_elastic2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
