#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_TIME_REVERSAL

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

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

      HyperbolicProblem<AcousticEquation<Dimension3> > Vars;
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" Acoustic Solver with Lobatto hexahedrals Qk"<<endl; 
      
      Vars.RunAll(file_name_data, type_element, type_equation);
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"time3D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
