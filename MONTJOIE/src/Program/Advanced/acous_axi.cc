#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#include "Instationary/MontjoieUnsteady.hxx"
#include "Hyperbolic/Acoustic/AxiSymAcoustic.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
#endif

#include "Hyperbolic/Acoustic/AxiSymAcoustic.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      HyperbolicProblem<AcousticEquationAxi> Vars;
      if (Vars.var_harmonic.GetRankProcMode() == 0)
	cout<<" Acoustic Solver with continuous formulation"<<endl; 
	  
      Vars.RunAll(file_name_data, type_element, type_equation);
      if (Vars.var_harmonic.GetRankProcMode() == 0)
        cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"acous_axi.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
