#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_HP_H1
#define MONTJOIE_WITH_NODAL_DG

#define SELDON_WITH_PRECONDITIONING

#define MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
#define MONTJOIE_WITH_THIN_SLOT_MODEL
#define MONTJOIE_WITH_TRANSMISSION

//#define MONTJOIE_WITH_WIRES

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

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
      
      if ((type_equation == "HELMHOLTZ_DG") || (type_equation == "HELMHOLTZ_HDG") || (type_equation == "ACOUSTIC_HDG"))
        {
          EllipticProblem<HelmholtzEquationDG<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            {
	      if (type_equation == "HELMHOLTZ_DG")
		cout << "Helmholtz Solver with LDG formulation" << endl;
	      else if (type_equation == "HELMHOLTZ_HDG")
		cout << "Helmholtz Solver with HDG formulation" << endl;
	    }
          
          Vars.RunAll(file_name_data, type_element, type_equation);
          
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
      else
        {
          EllipticProblem<HelmholtzEquation<Dimension2> > Vars;	      
          if (Vars.GetRankProcMode() == 0)
            {
	      if (type_equation == "HELMHOLTZ_SIPG")
		cout<<"Helmholtz Solver with SIPG formulation "<<endl;
	      else
		cout<<"Helmholtz Solver with continuous formulation "<<endl;
	    }
          
          Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
