#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_NODAL_HDIV
#define MONTJOIE_WITH_HP_HDIV

#define SELDON_WITH_PRECONDITIONING

#define MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
#define MONTJOIE_WITH_THIN_SLOT_MODEL

// #define MONTJOIE_WITH_WIRES

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Elliptic/Helmholtz/HelmholtzHdiv.hxx"
#include "Elliptic/Helmholtz/HelmholtzHdivInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/HelmholtzHdiv.cxx"
#endif

//#include "FiniteElement/Pyramid/PyramidHdivOther.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{

  InitMontjoie(argc, argv);
  
  if (argc>1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      //DISP(file_name_data);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      //cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      if (type_equation == "LAPLACE_DIV")
	{
	  EllipticProblem<LaplaceEquationHdiv<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Helmholtz Solver with Lobatto quadrilaterals "<<endl; 
      
	  Vars.RunAll(file_name_data, type_element, type_equation);

          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
	}
      else
	{
	  EllipticProblem<HelmholtzEquationHdiv<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Helmholtz Solver with Lobatto quadrilaterals "<<endl; 
      
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
