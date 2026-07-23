#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

#include "Hyperbolic/Acoustic/VarAdvection.hxx"
#include "Hyperbolic/Acoustic/VarAdvectionInline.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc>1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      // DISP(file_name_data);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      if (type_equation == "HARMONIC_ADVECTION_DG")
        {
          EllipticProblem<HarmonicAdvectionEquationDG<Dimension2> > Vars;

          Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;          
        }
      else if (type_equation == "HARMONIC_ADVECTION")
        {
          EllipticProblem<HarmonicAdvectionEquation<Dimension2> > Vars;

          Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;         
        }
      /*      else if (type_equation == "HARMONIC_DOUBLE_ADVECTION")
        {
          EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension2> > Vars;

          Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;         
        }
      else if (type_equation == "HARMONIC_DOUBLE_ADVECTION_DG")
        {
          EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension2> > Vars;

          Vars.RunAll(file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;         
            }*/
      else
        {
          cout << "Unknown equation"  << endl;          
        }
      
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"advec.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
}
