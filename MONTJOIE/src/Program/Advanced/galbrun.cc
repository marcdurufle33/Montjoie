#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Hyperbolic/Aeroacoustic/MontjoieTimeGalbrun.hxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
        
  if (argc>1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      if (rank_proc == 0)
        cout << "TypeEquation = " << type_equation << endl;

      if (!type_equation.compare("HARMONIC_GALBRUN_MODEL_DIV"))
        {
          EllipticProblem<HarmonicGalbrunEquation<Dimension2> > Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else if (!type_equation.compare("HARMONIC_GALBRUN"))
        {
          EllipticProblem<HarmonicGalbrunEquationDG<Dimension2> > Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else if (!type_equation.compare("HARMONIC_GALBRUN_SIPG"))
        {
          EllipticProblem<HarmonicGalbrunEquationSipg<Dimension2> > Vars;
	  //EllipticProblem<HarmonicGalbrunEquationSipg<Dimension3> > Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else if (!type_equation.compare("HARMONIC_GALBRUN_H1"))
        {
          EllipticProblem<HarmonicGalbrunEquationH1> Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
      else
        {
          HyperbolicProblem<TimeGalbrunEquation<Dimension2> > Vars;
          Vars.RunAll(file_name_data, type_element, type_equation);
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  if (rank_proc == 0)
    cout<<"End of the program"<<endl; 

  return FinalizeMontjoie();
}

