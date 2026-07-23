
#include "Montjoie.hxx"

#include "Elliptic/Thermo/ViscoThermalEquation.cxx"
#include "Elliptic/Thermo/AxiViscoThermal.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"

using namespace Montjoie;

int main(int argc, char **argv)
{
  InitMontjoie(argc,argv);
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  if (argc > 1)
    {
      string input_file(argv[1]), type_element, type_equation;
      getElement_Equation(input_file, type_element, type_equation);

      int dim_N = 3;
      if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
	dim_N = 2;

      if (type_equation == "VISCO_THERMAL_AXI")
        {
          EllipticProblem<ViscoThermalEquationAxi<Complex_wp> > Vars;
          Vars.RunAll(input_file, type_element, type_equation);
        }
      else
        {
          if (dim_N == 2)
            {
              EllipticProblem<ViscoThermalEquation<Dimension2, Complex_wp> > Vars;        
              Vars.RunAll(input_file, type_element, type_equation);
            }
          else
            {
              EllipticProblem<ViscoThermalEquation<Dimension3, Complex_wp> > Vars;        
              Vars.RunAll(input_file, type_element, type_equation);
            }
        }
    }
  else
    {
      std::cout << "Please provide a data file" << std::endl;
      abort();
    }

  return FinalizeMontjoie();
}
