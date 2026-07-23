
#include "Montjoie.hxx"

#include "Elliptic/Thermo/ViscousCellEquation.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"

using namespace Montjoie;


// Solve viscous cell problem (from homogenized vicsothermal acoustics)
// author : Alexis THIBAULT

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

      EllipticProblem<ViscousCellEquation<Complex_wp> > Vars;
      Vars.RunAll(input_file, type_element, type_equation);
    }
  else
    {
      std::cout << "Please provide a data file" << std::endl;
      abort();
    }

  return FinalizeMontjoie();
}
