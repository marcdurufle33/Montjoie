
#include "Elliptic/Aeroacoustic/MontjoieAeroacoustic.hxx"
#include "Elliptic/Aeroacoustic/HarmonicGoldsteinEquation.hxx"

#include "MontjoieFlag.hxx"
#include "MontjoieHeader.hxx"
#include "MontjoieInline.hxx"
//#include "Elliptic/Aeroacoustic/MontjoieAeroacousticHeader.hxx"
//#include "Elliptic/Aeroacoustic/MontjoieAeroacousticInline.hxx"

#include "Elliptic/Aeroacoustic/HarmonicGoldsteinEquation.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"


#include "Elliptic/Aeroacoustic/HarmonicGalbrun.cxx"

#include "Harmonic/TransparencyCondition.hxx"
#include "Harmonic/TransparencyCondition.cxx"
#include "Harmonic/TransparencyConditionInline.cxx"

using namespace Montjoie;

//SELDON_EXTERN template class VarHarmonic<GoldsteinEquation<Complex_wp, Dimension2> >;

int main(int argc, char **argv)
{
  InitMontjoie(argc,argv);
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  if(argc > 1)
    {
      string input_file(argv[1]), type_element, type_equation;
      getElement_Equation(input_file, type_element, type_equation);
      
      if (type_equation == "GOLDSTEIN_DG")
        {
          EllipticProblem<GoldsteinEquationDG<Complex_wp, Dimension2> > Vars;        
          Vars.RunAll(input_file, type_element, type_equation);
        }
      else
        {
          EllipticProblem<GoldsteinEquation<Complex_wp, Dimension2> > Vars;        
          Vars.RunAll(input_file, type_element, type_equation);
        }
    }
  else
    {
      std::cout << "Please provide a data file" << std::endl;
      abort();
    }

  return FinalizeMontjoie();
}
