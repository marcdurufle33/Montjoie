#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Aeroacoustic/MontjoieAeroacoustic.hxx"
#include "Nathan/MontjoiePoisson.hxx"


using namespace Montjoie;

int main(int argc, char **argv)
{
  InitMontjoie(argc,argv);

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0)
#endif

  if(argc>1){
    string input_file(argv[1]), type_element, type_equation;
    getElement_Equation(input_file, type_element, type_equation);

    //    PoissonEquation<double, double>  eq;

    EllipticProblem<PoissonEquation<Real_wp,Dimension2> > Vars;
    //    EllipticProblem<PoissonEquation> Vars;
    Vars.RunAll(input_file, type_element, type_equation); 
  }   
  else{
    std::cout << "You need to provide a data file" << std::endl;
    abort();
  }
}
