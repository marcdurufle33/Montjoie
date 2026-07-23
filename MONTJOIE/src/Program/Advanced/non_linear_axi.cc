#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Instationary/MontjoieTime.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"
#include "Elliptic/Maxwell/3D/NonLinearMaxwellProblem.hxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.cxx"
#include "Elliptic/Maxwell/3D/NonLinearMaxwellProblem.cxx"

using namespace Montjoie;

#include <xmmintrin.h>

// #define CSR_FLUSH_TO_ZERO         (1 << 15)

// probleme temporel
template<class TypeEquation>
void RunTimeProblem(NonLinearMaxwellProblem<TypeEquation>& var,
                    const string& input_file, const string& name_element,
		    const string& name_equation)
{
  typedef Complex_wp T;
  
  // reading parameters of the time scheme
  All_TimeScheme<T> scheme;
  ReadInputFile(input_file, scheme);
  
  // all what is needed to complete time iterations is constructed
  var.var_laplace.SetThresholdMatrix(1e-300);
  var.ConstructAll(input_file, name_element);
  
  // initialisation of time scheme
  var.InitTimeIterations();
  int N = var.GetNbDof();
  T zero; SetComplexZero(zero);
  Vector<T> E0(N), E0_prime(N);
  E0.Fill(zero);
  E0_prime.Reallocate(N);
  E0_prime.Fill(zero);
  
  Vector<T> E_current(var.var_laplace.GetNbDof());
  E_current.Fill(zero);
  
  scheme.SetInitialCondition(var.GetInitialTime(), var.GetTimeStep(), E0, E0_prime, var);
  
  // time iterations
  glob_chrono.Reset(VirtualTimer::ALL);
  int nb_iterations_time = toInteger(ceil((var.GetFinalTime()-var.GetInitialTime())/var.GetTimeStep()));
  Real_wp t;
  for (int nt = 0; nt <= nb_iterations_time; nt++)
    {
      glob_chrono.Start(VirtualTimer::ALL);
      
      // current time t^n
      t = var.GetInitialTime() + nt*var.GetTimeStep();
      
      // outputs are done if needed
      var.WriteSnapshot(nt, t, E_current);
      
      // time scheme is advanced
      var.AdvanceScheme(nt, t, scheme, E_current);
      
      if (var.var_laplace.GetRankProcMode() == 0)
        if (nt%100 == 0)
          cout << "t = " << t*var.t0_adim << endl;
      
      glob_chrono.Stop(VirtualTimer::ALL);
      
      if (var.var_laplace.GetRankProcMode() == 0)
        if (nt%100 == 0)
          cout << "real time elapsed " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
      
      //int test_input; cout << "we wait" << endl; cin >> test_input;
    }
  
  var.CloseBuffers();
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Donnez un fichier de donnees" << endl;
      abort();
    }
  
  InitMontjoie(argc, argv);
  cout.setf(ios::scientific);
  
  // in order to avoid denormal values
  //unsigned csr = __builtin_ia32_stmxcsr();
  //csr |= CSR_FLUSH_TO_ZERO;
  //__builtin_ia32_ldmxcsr(csr);
  
  _mm_setcsr( _mm_getcsr() | 0x8040 );
  
  string input_file(argv[1]);

  // on recupere le parametre d'adimensionalisation
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  // on lance la simulation (entre crochets, on met l'ordre d'approximation spatial)
  NonLinearMaxwellProblem<LaplaceEquationAxi> var;  

  // simulation temporelle
  RunTimeProblem(var, input_file, type_element, type_equation);
  
  return FinalizeMontjoie();
}
