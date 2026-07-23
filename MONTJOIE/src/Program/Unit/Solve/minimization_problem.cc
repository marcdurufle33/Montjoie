#include "Share/MontjoieCommon.hxx"
#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

class BasicExample : public VirtualMinimizedFunction<Real_wp>
{
  
public:
  BasicExample()
  {
    this->n = 2;
  }

  void FindInitGuess(Vector<Real_wp>& x)
  {
    x.Reallocate(this->n);
    x(0) = -1;
    x(1) = -2;
  }

  void EvaluateFunction(const Vector<Real_wp>& x, Real_wp& f)
  {
    f = 100*pow(x(0)+3, 4) + pow(x(1)-3, 4);
  }
  
  void EvaluateFunctionGradient(const Vector<Real_wp>& x,
				Real_wp& feval, Vector<Real_wp>& fjac)
  {
    feval = 100*pow(x(0)+3, 4) + pow(x(1)-3, 4);
    fjac(0) = 400*pow(x(0)+3, 3);
    fjac(1) = 4*pow(x(1)-3, 3);    
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  BasicExample fct;
  fct.SetGslAlgorithm(fct.CG);
  
  VectReal_wp xsol;
  MinimizeParametersAlglib(fct, xsol);
  DISP(xsol);
  
  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
