#include "Share/MontjoieCommon.hxx"
#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

class ElectrostaticProblem
  : public NonLinearEquations_Newton<Real_wp, R3, R3>
{
  int nb_points;
  TinyMatrix<Real_wp, General, 3, 3> DF;

public:
  ElectrostaticProblem(int n)
  {
    nb_points = n;
  }

  void ComputeScheme(const R3& x, R3& F)
  {
    F.Zero();
    Real_wp one(1);
    for (int i = 0; i < nb_points; i++)
      {
	F(i) = one / (one + x(i)) - one / (one - x(i));
	for (int j = 0; j < nb_points; j++)
	  if (j != i)
	    F(i) += (x(i) - x(j)) / square(x(i) - x(j));
      }
  }

  void ComputeAndFactoriseDiff(const R3& x, const R3& scale)
  {
    Real_wp one(1);
    // DF.Reallocate(nb_points, nb_points);
    DF.Zero();
    for (int i = 0; i < nb_points; i++)
      {
	DF(i, i) = -one / square(one + x(i)) - one / square(one - x(i));
	for (int j = 0; j < nb_points; j++)
	  if (j != i)
	    {
	      DF(i, i) -= one / square(x(i) - x(j));
	      DF(i, j) = one / square(x(i) - x(j));
	    }
      }
    
    GetInverse(DF);
  }

  void SolveDifferential(const R3& b, R3& x)
  {
    // x = b;
    // SolveLU(DF, pivot, x);
    Mlt(DF, b, x);
  }

  Real_wp GetNorm2Vector(const R3& x)
  {
    return Norm2(x);
  }

};


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  ElectrostaticProblem var(3);  
  NewtonSolver<Real_wp, R3, R3> solver;
  solver.ForceReevaluationJacobian();
  solver.SetMaxNumberOfIterations(8);
  solver.SetPrintLevel(2);
  
  TinyVector<Real_wp, 3> x, scale;
  x(0) = -0.4; x(1) = 0; x(2) = 0.4;
  solver.Init(var, x);
  solver.Solve(var, x);
  
  DISP(x);

  return FinalizeMontjoie();
}
