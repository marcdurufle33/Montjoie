#include "Montjoie.hxx"

using namespace Montjoie;

class NonLinearSplitFormulation : public NonLinearEigenProblem_Base<Complex_wp>
{
  Vector<Matrix<Complex_wp, General, ArrayRowSparse> > Ai;

public:
  NonLinearSplitFormulation(int n)
  {
    Ai.Reallocate(n);
  }

  void SetOperator(int i, const string& name)
  {
    ReadMatrixMarket(name, Ai(i));
    this->Init(Ai(i).GetM());
  }
  
  void ComputeOperatorSplitExplicit(int i, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A)
  {
    Copy(Ai(i), A);
  }
  
};
  
int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  Matrix<complex<double> > A(6, 4);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      A(i, j) = Complex_wp(rand(), rand()) / Real_wp(RAND_MAX);
  
  A.Write("A.dat");
  GetPseudoInverse(A, 1e-15); A.Write("Ainv.dat");
  exit(0);
  
  if (argc < 3)
    {
      cout << "Provide at least two arguments" << endl;
      cout << "Usage ./check.x nb_mat  file_coef.dat" << endl;
      abort();
    }

  int nb_mat = atoi(argv[1]);
  NonLinearSplitFormulation var(nb_mat);
  var.SetSplitMatrices(nb_mat);
  
  ifstream file_in(argv[2]);
  int deg_n; double nr, ni, omega, dr, di;
  string nom;
  for (int i = 0; i < nb_mat; i++)
    {
      file_in >> deg_n >> nr >> ni >> omega >> dr >> di >> nom;
      VectComplex_wp coef_num(deg_n+1), coef_denom(3);
      coef_num.Zero(); coef_denom.Zero();
      coef_num(0) = Complex_wp(nr, ni);
      if (omega == 0)
        {
          coef_denom.Clear();
        }
      else
        {
          coef_denom(0) = omega*omega;
          coef_denom(2) = Complex_wp(dr, di);
        }
      
      DISP(coef_num); DISP(coef_denom);

      var.SetOperator(i, nom);
      var.SetNumeratorSplitFct(i, coef_num);
      var.SetDenominatorSplitFct(i, coef_denom);
    }
  
  var.SetExplicitMatrix();
  var.SetStoppingCriterion(1e-12);
  var.SetNbMaximumIterations(1000);

  var.SetTypeSpectrum(var.CENTERED_EIGENVALUES, Complex_wp(0.58, 0.02));

  SlepcParamNep& param = var.GetSlepcParameters();
  param.SetEigensolverType(param.NLEIGS);
  param.SetIntervalRegion(0.1, 10.0, 0.01, 10.0);
  param.EnableCommandLineOptions();
  param.SetDefaultPetscSolver();

  VectComplex_wp eigen_val, lambda_imag;
  Matrix<Complex_wp, General, ColMajor> eigen_vec;
  FindEigenvaluesSlepc(var, eigen_val, lambda_imag, eigen_vec);
  DISP(eigen_val);
  
  return FinalizeMontjoie();
}
