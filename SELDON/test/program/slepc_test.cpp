#include "SeldonLib.hxx"

using namespace Seldon;

#include "scatter_matrix.cxx"

void ComputeEigenvalues(Vector<DistributedMatrix<complex<double>, Symmetric, ArrayRowSymSparse> >& vec_Ai,
                        Vector<Vector<complex<double> > >& coef_num,
                        Vector<Vector<complex<double> > >& coef_denom, int rank_proc)
{
  SplitSparseNonLinearEigenProblem<complex<double>, Symmetric, ArrayRowSymSparse> var_eig;

  var_eig.SetPrintLevel(4);
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(10);
  var_eig.InitMatrix(vec_Ai, coef_num, coef_denom);
  
  complex<double> center(0.4, -0.1);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, center);
  SlepcParamNep& param = var_eig.GetSlepcParameters();
  param.SetEigensolverType(param.NLEIGS);
  
  param.SetIntervalRegion(0.3, 0.6, -0.3, -1e-4);
  param.EnableCommandLineOptions();

  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

  if (rank_proc == 0)
    cout << "Eigenvalues = " << lambda << endl;
}

void CheckLinearEigenProblem(int rank_proc, int nb_proc)
{
  int nb_mat = 2;
  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse> > vec_Ai(nb_mat);
  vec_Ai(0).ReadText("test/matrix/KhLin.dat");
  vec_Ai(1).ReadText("test/matrix/MhLin.dat");
  
  DistributedMatrixIntegerArray info;
  ScatterMatrix(vec_Ai, info, nb_proc, rank_proc, MPI_COMM_WORLD);

  SparseEigenProblem<complex<double>, DistributedMatrix<complex<double>, General, ArrayRowSparse>,
                     DistributedMatrix<complex<double>, General, ArrayRowSparse> > var_eig;
  
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(20);
  var_eig.InitMatrix(vec_Ai(0), vec_Ai(1));
  var_eig.SetPrintLevel(4);
  
  complex<double> center(0.1, 2.6);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, center);
  var_eig.SetComputationalMode(var_eig.INVERT_MODE);
  
  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;

  // eigenvalues to find :
  // (0.0631078,2.63279)	(0.0351263,2.54875)	(0.0897294,2.68955)	(0.0348643,2.53424)	(0.0503969,2.7758)	(0.0515185,2.7926)	(0.0782493,2.37889)	(0.0631624,2.34089)	(0.022766,2.30919)	(0.0705122,2.90859)	(0.0244425,2.25785)	(0.104481,2.97512)	(0.0686265,2.99333)	(0.112796,2.99834)	(0.0722426,3.03467)	(0.0134148,2.05918)	(0.0884987,3.16839)	(0.0640078,2.03181)	(0.0895122,3.20105)	(0.119927,3.23721)
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);
  if (rank_proc == 0)
    cout << "Eigenvalues = " << lambda << endl;
}

void CheckPolynomialEigenProblem(int rank_proc, int nb_proc)
{
  int nb_mat = 3;
  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse> > vec_Ai(nb_mat);
  vec_Ai(0).ReadText("test/matrix/KhPol.dat");
  vec_Ai(1).ReadText("test/matrix/ShPol.dat");
  vec_Ai(2).ReadText("test/matrix/MhPol.dat");
  
  DistributedMatrixIntegerArray info;
  ScatterMatrix(vec_Ai, info, nb_proc, rank_proc, MPI_COMM_WORLD);

  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse>* > list_op(2);
  list_op(0) = &vec_Ai(0);
  list_op(1) = &vec_Ai(1);

  PolynomialSparseEigenProblem<complex<double>, DistributedMatrix<complex<double>, General, ArrayRowSparse>,
                               DistributedMatrix<complex<double>, General, ArrayRowSparse> > var_eig;
  
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(10);
  var_eig.InitMatrix(list_op, vec_Ai(2));
  var_eig.SetPrintLevel(4);
  
  complex<double> center(0.4, -0.1);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, center);
  var_eig.SetSpectralTransformation(true);
  
  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;

  // eigenvalues to find :
  // (0.417693,-0.171441)	(0.439179,-0.170513)	(0.381409,-0.182114)	(0.344035,-0.167351)	(0.427752,-0.0142588)	(0.378343,-0.0122269)	(0.41907,-0.00993428)	(0.451252,-0.178169)	(0.372822,-0.00970444)	(0.40334,-0.0055489)
  FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);
  if (rank_proc == 0)
    cout << "Eigenvalues = " << lambda << endl;
}

void CheckNonLinearEigenProblem(int rank_proc, int nb_proc)
{
  int nb_mat = 4;
  Vector<DistributedMatrix<complex<double>, Symmetric, ArrayRowSymSparse> > vec_Ai(nb_mat);
  Vector<Vector<complex<double> > > numer(nb_mat), denom(nb_mat);
  vec_Ai(0).ReadText("test/matrix/Ai0.dat");
  vec_Ai(1).ReadText("test/matrix/Ai1.dat");
  vec_Ai(2).ReadText("test/matrix/Ai2.dat");
  vec_Ai(3).ReadText("test/matrix/Ai3.dat");
  int N = vec_Ai(0).GetM();
  if ((vec_Ai(1).GetM() != N) || (vec_Ai(1).GetM() != N) || (vec_Ai(1).GetM() != N))
    {
      cout << "Matrices are not of the same size" << endl;
      abort();
    }
  
  numer(0).Reallocate(3); numer(0).Zero(); numer(0)(0) = complex<double>(1, 0);
  numer(1).Reallocate(1); numer(1)(0) = complex<double>(1, 0);
  numer(2).Reallocate(4); numer(2).Zero();
  numer(2)(0) = complex<double>(0,342.837335497338);
  numer(2)(1) = complex<double>(-507.774837563131,0);
  denom(2).Reallocate(3);
  denom(2)(0) = complex<double>(39.4784176043574,0);
  denom(2)(1) = complex<double>(0,32.4875550324677);
  denom(2)(2) = complex<double>(-38.4093314884154,0);
  numer(3).Reallocate(4); numer(3).Zero();
  numer(3)(0) = complex<double>(0,76.8074081861214);
  numer(3)(1) = complex<double>(-23.3687323245872,0);
  denom(3).Reallocate(3);
  denom(3)(0) = complex<double>(39.4784176043574,0);
  denom(3)(1) = complex<double>(0,6.26767912318158);
  denom(3)(2) = complex<double>(-28.0097475738552,0);

  // eigenvalues to find :
  // (0.417693,-0.171441)	(0.439179,-0.170513)	(0.381409,-0.182114)	(0.344035,-0.167351)	(0.427752,-0.0142588)	(0.378343,-0.0122269)	(0.41907,-0.00993428)	(0.451252,-0.178169)	(0.372822,-0.00970444)	(0.40334,-0.0055489)	(0.405646,-0.00559058)
  
  DistributedMatrixIntegerArray info;
  ScatterMatrix(vec_Ai, info, nb_proc, rank_proc, MPI_COMM_WORLD);
  
  ComputeEigenvalues(vec_Ai, numer, denom, rank_proc);
}

int main(int argc, char **argv)
{
  InitSeldon(argc, argv);
  
  int rank_proc = 0, nb_proc = 1;
  
#ifdef SELDON_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#endif

  srand(0);
  // testing linear eigenvalue problem
  CheckLinearEigenProblem(rank_proc, nb_proc);

  srand(0);
  // testing polynomial eigenvalue problem
  CheckPolynomialEigenProblem(rank_proc, nb_proc);

  srand(0);
  // testing non-linear eigenvalue problem
  CheckNonLinearEigenProblem(rank_proc, nb_proc);

  cout << endl << "Non-linear eigenvalue problem tested successfully" << endl;
  
  return FinalizeSeldon();
}
