#include "SeldonLib.hxx"

using namespace Seldon;

#include "scatter_matrix.cxx"

void CheckLinearEigenProblem(int rank_proc, int nb_proc)
{
  int nb_mat = 2, nL3 = 2;
  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse> > vec_Ai(nb_mat);
  vec_Ai(0).ReadText("test/matrix/KhLin.dat");
  vec_Ai(1).ReadText("test/matrix/MhLin.dat");

  int color = rank_proc / nL3;
  MPI_Comm comm_solver; int rank_procL3 = 0, nb_procL3 = 1;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank_proc, &comm_solver);
  MPI_Comm_rank(comm_solver, &rank_procL3);
  MPI_Comm_size(comm_solver, &nb_procL3);
  
  DistributedMatrixIntegerArray info;
  ScatterMatrix(vec_Ai, info, nb_procL3, rank_procL3, comm_solver);
  
  SparseEigenProblem<complex<double>, DistributedMatrix<complex<double>, General, ArrayRowSparse>,
                     DistributedMatrix<complex<double>, General, ArrayRowSparse> > var_eig;

  var_eig.SetGlobalCommunicator(MPI_COMM_WORLD);
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(20);
  var_eig.InitMatrix(vec_Ai(0), vec_Ai(1));
  var_eig.SetPrintLevel(1);
  
  complex<double> center(0.1, 2.6);
  FeastParam& param = var_eig.GetFeastParameters();
  param.SetCircleSpectrum(center, 0.2);
  
  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;

  // eigenvalues to find :
  // (0.0631078,2.63279)	(0.0351263,2.54875)	(0.0897294,2.68955)	(0.0348643,2.53424)	(0.0503969,2.7758)	(0.0515185,2.7926)	(0.0782493,2.37889)	(0.0631624,2.34089)	(0.022766,2.30919)	(0.0705122,2.90859)	(0.0244425,2.25785)	(0.104481,2.97512)	(0.0686265,2.99333)	(0.112796,2.99834)	(0.0722426,3.03467)	(0.0134148,2.05918)	(0.0884987,3.16839)	(0.0640078,2.03181)	(0.0895122,3.20105)	(0.119927,3.23721)
  FindEigenvaluesFeast(var_eig, lambda, lambda_imag, eigen_vec);
  if (rank_proc == 0)
    cout << "Eigenvalues = " << lambda << endl;
}

void CheckPolynomialEigenProblem(int rank_proc, int nb_proc)
{
  int nb_mat = 3, nL3 = 2;
  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse> > vec_Ai(nb_mat);
  vec_Ai(0).ReadText("test/matrix/KhPol.dat");
  vec_Ai(1).ReadText("test/matrix/ShPol.dat");
  vec_Ai(2).ReadText("test/matrix/MhPol.dat");

  int color = rank_proc / nL3;
  MPI_Comm comm_solver; int rank_procL3 = 0, nb_procL3 = 1;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank_proc, &comm_solver);
  MPI_Comm_rank(comm_solver, &rank_procL3);
  MPI_Comm_size(comm_solver, &nb_procL3);

  DistributedMatrixIntegerArray info;
  ScatterMatrix(vec_Ai, info, nb_procL3, rank_procL3, comm_solver);

  Vector<DistributedMatrix<complex<double>, General, ArrayRowSparse>* > list_op(2);
  list_op(0) = &vec_Ai(0);
  list_op(1) = &vec_Ai(1);

  PolynomialSparseEigenProblem<complex<double>, DistributedMatrix<complex<double>, General, ArrayRowSparse>,
                               DistributedMatrix<complex<double>, General, ArrayRowSparse> > var_eig;

  var_eig.SetGlobalCommunicator(MPI_COMM_WORLD);
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(10);
  var_eig.InitMatrix(list_op, vec_Ai(2));
  var_eig.SetPrintLevel(1);
  
  complex<double> center(0.4, -0.15);
  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, center);

  FeastParam& param = var_eig.GetFeastParameters();
  param.SetCircleSpectrum(center, 0.08);
  
  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;

  // eigenvalues to find :
  // (0.417693,-0.171441)	(0.439179,-0.170513)	(0.381409,-0.182114)	(0.344035,-0.167351)	(0.427752,-0.0142588)	(0.378343,-0.0122269)	(0.41907,-0.00993428)	(0.451252,-0.178169)	(0.372822,-0.00970444)	(0.40334,-0.0055489)
  FindEigenvaluesFeast(var_eig, lambda, lambda_imag, eigen_vec);
  if (rank_proc == 0)
    cout << "Eigenvalues = " << lambda << endl;
}


int main(int argc, char **argv)
{
  int rank_proc = 0, nb_proc = 1;
  
#ifdef SELDON_WITH_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#endif
  
  srand(0);
  // testing linear eigenvalue problem
  CheckLinearEigenProblem(rank_proc, nb_proc);

  srand(0);
  // testing polynomial eigenvalue problem
  CheckPolynomialEigenProblem(rank_proc, nb_proc);

  if (rank_proc == 0)
    cout << endl << "Feast tested successfully" << endl;
  
  SlepcFinalize();
      
#ifdef SELDON_WITH_MPI
  MPI_Finalize();
#endif

  return 0;
}
