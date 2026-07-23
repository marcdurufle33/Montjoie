#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

template<class MatrixSparse>
void ComputeLaplaceMatrix(MatrixSparse& A, double dx, int Nx, int Ny, Vector<double>& mass)
{  
  for (int i = 0; i < Nx-1; i++)
    for (int j = 0; j < Ny-1; j++)
      {
	int num_row = j*Nx + i;
	// interaction between node (i, j) and other nodes
        mass(num_row) += 0.25*dx*dx;
        A.AddInteraction(num_row, num_row, 1.0);
        A.AddInteraction(num_row, num_row+1, -0.5);
        A.AddInteraction(num_row, num_row+Nx, -0.5);
        	
	// interaction between node (i+1, j) and other nodes
        mass(num_row+1) += 0.25*dx*dx;
        A.AddInteraction(num_row+1, num_row+1, 1.0);
        A.AddInteraction(num_row+1, num_row, -0.5);
        A.AddInteraction(num_row+1, num_row+Nx+1, -0.5);

	// interaction between node (i+1, j+1) and other nodes
        mass(num_row+Nx+1) += 0.25*dx*dx;
        A.AddInteraction(num_row+Nx+1, num_row+Nx+1, 1.0);
        A.AddInteraction(num_row+Nx+1, num_row+1, -0.5);
        A.AddInteraction(num_row+Nx+1, num_row+Nx, -0.5);
        
	// interaction between node (i, j+1) and other nodes
        mass(num_row+Nx) += 0.25*dx*dx;
        A.AddInteraction(num_row+Nx, num_row+Nx, 1.0);
        A.AddInteraction(num_row+Nx, num_row, -0.5);
        A.AddInteraction(num_row+Nx, num_row+Nx+1, -0.5);
      }
}

void CreateLaplacianMatrix(DistributedMatrix<complex<double>, General, ArrayRowSparse>& A,
                           DistributedMatrixIntegerArray& info)
{
  int nb_processors, rank_processor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);

  double L = 10.0;
  int nb_points_glob = 10;
  int Nglob = nb_points_glob*nb_points_glob;
  double dx = L/(nb_points_glob - 1);
  
  // splitting the 2-D domain between processors
  int nb_proc_row = 2, nb_proc_col = 2;
  if (nb_proc_row*nb_proc_col != nb_processors)
    {
      cout << "Run it with the correct number of processors" << endl;
      abort();
    }
  
  int num_proc_row = rank_processor/nb_proc_col;
  int num_proc_col = rank_processor%nb_proc_col;

  int nb_row = nb_points_glob/nb_proc_row;
  int i0_proc = num_proc_row*nb_row;
  int i1_proc = (num_proc_row+1)*nb_row+1;
  if (num_proc_row == nb_proc_row-1)
    i1_proc = nb_points_glob;

  int j0_proc = num_proc_col*nb_row;
  int j1_proc = (num_proc_col+1)*nb_row+1;
  if (num_proc_col == nb_proc_col-1)
    j1_proc = nb_points_glob;

  int N = (i1_proc-i0_proc)*(j1_proc-j0_proc);

  // global row numbers
  Vector<int>& GlobalRowNumbers = info.GlobalRowNumbers;
  GlobalRowNumbers.Reallocate(N);
  int Nx = i1_proc-i0_proc; // Ny = j1_proc-j0_proc;
  for (int i = i0_proc; i < i1_proc; i++)
    for (int j = j0_proc; j < j1_proc; j++)
      GlobalRowNumbers((j-j0_proc)*Nx + i-i0_proc) = j*nb_points_glob + i;

  IVect& ListProcInterac = info.ProcSharingRows;
  Vector<IVect>& ListRowInterac = info.SharingRowNumbers;
  IVect& OverlapRowNumbers = info.OverlapRowNumbers;
  IVect& OriginalProcNumbers = info.OverlapProcNumbers;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // the arrays ListProcInterac, ListRowInterac, 
  // OverlapRowNumbers and OriginalProcNumbers
  // are recovered from the global row numbers
  
  // GlobalRowNumbers should be sorted before calling this function
  DistributedMatrix<complex<double>, General, ArrayRowSparse> Atest;
  Atest.Init(GlobalRowNumbers, OverlapRowNumbers, OriginalProcNumbers,
                 ListProcInterac, ListRowInterac, comm);
  
  // computes Laplace matrix
  A.Reallocate(N, N);
  A.Init(Nglob, &GlobalRowNumbers, &OverlapRowNumbers, &OriginalProcNumbers,
	 N, 1, &ListProcInterac, &ListRowInterac, comm);

  Vector<double> mass(N); mass.Zero();
  ComputeLaplaceMatrix(A, dx, i1_proc-i0_proc, j1_proc-j0_proc, mass);

  AssembleVector(mass, MPI_SUM, ListProcInterac, ListRowInterac, MPI_COMM_WORLD, N, 1, 1);

  for (int i = 0; i < N; i++)
    mass(i) = 1.0/mass(i);

  ScaleLeftMatrix(A, mass);
}

void CreateSmallMatrix(DistributedMatrix<complex<double>, General, ArrayRowSparse>& A,
                       DistributedMatrixIntegerArray& info)
{
  int nb_processors, rank_processor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);

  int Nglob = 6;
  int N = 3;
  if (nb_processors == 1)
    N = 6;
  
  Vector<int>& GlobalRowNumbers = info.GlobalRowNumbers;
  GlobalRowNumbers.Reallocate(N);

  // tableaux vides car pas de ddls partages
  IVect& ListProcInterac = info.ProcSharingRows;
  Vector<IVect>& ListRowInterac = info.SharingRowNumbers;
  IVect& OverlapRowNumbers = info.OverlapRowNumbers;
  IVect& OriginalProcNumbers = info.OverlapProcNumbers;
  MPI_Comm comm = MPI_COMM_WORLD;

  GlobalRowNumbers.Fill();
  if (nb_processors == 2)
    {
      if (rank_processor == 1)
        {
          GlobalRowNumbers(0) = 3;
          GlobalRowNumbers(1) = 4;
          GlobalRowNumbers(2) = 5;
          //GlobalRowNumbers(3) = 7;
        }
      //else
        //GlobalRowNumbers(3) = 6;
    }
  
  // on alloue la matrice
  A.Reallocate(N, N);
  A.Init(Nglob, &GlobalRowNumbers, &OverlapRowNumbers, &OriginalProcNumbers,
	 N, 1, &ListProcInterac, &ListRowInterac, comm);
  
  if (nb_processors == 1)
    {
      A.AddInteraction(0, 0, complex<double>(1, 1));
      A.AddInteraction(0, 3, complex<double>(2, 0));
      A.AddInteraction(0, 4, complex<double>(0, 1));
      A.AddInteraction(1, 0, complex<double>(2, 0));
      A.AddInteraction(1, 1, complex<double>(1, 0));
      A.AddInteraction(1, 2, complex<double>(2, 0));
      A.AddInteraction(1, 4, complex<double>(4, 0));
      A.AddInteraction(2, 3, complex<double>(2, 1));
      A.AddInteraction(2, 4, complex<double>(0, 1));

      A.AddInteraction(3, 1, complex<double>(3, 0));
      A.AddInteraction(3, 2, complex<double>(4, 0));
      A.AddInteraction(3, 3, complex<double>(2, 0));

      A.AddInteraction(4, 0, complex<double>(0, 1));
      A.AddInteraction(4, 1, complex<double>(2, 0));
      A.AddInteraction(4, 3, complex<double>(2, 0));

      A.AddInteraction(5, 4, complex<double>(1, 0));
      A.AddInteraction(5, 5, complex<double>(2, 0));
      
      //A.AddInteraction(6, 6, complex<double>(1, 0));
      //A.AddInteraction(7, 7, complex<double>(1, 0));
    }
  else
    {
      if (rank_processor == 0)
        {
          A.AddInteraction(0, 0, complex<double>(1, 1));
          A.AddDistantInteraction(0, 3, 1, complex<double>(2, 0));
          A.AddDistantInteraction(0, 4, 1, complex<double>(0, 1));
          A.AddInteraction(1, 0, complex<double>(2, 0));
          A.AddInteraction(1, 1, complex<double>(1, 0));
          A.AddInteraction(1, 2, complex<double>(2, 0));
          A.AddDistantInteraction(1, 4, 1, complex<double>(4, 0));
          A.AddDistantInteraction(2, 3, 1, complex<double>(2, 1));
          A.AddDistantInteraction(2, 4, 1, complex<double>(0, 1));
          //A.AddInteraction(3, 3, complex<double>(1, 0));
        }
      else
        {
          A.AddDistantInteraction(0, 1, 0, complex<double>(3, 0));
          A.AddDistantInteraction(0, 2, 0, complex<double>(4, 0));
          A.AddInteraction(0, 0, complex<double>(2, 0));
          
          A.AddDistantInteraction(1, 0, 0, complex<double>(0, 1));
          A.AddDistantInteraction(1, 1, 0, complex<double>(2, 0));
          A.AddInteraction(1, 0, complex<double>(2, 0));
          
          A.AddInteraction(2, 1, complex<double>(1, 0));
          A.AddInteraction(2, 2, complex<double>(2, 0));
          //A.AddInteraction(3, 3, complex<double>(1, 0));
        }
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  // on cree la matrice du laplacien
  DistributedMatrixIntegerArray info;
  DistributedMatrix<complex<double>, General, ArrayRowSparse> A;
  
  //CreateLaplacianMatrix(A, info); int nev = 10;
  CreateSmallMatrix(A, info); int nev = 2;
  //DISP(A.GetM()); DISP(A.GetGlobalRowNumber());
  A.WriteText("Lap.dat");
  
  SparseEigenProblem<complex<double>, DistributedMatrix<complex<double>, General, ArrayRowSparse>,
                     DistributedMatrix<double, Symmetric, ArrayRowSymSparse> > var_eig;

  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(nev);
  var_eig.SetNbArnoldiVectors(nev+2);
  var_eig.SetComputationalMode(var_eig.REGULAR_MODE);

  var_eig.InitMatrix(A);

  //var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.1, var_eig.SORTED_MODULUS);
  var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0.0, var_eig.SORTED_MODULUS);

  DISP(var_eig.GetM()); DISP(A.GetM());
  Vector<complex<double> > lambda, lambda_imag;
  Matrix<complex<double>, General, ColMajor> eigen_vec;
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
  DISP(lambda); DISP(lambda_imag);

  
  return FinalizeMontjoie();
}
