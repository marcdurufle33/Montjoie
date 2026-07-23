#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

template<class MatrixSparse>
void ComputeHelmholtzMatrix(MatrixSparse& A, double dx, int Nx, int Ny, double w2)
{
  double coef = 0.25*w2*dx*dx;
  for (int i = 0; i < Nx-1; i++)
    for (int j = 0; j < Ny-1; j++)
      {
	int num_row = j*Nx + i;
	// interaction between node (i, j) and other nodes
        A.AddInteraction(num_row, num_row, 1.0-coef);
        A.AddInteraction(num_row, num_row+1, -0.5);
        A.AddInteraction(num_row, num_row+Nx, -0.5);
        	
	// interaction between node (i+1, j) and other nodes
        A.AddInteraction(num_row+1, num_row+1, 1.0-coef);
        A.AddInteraction(num_row+1, num_row, -0.5);
        A.AddInteraction(num_row+1, num_row+Nx+1, -0.5);

	// interaction between node (i+1, j+1) and other nodes
        A.AddInteraction(num_row+Nx+1, num_row+Nx+1, 1.0-coef);
        A.AddInteraction(num_row+Nx+1, num_row+1, -0.5);
        A.AddInteraction(num_row+Nx+1, num_row+Nx, -0.5);
        
	// interaction between node (i, j+1) and other nodes
        A.AddInteraction(num_row+Nx, num_row+Nx, 1.0-coef);
        A.AddInteraction(num_row+Nx, num_row, -0.5);
        A.AddInteraction(num_row+Nx, num_row+Nx+1, -0.5);
      }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  {
  int nb_processors, rank_processor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);
  
  double w = 0.31*M_PI, w2 = w*w;
  double L = 10.0;
  int nb_points_glob = 5;
  int Nglob = nb_points_glob*nb_points_glob;
  double dx = L/(nb_points_glob + 1);
  


  // splitting the 2-D domain between processors
  int nb_proc_row = 2, nb_proc_col = 2;
  if (nb_proc_row*nb_proc_col != nb_processors)
    {
      cout << "Run it with the correct number of processors" << endl;
      abort();
    }
  
  int num_proc_row = rank_processor/nb_proc_col;
  int num_proc_col = rank_processor%nb_proc_col;
  //DISP(num_proc_row); DISP(num_proc_col);
  int nb_row = nb_points_glob/nb_proc_row;
  int i0_proc = num_proc_row*nb_row;
  int i1_proc = (num_proc_row+1)*nb_row+1;
  if (num_proc_row == nb_proc_row-1)
    i1_proc = nb_points_glob;

  int j0_proc = num_proc_col*nb_row;
  int j1_proc = (num_proc_col+1)*nb_row+1;
  if (num_proc_col == nb_proc_col-1)
    j1_proc = nb_points_glob;

  //DISP(i0_proc); DISP(i1_proc); DISP(j0_proc); DISP(j1_proc);
  
  int N = (i1_proc-i0_proc)*(j1_proc-j0_proc);
  //DISP(N); DISP(dx);
  


  // global row numbers
  IVect GlobalRowNumbers(N);
  int Nx = i1_proc-i0_proc, Ny = j1_proc-j0_proc;
  for (int i = i0_proc; i < i1_proc; i++)
    for (int j = j0_proc; j < j1_proc; j++)
      GlobalRowNumbers((j-j0_proc)*Nx + i-i0_proc) = j*nb_points_glob + i;
  
  //DISP(GlobalRowNumbers);
  IVect ListProcInterac;
  Vector<IVect> ListRowInterac;
  IVect OverlapRowNumbers, OriginalProcNumbers;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  if (true)
    {
      // first choice : the arrays ListProcInterac, ListRowInterac, 
      // OverlapRowNumbers and OriginalProcNumbers
      // are recovered from the global row numbers
      
      // with this choice, you only need to give the global row numbers for each processor
      // however, more communications and memory will be needed to construct the arrays ListProcInterac, etc
      
      // GlobalRowNumbers should be sorted before calling this function
      DistributedMatrix<double, General, ArrayRowSparse> Atest;
      Atest.Init(GlobalRowNumbers, OverlapRowNumbers, OriginalProcNumbers,
                 ListProcInterac, ListRowInterac, comm);
    }
  else
    {
      // second choice : the arrays ListProcInterac, ListRowInterac, 
      // OverlapRowNumbers and OriginalProcNumbers
      // are constructed manually with the knowledge of the geometry
      
      // Here, since we are using finite difference discretization
      // it is easy to know the numbers of overlapped rows
      
      // row numbers that are "overlapped" by another processor
      int Ns = 0;
      if (num_proc_row > 0)
        {
          Ns = Ny;
          if (num_proc_col > 0)
            Ns += Nx-1;
        }
      else if (num_proc_col > 0)
        Ns = Nx;
      
      OverlapRowNumbers.Reallocate(Ns); OriginalProcNumbers.Reallocate(Ns);
      OverlapRowNumbers.Fill(-1); OriginalProcNumbers.Fill(-1);
      if (num_proc_col > 0)
        {
          for (int i = 0; i < Nx; i++)
            {
              OverlapRowNumbers(i) = i;
              OriginalProcNumbers(i) = rank_processor - 1;
            }
          
          if (num_proc_row > 0)
            {
              OriginalProcNumbers(0) = rank_processor - 1 - nb_proc_col;
              for (int i = 1; i < Ny; i++)
                {
                  OverlapRowNumbers(Nx+i-1) = i*Nx;
                  OriginalProcNumbers(Nx+i-1) = rank_processor - nb_proc_col;
                }
            }
        }
      else if (num_proc_row > 0)
        {
          for (int i = 0; i < Ny; i++)
            {
              OverlapRowNumbers(i) = i*Nx;
              OriginalProcNumbers(i) = rank_processor - nb_proc_col;
            }      
        }
      
      //DISP(OverlapRowNumbers); DISP(OriginalProcNumbers);
      
      // dofs shared with other processors
      int nb_proc_interac = 0;
      if (num_proc_col > 0)
        {
          nb_proc_interac++;
          if (num_proc_row > 0)
            nb_proc_interac++;
          
          if (num_proc_row < nb_proc_row-1)
            nb_proc_interac++;      
        }
      
      if (num_proc_col < nb_proc_col-1)
        {
          nb_proc_interac++;
          if (num_proc_row > 0)
            nb_proc_interac++;
          
          if (num_proc_row < nb_proc_row-1)
            nb_proc_interac++;      
        }
      
      if (num_proc_row > 0)
        nb_proc_interac++;
      
      if (num_proc_row < nb_proc_row-1)
        nb_proc_interac++;
      
      ListProcInterac.Reallocate(nb_proc_interac);
      ListRowInterac.Reallocate(nb_proc_interac);
      
      nb_proc_interac = 0;
      if (num_proc_col > 0)
        {
          ListProcInterac(nb_proc_interac) = rank_processor-1;
          ListRowInterac(nb_proc_interac).Reallocate(Nx);
          for (int i = 0; i < Nx; i++)
            ListRowInterac(nb_proc_interac)(i) = i;
          
          nb_proc_interac++;
          
          if (num_proc_row > 0)
            {
              ListProcInterac(nb_proc_interac) = rank_processor-1-nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = 0;
              nb_proc_interac++;
            }
          
          if (num_proc_row < nb_proc_row-1)
            {
              ListProcInterac(nb_proc_interac) = rank_processor-1+nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx-1;
              nb_proc_interac++;      
            }
        }
      
      if (num_proc_col < nb_proc_col-1)
        {
          ListProcInterac(nb_proc_interac) = rank_processor+1;
          ListRowInterac(nb_proc_interac).Reallocate(Nx);
          for (int i = 0; i < Nx; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*(Ny-1) + i;
          
          nb_proc_interac++;
          
          if (num_proc_row > 0)
            {
              ListProcInterac(nb_proc_interac) = rank_processor+1-nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx*(Ny-1);
              nb_proc_interac++;
            }
          
          if (num_proc_row < nb_proc_row-1)
            {
              ListProcInterac(nb_proc_interac) = rank_processor+1+nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx*Ny-1;
              nb_proc_interac++;      
            }
        }
      
      if (num_proc_row > 0)
        {
          ListProcInterac(nb_proc_interac) = rank_processor-nb_proc_col;
          ListRowInterac(nb_proc_interac).Reallocate(Ny);
          for (int i = 0; i < Ny; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*i;
          
          nb_proc_interac++;
        }
      
      if (num_proc_row < nb_proc_row-1)
        {
          ListProcInterac(nb_proc_interac) = rank_processor+nb_proc_col;
          ListRowInterac(nb_proc_interac).Reallocate(Ny);
          for (int i = 0; i < Ny; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*i+Nx-1;
          
          nb_proc_interac++;
        }
      
      for (int j = 0; j < nb_proc_interac; j++)
        {
          //DISP(j);
          //DISP(ListProcInterac(j));
          //DISP(ListRowInterac(j));
        }
    }
  
  // computes Helmholtz matrix
  DistributedMatrix<double, General, ArrayRowSparse> A;
  A.Reallocate(N, N);
  A.Init(Nglob, &GlobalRowNumbers, &OverlapRowNumbers, &OriginalProcNumbers,
	 N, 1, &ListProcInterac, &ListRowInterac, comm);
  
  ComputeHelmholtzMatrix(A, dx, i1_proc-i0_proc, j1_proc-j0_proc, w2);
  A.WriteText("mat.dat");
  
  // factorizes matrix
  bool iterative = false;
  SparseDistributedSolver<double> mat_lu;
  mat_lu.ShowFullHistory();
  //mat_lu.SelectOrdering(SparseMatrixOrdering::METIS);
  
  MPI_Barrier(comm);
  if (!iterative)
    mat_lu.Factorize(A);
  
  // computes right hand side (dirac at the center)
  Vector<double> x_sol(N);
  x_sol.Fill(0.0);
  int iglob0 = nb_points_glob/2;
  if ((iglob0 >= i0_proc) && (iglob0 < i1_proc))
    if ((iglob0 >= j0_proc) && (iglob0 < j1_proc))
      {
	int i = iglob0 - i0_proc;
	int j = iglob0 - j0_proc;
	x_sol(j*Nx + i) = 1.0;
      }
  
  // solves linear sytem  
  if (!iterative)
    mat_lu.Solve(x_sol);
  else
    {
      DistributedVector<double> source(OverlapRowNumbers, comm), sol(OverlapRowNumbers, comm);
      
      source.SetData(x_sol);
      sol.Reallocate(N);

      Preconditioner_Base<double> prec;
      Iteration<double> iter(10000, 1e-6);
      iter.HideMessages();
      if (rank_processor == 0)
        iter.ShowFullHistory();
      
      Cg(A, sol, source, prec, iter);
      for (int i = 0; i < N; i++)
        x_sol(i) = sol(i);
    }
  
  if (rank_processor == 0)
    cout << "Solution computed " << endl;
  
  // writes solution
  if (rank_processor != 0)
    {
      IVect num(4);
      num(0) = i0_proc; num(1) = i1_proc; num(2) = j0_proc; num(3) = j1_proc;
      MPI_Ssend(num.GetData(), 4, MPI_INTEGER, 0, 2, comm);
      
      MPI_Ssend(x_sol.GetData(), x_sol.GetM(), MPI_DOUBLE, 0, 3, comm);
    }
  else
    {
      Matrix<double> x_mat(nb_points_glob, nb_points_glob);
      for (int i = i0_proc; i < i1_proc; i++)
	for (int j = j0_proc; j < j1_proc; j++)
	  x_mat(i, j) = x_sol((j-j0_proc)*Nx + i-i0_proc);
      
      IVect num(4);
      Vector<double> xtmp; MPI_Status status;
      for (int rank = 1; rank < nb_processors; rank++)
	{
	  MPI_Recv(num.GetData(), 4, MPI_INTEGER, rank, 2, comm, &status);
	  int i0 = num(0), i1 = num(1), j0 = num(2), j1 = num(3);
	  xtmp.Reallocate((i1-i0)*(j1-j0));
	  MPI_Recv(xtmp.GetData(), xtmp.GetM(), MPI_DOUBLE, rank, 3, comm, &status);
	  
	  for (int i = i0; i < i1; i++)
	    for (int j = j0; j < j1; j++)
	      x_mat(i, j) = xtmp((j-j0)*(i1-i0) + i-i0);
	}
      
  
      x_mat.WriteText("Sol.dat");
    }
  
  MPI_Barrier(comm);
  cout << "Destroying variables" << endl;
  }
  
  return FinalizeMontjoie();
}
