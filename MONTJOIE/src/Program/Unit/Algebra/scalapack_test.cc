#define SELDON_WITH_SCALAPACK

#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

Real_wp threshold;
int rank_processor, nb_processors;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T, class Prop, class Storage>
void GenerateRandomMatrix(Matrix<T, Prop, Storage>& A, int m, int n)
{
  A.Reallocate(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      GetRandNumber(A(i, j));
}

template<class T>
bool CheckVector(Vector<T>& x)
{
  bool test = true;
  T icplx, coef;
  SetComplexReal(x.GetM(), icplx);
  SetComplexOne(coef); coef /= icplx;
  for (int i = 0; i < x.GetM(); i++)
    {
      SetComplexReal(i, icplx);
      if ((abs(x(i) - icplx)*coef > threshold) || isnan(x(i)))
        test = false;
    }

  return test;
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(x(i)) || isnan(y(i)))
      return false;
  
  return true;
}

template<class T, class Prop, class Storage,
         class T2, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop, Storage>& x,
                 const Matrix<T2, Prop2, Storage2>& y,
		 Real_wp eps = threshold)
{
  if ((x.GetM() != y.GetM()) || (x.GetN() != y.GetN()))
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    for (int j = 0; j < x.GetN(); j++)
      if ((abs(x(i, j) - y(i, j)) > eps) || isnan(x(i, j)) || isnan(y(i, j)))
        return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckDenseCholesky(Matrix<T, Prop, Storage, Allocator>& A)
{
  Matrix<T> Bdense, Cdense;
  Matrix<T, Prop, Storage, Allocator> Adense;
  Vector<T> xglob, bglob;
  
  T zero, one; SetComplexZero(zero); SetComplexOne(one);
  
  int n = 50;
  GenerateRandomMatrix(Cdense, n, n);
  GenerateRandomVector(xglob, n);
  bglob.Reallocate(n);
  
  Bdense.Reallocate(n, n);
  Bdense.Zero();
  
  MltAdd(one, SeldonNoTrans, Cdense, SeldonTrans, Cdense, zero, Bdense);

  Adense.Reallocate(n, n);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      Adense(i, j) = Bdense(i, j);
  
  DistributedMatrix<double, General, ColMajor> P, x, y, b;
  
  int nb = 5;
  P.Init(global_blacs_handle, n, n, nb, nb);
  x.Init(global_blacs_handle, n, 1, nb, nb);
  y.Init(global_blacs_handle, n, 1, nb, nb);
  b.Init(global_blacs_handle, n, 1, nb, nb);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      P.SetGlobal(i, j, Adense(i, j));

  VectReal_wp wr, wi;
  GetEigenvalues(P, wr, wi);

  if (rank_processor == 0)
    { DISP(wr); DISP(wi); }
  /* 
  GetCholesky(P);

  Mlt(Bdense, xglob, bglob);
  
  DistributeVector(bglob, b);
  DistributeVector(xglob, y);

  DISP(b);

  x = b;
  SolveCholesky(P, x);
  DISP(x); DISP(y);
  if (!EqualMatrix(x, y, 1e2*threshold))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }
  */
}

template<class T>
void CheckDenseRowMajor(DistributedMatrix<T, General, RowMajor>& A)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  
  int rank_proc;
  MPI_Comm_rank(comm, &rank_proc);

  srand(1+3*rank_proc);

  int n = 10;
  int m = 2+rank_proc;
  
  A.SetCommunicator(comm);
  A.Reallocate(m, n);
  A.FillRand(); A *= 1e-9;
  
  A.WriteText("mat" + to_str(rank_proc) + ".dat");

  Vector<T> X(m), Y(n);
  X.FillRand(); X *= 1e-9;
  
  X.WriteText("vec" + to_str(rank_proc) + ".dat");

  Mlt(SeldonTrans, A, X, Y);
  
  Y.WriteText("res" + to_str(rank_proc) + ".dat");

  T alpha(0.23), beta(0.082);
  MltAdd(alpha, SeldonTrans, A, X, beta, Y);

  Y.WriteText("resA" + to_str(rank_proc) + ".dat");

  Matrix<T, Symmetric, RowSymPacked> C(n, n);
  
  T one, zero; SetComplexOne(one); SetComplexZero(zero);
  MltAdd(one, SeldonTrans, A, SeldonNoTrans, A, zero, C);

  C.WriteText("Ca" + to_str(rank_proc) + ".dat");
  
  C.Zero();
  MltAdd(alpha, SeldonTrans, A, SeldonNoTrans, A, zero, C);
  
  C.WriteText("Cb" + to_str(rank_proc) + ".dat");

  C.SetIdentity();

  MltAdd(alpha, SeldonTrans, A, SeldonNoTrans, A, beta, C);
  
  C.WriteText("Cc" + to_str(rank_proc) + ".dat");

}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 5e3*epsilon_machine;
  DISP(threshold);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);
  
  {
    //DistributedMatrix<Real_wp, General, RowMajor> A;
    //CheckDenseRowMajor(A);
  }

  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckDenseCholesky(A);
  }

}
