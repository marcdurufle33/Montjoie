#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

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
    x = complex<T>(0, rand())/T(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/T(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T>
void GenerateRandomVector(Vector<T, VectSparse>& x, int n, int nnz)
{
  x.Clear();
  for (int i = 0; i < nnz; i++)
    {
      int j = rand()%n;
      GetRandNumber(x.Get(j));
    }
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n)
{
  A.Reallocate(m, n);
  T x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(x);
        A(i, j) = x;
      }
}

template<class T, class Allocator>
void GenerateRandomMatrix(Matrix<T, General, RowSparse, Allocator>& Af,
                          int m, int n, int nnz)
{
  T x;
  Matrix<T, General, ArrayRowSparse, Allocator> A;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      A.Set(i, j, x);
    }
  
  Copy(A, Af);
}

template<class T, class Allocator>
void GenerateRandomMatrix(Matrix<T, General, ArrayRowSparse, Allocator>& A,
                          int m, int n, int nnz)
{
  T x;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      A.Set(i, j, x);
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      abort();
    }
  
  int nb_iter = atoi(argv[1]);
  int n = atoi(argv[2]);
  DISP(n);
  
  // testing DotProd 
  // multi-threaded with n large enough
  // more efficient with mkl
  if (false)
  {
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    double sum(0); 
    for (int i = 0; i < nb_iter; i++)
      sum += DotProd(x, y);
    
    DISP(sum);
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
  }
  
  // ApplyRot
  // multithreaded with mkl (for n large enough)
  // more efficient with mkl
  if (false)
  {
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    double c, s, theta;
    GetRandNumber(theta);
    c = cos(theta); s = sin(theta);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      ApplyRot(x, y, c, s);
    
    DISP(Norm2(x)); DISP(Norm2(y));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing Mlt with a vector
  // mono-threaded with mkl
  // similar efficiency
  if (false)
  {
    Vector<double> x;
    
    GenerateRandomVector(x, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        Mlt(1.0/Norm2(x), x);
      }
    
    DISP(Norm2(x));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    

  }
  
  // testing Copy
  // multi-threaded with n large enough
  // more efficient with mkl
  if (false)
  {
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      Copy(x, y);
    
    DISP(Norm2(x)); DISP(Norm2(y));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing Add
  // multi-threaded for n large enough
  // more efficient with mkl
  if (false)
  {
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    double alpha;
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        Add(alpha, x, y);
      }
    
    DISP(Norm2(x)); DISP(Norm2(y));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing product with triangular matrix (Blas level 2)
  // multi-threaded with mkl
  // more efficient with mkl
  if (false)
  {
    Matrix<double, General, RowUpTriang> T;
    T.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
        GetRandNumber(T.Get(i, j));
    
    Vector<double> x;    
    GenerateRandomVector(x, n);
    
    double alpha;
    GetRandNumber(alpha);
    
    glob_chrono.Start(VirtualTimer::ALL);
 
    for (int i = 0; i < nb_iter; i++)
      {
        if (i%50 == 0)
          GenerateRandomVector(x, n);
        
        Mlt(T, x);
      }
    
    DISP(Norm2(x));
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing matrix vector product (Blas level 2)
  // multi-threaded with mkl
  // more efficient with mkl
  if (false)
  {
    Matrix<double> A;
    GenerateRandomMatrix(A, n, n);
    
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    double alpha, beta;
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        GetRandNumber(beta);
        MltAdd(alpha, A, x, beta, y);
      }
    
    DISP(Norm2(y));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
        
  }

  // testing Rank1Update
  // mono-threaded with mkl
  // less efficient with mkl
  if (false)
  {
    Matrix<double> A;
    GenerateRandomMatrix(A, n, n);
    
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    double alpha;
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        Rank1Update(alpha, x, y, A);
      }
    
    DISP(MaxAbs(A));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
    
  }

  // testing Rank2Update
  // mono-threaded with mkl
  // less efficient with mkl
  if (false)
  {
    Matrix<double, Symmetric, RowSymPacked> A;
    GenerateRandomMatrix(A, n, n);
    
    Vector<double> x, y;
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    double alpha;
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        Rank2Update(alpha, x, y, A);
      }
    
    DISP(MaxAbs(A));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
    
  }

  // testing solve with triangular matrix (Blas level 2)
  // mono-threaded with mkl
  // more efficient with mkl
  if (false)
  {
    Matrix<double, General, RowUpTriang> T;
    T.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
          GetRandNumber(T.Get(i, j));
        
        T.Get(i, i) += 1.0;
      }
    
    Vector<double> x;    
    GenerateRandomVector(x, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
 
    for (int i = 0; i < nb_iter; i++)
      {
        if (i%50 == 0)
          GenerateRandomVector(x, n);
        
        Solve(T, x);
      }
    
    DISP(Norm2(x));
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing Mlt between matrices (Blas Level 3)
  // multi-threaded with mkl
  // much more efficient with mkl
  if (false)
  {
    Matrix<double> A, B, C;
    GenerateRandomMatrix(A, n, n);
    GenerateRandomMatrix(B, n, n);
    GenerateRandomMatrix(C, n, n);
    double alpha, beta;

    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        GetRandNumber(beta);
        MltAdd(alpha, A, B, beta, C);
      }
    
    DISP(MaxAbs(C));

    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
  }

  // testing Mlt with triangular matrices
  // multi-threaded with mkl
  // much more efficient with mkl
  if (false)
  {
    Matrix<double, General, RowUpTriang> T;
    T.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
        GetRandNumber(T.Get(i, j));
    
    Matrix<double> A;    
    GenerateRandomMatrix(A, n, n);
    
    double alpha;
    GetRandNumber(alpha);
    
    glob_chrono.Start(VirtualTimer::ALL);
 
    for (int i = 0; i < nb_iter; i++)
      {
        GenerateRandomMatrix(A, n, n);
        Mlt(SeldonLeft, alpha, T, A);
      }
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }


  // testing Solve with triangular matrices
  // multi-threaded with mkl  
  // much more efficient with mkl
  if (false)
  {
    Matrix<double, General, RowUpTriang> T;
    T.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
          GetRandNumber(T.Get(i, j));
        
        T.Get(i, i) += 1.0;
      }
    
    Matrix<double> A;    
    GenerateRandomMatrix(A, n, n);
    
    double alpha;
    GetRandNumber(alpha);
    
    glob_chrono.Start(VirtualTimer::ALL);
 
    for (int i = 0; i < nb_iter; i++)
      {
        if (i%10 == 0)
          GenerateRandomMatrix(A, n, n);
        
        Solve(SeldonLeft, alpha, T, A);
      }
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;    
  }
  
  // testing Blas Sparse Level 1
  // mono-threaded with mkl
  // same efficiency MKL or Seldon
  if (false)  
  {
    Vector<double, VectSparse> x;
    Vector<double> y;
    
    GenerateRandomVector(y, n);
    GenerateRandomVector(x, n, n/10);
    
    double alpha;
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        Add(double(2.0*alpha-1.0), x, y);
      }
    
    DISP(Norm2(y));
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
  }
  
  // testing Blas Sparse Level 2
  // multi-threaded with mkl
  // RowSparse is more efficient than MKL which is more efficient than ArrayRowSparse
  if (false)
  {
    Matrix<double, General, RowSparse> A;
    Vector<double> x(n), y(n);
    
    GenerateRandomMatrix(A, n, n, 20*n);
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      Mlt(A, x, y);
    
    DISP(Norm2(y));
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
    
  }
  
  // MKL and RowSparse have same efficiency
  if (false)
  {
    Matrix<double, General, RowSparse> A;
    Vector<double> x(n), y(n);
    
    GenerateRandomMatrix(A, n, n, 20*n);
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    double alpha, beta;
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        GetRandNumber(beta);
        MltAdd(alpha, A, x, beta, y);
      }
    
    DISP(Norm2(y));
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
    
  }
 
  // MKL and RowSparse have same efficiency
  if (false)
  {
    Matrix<double, General, RowSparse> A;
    Vector<double> x(n), y(n);
    
    GenerateRandomMatrix(A, n, n, 20*n);
    
    GenerateRandomVector(x, n);
    GenerateRandomVector(y, n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      Solve(SeldonLower, SeldonNoTrans, SeldonUnit, A, x, y);
    
    DISP(Norm2(y));
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
    
  }
  
  
  // testing Sparse Blas Level 3
  // multi-threaded with MKL
  // faster with RowSparse than with MKL
  if (false)
  {
    Matrix<double, General, RowSparse> A;
    Matrix<double> B, C;
    
    GenerateRandomMatrix(A, n, n, 20*n);
    GenerateRandomMatrix(B, n, 100);
    GenerateRandomMatrix(C, n, 100);
    
    double alpha, beta;
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        GetRandNumber(beta);
        MltAdd(alpha, SeldonNoTrans, A, SeldonNoTrans, B, beta, C);
      }
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    DISP(MaxAbs(C));
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;            
  }

#ifdef SELDON_WITH_MKL
  // monothreaded
  if (false)
    {
    Matrix<double, General, RowSparse> A;
    Matrix<double> B, C;
    
    GenerateRandomMatrix(A, n, n, 20*n);
    GenerateRandomMatrix(B, n, 100);
    GenerateRandomMatrix(C, n, 100);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      Solve(SeldonLower, SeldonNoTrans, SeldonUnit, A, B, C);
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    DISP(MaxAbs(C));
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;            
    }
#endif
  
  // testing Add between sparse matrices
  // mono-threaded with MKL
  // faster with RowSparse, then MKL, then ArrayRowSparse is the slowest
  if (false)
  {
    Matrix<double, General, ArrayRowSparse> A, B, C;
    
    GenerateRandomMatrix(A, n, n, 20*n);
    GenerateRandomMatrix(B, n, n, 20*n);
    
    double alpha;
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      {
        GetRandNumber(alpha);
        Add(alpha, A, B);
      }
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    DISP(MaxAbs(B));
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
  }

  // testing Mlt between sparse matrices
  // multi-threaded with MKL
  // faster with ArrayRowSparse than by using MKL
  //if (false)
  {
    Matrix<double, General, RowSparse> A, B, C;
    //Matrix<double, General, ArrayRowSparse> A, B, C;
    
    GenerateRandomMatrix(A, n, n, 20*n);
    GenerateRandomMatrix(B, n, n, 20*n);
    
    glob_chrono.Start(VirtualTimer::ALL);
    
    for (int i = 0; i < nb_iter; i++)
      Mlt(A, B, C);
    
    glob_chrono.Stop(VirtualTimer::ALL);
    
    DISP(MaxAbs(C)); DISP(C.GetDataSize());
    
    cout << "Temps = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;        
  }
  
  return FinalizeMontjoie();
}
