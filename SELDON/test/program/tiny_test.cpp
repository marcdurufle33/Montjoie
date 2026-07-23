#include "SeldonLib.hxx"
#include <chrono>

using namespace Seldon;

template<int m>
void RunTestGetInverse(TinyMatrix<double, General, m, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les matrices Seldon
  Matrix<double> A(m, m);
  auto start = std::chrono::steady_clock::now();
  double vloc = 0; srand(0);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      A.FillRand(); A *= 1e-9;
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(A);

      
      vloc += A(0, 0);
    }

  DISP(vloc);
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> sec = end-start;
  double time_std = sec.count();
  cout<<"Time with dense matrices "<< sec.count() << endl;

  vloc = 0; srand(0);
  TinyMatrix<double, General, m, m> a;
  start = std::chrono::steady_clock::now();  
  for (int n = 0; n < nb_outer_iter; n++)
    {
      a.FillRand(); a *= 1e-9;
      
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(a);

      vloc += a(0, 0);
    }

  DISP(vloc);
  
  end = std::chrono::steady_clock::now();
  sec = end-start;
  double time_tiny = sec.count();
  cout<<"Time with tiny matrices "<<sec.count()<<endl;
  cout << "Ratio = " << time_std / time_tiny << endl;
}

template<int m>
void RunTestGetInverse(TinyMatrix<double, Symmetric, m, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les matrices Seldon
  double vloc = 0; srand(0);
  Matrix<double, Symmetric, RowSymPacked> A(m, m);
  auto start = std::chrono::steady_clock::now();
  for (int n = 0; n < nb_outer_iter; n++)
    {
      A.FillRand(); A *= 1e-9;

      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(A);

      vloc += A(0, 0);
    }

  DISP(vloc);
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> sec = end-start;
  double time_std = sec.count();
  cout<<"Time with dense matrices "<< sec.count() << endl;

  vloc = 0; srand(0);
  TinyMatrix<double, Symmetric, m, m> a;
  start = std::chrono::steady_clock::now();  
  for (int n = 0; n < nb_outer_iter; n++)
    {
      a.FillRand(); a *= 1e-9;

      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(a);

      vloc += a(0, 0);
    }

  DISP(vloc);
  
  end = std::chrono::steady_clock::now();
  sec = end-start;
  double time_tiny = sec.count();
  cout<<"Time with tiny matrices "<<sec.count()<<endl;
  cout << "Ratio = " << time_std / time_tiny << endl;
}

// testing class TinyVector/TinyMatrix
int main(int argc, char** argv)
{
  cout.precision(15);

  
  int nb_inner = atoi(argv[1]);
  int nb_outer = atoi(argv[2]);
  {
    TinyMatrix<double, General, 3, 3> A;
  
    RunTestGetInverse(A, nb_inner, nb_outer);
  }
  
  {
    TinyMatrix<double, Symmetric, 3, 3> A;
    
    RunTestGetInverse(A, nb_inner, nb_outer);
  }
  
  /*
  TinyMatrix<double, Symmetric, 8, 8> A;
  A(0, 0) = 2.0; A(0, 1) = 4.0; A(0, 2) = 5.0;
  A(1, 1) = 1.0; A(1, 2) = -2.0; A(2, 2) = 0.0;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = i; j < A.GetM(); j++)
      A(i, j) = Real_wp(rand()) / RAND_MAX;
  
  Matrix<double, Symmetric, RowSymPacked> B;
  B.Reallocate(A.GetM(), A.GetM());
  for (int i = 0; i < A.GetM(); i++)
    for (int j = i; j < A.GetM(); j++)
      B(i, j) = A(i, j);
  
  DISP(A);
  GetInverse(A);  

  cout << endl << "Appel Fortran" << endl;
  
  Vector<int> pivot(B.GetM());
  char uplo('L'); int n = B.GetM(), info;
  dsptrf_(&uplo, &n, B.GetData(), pivot.GetData(), &info);
  DISP(B); DISP(pivot);

  Vector<double> work(B.GetM());
  dsptri_(&uplo, &n, B.GetData(), pivot.GetData(), work.GetData(), &info);
  DISP(B);

  double err(0);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = i; j < A.GetM(); j++)
      err = max(err, abs(A(i, j) - B(i, j)));

  DISP(err);
  */
  return 0;
}
