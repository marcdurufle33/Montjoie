#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

template<class T, int d>
void CheckBandMatrix(Matrix<T, General, BandedCol>& A,
                     TinyBandMatrix<T, d>& B, int n, int nb_iter)
{
  A.Reallocate(n, n, d, d);
  
  glob_chrono.Start(VirtualTimer::MASS);
  for (int k = 0; k < nb_iter; k++)
    {
      A.FillRand();
      Mlt(1e-9, A);
      
      GetLU(A);
    }
  
  glob_chrono.Stop(VirtualTimer::MASS);
  
  cout << "Time to factorize general banded matrices : "
       << glob_chrono.GetSeconds(VirtualTimer::MASS) << endl;
  
  // Lapack interface
  IVect pivot(n);
  Matrix<T, General, BandedCol> C;
  C.Reallocate(n, n, d, d);
  
  glob_chrono.Start(VirtualTimer::STIFFNESS);
  for (int k = 0; k < nb_iter; k++)
    {
      C.FillRand();
      Mlt(1e-9, C);
      
      GetLU(C, pivot);
    }
  
  glob_chrono.Stop(VirtualTimer::STIFFNESS);

  cout << "Time to factorize general banded matrices with Lapack : "
       << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
  
  // optimized class
  B.Reallocate(n, n);
  
  glob_chrono.Start(VirtualTimer::FACTO);
  for (int k = 0; k < nb_iter; k++)
    {
      B.FillRand();
      Mlt(1e-9, B);
      
      GetLU(B);
    }
  
  glob_chrono.Stop(VirtualTimer::FACTO);
  
  cout << "Time to factorize tiny banded matrices : "
       << glob_chrono.GetSeconds(VirtualTimer::FACTO) << endl;  

  glob_chrono.Start(VirtualTimer::MASS);  
  Vector<T> x(n);
  x.FillRand(); Mlt(1e-9, x);
  glob_chrono.Reset(VirtualTimer::MASS);
  glob_chrono.Start(VirtualTimer::MASS);
  for (int k = 0; k < nb_iter; k++)
    {
      SolveLU(A, x);
      Mlt(1.0/Norm2(x), x);
    }

  glob_chrono.Stop(VirtualTimer::MASS);
  
  cout << "Time to solve general banded matrices : "
       << glob_chrono.GetSeconds(VirtualTimer::MASS) << endl;
  
  glob_chrono.Reset(VirtualTimer::STIFFNESS);
  glob_chrono.Start(VirtualTimer::STIFFNESS);
  for (int k = 0; k < nb_iter; k++)
    {
      SolveLU(C, pivot, x);
      Mlt(1.0/Norm2(x), x);
    }

  glob_chrono.Stop(VirtualTimer::STIFFNESS);

  cout << "Time to solve general banded matrices with Lapack : "
       << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
  
  glob_chrono.Reset(VirtualTimer::FACTO);
  glob_chrono.Start(VirtualTimer::FACTO);
  for (int k = 0; k < nb_iter; k++)
    {
      SolveLU(B, x);
      Mlt(1.0/Norm2(x), x);
    }

  glob_chrono.Stop(VirtualTimer::FACTO);
  
  cout << "Time to solve tiny banded matrices : "
       << glob_chrono.GetSeconds(VirtualTimer::FACTO) << endl;  


}

int main(int argc, char** argv)
{
  cout.precision(15);
  bool overall_success = true;
  srand(0);
  
  int nb_iter = atoi(argv[1]);
  
  {
    // testing band-matrices
    Matrix<double, General, BandedCol> A;
    TinyBandMatrix<double, 20> B;
    int n = 50000;
    
    CheckBandMatrix(A, B, n, nb_iter);
  }
  
  if (overall_success)
    cout << "All tests passed successfully" << endl;
  
  return 0;
}
