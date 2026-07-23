#define MONTJOIE_WITH_TWO_DIM

#include "Harmonic/MontjoieHarmonic.hxx"

using namespace Montjoie;

template<int m>
void RunTestGetInverse(TinyMatrix<double, General, m, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les matrices Seldon
  int jeton = glob_chrono.GetNumber();
  Matrix<double> A(m, m);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      A.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(A);
    }

  glob_chrono.Stop(jeton);
  cout<<"Big vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
    
  TinyMatrix<double, General, m, m> a;
  glob_chrono.Reset(jeton);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      a.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(a);
    }

  glob_chrono.Stop(jeton);
  cout<<"Small vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
}

template<int m>
void RunTestGetInverse(TinyMatrix<double, Symmetric, m, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les matrices Seldon
  int jeton = glob_chrono.GetNumber();
  Matrix<double, Symmetric, RowSymPacked> A(m, m);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      A.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(A);
    }

  glob_chrono.Stop(jeton);
  cout<<"Big vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
    
  TinyMatrix<double, General, m, m> a;
  glob_chrono.Reset(jeton);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      a.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	GetInverse(a);
    }

  glob_chrono.Stop(jeton);
  cout<<"Small vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
}

template<class T>
void GetRand(TinyVector<TinyVector<T, 2>, 2>& du)
{
  du(0)(0) = 1e-9*rand();
  du(1)(0) = 1e-9*rand();
  du(0)(1) = 1e-9*rand();
  du(1)(1) = 1e-9*rand();
}

template<class T>
void GetRand(TinyVector<TinyVector<T, 3>, 3>& du)
{
  du(0)(0) = 1e-9*rand();
  du(1)(0) = 1e-9*rand();
  du(0)(1) = 1e-9*rand();
  du(1)(1) = 1e-9*rand();
  du(0)(2) = 1e-9*rand();
  du(1)(2) = 1e-9*rand();
  du(2)(0) = 1e-9*rand();
  du(2)(1) = 1e-9*rand();
  du(2)(2) = 1e-9*rand();
}

template<class T, int m>
void RunTestSymTensor(TinySymmetricTensor<T, m>& C, int nb_inner_iter, int nb_outer_iter)
{
  for (int i = 0; i < m*(m+1)/2; i++)
    for (int j = 0; j < m*(m+1)/2; j++)
      C(i, j) = 1e-10*rand();

  for (int i = 0; i < m; i++)
    for (int j = i; j < m; j++)
      for (int k = 0; k < m; k++)
        for (int l = k+1; l < m; l++)
          if ((i != k) || (j != l))
            C(i, j, k, l) = 0.0;
  
  VarPhysicalProblem var;
  ElasticPhysicalIndice<Dimension2, m, T> index;
  index.SetConstant(C);
  DISP(index.GetAnisotropyType()); DISP(nb_inner_iter);
  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);
  TinyVector<TinyVector<T, m>, m> du, dv;
  T sum = 0; DISP(m);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      GetRand(du);
      for (int p = 0; p < nb_inner_iter; p++)
        {
          index.MltMatrix(var, 0, 0, du, dv);
          //C.MltOrthotrope(du, dv);
          //C.Mlt(du, dv);
          du = dv;
        }
      
      sum += DotProd(dv(0), dv(1));
    }
  DISP(sum);
  
  glob_chrono.Stop(VirtualTimer::ALL);
  DISP(glob_chrono.GetSeconds(VirtualTimer::ALL));
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  cout.precision(7);
  
  if (argc < 3)
    {
      cout<<"Enter more arguments"<<endl;
      return -1;
    }

  int nb_inner_iter = atoi(argv[1]);
  int nb_outer_iter = atoi(argv[2]);

  TinyMatrix<double, Symmetric, 10, 10> A;
  //  RunTestGetInverse(A, nb_inner_iter, nb_outer_iter);
  
  TinySymmetricTensor<double, 2> C;
  RunTestSymTensor(C, nb_inner_iter, nb_outer_iter);

  TinySymmetricTensor<double, 3> C3;
  RunTestSymTensor(C3, nb_inner_iter, nb_outer_iter);
  
  return FinalizeMontjoie();
}
