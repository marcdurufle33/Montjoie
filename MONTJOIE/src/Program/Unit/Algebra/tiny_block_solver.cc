#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold;

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

template<class T>
void GenerateRandomMatrix(Matrix<T>& A, int m, int n)
{
  A.Reallocate(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      GetRandNumber(A(i, j));
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;

  if (Norm2(y) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(abs(x(i)-y(i))))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, int r>
void CheckSolver(TinyBlockSolver1D<T, r>& A)
{
  Matrix<T, General, BandedCol> B;

  int nb_elem = 10;
  int nb_dof = nb_elem*(r+1) + 1;
  B.Reallocate(nb_dof, nb_dof, r+1, r+1);
  
  Matrix<T> block;
  int offset = 0;
  for (int i = 0; i < nb_elem; i++)
    {
      if (i == 0)
	GenerateRandomMatrix(block, r+2, r+2);
      
      for (int i = 0; i <= r+1; i++)
	for (int j = 0; j <= r+1; j++)
	  B.AddInteraction(offset+i, offset+j, block(i, j));

      offset += r+1;
    }

  Vector<T> x, b(nb_dof);
  GenerateRandomVector(x, nb_dof);

  Mlt(B, x, b);
  
  A.Factorize(B, true);
  A.Solve(b);

  if (!EqualVector(x, b))
    {
      cout << "Factorize/Solve incorrect" << endl;
      abort();
    }

  // non homogeneous case
  B.Reallocate(nb_dof, nb_dof, r+1, r+1);
  B.Zero();
  offset = 0;
  for (int i = 0; i < nb_elem; i++)
    {
      GenerateRandomMatrix(block, r+2, r+2);
      
      for (int i = 0; i <= r+1; i++)
	for (int j = 0; j <= r+1; j++)
	  B.AddInteraction(offset+i, offset+j, block(i, j));

      offset += r+1;
    }

  Mlt(B, x, b);
  
  A.Factorize(B);
  A.Solve(b);
  
  if (!EqualVector(x, b))
    {
      cout << "Factorize/Solve incorrect" << endl;
      abort();
    }

}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e4*epsilon_machine;
  
  {
    // testing band-matrices
    TinyBlockSolver1D<Real_wp, 9> A;
    CheckSolver(A);
  }

  cout << "All tests passed successfully" << endl;

  FinalizeMontjoie();
}

