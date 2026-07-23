#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold;

template<class T>
void GetRandNumber(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T> & x)
{
  x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
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

template<class T, class Prop, class T2>
bool EqualVector(const Matrix<T, Prop, DiagonalRow>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;

  if (Norm2(y) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i, i) - y(i)) > eps) || isnan(abs(x(i, i)-y(i))))
      {
        DISP(x(i, i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class Prop>
void CheckDiagonalMatrix(Matrix<T, Prop, DiagonalRow>&)
{
  int m = 20;
  
  Matrix<T, Prop, DiagonalRow> A;
  Matrix<T, Prop, DiagonalRow> B(m, m);

  if ((A.GetM() != 0) || (A.GetN() != 0))
    {
      cout << "Constructor incorrect" << endl;
      abort();
    }

  if ((B.GetM() != m) || (B.GetN() != m))
    {
      cout << "Constructor incorrect" << endl;
      abort();
    }

  if (B.GetMemorySize() < m*sizeof(T))
    {
      cout << "GetMemorySize incorrect" << endl;
      abort();
    }

  A.Reallocate(m, m);

  if ((A.GetM() != m) || (A.GetN() != m))
    {
      cout << "Reallocate incorrect" << endl;
      abort();
    }

  Vector<T> x(m);
  GenerateRandomVector(x, m);

  for (int i = 0; i < m; i++)
    A.Get(i, i) = x(i);

  if (!EqualVector(A, x))
    {
      cout << "Get incorrect" << endl;
      abort();
    }

  Vector<T> y(m), z(m);
  GenerateRandomVector(y, m);

  for (int i = 0; i < m; i++)
    A.AddInteraction(i, i, y(i));

  z = x + y;
  if (!EqualVector(A, z))
    {
      cout << "AddInteraction incorrect" << endl;
      abort();
    }

  for (int i = 0; i < m; i++)
    A.Get(i, i) = x(i);

  GenerateRandomVector(y, m);

  for (int i = 0; i < m; i++)
    {
      IVect col; Vector<T> val;
      int p = rand()%10; T vloc;
      col.PushBack(i); val.PushBack(y(i));
      for (int k = 0; k < p; k++)
	{
	  int j = rand()%m;
	  if (j != i)
	    {
	      GetRandNumber(vloc);
	      col.PushBack(j);
	      val.PushBack(vloc);
	    }
	}
      
      A.AddInteractionRow(i, col.GetM(), col, val);
    }

  z = x + y;
  if (!EqualVector(A, z))
    {
      cout << "AddInteractionRow incorrect" << endl;
      abort();
    }

  for (int i = 0; i < m; i++)
    if (A.GetRowSize(i) != 1)
      {
	cout << "GetRowSize incorrect" <<endl;
	abort();
      }

  B = A;
  if (!EqualVector(B, z))
    {
      cout << "AddInteractionRow incorrect" << endl;
      abort();
    }

  B.Resize(m+10, m+10);
  z.Resize(m+10);
  for (int i = m; i < B.GetM(); i++)
    {
      GetRandNumber(z(i));
      B.Get(i, i) = z(i);
    }
  
  if (!EqualVector(B, z))
    {
      cout << "Resize incorrect" << endl;
      abort();
    }

  B.Clear();
  if ((B.GetM() != 0) || (B.GetN() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  z.Resize(m);
  B = A;
  int num_row = rand()%m;
  B.ClearRow(num_row); SetComplexZero(z(num_row));

  if (!EqualVector(B, z))
    {
      cout << "ClearRow incorrect" << endl;
      abort();
    }

  A.SetData(x.GetM(), x.GetData());
  if (!EqualVector(A, x))
    {
      cout << "SetData incorrect" << endl;
      abort();
    }

  z = x;
  Mlt(A, y, z);
  Vector<T> res; res = x*y;
  if (!EqualVector(res, z))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  A.Nullify();
  if ((A.GetM() != 0) || (A.GetN() != 0))
    {
      cout << "Nullify incorrect" << endl;
      abort();
    }

  B.WriteText("B.dat");
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  threshold = 1e4*epsilon_machine;
  
  {
    Matrix<Real_wp, Symmetric, DiagonalRow> A;
    CheckDiagonalMatrix(A);
  }

  {
    Matrix<Complex_wp, Symmetric, DiagonalRow> A;
    CheckDiagonalMatrix(A);
  }

  cout << "All tests successfully passed"  << endl;
  
  FinalizeMontjoie();
}
