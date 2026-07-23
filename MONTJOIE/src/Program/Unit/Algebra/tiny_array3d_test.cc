#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold;

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  int type = rand()%3;
  //int type = 2;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T, int m, int n, int p>
void GenerateRandomArray(TinyArray3D<T, m, n, p>& A)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	GetRand(A(i, j, k));
}

template<class T, int m, int n, int p>
void CheckTinyArray3D(TinyArray3D<T, m, n, p>& A)
{
  TinyVector<T, m>::threshold = threshold;
  TinyVector<T, n>::threshold = threshold;
  TinyVector<T, p>::threshold = threshold;
  
  A.FillRand(); A *= 1e-9;
  
  if (A.GetLength1() != m)
    {
      cout << "GetLength1 incorrect" << endl;
      abort();
    }

  if (A.GetLength2() != n)
    {
      cout << "GetLength2 incorrect" << endl;
      abort();
    }

  if (A.GetLength3() != p)
    {
      cout << "GetLength3 incorrect" << endl;
      abort();
    }

  if (A.GetSize() != m*n*p)
    {
      cout << "GetSize incorrect" << endl;
      abort();
    }
  
  // checking Zero, Fill() and operator *=
  A.Fill(2.3);  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        if ((abs(A(i, j, k) - 2.3) > threshold) || isnan(abs(A(i, j, k) - 2.3)))
          {
            cout << "Fill incorrect" << endl;
            abort();
          }

  A.Zero();
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        if ((abs(A(i, j, k)) > threshold) || isnan(A(i, j, k)))
          {
            cout << "Zero incorrect" << endl;
            abort();
          }
  
  A.FillRand(); A *= 1e-9;
  TinyArray3D<T, m, n, p> B; B = A;
  B *= 1.7;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        if ((abs(B(i, j, k) - 1.7*A(i, j, k)) > threshold) || isnan(abs(B(i, j, k) - 1.7*A(i, j, k))))
          {
            cout << "Operator *= incorrect" << endl;
            abort();
          }
  
  // checking operator *
  TinyVector<T, n> x;
  TinyMatrix<T, General, m, p> y, z;
  
  x.FillRand(); x *= 1e-9;
  y = dot(A, x);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < p; j++)
      {
        T vloc; SetComplexZero(vloc);
        for (int k = 0; k < n; k++)
          vloc += A(i, k, j)*x(k);
        
        if ((abs(vloc - y(i, j)) > threshold) || isnan(abs(vloc - y(i, j))))
          {
            cout << "dot incorrect" << endl;
            abort();
          }
      }
      
  // checking ExtractMatrix
  ExtractMatrix(A, 1, y);
  TinyVector<int, p> theta;
  for (int i = 0; i < p; i++)
    theta(i) = (i+2)%p;
  
  ExtractMatrix(A, theta, z);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < p; j++)
      if ( (abs(y(i, j) - A(i, 1, j)) >threshold) || 
           (abs(z(i, j) - A(i, theta(j), j)) >threshold)
           || isnan(abs(y(i, j) - A(i, 1, j))) || isnan(abs(z(i, j) - A(i, theta(j), j))) )
        {
          cout << "ExtractMatrix incorrect" << endl;
          abort();
        }

  // checking operator = expression
  typedef typename ClassComplexType<T>::Tcplx Tcplx;
  TinyArray3D<T, m, n, p> C, D;
  TinyArray3D<Tcplx, m, n, p> Ac, Bc, Cc, Dc;
  GenerateRandomArray(A); GenerateRandomArray(Ac);
  GenerateRandomArray(B); GenerateRandomArray(Bc);
  GenerateRandomArray(C); GenerateRandomArray(Cc);

  T alpha, beta;
  GetRand(alpha); GetRand(beta);
  
  Tcplx alpha_c;
  GetRand(alpha_c);
  
  // testing operator +, - with expressions
  D = alpha*(A+B) - beta*C;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha*(A(i, j, k) + B(i, j, k)) - beta*C(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator +,-,= incorrect" << endl;
	      abort();
	    }
	}

  // testing operator +=
  D.Zero();
  D = -beta*C;
  D += alpha*(A+B);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha*(A(i, j, k) + B(i, j, k)) - beta*C(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator += incorrect" << endl;
	      abort();
	    }
	}

  // testing operator -=
  D.Zero();
  D = alpha*(A+B);
  D -= beta*C; 
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha*(A(i, j, k) + B(i, j, k)) - beta*C(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator -= incorrect" << endl;
	      abort();
	    }
	}

  // testing operator +
  Dc = Ac + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) + B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) + Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) - B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) - Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  
  // testing operator *, / with expressions
  D = alpha*(A*B) / (beta*C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha*(A(i, j, k) * B(i, j, k)) / (beta*C(i, j, k));
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator *, / incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) * B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) * Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) / B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) / Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  // testing operator + with scalar
  D = alpha + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha + B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  D = B + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha + B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha_c + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha_c + B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A + alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) + alpha_c;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) + alpha;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha + Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator + incorrect" << endl;
	      abort();
	    }
	}
  
  // testing operator - with scalar
  D = alpha - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha - B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  D = B - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = -alpha + B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha_c - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha_c - B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A - alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) - alpha_c;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) - alpha;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha - Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
  
  // testing operator * with scalar
  D = alpha * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha * B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  D = B * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = B(i, j, k) * alpha;
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha_c * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha_c * B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A * alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) * alpha_c;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) * alpha;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha * Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator * incorrect" << endl;
	      abort();
	    }
	}
  
  // testing operator / with scalar
  D = alpha / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = alpha / B(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  D = B / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = B(i, j, k) / alpha;
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha_c / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha_c / B(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = A / alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = A(i, j, k) / alpha_c;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = Ac / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = Ac(i, j, k) / alpha;
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  Dc = alpha / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  Tcplx vexact = alpha / Bc(i, j, k);
	  if ( (abs(vexact-Dc(i, j, k)) > threshold) || isnan(abs(vexact-Dc(i, j, k))))
	    {
	      cout << "Operator / incorrect" << endl;
	      abort();
	    }
	}
  
  // testing opposite
  D = -A;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
	{
	  T vexact = -A(i, j, k);
	  if ( (abs(vexact-D(i, j, k)) > threshold) || isnan(abs(vexact-D(i, j, k))))
	    {
	      cout << "Operator - incorrect" << endl;
	      abort();
	    }
	}
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  bool overall_success = true;  
  threshold = 1e4*epsilon_machine;  
  cout.precision(15);
  
  {
    TinyArray3D<Real_wp, 5, 4, 3> A;
    CheckTinyArray3D(A);  
  }

  {
    TinyArray3D<Complex_wp, 3, 5, 4> A;
    CheckTinyArray3D(A);  
  }

  if (overall_success)
    cout << "All tests passed successfully " << endl;
  
  return FinalizeMontjoie();
}
