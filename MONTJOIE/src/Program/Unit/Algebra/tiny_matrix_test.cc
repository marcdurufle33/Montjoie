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

template<class T, class Prop, int m, int n>
void GenerateRandomMatrix(TinyMatrix<T, Prop, m, n>& A)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      GetRand(A(i, j));
}


template<class T, int m, int n>
void CheckTinyGeneralMatrix(TinyMatrix<T, General, m, n>& A)
{
  TinyVector<T, m>::threshold = threshold;
  TinyVector<T, n>::threshold = threshold;
  
  // checking constructors
  TinyVector<TinyVector<T, m>, n> cols;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      SetComplexReal(Real_wp(i) + 0.5*Real_wp(j*j), cols(j)(i));
  
  A.FillRand(); Mlt(1e-8, A);
  TinyMatrix<T, General, m, n> B(A), C, D(cols);
  
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      {
        if ((abs(B(i, j) - A(i, j)) > threshold) || isnan(abs(B(i, j) - A(i, j))))
          {
            cout << "Constructor incorrect " << endl;
            abort();
          }

        if ((abs(C(i, j)) > threshold) || isnan(C(i, j)))
          {
            cout << "Constructor incorrect " << endl;
            abort();
          }

        if ((abs(D(i, j) - cols(j)(i)) > threshold) || isnan(abs(D(i, j) - cols(j)(i))))
          {
            cout << "Constructor incorrect " << endl;
            abort();
          }
      }
  
  // checking GetM, GetN, GetSize
  if (A.GetM() != m)
    {
      cout << "GetM incorrect" << endl;
      abort();
    }

  if (A.GetN() != n)
    {
      cout << "GetN incorrect" << endl;
      abort();
    }

  if (A.GetSize() != m*n)
    {
      cout << "GetSize incorrect" << endl;
      abort();
    }
  
  // checking operator *
  T alpha, beta;
  SetComplexReal(0.4, alpha); 
  SetComplexReal(1.5, beta);
  B = alpha*D;
  C = D*beta;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      if ((abs(B(i, j) - alpha*D(i, j)) > threshold) || (abs(C(i, j) - beta*D(i, j)) > threshold) ||
          isnan(abs(B(i, j) - alpha*D(i, j))) || isnan(abs(C(i, j) - beta*D(i, j))))
        {
          cout << "Operator * incorrect" << endl;
          abort();
        }
  
  TinyVector<T, n> x; TinyVector<T, m> y;
  x.FillRand(); Mlt(Real_wp(1e-10), x);
  y = dot(A, x);
  for (int i = 0; i < m; i++)
    {
      T vloc; SetComplexZero(vloc);
      for (int j = 0; j < n; j++)
        vloc += A(i, j)*x(j);
      
      if ((abs(vloc - y(i)) > threshold) || isnan(abs(vloc - y(i))))
        {
          cout << "Operator * incorrect " << endl;
          abort();
        }
    }
  
  TinyMatrix<T, General, m, 6> Ae;
  TinyMatrix<T, General, n, 6> E;
  E.FillRand(); Mlt(Real_wp(1e-10), E);
  Ae = dot(A, E);
  for (int k = 0; k < 6; k++)
    for (int i = 0; i < m; i++)
      {
        T vloc; SetComplexZero(vloc);
        for (int j = 0; j < n; j++)
          vloc += A(i, j)*E(j, k);
        
        if ((abs(vloc - Ae(i, k)) > threshold) || isnan(abs(vloc - Ae(i, k))))
          {
            cout << "Operator * incorrect " << endl;
            abort();
          }
      }
  
  // checking operators +, -
  C = A+B;
  D = A-B;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      if (  (abs(C(i, j) - A(i, j) - B(i, j)) > threshold) ||
            (abs(D(i, j) - A(i, j) + B(i, j)) > threshold) ||
            isnan(abs(C(i, j) - A(i, j) - B(i, j))) ||
            isnan(abs(D(i, j) - A(i, j) + B(i, j))) )
        {
          cout << "Operator +, - incorrect " << endl;
          abort();
        }
  
  // checking operator ==
  if ((C == A+B) && (C != A-B))
    {
    }
  else
    {
      cout << "Operator ==, != incorrect" << endl;
      abort();
    }
  
  if ((D==A+B) || (D != A-B))
    {
      cout << "Operator ==, != incorrect" << endl;
      abort();
    }

  // checking operator *=, +=, and -=
  D = A; D += B;
  if (D != C)
    {
      cout << "Operator += incorrect " << endl;
      abort();
    }

  C = A-B;
  D = A; D -= B;
  if (D != C)
    {
      cout << "Operator -= incorrect " << endl;
      abort();
    }

  C = T(1.5)*A;
  D = A; D *= 1.5;
  if (D != C)
    {
      cout << "Operator *= incorrect " << endl;
      abort();
    }
  
  // checking operator =
  C = A;
  if (A != C)
    {
      cout << "operator = incorrect " << endl;
      abort();
    }
  
  C = 0.5;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      if ((abs(C(i, j) - 0.5) > threshold) || isnan(C(i, j)))
        {
          cout << "operator = incorrect " << endl;
          abort();
        }
  
  // checking matrix-vector operations
  y.Fill(0);
  Mlt(A, x, y);
  TinyVector<T, n> x2;
  TinyVector<T, m> y2;
  y2 = dot(A, x);
  if (y != y2)
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  y.FillRand(); Mlt(Real_wp(1e-10), y);
  y2 = y; y2 += Real_wp(2.1)*dot(A, x);
  MltAdd(Real_wp(2.1), A, x, y);
  if (y != y2)
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }
  
  MltTrans(A, y, x);
  for (int i = 0; i < n; i++)
    {
      T vloc; SetComplexZero(vloc);
      for (int j = 0; j < m; j++)
        vloc += A(j, i)*y(j);
      
      if ((abs(vloc-x(i)) > threshold) || isnan(abs(vloc-x(i))))
        {
          cout << "MltTrans incorrect" << endl;
          abort();
        }
    }
  
  C = A;
  Rank1Update(Real_wp(1.3), y, x, C);
  D = A;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      D(i, j) += Real_wp(1.3)*y(i)*x(j);
  
  if (C != D)
    {
      cout << "Rank1Update incorrect" << endl;
      abort();
    }

  Rank1Matrix(y, x, C);
  D.Fill(0);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      D(i, j) = y(i)*x(j);
  
  if (C != D)
    {
      cout << "Rank1Matrix incorrect" << endl;
      abort();
    }
  
  x2.FillRand(); Mlt(Real_wp(1e-10), x2);
  TinyMatrix<T, General, m, n> Acopy(A);
  for (int i = 0; i < m; i++)
    {
      GetRow(A, i, x2);
      for (int j = 0; j < n; j++)
        if ((abs(x2(j) - A(i, j)) > threshold) || isnan(abs(x2(j) - A(i, j))))
          {
            cout << "GetRow incorrect" << endl;
            abort();
          }
      
      SetRow(x, i, A);
      for (int j = 0; j < n; j++)
        if ((abs(x(j) - A(i, j)) > threshold) || isnan(abs(x(j) - A(i, j))))
          {
            cout << "SetRow incorrect" << endl;
            abort();
          }      
    }
  
  for (int j = 0; j < n; j++)
    {
      GetCol(A, j, y2);
      for (int i = 0; i < m; i++)
        if ((abs(y2(i) - A(i, j)) > threshold) || isnan(abs(y2(i) - A(i, j))))
          {
            cout << "GetCol incorrect" << endl;
            abort();
          }
      
      SetCol(y, j, A);
      for (int i = 0; i < m; i++)
        if ((abs(y(i) - A(i, j)) > threshold) || isnan(abs(y(i) - A(i, j))))
          {
            cout << "SetCol incorrect" << endl;
            abort();
          }      
    }
  
  A = Acopy;
  
  // matrix-matrix operations
  D = C; 
  Mlt(Real_wp(2.5), C);
  D *= Real_wp(2.5);
  if (C != D)
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Mlt(A, E, Ae);
  TinyMatrix<T, General, m, 6> Ae2;
  Ae2 = dot(A, E);
  if (Ae2 != Ae)
    {
      cout << "Mlt incorrect " << endl;
      abort();
    }

  D = A+B;
  Add(A, B, C);
  if (C != D)
    {
      cout << "Add incorrect" << endl;
      abort();
    }
  
  // check Zero, SetIdentity, Fill
  C.Zero();
  D.FillRand(); Mlt(Real_wp(1e-10), D);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j)) > threshold) || isnan(abs(C(i, j))))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  if ((!C.IsZero()) || (D.IsZero()))
    {
      cout << "IsZero incorrect" << endl;
      abort();
    }
  
  
  TinyMatrix<T, General, m, m> As, Bs, Cs;
  As.SetIdentity();
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
          SetComplexOne(vexact);
        
        if ((abs(vexact - As(i, j)) > threshold) || isnan(abs(vexact - As(i, j))))
          {
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  C.Fill(Real_wp(-1.6));
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(C(i, j) + Real_wp(1.6)) > threshold) || isnan(abs(C(i, j) + Real_wp(1.6))))
        {
          cout << "Fill incorrect" << endl;
          abort();
        }
  
  C.Fill();
  int nb = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {        
        if ((abs(C(i, j) - nb) > threshold) || isnan(abs(C(i, j) - nb)))
          {
            cout << "Fill incorrect" << endl;
            abort();
          }
        
        nb++;
      }
  
  // checking Transpose and MltTrans
  TinyMatrix<T, General, 6, n> Et;
  Transpose(E, Et);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 6; j++)
      if ((abs(E(i, j) - Et(j, i)) > threshold) || isnan(abs(E(i, j) - Et(j, i))))
        {
          cout << "Transpose incorrect" << endl;
          abort();
        }
  
  Ae2 = dot(A, E);
  MltTrans(A, Et, Ae);
  if (Ae != Ae2)
    {
      cout << "MltTrans incorrect" << endl;
      abort();  
    }

  // checking MaxAbs
  typename ClassComplexType<T>::Treal Emax, Etest(0);
  Emax = MaxAbs(A);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Etest = max(Etest, abs(A(i, j)));
  
  if ((abs(Etest-Emax) > threshold) || isnan(abs(Etest-Emax)))
    {
      cout << "MaxAbs incorrect" << endl;
      abort();
    }
  
  // checking GetInverse
  As.FillRand(); Mlt(Real_wp(1e-10), As);
  GetInverse(As, Bs);
  Mlt(As, Bs, Cs);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
          SetComplexOne(vexact);
        
        if ((abs(vexact - Cs(i, j)) > threshold) || isnan(abs(vexact - Cs(i, j))))
          {
            cout << "GetInverse incorrect" << endl;
            abort();
          }
      }

  // checking SetDiagonal
  GetRand(alpha);   GetRand(beta);
  As.SetDiagonal(alpha);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
	  vexact = alpha;
        
        if ((abs(vexact - As(i, j)) > threshold) || isnan(abs(vexact - As(i, j))))
          {
            cout << "SetDiagonal incorrect" << endl;
            abort();
          }
      }

  typedef typename ClassComplexType<T>::Tcplx Tcplx;
  TinyMatrix<Tcplx, General, m, n> Ac, Bc, Cc, Dc;
  GenerateRandomMatrix(A); GenerateRandomMatrix(Ac);
  GenerateRandomMatrix(B); GenerateRandomMatrix(Bc);
  GenerateRandomMatrix(C); GenerateRandomMatrix(Cc);

  Tcplx alpha_c;
  GetRand(alpha_c);
  
  // testing operator +, - with expressions
  D = alpha*(A+B) - beta*C;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha*(A(i, j) + B(i, j)) - beta*C(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator +/- incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) + B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = A + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) + Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) - B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) - Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }


  // testing operator *, / with expressions
  D = alpha*(A*B) / (beta*C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha*(A(i, j) * B(i, j)) / (beta*C(i, j));
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator *, / incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) * B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = A * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) * Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) / B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) / Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  // testing operator + with scalar
  D = alpha + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  D = B + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha_c + B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = A + alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) + alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) + alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha + Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  // testing operator - with scalar
  D = alpha - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha - B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  D = B - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = -alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha_c - B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A - alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) - alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) - alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha - Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }


  // testing operator * with scalar
  D = alpha * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha * B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  D = B * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = B(i, j) * alpha;
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha_c * B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = A * alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) * alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) * alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha * Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  // testing operator / with scalar
  D = alpha / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = alpha / B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  D = B / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = B(i, j) / alpha;
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha_c / B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = A / alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = A(i, j) / alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = Ac(i, j) / alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	Tcplx vexact = alpha / Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  // testing opposite
  D = -A;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
	T vexact = -A(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  // testing transpose
  TinyMatrix<T, General, n, m> Dt;
  Dt = transpose(A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = A(j, i);
	if ( (abs(vexact-Dt(i, j)) > threshold) || isnan(abs(vexact-Dt(i, j))))
	  {
	    cout << "transpose incorrect" << endl;
	    abort();
	  }
      }

  // testing DotProdCol
  x2 = x;
  DotProdCol(A, B, x);
  for (int i = 0; i < n; i++)
    {
      T vexact = x2(i);
      for (int j = 0; j < m; j++)
	vexact += A(j, i)*B(j, i);

      if ((abs(vexact-x(i)) > threshold) || isnan(abs(vexact-x(i))))
	{
	  cout << "DotProdCol incorrect" << endl;
	  abort();
	}
    }

  // testing Norm2_column
  typename ClassComplexType<T>::Treal am, aref;
  int first_row = 1;
  for (int i = 0; i < n; i++)
    {
      am = Norm2_Column(A, first_row, i);
      SetComplexZero(aref);
      for (int j = first_row; j < m; j++)
	aref += absSquare(A(j, i));

      aref = sqrt(aref);
      if ( (abs(am-aref) > threshold) || isnan(abs(am-aref)))
	{
	  cout << "Norm2_column incorrect" << endl;
	  abort();
	}
    }
}


template<class T, int m>
void CheckTinySymmetricMatrix(TinyMatrix<T, Symmetric, m, m>& A)
{
  TinyVector<T, m>::threshold = threshold;
  
  // checking constructors
  A.FillRand(); Mlt(Real_wp(1e-8), A);
  TinyMatrix<T, Symmetric, m, m> B(A), C, D;
  
  for (int j = 0; j < m; j++)
    for (int i = 0; i < m; i++)
      {
        if ((abs(B(i, j) - A(i, j)) > threshold) || isnan(abs(B(i, j) - A(i, j))))
          {
            cout << "Constructor incorrect " << endl;
            abort();
          }

        if ((abs(C(i, j)) > threshold) || isnan(C(i, j)))
          {
            cout << "Constructor incorrect " << endl;
            abort();
          }
      }
  
  // checking GetM, GetN, GetSize
  if (A.GetM() != m)
    {
      cout << "GetM incorrect" << endl;
      abort();
    }

  if (A.GetN() != m)
    {
      cout << "GetN incorrect" << endl;
      abort();
    }

  if (A.GetSize() != m*(m+1)/2)
    {
      cout << "GetSize incorrect" << endl;
      abort();
    }
  
  // checking operator *
  T alpha, beta; SetComplexReal(0.4, alpha);
  SetComplexReal(1.5, beta);
  B = alpha*A;
  C = A*beta;
  for (int j = 0; j < m; j++)
    for (int i = 0; i < m; i++)
      if ((abs(B(i, j) - alpha*A(i, j)) > threshold)
          || (abs(C(i, j) - beta*A(i, j)) > threshold)
          || isnan(abs(B(i, j) - alpha*A(i, j))) || isnan(abs(C(i, j) - beta*A(i, j))))
        {
          cout << "Operator * incorrect" << endl;
          abort();
        }
  
  TinyVector<T, m> x; TinyVector<T, m> y;
  x.FillRand(); Mlt(Real_wp(1e-10), x);
  y = dot(A, x);
  for (int i = 0; i < m; i++)
    {
      T vloc; SetComplexZero(vloc);
      for (int j = 0; j < m; j++)
        vloc += A(i, j)*x(j);
      
      if ((abs(vloc - y(i)) > threshold) || isnan(abs(vloc - y(i))))
        {
          cout << "Operator * incorrect " << endl;
          abort();
        }
    }
  
  TinyMatrix<T, General, m, m> Cs;
  Cs = dot(A, B);
  for (int k = 0; k < m; k++)
    for (int i = 0; i < m; i++)
      {
        T vloc; SetComplexZero(vloc);
        for (int j = 0; j < m; j++)
          vloc += A(i, j)*B(j, k);
        
        if ((abs(vloc - Cs(i, k)) > threshold) || isnan(abs(vloc - Cs(i, k))))
          {
            cout << "Operator * incorrect " << endl;
            abort();
          }
      }
  
  // checking operators +, -
  C = A+B;
  D = A-B;
  for (int j = 0; j < m; j++)
    for (int i = 0; i < m; i++)
      if (  (abs(C(i, j) - A(i, j) - B(i, j)) > threshold) ||
            (abs(D(i, j) - A(i, j) + B(i, j)) > threshold) ||
            isnan(abs(C(i, j) - A(i, j) - B(i, j))) ||
            isnan(abs(D(i, j) - A(i, j) + B(i, j))) )
        {
          cout << "Operator +, - incorrect " << endl;
          abort();
        }
  
  // checking operator ==
  if ((C == A+B) && (C != A-B))
    {
    }
  else
    {
      cout << "Operator ==, != incorrect" << endl;
      abort();
    }
  
  if ((D==A+B) || (D != A-B))
    {
      cout << "Operator ==, != incorrect" << endl;
      abort();
    }

  // checking operator *=, +=, and -=
  D = A; D += B;
  if (D != C)
    {
      cout << "Operator += incorrect " << endl;
      abort();
    }

  C = A-B;
  D = A; D -= B;
  if (D != C)
    {
      cout << "Operator -= incorrect " << endl;
      abort();
    }

  C = Real_wp(1.5)*A;
  D = A; D *= Real_wp(1.5);
  if (D != C)
    {
      cout << "Operator *= incorrect " << endl;
      abort();
    }
  
  // checking operator =
  C = A;
  if (A != C)
    {
      cout << "operator = incorrect " << endl;
      abort();
    }
  
  C = 0.5;
  for (int j = 0; j < m; j++)
    for (int i = 0; i < m; i++)
      if ((abs(C(i, j) - 0.5) > threshold) || isnan(abs(C(i, j) - 0.5)))
        {
          cout << "operator = incorrect " << endl;
          abort();
        }
  
  // checking matrix-vector operations
  y.Fill(0);
  Mlt(A, x, y);
  TinyVector<T, m> x2;
  TinyVector<T, m> y2;
  y2 = dot(A, x);
  if (y != y2)
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  y.FillRand(); Mlt(Real_wp(1e-10), y);
  y2 = y; y2 += Real_wp(2.1)*dot(A, x);
  MltAdd(Real_wp(2.1), A, x, y);
  if (y != y2)
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }
  
  TinyMatrix<T, Symmetric, m, m> Acopy(A);
  x2.FillRand(); Mlt(Real_wp(1e-10), x2);
  for (int i = 0; i < m; i++)
    {
      GetRow(A, i, x2);
      for (int j = 0; j < m; j++)
        if ((abs(x2(j) - A(i, j)) > threshold) || isnan(abs(x2(j) - A(i, j))))
          {
            cout << "GetRow incorrect" << endl;
            abort();
          }
      
      SetRow(x, i, A);
      for (int j = 0; j < m; j++)
        if ((abs(x(j) - A(i, j)) > threshold) || isnan(abs(x(j) - A(i, j))))
          {
            cout << "SetRow incorrect" << endl;
            abort();
          }      
    }
  
  for (int j = 0; j < m; j++)
    {
      GetCol(A, j, y2);
      for (int i = 0; i < m; i++)
        if ((abs(y2(i) - A(i, j)) > threshold) || isnan(abs(y2(i) - A(i, j))))
          {
            cout << "GetCol incorrect" << endl;
            abort();
          }
      
      SetCol(y, j, A);
      for (int i = 0; i < m; i++)
        if ((abs(y(i) - A(i, j)) > threshold) || isnan(abs(y(i) - A(i, j))))
          {
            cout << "SetCol incorrect" << endl;
            abort();
          }      
    }
  
  A = Acopy;
  
  // matrix-matrix operations
  D = C; 
  Mlt(Real_wp(2.5), C);
  D *= Real_wp(2.5);
  if (C != D)
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Mlt(A, B, Cs);
  TinyMatrix<T, General, m, m> Ds;
  Ds = dot(A, B);
  if (Cs != Ds)
    {
      cout << "Mlt incorrect " << endl;
      abort();
    }
  
  D = A+B;
  Add(A, B, C);
  if (C != D)
    {
      cout << "Add incorrect" << endl;
      abort();
    }
  
  // check Zero, SetIdentity, Fill
  C.Zero();
  D.FillRand(); Mlt(Real_wp(1e-10), D);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(C(i, j)) > threshold) || isnan(C(i, j)))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  if ((!C.IsZero()) || (D.IsZero()))
    {
      cout << "IsZero incorrect" << endl;
      abort();
    }
  
  
  Cs.SetIdentity();
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
          SetComplexOne(vexact);
        
        if ((abs(vexact - Cs(i, j)) > threshold) || isnan(abs(vexact - Cs(i, j))))
          {
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  C.Fill(Real_wp(-1.6));
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if ((abs(C(i, j) + Real_wp(1.6)) > threshold) || isnan(abs(C(i, j) + Real_wp(1.6))))
        {
          cout << "Fill incorrect" << endl;
          abort();
        }
  
  C.Fill();
  int nb = 0;
  for (int i = 0; i < m; i++)
    for (int j = i; j < m; j++)
      {        
        if ((abs(C(i, j) - nb) > threshold) || isnan(abs(C(i, j) - nb)))
          {
            DISP(C(i, j)); DISP(nb); DISP(C(i, j)-nb);
            cout << "Fill incorrect" << endl;
            abort();
          }
        
        nb++;
      }

  // checking MaxAbs
  typename ClassComplexType<T>::Treal Emax, Etest(0);
  Emax = MaxAbs(A);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      Etest = max(Etest, abs(A(i, j)));
  
  if ((abs(Etest-Emax) > threshold) || isnan(abs(Etest-Emax)))
    {
      cout << "MaxAbs incorrect" << endl;
      abort();
    }
  
  // checking GetInverse  
  GetInverse(A, B);
  Mlt(A, B, Cs);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
          SetComplexOne(vexact);
        
        if ((abs(vexact - Cs(i, j)) > threshold) || isnan(abs(vexact - Cs(i, j))))
          {
            cout << "GetInverse incorrect" << endl;
            abort();
          }
      }  
  
  // checking GetEigenvaluesEigenvectors
#ifndef MONTJOIE_WITH_MPFR
  TinyVector<T, m> lambda;
  B = A;
  GetEigenvaluesEigenvectors(B, lambda, Cs);
  for (int j = 0; j < m; j++)
    {
      GetCol(Cs, j, x);
      Mlt(A, x, y);
      Real_wp vmax; T Li;
      SetComplexZero(vmax); SetComplexZero(Li);
      for (int i = 0; i < m; i++)
        {
          if (abs(x(i)) > vmax)
            {
              vmax = abs(x(i));
              Li = y(i)/x(i);
            }
        }
      
      if ((abs(Li - lambda(j)) > threshold) || isnan(abs(Li - lambda(j))))
        {
          cout << "GetEigenvalueEigenvectors incorrect" << endl;
          abort();
        }
      
      Mlt(1.0/Li, y);
      if (x != y)
        {
          cout << "GetEigenvalueEigenvectors incorrect" << endl;
          abort();
        }
    }
#endif
  
  // checking GetSquareRoot/GetAbsoluteValue
  A.FillRand(); Mlt(Real_wp(2e-10), A);
  B = dot(A, A); D = B;
  B.WriteText("AAh.dat");

#ifndef MONTJOIE_WITH_MPFR
  GetSquareRoot(B);
  A.WriteText("Ah.dat");
  B.WriteText("Bh.dat");
  C = dot(B, B);
  if (C != D)
    {
      cout << "GetSquareRoot incorrect" << endl;
      abort();
    }

  if (!IsComplexNumber(alpha))
    {
      C = A;
      GetAbsoluteValue(C);
      C.WriteText("Ch.dat");
      
      if (B != C)
	{
	  cout << "GetAbsoluteValue incorrect" << endl;
	  abort();
	}
    }
#endif

  // checking GetTangentialProjector/GetNormalProjector
  TinyVector<Real_wp, m> normale;
  normale.FillRand(); Mlt(1.0/Norm2(normale), normale);
  
  TinyMatrix<Real_wp, Symmetric, m, m> Proj;
  GetTangentialProjector(normale, Proj);
  
  TinyMatrix<Real_wp, Symmetric, m, m> Cmat;
  
  Cmat.SetIdentity();
  for (int i = 0; i < m; i++)
    for (int j = i; j < m; j++)
      Cmat(i, j) -= normale(i)*normale(j);
    
  if (Cmat != Proj)
    {
      DISP(Cmat); DISP(Proj);
      cout << "GetTangentialProjector incorrect" << endl;
      abort();
    }

  GetNormalProjector(normale, Proj);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      Cmat(i, j) = normale(i)*normale(j);
    
  if (Cmat != Proj)
    {
      cout << "GetNormalProjector incorrect" << endl;
      abort();
    }

  // checking SetDiagonal
  GetRand(alpha);   GetRand(beta);
  A.SetDiagonal(alpha);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        T vexact; SetComplexZero(vexact);
        if (i == j)
	  vexact = alpha;
        
        if ((abs(vexact - A(i, j)) > threshold) || isnan(abs(vexact - A(i, j))))
          {
            cout << "SetDiagonal incorrect" << endl;
            abort();
          }
      }

  typedef typename ClassComplexType<T>::Tcplx Tcplx;
  TinyMatrix<Tcplx, Symmetric, m, m> Ac, Bc, Cc, Dc;
  GenerateRandomMatrix(A); GenerateRandomMatrix(Ac);
  GenerateRandomMatrix(B); GenerateRandomMatrix(Bc);
  GenerateRandomMatrix(C); GenerateRandomMatrix(Cc);

  Tcplx alpha_c;
  GetRand(alpha_c);
  
  // testing operator +, - with expressions
  D = alpha*(A+B) - beta*C;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha*(A(i, j) + B(i, j)) - beta*C(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator +/- incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) + B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = A + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) + Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) - B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) - Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }


  // testing operator *, / with expressions
  D = alpha*(A*B) / (beta*C);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha*(A(i, j) * B(i, j)) / (beta*C(i, j));
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator *, / incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) * B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = A * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) * Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) / B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) / Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  // testing operator + with scalar
  D = alpha + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  D = B + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c + B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha_c + B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = A + alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) + alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac + alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) + alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha + Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha + Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator + incorrect" << endl;
	    abort();
	  }
      }

  // testing operator - with scalar
  D = alpha - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha - B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  D = B - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = -alpha + B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c - B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha_c - B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = A - alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) - alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac - alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) - alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha - Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha - Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }


  // testing operator * with scalar
  D = alpha * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha * B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  D = B * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = B(i, j) * alpha;
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c * B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha_c * B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = A * alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) * alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac * alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) * alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha * Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha * Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator * incorrect" << endl;
	    abort();
	  }
      }

  // testing operator / with scalar
  D = alpha / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = alpha / B(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  D = B / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = B(i, j) / alpha;
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha_c / B;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha_c / B(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = A / alpha_c;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = A(i, j) / alpha_c;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = Ac / alpha;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = Ac(i, j) / alpha;
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  Dc = alpha / Bc;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	Tcplx vexact = alpha / Bc(i, j);
	if ( (abs(vexact-Dc(i, j)) > threshold) || isnan(abs(vexact-Dc(i, j))))
	  {
	    cout << "Operator / incorrect" << endl;
	    abort();
	  }
      }

  // testing opposite
  D = -A;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
	T vexact = -A(i, j);
	if ( (abs(vexact-D(i, j)) > threshold) || isnan(abs(vexact-D(i, j))))
	  {
	    cout << "Operator - incorrect" << endl;
	    abort();
	  }
      }

  // testing Norm2_column
  typename ClassComplexType<T>::Treal am, aref;
  int first_row = 1;
  for (int i = 0; i < m; i++)
    {
      am = Norm2_Column(A, first_row, i);
      SetComplexZero(aref);
      for (int j = first_row; j < m; j++)
	aref += absSquare(A(j, i));

      aref = sqrt(aref);
      if ( (abs(am-aref) > threshold) || isnan(abs(am-aref)))
	{
	  cout << "Norm2_column incorrect" << endl;
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
    TinyMatrix<Real_wp, General, 5, 4> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Real_wp, General, 7, 7> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Complex_wp, General, 5, 4> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Complex_wp, General, 7, 7> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Real_wp, General, 2, 2> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Real_wp, General, 3, 3> A;
    CheckTinyGeneralMatrix(A);
  }

  {
    TinyMatrix<Real_wp, Symmetric, 5, 5> A;
    CheckTinySymmetricMatrix(A);
  }

  {
    TinyMatrix<Real_wp, Symmetric, 2, 2> A;
    CheckTinySymmetricMatrix(A);
  }

  {
    TinyMatrix<Real_wp, Symmetric, 3, 3> A;
    CheckTinySymmetricMatrix(A);
  }

  {
    TinyMatrix<Complex_wp, Symmetric, 5, 5> A;
    CheckTinySymmetricMatrix(A);
  }

  {
    TinyMatrix<Complex_wp, Symmetric, 2, 2> A;
    CheckTinySymmetricMatrix(A);
  }

  {
    TinyMatrix<Complex_wp, Symmetric, 3, 3> A;
    CheckTinySymmetricMatrix(A);
  }
  
  if (overall_success)
    cout << "All tests passed successfully " << endl;
  
  return 0;
}
