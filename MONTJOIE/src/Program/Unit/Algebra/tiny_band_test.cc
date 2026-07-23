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

template<class T, int d>
void CheckTinyBandMatrix(TinyBandMatrix<T, d>& A, int n)
{
  TinyVector<T, d>::threshold = threshold;
  
  A.Reallocate(n, n);
  
  // checking GetM/GetN
  if ( (A.GetM() != n) || (A.GetN() != n) )
    {
      cout << "GetM/GetN incorrect" << endl;
      abort();
    }
  
  // checking Clear
  A.Clear();
  
  if (A.GetM() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // checking Zero
  A.Reallocate(n, n);
  
#ifndef MONTJOIE_WITH_MPFR
  A.Zero();
#else
  A.Fill(0);
#endif
  
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      if ((abs(A(i, j)) > threshold) || isnan(A(i, j)))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  // checking SetIdentity
  A.SetIdentity();
  
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val_ex; SetComplexZero(val_ex);
        if (i == j)
          SetComplexOne(val_ex);
        
        if ((abs(A(i, j)-val_ex) > threshold) || isnan(abs(A(i, j)-val_ex)))
          {
            DISP(val_ex); DISP(A(i, j)); DISP(abs(A(i, j)-val_ex));
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  // checking conversion from sparse matrices to band matrices
  Matrix<T, General, ArrayRowSparse> Acsr(n, n);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val;
        SetComplexReal(0.2*abs(j-i) + 0.01*(j+1), val);
        if (j == i)
          val += 3.5;
        
        Acsr.AddInteraction(i, j, val);
      }
  
  Copy(Acsr, A);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      if ((abs(Acsr(i, j) - A(i, j)) > threshold) || isnan(abs(Acsr(i, j) - A(i, j))))
        {
          cout << "Copy incorrect" << endl;
          abort();
        }
  
  // checking Add
  TinyBandMatrix<T, d> B;
  B.Reallocate(n, n);
  B.SetIdentity();
  Add(-1.8, A, B);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val_ex = -1.8*A(i, j);
        if (i == j)
          val_ex += 1.0;
        
        if ((abs(B(i, j) - val_ex) > threshold) || isnan(abs(B(i, j) - val_ex)))
          {
            cout << "Add incorrect" << endl;
            abort();
          }
      }
  
  // checking MltAdd
  Vector<T> x(n), y(n), y2(n);
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, n); y2 = y;
  MltAdd(Real_wp(2.2), A, x, Real_wp(1.2), y);
  MltAdd(Real_wp(2.2), Acsr, x, Real_wp(1.2), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonTrans, A, x, Real_wp(0.8), y);
  MltAdd(Real_wp(1.3), SeldonTrans, Acsr, x, Real_wp(0.8), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        cout << "MltAdd with transpose incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonConjTrans, A, x, Real_wp(0.8), y);
  MltAdd(Real_wp(1.3), SeldonConjTrans, Acsr, x, Real_wp(0.8), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        DISP(y(i)); DISP(y2(i)); DISP(abs(y(i)-y2(i)));
        cout << "MltAdd with conjugate transpose incorrect" << endl;
        abort();
      }
  
  // with a small matrix and vectors
  TinyBandMatrix<T, d> Cs;
  Cs.Reallocate(d+2, d+2);  
  Vector<T> xs, ys, ys2; Matrix<T> Cf(d+2, d+2);
  Cs.FillRand(); T val;
  for (int i = 0; i < Cs.GetM(); i++)
    for (int j = 0; j < Cs.GetN(); j++)
      {
        GetRandNumber(val);
        if (Cs(i, j) != 0)
          Cs.Get(i, j) = val;
        else
          SetComplexZero(val);
        
        Cf(i, j) = val;
      }
  
  GenerateRandomVector(xs, Cs.GetM());
  GenerateRandomVector(ys, Cs.GetM()); ys2 = ys;
  MltAdd(Real_wp(2.7), Cs, xs, Real_wp(1.6), ys);
  MltAdd(Real_wp(2.7), Cf, xs, Real_wp(1.6), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonTrans, Cs, xs, Real_wp(0.8), ys);
  MltAdd(Real_wp(1.3), SeldonTrans, Cf, xs, Real_wp(0.8), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        cout << "MltAdd with transpose incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonConjTrans, Cs, xs, Real_wp(0.8), ys);
  MltAdd(Real_wp(1.3), SeldonConjTrans, Cf, xs, Real_wp(0.8), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        cout << "MltAdd with conjugate transpose incorrect" << endl;
        abort();
      }  
  
  // checking GetLU/ SolveLU
  TinyBandMatrix<T, d> C(A);
  MltAdd(Real_wp(1.0), A, x, Real_wp(0.0), y);
  GetLU(A, B);
  SolveLU(B, y);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - x(i)) > threshold) || isnan(abs(y(i) - x(i))))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }

  // scaling matrix
  B = C;
  Vector<typename ClassComplexType<T>::Treal> scale_row(n);
  scale_row.FillRand(); scale_row *= 1e-9;
  ScaleLeftMatrix(B, scale_row);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(B(i, j) - scale_row(i)*C(i, j)) > threshold) || isnan(abs(B(i, j) - scale_row(i)*C(i, j))))
        {
          DISP(i); DISP(j); DISP(B(i, j));
          DISP(scale_row(i)*C(i, j));
          cout << "ScaleLeftMatrix incorrect" << endl;
          abort();
        }
  
  // computing the sum
  GetRowSum(scale_row, C);
  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 0.0;
      for (int j = 0; j < n; j++)
        sum += abs(C(i, j));
      
      if ((abs(sum - scale_row(i)) > threshold) || isnan(abs(sum - scale_row(i))))
        {
          cout << "GetRowSum incorrect" << endl;
          abort();
        }
    }
  
}

template<class T, int d, int m>
void CheckTinyArrowMatrix(TinyArrowMatrix<T, d, m>& A, int n)
{
  TinyVector<T, d>::threshold = threshold;
  
  A.Reallocate(n, n);
  
  // checking GetM/GetN
  if ( (A.GetM() != n) || (A.GetN() != n) )
    {
      cout << "GetM/GetN incorrect" << endl;
      abort();
    }
  
  // checking Clear
  A.Clear();
  
  if (A.GetM() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // checking Zero
  A.Reallocate(n, n);
  
#ifndef MONTJOIE_WITH_MPFR
  A.Zero();
#else
  A.Fill(0);
#endif
  
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      if ((abs(A(i, j)) > threshold) || isnan(A(i, j)))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  // checking SetIdentity
  A.SetIdentity();
  
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val_ex; SetComplexZero(val_ex);
        if (i == j)
          SetComplexOne(val_ex);
        
        if ((abs(A(i, j)-val_ex) > threshold) || isnan(abs(A(i, j)-val_ex)))
          {
            DISP(i); DISP(j); DISP(A(i, j)); DISP(val_ex); DISP(abs(A(i, j)-val_ex));
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  // checking conversion from sparse matrices to band matrices
  Matrix<T, General, ArrayRowSparse> Acsr(n, n);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val;
        SetComplexReal(0.2*abs(j-i) + 0.01*(j+1), val);
        if (j == i)
          val += 3.5;
        
        Acsr.AddInteraction(i, j, val);
      }
  
  for (int j = 0; j < m; j++)
    for (int i = 0; i < n-m; i++)
      {
        T val;
        SetComplexReal(-0.02*Real_wp(i+1) + 0.1*Real_wp(j*j), val);
        Acsr.AddInteraction(i, n-1-j, val);
        Acsr.AddInteraction(n-1-j, i, -val);
      }
  
  Copy(Acsr, A);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      if ((abs(Acsr(i, j) - A(i, j)) > threshold) || isnan(abs(Acsr(i, j) - A(i, j))))
        {
          cout << "Copy incorrect" << endl;
          abort();
        }

  for (int j = 0; j < m; j++)
    for (int i = 0; i < n-m; i++)
      if ((abs(Acsr(i, n-1-j) - A(i, n-1-j)) > threshold) || isnan(abs(Acsr(i, n-1-j) - A(i, n-1-j))))
        {
          cout << "Copy incorrect" << endl;
          abort();
        }
  
  // checking Add
  TinyArrowMatrix<T, d, m> B;
  B.Reallocate(n, n);
  B.SetIdentity();
  Add(-1.8, A, B);
  for (int i = 0; i < n; i++)
    for (int j = max(0, i-d); j <= min(n-1, i+d); j++)
      {
        T val_ex = -1.8*A(i, j);
        if (i == j)
          val_ex += 1.0;
        
        if ((abs(B(i, j) - val_ex) > threshold) || isnan(abs(B(i, j) - val_ex)))
          {
            cout << "Add incorrect" << endl;
            abort();
          }
      }
  
  // checking MltAdd
  Vector<T> x(n), y(n), y2(n);
  x.Fill(); y.Fill(1); y2.Fill(1);
  MltAdd(Real_wp(2.2), A, x, Real_wp(1.2), y);
  MltAdd(Real_wp(2.2), Acsr, x, Real_wp(1.2), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonTrans, A, x, Real_wp(0.8), y);
  MltAdd(Real_wp(1.3), SeldonTrans, Acsr, x, Real_wp(0.8), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        cout << "MltAdd with transpose incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonConjTrans, A, x, Real_wp(0.8), y);
  MltAdd(Real_wp(1.3), SeldonConjTrans, Acsr, x, Real_wp(0.8), y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(abs(y(i) - y2(i))))
      {
        cout << "MltAdd with conjugate transpose incorrect" << endl;
        abort();
      }
  
  // with a small matrix and vectors
  TinyArrowMatrix<T, d, m> Cs;
  Cs.Reallocate(d+m+2, d+m+2);  
  Vector<T> xs, ys, ys2; Matrix<T> Cf(d+m+2, d+m+2);
  Cs.FillRand(); T val;
  for (int i = 0; i < Cs.GetM(); i++)
    for (int j = 0; j < Cs.GetN(); j++)
      {
        GetRandNumber(val);
        if (Cs(i, j) != 0)
          Cs.Get(i, j) = val;
        else
          SetComplexZero(val);
        
        Cf(i, j) = val;
      }
  
  GenerateRandomVector(xs, Cs.GetM());
  GenerateRandomVector(ys, Cs.GetM()); ys2 = ys;
  MltAdd(Real_wp(2.7), Cs, xs, Real_wp(1.6), ys);
  MltAdd(Real_wp(2.7), Cf, xs, Real_wp(1.6), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        DISP(ys(i)); DISP(ys2(i)); DISP(abs(ys(i)-ys2(i)));
        cout << "MltAdd incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonTrans, Cs, xs, Real_wp(0.8), ys);
  MltAdd(Real_wp(1.3), SeldonTrans, Cf, xs, Real_wp(0.8), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        cout << "MltAdd with transpose incorrect" << endl;
        abort();
      }

  MltAdd(Real_wp(1.3), SeldonConjTrans, Cs, xs, Real_wp(0.8), ys);
  MltAdd(Real_wp(1.3), SeldonConjTrans, Cf, xs, Real_wp(0.8), ys2);
  for (int i = 0; i < Cs.GetM(); i++)
    if ((abs(ys(i) - ys2(i)) > threshold) || isnan(abs(ys(i) - ys2(i))))
      {
        cout << "MltAdd with conjugate transpose incorrect" << endl;
        abort();
      }  
  
  // checking GetLU/ SolveLU
  MltAdd(Real_wp(1.0), A, x, Real_wp(0.0), y);
  GetLU(A, B, true);
  SolveLU(B, y);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - x(i)) > threshold) || isnan(abs(y(i) - x(i))))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }
  
  // checking operator =
  TinyArrowMatrix<T, 2, 1> C;

  C.Reallocate(n, n);
  C.FillRand();
  C *= Real_wp(1e-9);
  
  B = C;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((B(i, j) != C(i, j)) || isnan(B(i, j)) || isnan(C(i, j)) )
        {
          DISP(i); DISP(j);
          DISP(B(i, j)); DISP(C(i, j));
          abort();
        }
  
  // and operation Add
  A.FillRand(); A *= Real_wp(1e-9);
  B = A;
  Add(Real_wp(1.34), C, B);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(B(i, j) - A(i, j) - 1.34*C(i, j)) > threshold) || isnan(abs(B(i, j) - A(i, j) - 1.34*C(i, j))))
        {
          DISP(i); DISP(j);
          DISP(B(i, j)); DISP(A(i, j) + 1.34*C(i, j));
          abort();
        }

  // scaling matrix
  B = C;
  Vector<typename ClassComplexType<T>::Treal> scale_row(n);
  scale_row.FillRand(); scale_row *= Real_wp(1e-9);
  ScaleLeftMatrix(B, scale_row);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(B(i, j) - scale_row(i)*C(i, j)) > threshold) || isnan(abs(B(i, j) - scale_row(i)*C(i, j))))
        {
          DISP(i); DISP(j); DISP(B(i, j));
          DISP(scale_row(i)*C(i, j));
          cout << "ScaleLeftMatrix incorrect" << endl;
          abort();
        }
  
  // computing the sum
  GetRowSum(scale_row, C);
  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 0.0;
      for (int j = 0; j < n; j++)
        sum += abs(C(i, j));
      
      if (abs(sum - scale_row(i)) > threshold)
        {
          cout << "GetRowSum incorrect" << endl;
          abort();
        }
    }
  
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  bool overall_success = true;
  threshold = 1e4*epsilon_machine;
  
  {
    // testing band-matrices
    TinyBandMatrix<Real_wp, 3> A;
    int n = 12;
    
    CheckTinyBandMatrix(A, n);
  }

  {
    // testing band-matrices
    TinyBandMatrix<Real_wp, 2> A;
    int n = 9;
    
    CheckTinyBandMatrix(A, n);
  }
  
  /*{
    // testing arrow matrices
    TinyArrowMatrix<Real_wp, 3, 5> A;
    int n = 15;
    
    CheckTinyArrowMatrix(A, n);
    }*/

  
  {
    // testing band-matrices
    TinyBandMatrix<Complex_wp, 3> A;
    int n = 12;
    
    CheckTinyBandMatrix(A, n);
  }

  {
    // testing band-matrices
    TinyBandMatrix<Complex_wp, 2> A;
    int n = 9;
    
    CheckTinyBandMatrix(A, n);
  }

  /*
  {
    // testing arrow matrices
    TinyArrowMatrix<Complex_wp, 3, 5> A;
    int n = 15;
    
    CheckTinyArrowMatrix(A, n);
  }
  */
  
  if (overall_success)
    cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
