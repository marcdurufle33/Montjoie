#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

template<class T, class Allocator>
int TestSkyLineMatrix(Matrix<T, Symmetric, SymColSkyLine, Allocator>& A)
{
  Real_wp threshold = 1e4*epsilon_machine;
  
  // basic matrix
  // A = | 2    4.5   0     0.5    0    0 |
  //     | 4.5  3.0  -1.2   0.8    0    0 |
  //     | 0   -1.2   5.2    0    1.3   0 |
  //     | 0.5  0.8    0    2.7    0   1.9|
  //     | 0     0    1.3    0    3.3  1.2|
  //     | 0     0     0    1.9   1.2   2.4|
  
  // Ind(i+1) - Ind(i) is the number of non-zero entries in the column i
  // (upper part of the matrix)
  IVect Ind(7);
  Ind(0) = 0;
  Ind(1) = 1;
  Ind(2) = 3;
  Ind(3) = 5;
  Ind(4) = 9;
  Ind(5) = 12;
  Ind(6) = 15;
  
  // then all the values of the upper part (by columns)
  Vector<T> Val;
  string data("2 4.5 3.0 -1.2 5.2 0.5 0.8 0 2.7 1.3 0 3.3 1.9 1.2 2.4 ");
  istringstream stream_data(data);
  Val.ReadText(stream_data);  
  
  // after SetData, Ind and Val will be empty
  A.SetData(Ind, Val);
  
  if (A.GetM() != 6)
    {
      cout << "GetM incorrect " << endl;
      abort();      
   }

  if (A.GetDataSize() != 15)
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  Matrix<T, Symmetric, SymColSkyLine, Allocator> B(A);
  //B.WriteText("Ah,dat");
  
  // testing Fill
  A.Fill(2.0);
  Real_wp* a_ = A.GetData();
  for (int i = 0; i < A.GetDataSize(); i++)
    if ((abs(a_[i] - 2.0) > threshold) || isnan(a_[i]))
      {
        cout << "Fill incorrect" << endl;
        abort();
      }

#ifndef MONTJOIE_WITH_MPFR
  A.Zero();
#else
  A.Fill(0);
#endif
  
  for (int i = 0; i < A.GetDataSize(); i++)
    if ((abs(a_[i]) > threshold) || isnan(a_[i]))
      {
        cout << "Zero incorrect" << endl;
        abort();
      }
  
  A.FillRand();
  A.Clear();
  if ((A.GetM() != 0) || (A.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  Matrix<T, Symmetric, ArrayRowSymSparse> Acsr;
  Acsr.ReadText("example/unit/MhSparse.dat");
  
  // conversion to skyline matrix
  Copy(Acsr, A);
  
  if ((A.GetM() != 30) || (A.GetN() != 30))
    {
      cout << "GetM incorrect " << endl;
      abort();
    }
  
  if (A.GetDataSize() != 315)
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  // testing operator()
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      if ((abs(A(i, j) - Acsr(i, j)) > threshold) || isnan(abs(A(i, j) - Acsr(i, j))))
        {
          cout << "operator() incorrect" << endl;
          abort();
        }
  
  // testing WriteText
  A.WriteText("mat.dat");
  Acsr.ReadText("mat.dat");
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      if ((abs(A(i, j) - Acsr(i, j)) > threshold) || isnan(abs(A(i, j) - Acsr(i, j))))
        {
          cout << "WriteText incorrect" << endl;
          abort();
        }
  
  // testing GetCholesky/SolveCholesky
  GetCholesky(A);
  // A.WriteText("chol.dat");
  Vector<Real_wp> X(A.GetM()), Y(A.GetM()), Yref;
  X.Fill(); Y.Fill(0);
  Mlt(Acsr, X, Y);
  
  X = Y;
  SolveCholesky(SeldonNoTrans, A, X);
  SolveCholesky(SeldonTrans, A, X);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(X(i) - Real_wp(i)) > threshold) || isnan(abs(X(i) - Real_wp(i))))
      {
        cout << "GetCholesky/SolveCholesky incorrect" << endl;
        abort();
      }
  
  // testing MltCholesky
  MltCholesky(SeldonTrans, A, X);
  MltCholesky(SeldonNoTrans, A, X);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(X(i) - Y(i)) > threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "MltCholesky incorrect" << endl;
        abort();
      }  
  
  // testing MltAdd
  Copy(Acsr, A);
  
  Yref = Y;
  X.Fill();
  MltAdd(Real_wp(1), A, X, Real_wp(0), Y);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(Yref(i) - Y(i)) > threshold) || isnan(abs(Yref(i) - Y(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }  
  
  // testing GetLU/SolveLU
  X.Reallocate(B.GetM()); Y.Reallocate(B.GetM());
  X.Fill(); Y.Fill(0);
  Mlt(B, X, Y);
  
  A = B;
  GetLU(B);
  // B.WriteText("ldlt.dat");
  X = Y;
  SolveLU(SeldonNoTrans, B, X);
  for (int i = 0; i < B.GetM(); i++)
    if ((abs(X(i) - Real_wp(i)) > threshold) || isnan(abs(X(i) - Real_wp(i))))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }
  
  return 0;
}


template<class T, class Allocator>
int TestSkyLineMatrix(Matrix<T, General, RowSkyLine, Allocator>& A)
{
  Real_wp threshold = 1e4*epsilon_machine;
  
  // basic matrix
  // A = | 2    4.5   0     0.5    0    0 |
  //     | 4.5  3.0  -1.2   0.8    0    0 |
  //     | 0   -1.2   5.2    0    1.3   0 |
  //     | 0.5  0.8    0    2.7    0   1.9|
  //     | 0     0    1.3    0    3.3  1.2|
  //     | 0     0     0    1.9   1.2   2.4|
  A.Reallocate(6, 6);
  A.ReallocateRow(0, 4, 0);
  A.Value(0, 0) = 2.0; A.Value(0, 1) = 4.5; A.Value(0, 2) = 0; A.Value(0, 3) = 0.5;
  
  A.ReallocateRow(1, 4, 0);
  A.Value(1, 0) = 4.5; A.Value(1, 1) = 3.0; A.Value(1, 2) = -1.2; A.Value(1, 3) = 0.8;
  
  A.ReallocateRow(2, 4, 1);
  A.Value(2, 0) = -1.2; A.Value(2, 1) = 5.2; A.Value(2, 2) = 0; A.Value(2, 3) = 1.3;
  
  A.ReallocateRow(3, 6, 0);
  A.Value(3, 0) = 0.5; A.Value(3, 1) = 0.8; A.Value(3, 2) = 0; A.Value(3, 3) = 2.7;
  A.Value(3, 4) = 0; A.Value(3, 5) = 1.9;
  
  A.ReallocateRow(4, 4, 2);
  A.Value(4, 0) = 1.3; A.Value(4, 1) = 0; A.Value(4, 2) = 3.3; A.Value(4, 3) = 1.2;
  
  A.ReallocateRow(5, 3, 3);
  A.Value(5, 0) = 1.9; A.Value(5, 1) = 1.2; A.Value(5, 2) = 2.4;
  
  if (A.GetM() != 6)
    {
      cout << "GetM incorrect " << endl;
      abort();      
   }

  if (A.GetDataSize() != 25)
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  Matrix<T, General, RowSkyLine> B(A);
  
  // testing Fill
  A.Fill(2.0);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetRowSize(i); j++)
      if ((abs(A.Value(i, j) - 2.0) > threshold) || isnan(abs(A.Value(i, j) - 2.0)))
        {
          cout << "Fill incorrect" << endl;
          abort();
        }
  
#ifndef MONTJOIE_WITH_MPFR
  A.Zero();
#else
  A.Fill(0);
#endif
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetRowSize(i); j++)
      if ((abs(A.Value(i, j)) > threshold) || isnan(abs(A.Value(i, j))))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  A.FillRand();
  A.Clear();
  if ((A.GetM() != 0) || (A.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  Matrix<T, General, ArrayRowSparse> Acsr;
  Acsr.ReadText("example/unit/MhSparse.dat");
  
  // conversion to skyline matrix
  Copy(Acsr, A);
  
  if ((A.GetM() != 30) || (A.GetN() != 30))
    {
      cout << "GetM incorrect " << endl;
      abort();
    }
  
  // testing operator()
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      if ((abs(A(i, j) - Acsr(i, j)) > threshold) || isnan(abs(A(i, j) - Acsr(i, j))))
        {
          DISP(A(i, j));
          DISP(Acsr(i, j));
          cout << "operator() incorrect" << endl;
          abort();
        }
  
  // testing WriteText
  A.WriteText("mat.dat");
  Acsr.ReadText("mat.dat");
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      if ((abs(A(i, j) - Acsr(i, j)) > threshold) || isnan(abs(A(i, j) - Acsr(i, j))))
        {
          cout << "WriteText incorrect" << endl;
          abort();
        }
  
  // testing GetLU/SolveLU
  GetLU(A);
  //A.WriteText("lu.dat");
  Vector<Real_wp> X(A.GetM()), Y(A.GetM()), Yt(A.GetM()), Yref;
  X.Fill(); Y.Fill(0); Yt.Fill(0);
  Mlt(Acsr, X, Y);
  Mlt(SeldonTrans, Acsr, X, Yt);
  
  X = Y;
  SolveLU(SeldonNoTrans, A, X);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(X(i) - Real_wp(i)) > threshold) || isnan(abs(X(i) - Real_wp(i))))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }

  X = Yt;
  SolveLU(SeldonTrans, A, X);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(X(i) - Real_wp(i)) > threshold) || isnan(abs(X(i) - Real_wp(i))))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }

  
  // testing MltAdd
  Copy(Acsr, A);
  
  Yref = Y;
  X.Fill();
  MltAdd(Real_wp(1), A, X, Real_wp(0), Y);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(Yref(i) - Y(i)) > threshold) || isnan(abs(Yref(i) - Y(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }  

  MltAdd(Real_wp(1), SeldonTrans, A, X, Real_wp(0), Y);
  for (int i = 0; i < A.GetM(); i++)
    if ((abs(Yt(i) - Y(i)) > threshold) || isnan(abs(Yt(i) - Y(i))))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }  
  
  A = B;
  // other way to construct matrix, with AddInteraction
  // A = | 2    4.5   0     0.5    0    0 |
  //     | 4.5  3.0  -1.2   0.8    0    0 |
  //     | 0   -1.2   5.2    0    1.3   0 |
  //     | 0.5  0.8    0    2.7    0   1.9|
  //     | 0     0    1.3    0    3.3  1.2|
  //     | 0     0     0    1.9   1.2   2.4|
  B.Reallocate(6, 6);
  B.AddInteraction(0, 0, 2.0);
  B.AddInteraction(0, 1, 4.5);
  B.AddInteraction(0, 3, 0.5);
  
  B.AddInteraction(1, 1, 3.0);
  B.AddInteraction(1, 2, -1.2);
  B.AddInteraction(1, 0, 4.5);
  B.AddInteraction(1, 3, 0.8);
  
  B.AddInteraction(2, 2, 5.2);
  IVect col(4); Vector<T> val(4);
  col(0) = 1; val(0) = -1.2;
  col(1) = 4; val(1) = 1.3;
  B.AddInteractionRow(2, 2, col, val);
  
  col(0) = 0; val(0) = 0.5;
  col(1) = 1; val(1) = 0.8;
  col(2) = 3; val(2) = 2.7;
  col(3) = 5; val(3) = 1.9;
  B.AddInteractionRow(3, 4, col, val);
  
  B.AddInteraction(4, 5, 1.2);
  col(0) = 2; val(0) = 1.3;
  col(1) = 4; val(1) = 3.3;
  B.AddInteractionRow(4, 2, col, val);
  
  B.AddInteraction(5, 3, 1.9);
  col(0) = 4; val(0) = 1.2;
  col(1) = 5; val(1) = 2.4;
  B.AddInteractionRow(5, 2, col, val);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetM(); j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
        {
          cout << "AddInteraction/AddInteractionRow incorrect" << endl;
          abort();
        }
 
  return 0;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  bool overall_success = true;
  
  int n = toInteger(-log10(epsilon_machine));
  cout.precision(n);
  
  {
    Matrix<Real_wp, Symmetric, SymColSkyLine> A;
    int success = TestSkyLineMatrix(A);
    if (success != 0)
      {
        overall_success = false;
        cout << "Problem on symmetric SkyLine Matrix" << endl;
      }
  }

  {
    Matrix<Real_wp, General, RowSkyLine> A;
    int success = TestSkyLineMatrix(A);
    if (success != 0)
      {
        overall_success = false;
        cout << "Problem on general SkyLine Matrix" << endl;
      }
  }
  
  if (overall_success)
    cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}

