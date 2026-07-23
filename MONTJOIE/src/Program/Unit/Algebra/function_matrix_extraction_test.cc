#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-12;

namespace Montjoie
{
  template<class MatrixSparse, class Vector1>
  void ScaleLeftMatrix(MatrixSparse& A, const Vector1& x) {}
}

using namespace Montjoie;

void GenerateRandomPermutation(int n, IVect& permut)
{
  Vector<bool> NumUsed(n);
  NumUsed.Fill(false);
  permut.Reallocate(n);
  permut.Fill(-1);
  int nb = 0;
  // premiere iteration
  for (int i = 0; i < n; i++)
    {
      int i2 = rand()%n;
      if (!NumUsed(i2))
        {
          NumUsed(i2) = true;
          permut(i) = i2;
          nb++;
        }
    }
  
  while (nb < n)
    {
      // on recupere les numeros non-selectionnes
      IVect non_selec(n-nb);
      int k = 0;
      for (int i = 0; i < n; i++)
        if (!NumUsed(i))
          non_selec(k++) = i;
      
      // iteration suivante
      for (int i = 0; i < n; i++)
        if (permut(i) == -1)
          {
            int i2 = rand()%(n-nb);
            int j = non_selec(i2);
            if (!NumUsed(j))
              {
                NumUsed(j) = true;
                permut(i) = j;
                nb++;
              }
          }
    }
}

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
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
void GenerateRandomVector(Vector<T>& x, int m)
{
  x.Reallocate(m);
  for (int i = 0; i < m; i++)
    GetRand(x(i));  
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRand(x);
      A.Set(i, j, x);
    }
}

template<class T>
T AbsCplx(const T& x)
{
  return abs(x);
}

template<class T>
T AbsCplx(const complex<T>& x)
{
  return abs(real(x)) + abs(imag(x));
}

template<class MatrixSparse>
void TestFunctionExtract(MatrixSparse& A, bool norm_cplx = true)
{
  // testing functions in FunctionMatrixExtraction.cxx
  int m = 30, n = 20, nnz = 100;
  if (IsSymmetricMatrix(A))
    n = m;
  
  IVect permut_row, permut_col;
  GenerateRandomPermutation(n, permut_col);
  GenerateRandomPermutation(m, permut_row);
  
  IVect col_num(3), row_num(4);
  for (int k = 0; k < col_num.GetM(); k++)
    col_num(k) = permut_col(k);

  for (int k = 0; k < row_num.GetM(); k++)
    row_num(k) = permut_row(k);
  
  GenerateRandomMatrix(A, m, n, nnz);
  
  Vector<bool> RowUsed(m), ColUsed(n);
  RowUsed.Fill(false); ColUsed.Fill(false);
  for (int i = 0; i < row_num.GetM(); i++)
    RowUsed(row_num(i)) = true;

  for (int i = 0; i < col_num.GetM(); i++)
    ColUsed(col_num(i)) = true;

  // testing GetCol
  typedef typename MatrixSparse::entry_type Complexe;
  Vector<Vector<Complexe, VectSparse>, VectSparse> V;
  GetCol(A, col_num, V);  
  
  for (int jloc = 0; jloc < col_num.GetM(); jloc++)
    {
      int j = col_num(jloc);
      for (int i = 0; i < m; i++)
        if ((A(i, j) != V.Value(jloc)(i)) || isnan(A(i, j)) || isnan(V.Value(jloc)(i)))
          {
            cout << "GetCol incorrect" << endl;
            DISP(i); DISP(j); DISP(A(i, j)); DISP(V.Value(jloc)(i));
            abort();
          }
    }

  // testing EraseCol
  MatrixSparse B(A);
  EraseCol(col_num, A);
  
  Complexe zero; SetComplexZero(zero);
  for (int j = 0; j < n; j++)
    {
      if (ColUsed(j))
        {
          for (int i = 0; i < m; i++)
            if ((A(i, j) != zero) || isnan(A(i, j)))
              {
                cout << "EraseCol incorrect" << endl;
                DISP(i); DISP(j); DISP(A(i, j));
                abort();
              }
        }
      else
        {
          for (int i = 0; i < m; i++)
            {              
              if (IsSymmetricMatrix(A) && ColUsed(i))
                {
                  if ((A(i, j) != zero) || isnan(A(i, j)))
                    {
                      cout << "EraseCol incorrect" << endl;
                      DISP(i); DISP(j); DISP(A(i, j));
                      abort();
                    }
                }
              else if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
                {
                  cout << "EraseCol incorrect" << endl;
                  DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
                  abort();
                }
            }
        }
    }
  
  // testing EraseRow
  A = B;
  EraseRow(row_num, A);
  for (int i = 0; i < m; i++)
    {
      if (RowUsed(i))
        {
          for (int j = 0; j < n; j++)
            if ((A(i, j) != zero) || isnan(A(i, j)))
              {
                cout << "EraseRow incorrect" << endl;
                DISP(i); DISP(j); DISP(A(i, j));
                abort();
              }
        }
      else
        {
          for (int j = 0; j < n; j++)
            {
              if (IsSymmetricMatrix(A) && RowUsed(j))
                {
                  if ((A(i, j) != zero) || isnan(A(i, j)))
                    {
                      cout << "EraseRow incorrect" << endl;
                      DISP(i); DISP(j); DISP(A(i, j));
                      abort();
                    }
                }
              else if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
                {
                  cout << "EraseRow incorrect" << endl;
                  DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
                  abort();
                }
            }
        }
    }
  
  // testing GetRowSum
  Vector<Real_wp> sum_row, sum_col;
  GetRowSum(sum_row, B);
  
  for (int i = 0; i < m; i++)
    {
      Real_wp val = 0;
      if (norm_cplx)
        for (int j = 0; j < n; j++)
          val += abs(B(i, j));
      else
        for (int j = 0; j < n; j++)
          val += AbsCplx(B(i, j));
      
      if ((abs(sum_row(i) - val) > threshold) || isnan(sum_row(i)) || isnan(val))
        {
          cout << "GetRowSum incorrect" << endl;
          DISP(i); DISP(sum_row(i)); DISP(val);
          abort();
        }
    }

  GetRowColSum(sum_row, sum_col, B);
  
  for (int i = 0; i < m; i++)
    {
      Real_wp val = 0;
      if (norm_cplx)
        for (int j = 0; j < n; j++)
          val += abs(B(i, j));
      else
        for (int j = 0; j < n; j++)
          val += AbsCplx(B(i, j));
      
      if ((abs(sum_row(i) - val) > threshold) || isnan(abs(sum_row(i) - val)))
        {
          cout << "GetRowColSum incorrect" << endl;
          DISP(i); DISP(sum_row(i)); DISP(val);
          abort();
        }
    }
  
  for (int j = 0; j < n; j++)
    {
      Real_wp val = 0;
      if (norm_cplx)
        for (int i = 0; i < m; i++)
          val += abs(B(i, j));
      else
        for (int i = 0; i < m; i++)
          val += AbsCplx(B(i, j));
      
      if ((abs(sum_col(j) - val) > threshold) || isnan(abs(sum_col(j) - val)))
        {
          cout << "GetRowColSum incorrect" << endl;
          DISP(j); DISP(sum_col(j)); DISP(val);
          abort();
        }
    }

  // checking NormFro
  {
    Real_wp val_ex = NormFro(B);
    Real_wp val = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        val += absSquare(B(i, j));
    
    val = sqrt(val);
    if ((abs(val - val_ex) > threshold) || isnan(abs(val - val_ex)))
      {
        cout << "NormFro incorrect" << endl;
        DISP(val); DISP(val_ex);
        abort();
      }
  }
  
  A = B;

  // ScaleMatrix
  GenerateRandomVector(sum_row, m);
  GenerateRandomVector(sum_col, n);
  if (IsSymmetricMatrix(A))
    sum_col = sum_row;

  ScaleMatrix(B, sum_row, sum_col);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(B(i, j) - A(i, j)*sum_row(i)*sum_col(j)) > threshold)
          || isnan(abs(B(i, j) - A(i, j)*sum_row(i)*sum_col(j))))
        {
          cout << "ScaleMatrix incorrect" << endl;
          DISP(i); DISP(j);
          abort();
        }

  if (!IsSymmetricMatrix(A))
    {
      B = A;
      ScaleLeftMatrix(B, sum_row);
      
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          if ((abs(B(i, j) - A(i, j)*sum_row(i)) > threshold) || isnan(abs(B(i, j) - A(i, j)*sum_row(i))))
            {
              cout << "ScaleLeftMatrix incorrect" << endl;
              DISP(i); DISP(j);
              abort();
            }
    }
  
}


template<class MatrixSparse, class MatrixSparseReal>
void CheckCopyReal(const MatrixSparse& A, MatrixSparseReal& B)
{
  CopyReal(A, B);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(B(i, j) - real(A(i, j))) > threshold) || isnan(abs(B(i, j) - real(A(i, j)))))
        {
          cout << "CopyReal incorrect" << endl;
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          abort();
        }
}

template<class MatrixSparse>
void CheckSubMatrix(MatrixSparse& A)
{
  int m = 1, n = 8;
  MatrixSparse B;
  GetSubMatrix(A, m, n, B);
  
  if ( (B.GetM() != n-m) || B.GetN() != n-m)
    {
      cout << "GetSubMatrix incorrect" << endl;
      abort();
    }
  
  for (int i = m; i < n; i++)
    for (int j = m; j < n; j++)
      if ((A(i, j) != B(i-m, j-m)) || isnan(A(i, j)) || isnan(B(i-m, j-m)))
        {
          cout << "GetSubMatrix incorrect" << endl;
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i-m, j-m));
          abort();
        }
  
  int m1 = 2, m2 = 10, n1 = 3, n2 = 7;
  if (IsSymmetricMatrix(A))
    {
      n1 = m1; n2 = m2;
    }
  
  GetSubMatrix(A, m1, m2, n1, n2, B);
  
  if ( (B.GetM() != m2-m1) || B.GetN() != n2-n1)
    {
      cout << "GetSubMatrix incorrect" << endl;
      abort();
    }
  
  for (int i = m1; i < m2; i++)
    for (int j = n1; j < n2; j++)
      if ((A(i, j) != B(i-m1, j-n1)) || isnan(A(i, j)) || isnan(B(i-m1, j-n1)))
        {
          cout << "GetSubMatrix incorrect" << endl;
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i-m1, j-n1));
          abort();
        }  
}


template<class MatrixSparse, class MatrixSparseTarget>
void CheckSubMatrix(MatrixSparse& A, MatrixSparseTarget& B)
{
  A.Clear(); B.Clear();
  
  int m = 40, n = 32, nnz = 300;
  if (IsSymmetricMatrix(A))
    n = m;
  
  GenerateRandomMatrix(A, m, n, nnz);
  
  int nrow = rand()%m + 1;
  IVect row_number_ref(nrow), permut(m);
  GenerateRandomPermutation(m, permut);
  for (int i = 0; i < nrow; i++)
    row_number_ref(i) = permut(i);

  int ncol = rand()%n + 1;
  IVect col_number_ref(ncol);
  GenerateRandomPermutation(n, permut);
  for (int i = 0; i < ncol; i++)
    col_number_ref(i) = permut(i);

  CopySubMatrix(A, row_number_ref, col_number_ref, B);
  
  Vector<bool> RowKept(m), ColKept(n);
  RowKept.Fill(false); ColKept.Fill(false);
  for (int i = 0; i < nrow; i++)
    RowKept(row_number_ref(i)) = true;
  
  for (int i = 0; i < ncol; i++)
    ColKept(col_number_ref(i)) = true;
  
  typename MatrixSparse::entry_type val_exact;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        if (RowKept(i) && ColKept(j))
          val_exact = A(i, j);
        else
          SetComplexZero(val_exact);
        
        if ((B(i, j) != val_exact) || isnan(B(i, j)) || isnan(val_exact))
          {
            cout << "CopySubMatrix incorrect" << endl;
            abort();
          }
      }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e4*epsilon_machine;
  
  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    TestFunctionExtract(A);
    
    CheckSubMatrix(A);
    
    Matrix<Real_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    TestFunctionExtract(A);

    CheckSubMatrix(A);

    Matrix<Real_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    TestFunctionExtract(A);

    CheckSubMatrix(A);

    Matrix<Complex_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    TestFunctionExtract(A);
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> B;
    CheckCopyReal(A, B);

    CheckSubMatrix(A);

    Matrix<Complex_wp, General, RowSparse> C;
    CheckSubMatrix(A, C);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    TestFunctionExtract(A, false);
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> B;
    CheckCopyReal(A, B);

    Matrix<Complex_wp, General, RowSparse> C;
    CheckSubMatrix(A, C);
  }

  {
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    TestFunctionExtract(A, false);

    Matrix<Complex_wp, General, RowSparse> C;
    CheckSubMatrix(A, C);
  }

  {
    Matrix<Real_wp, General, RowSparse> A;
    Matrix<Real_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    Matrix<Real_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Complex_wp, General, RowComplexSparse> A;
    Matrix<Complex_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    Matrix<Complex_wp, General, RowSparse> B;
    CheckSubMatrix(A, B);
  }

  {
    //Matrix<double, General, BlockDiagRow> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<double, Symmetric, BlockDiagRowSym> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<double, Symmetric, DiagonalRow> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<complex<double>, General, BlockDiagRow> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<complex<double>, Symmetric, BlockDiagRowSym> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<complex<double>, General, DiagonalRow> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<double, General, BandedCol> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<double, General, ArrowCol> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<complex<double>, General, BandedCol> A;
    //TestFunctionExtract(A);
  }

  {
    //Matrix<complex<double>, General, ArrowCol> A;
    //TestFunctionExtract(A);
  }
  
  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
