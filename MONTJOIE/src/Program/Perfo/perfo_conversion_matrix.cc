#include "Share/MontjoieCommon.hxx"

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

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  T x;
  Matrix<T, General, ArrayRowSparse> B(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      B.AddInteraction(i, j, x);
    }
  
  Vector<int> IndRow, IndCol;
  Vector<T> Val;
  ConvertMatrix_to_Coordinates(B, IndRow, IndCol, Val);
  A.Clear();
  A.Reallocate(m, n);
  ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, A);
}


template<class T>
void GenerateRandomTriplet(Vector<int>& IndRow, Vector<int>& IndCol, Vector<T>& Val,
			   int m, int n, int nnz, bool sym)
{
  T zero, x;
  SetComplexZero(zero);
  
  IndRow.Reallocate(nnz);
  IndCol.Reallocate(nnz);
  Val.Reallocate(nnz);
  Matrix<T, General, ArrayRowSparse> B(m, n);
  int nb = 0;
  while (nb < nnz)
    {
      int i = rand()%m;
      int j = rand()%n;
      if (sym)
	{
	  if (i > j)
	    {
	      int k = i;
	      i = j;
	      j = k;
	    }
	}
      
      GetRandNumber(x);
      if (B(i, j) == zero)
	{
	  IndRow(nb) = i;
	  IndCol(nb) = j;
	  Val(nb) = x;
	  B.AddInteraction(i, j, x);
	  nb++;
	}
    }
}


template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2>
bool EqualMatrix(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
                 const Matrix<T2, Prop2, Storage2, Allocator2>& B)
{
  if ( (A.GetM() != B.GetM())  || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          return false;
        }
  
  return true;
}

template<class T1, class Prop1, class Storage1>
void TestConversionAlone(Matrix<T1, Prop1, Storage1>& A,
			 int m, int nb_iter, bool sym_pattern = true)
{
  int n = m;
  int nnz = 10*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  Vector<int> IndRow, IndCol, Ptr;
  Vector<T1> Val;
  for (int i = 0; i < nb_iter; i++)
    {      
      ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    }

  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertMatrix_to_Coordinates = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
  
  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      General prop;
      ConvertToCSC(A, prop, Ptr, IndRow, Val, sym_pattern);
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertToCSC = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;


  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      General prop;
      ConvertToCSR(A, prop, Ptr, IndRow, Val);
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertToCSR = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
  
  bool sym = IsSymmetricMatrix(A);
  GenerateRandomTriplet(IndRow, IndCol, Val, m, n, nnz, sym);

  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      Vector<int> IndRow_(IndRow), IndCol_(IndCol);
      Vector<T1> Val_(Val);
      ConvertMatrix_from_Coordinates(IndRow_, IndCol_, Val_, A, 0);
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertMatrix_from_Coordinates = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
    
}

template<class T1, class Prop1, class Storage1>
void TestConversionSymmetric(Matrix<T1, Prop1, Storage1>& A, int m, int nb_iter)
{
  int n = m;
  int nnz = 10*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  Vector<int> IndRow, Ptr;
  Vector<T1> Val;
  Symmetric prop;

  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      ConvertToCSC(A, prop, Ptr, IndRow, Val);
    }

  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertToCSC (symmetric) = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;

  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      ConvertToCSR(A, prop, Ptr, IndRow, Val);
    }

  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for ConvertToCSR (symmetric) = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
}

template<class T1, class Prop1, class Storage1,
	 class T2, class Prop2, class Storage2>
void TestConversion(Matrix<T1, Prop1, Storage1>& A, Matrix<T2, Prop2, Storage2>& B,
		    int m, int nb_iter)
{
  int n = m;
  int nnz = 10*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);

  for (int i = 0; i < nb_iter; i++)
    {      
      Copy(A, B);
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Time for Copy = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
  
}


template<class T1, class Prop1, class Storage1>
void CheckConversionAlone(Matrix<T1, Prop1, Storage1>& A, bool sorted_by_rows,
			  bool sym_pattern = true)
{
  int m = 40, n = 30;
  if (IsSymmetricMatrix(A))
    n = m;
  
  int nnz = 4*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  Vector<int> IndRow, IndCol, Ptr;
  Vector<T1> Val;
  
  ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 1, true);
  
  if (sorted_by_rows)
    {
      for (int i = 0; i < IndRow.GetM()-1; i++)
	{
	  if (IndRow(i+1) < IndRow(i))
	    {
	      cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	      DISP(i); DISP(IndRow(i)); DISP(IndRow(i+1));
	      abort();
	    }
	  else if (IndRow(i+1) == IndRow(i))
	    {
	      if (IndCol(i+1) <= IndCol(i))
		{
		  cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
		  DISP(i); DISP(IndCol(i)); DISP(IndCol(i+1));
		  abort();
		}	      
	    }
	}
    }
  else
    {
      for (int i = 0; i < IndRow.GetM()-1; i++)
	{
	  if (IndCol(i+1) < IndCol(i))
	    {
	      cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	      DISP(i); DISP(IndCol(i)); DISP(IndCol(i+1));
	      abort();
	    }
	  else if (IndCol(i+1) == IndCol(i))
	    {
	      if (IndRow(i+1) <= IndRow(i))
		{
		  cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
		  DISP(i); DISP(IndRow(i)); DISP(IndRow(i+1));
		  abort();
		}	      
	    }
	}
    }

  //DISP(IndRow); DISP(IndCol); DISP(Val);
  Matrix<T1, General, ArrayRowSparse> C(m, n);
  for (int i = 0; i < IndRow.GetM(); i++)
    C.AddInteraction(IndRow(i)-1, IndCol(i)-1, Val(i));
      
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
      abort();
    }

  // testing ConvertToCSC with m = n and symmetrization of pattern
  GenerateRandomMatrix(A, m, n, nnz);

  General prop;
  ConvertToCSC(A, prop, Ptr, IndRow, Val, sym_pattern);
  
  C.Clear();
  C.Reallocate(m, n);
  for (int i = 0; i < n; i++)
    for (int j = Ptr(i); j < Ptr(i+1); j++)
      {
	if (j < Ptr(i+1)-1)
	  {
	    if (IndRow(j+1) < IndRow(j))
	      {
		cout << "ConvertToCSC incorrect " << endl;
		abort();
	      }
	  }

	C.AddInteraction(IndRow(j), i, Val(j));
      }
  
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertToCSC incorrect" << endl;
      abort();
    }

  if (sym_pattern)
    {
      C.Fill(2.0);
      for (int i = 0; i < n; i++)
	for (int j = i; j < n; j++)
	  if (C(i, j) != C(j, i))
	    {
	      cout << "ConvertToCSC incorrect" << endl;
	      abort();
	    }
    }

  // testing ConvertMatrix_from_Coordinates
  bool sym = IsSymmetricMatrix(A);
  GenerateRandomTriplet(IndRow, IndCol, Val, m, n, nnz, sym);
  for (int i = 0; i < IndRow.GetM(); i++)
    {
      IndRow(i)++;
      IndCol(i)++;
    }

  C.Clear();
  C.Reallocate(m, n);
  for (int i = 0; i < IndRow.GetM(); i++)
    {
      C.AddInteraction(IndRow(i)-1, IndCol(i)-1, Val(i));
      if (sym && (IndRow(i) != IndCol(i)))
	C.AddInteraction(IndCol(i)-1, IndRow(i)-1, Val(i));
    }

  A.Reallocate(m, n);
  ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, A, 1);
  //cout << "A = " << endl; A.WriteText(cout); cout << endl;
  
  //cout << "C = " << endl; C.WriteText(cout); cout << endl;
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertMatrix_from_Coordinates incorrect" << endl;
      abort();
    }

  // testing ConvertToCSR
  GenerateRandomMatrix(A, m, n, nnz);

  ConvertToCSR(A, prop, Ptr, IndRow, Val);
  
  C.Clear();
  C.Reallocate(m, n);
  for (int i = 0; i < m; i++)
    for (int j = Ptr(i); j < Ptr(i+1); j++)
      {
	if (j < Ptr(i+1)-1)
	  {
	    if (IndRow(j+1) < IndRow(j))
	      {
		cout << "ConvertToCSR incorrect " << endl;
		abort();
	      }
	  }

	C.AddInteraction(i, IndRow(j), Val(j));
      }
  
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertToCSR incorrect" << endl;
      abort();
    }

}

template<class T1, class Prop1, class Storage1,
	 class T2, class Prop2, class Storage2>
void CheckConversion(Matrix<T1, Prop1, Storage1>& A,
		     Matrix<T2, Prop2, Storage2>& B)
{
  int m = 40, n = 30;
  if (IsSymmetricMatrix(A))
    n = m;
  
  int nnz = 4*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  Copy(A, B);

  if (!EqualMatrix(A, B))
    {
      cout << "Copy incorrect" << endl;
      abort();
    }
}

template<class T1, class Prop1, class Storage1>
void CheckConversionSymmetric(Matrix<T1, Prop1, Storage1>& A)
{
  int m = 10, n = m;
  int nnz = 4*n;
  GenerateRandomMatrix(A, m, n, nnz);
  
  Vector<int> IndRow, Ptr;
  Vector<T1> Val;

  Symmetric prop;
  ConvertToCSC(A, prop, Ptr, IndRow, Val);
  
  Matrix<T1, Symmetric, ArrayRowSymSparse> C;
  
  C.Clear();
  C.Reallocate(m, n);
  for (int i = 0; i < n; i++)
    for (int j = Ptr(i); j < Ptr(i+1); j++)
      {
	if (j < Ptr(i+1)-1)
	  {
	    if (IndRow(j+1) < IndRow(j))
	      {
		cout << "ConvertToCSC incorrect " << endl;
		abort();
	      }
	  }

	C.AddInteraction(IndRow(j), i, Val(j));
      }
  
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertToCSC incorrect" << endl;
      abort();
    }

  // testing ConvertToCSR
  ConvertToCSR(A, prop, Ptr, IndRow, Val);
  
  C.Clear();
  C.Reallocate(m, n);
  for (int i = 0; i < n; i++)
    for (int j = Ptr(i); j < Ptr(i+1); j++)
      {
	if (j < Ptr(i+1)-1)
	  {
	    if (IndRow(j+1) < IndRow(j))
	      {
		cout << "ConvertToCSR incorrect " << endl;
		abort();
	      }
	  }

	C.AddInteraction(i, IndRow(j), Val(j));
      }
  
  if (!EqualMatrix(A, C))
    {
      cout << "ConvertToCSR incorrect" << endl;
      abort();
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Entrez la taille de la matrice et le nombre d'iterations" << endl;
      abort();
    }

  threshold = epsilon_machine;
  int m = atoi(argv[1]);
  int nb_iter = atoi(argv[2]);
  
  {
    cout << "Testing RowSparse matrices" << endl;
    Matrix<Real_wp, General, RowSparse> A;
    CheckConversionAlone(A, true);
    TestConversionAlone(A, m, nb_iter);

    Matrix<Real_wp, General, ColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, General, ArrayColSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, General, ArrayRowSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ColSymSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> F;
    CheckConversion(F, A);
    TestConversion(F, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, RowSymSparse> G;
    CheckConversion(G, A);
    TestConversion(G, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> H;
    CheckConversion(H, A);
    TestConversion(H, A, m, nb_iter);

  }
  
  {
    cout << "Testing ColSparse matrices" << endl;
    Matrix<Real_wp, General, ColSparse> A;
    CheckConversionAlone(A, false);
    TestConversionAlone(A, m, nb_iter);

    Matrix<Real_wp, General, ArrayColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, General, RowSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, General, ArrayRowSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, RowSymSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> F;
    CheckConversion(F, A);
    TestConversion(F, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ColSymSparse> G;
    CheckConversion(G, A);
    TestConversion(G, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> H;
    CheckConversion(H, A);
    TestConversion(H, A, m, nb_iter);

  }

  {
    cout << "Testing RowSymSparse matrices" << endl;
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, true);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    //Matrix<Real_wp, General, ArrayRowSparse> B;
    //CheckConversion(B, A);
    //TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ColSymSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);
  }

  {
    cout << "Testing ColSymSparse matrices" << endl;
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Real_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayRowSparse matrices" << endl;
    Matrix<Real_wp, General, ArrayRowSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, true);

    Matrix<Real_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, General, RowSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, General, ColSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);    

    Matrix<Real_wp, General, ArrayColSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);    

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> F;
    CheckConversion(F, A);
    TestConversion(F, A, m, nb_iter);    

  }

  {
    cout << "Testing ArrayColSparse matrices" << endl;
    Matrix<Real_wp, General, ArrayColSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, false);

    Matrix<Real_wp, General, RowSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, General, ColSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);    

    Matrix<Real_wp, General, ArrayRowSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);    

  }

  {
    cout << "Testing ArrayRowSymSparse matrices" << endl;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, true);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Real_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ColSymSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

  }

  {
    cout << "Testing ArrayColSymSparse matrices" << endl;
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    TestConversionAlone(A, m, nb_iter);
    CheckConversionAlone(A, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Real_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

  }


  {
    cout << "Testing RowComplexSparse matrices" << endl;
    Matrix<Complex_wp, General, RowComplexSparse> A;
    CheckConversionAlone(A, true, false);
    TestConversionAlone(A, m, nb_iter, false);

    Matrix<Complex_wp, General, ColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayRowComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

  }

  {
    cout << "Testing ColComplexSparse matrices" << endl;
    Matrix<Complex_wp, General, ColComplexSparse> A;
    CheckConversionAlone(A, false, false);
    TestConversionAlone(A, m, nb_iter, false);

    Matrix<Complex_wp, General, ColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing RowSymComplexSparse matrices" << endl;
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    CheckConversionAlone(A, true, false);
    TestConversionAlone(A, m, nb_iter, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Complex_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);
  }

  {
    cout << "Testing ColSymComplexSparse matrices" << endl;
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    CheckConversionAlone(A, false, false);
    TestConversionAlone(A, m, nb_iter, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Complex_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayRowComplexSparse matrices" << endl;
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    CheckConversionAlone(A, true, false);
    TestConversionAlone(A, m, nb_iter, false);

    Matrix<Complex_wp, General, ColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayColComplexSparse matrices" << endl;
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    CheckConversionAlone(A, false, false);
    TestConversionAlone(A, m, nb_iter, false);

    Matrix<Complex_wp, General, ColSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayRowSymComplexSparse matrices" << endl;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    CheckConversionAlone(A, true, false);
    TestConversionAlone(A, m, nb_iter, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Complex_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayColSymComplexSparse matrices" << endl;
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    CheckConversionAlone(A, false, false);
    TestConversionAlone(A, m, nb_iter, false);
    TestConversionSymmetric(A, m, nb_iter);
    CheckConversionSymmetric(A);

    Matrix<Complex_wp, Symmetric, RowSymSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);
  }

  {
    cout << "Testing ArrayRowSparse conversion" << endl;
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayRowComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Complex_wp, General, RowComplexSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

  }

  {
    cout << "Testing ArrayRowSymSparse conversion" << endl;
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);
  }

  {
    cout << "Testing RowSparse conversion" << endl;
    Matrix<Complex_wp, General, RowSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayRowComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayColComplexSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> F;
    CheckConversion(F, A);
    TestConversion(F, A, m, nb_iter);

    Matrix<Complex_wp, General, RowComplexSparse> G;
    CheckConversion(G, A);
    TestConversion(G, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> H;
    CheckConversion(H, A);
    TestConversion(H, A, m, nb_iter);

    Matrix<Complex_wp, General, ColComplexSparse> J;
    CheckConversion(J, A);
    TestConversion(J, A, m, nb_iter);

  }

  {
    cout << "Testing RowSymSparse conversion" << endl;
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

  }


  {
    cout << "Testing ColSparse conversion" << endl;
    Matrix<Complex_wp, General, ColSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayRowComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Complex_wp, General, ArrayColComplexSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> F;
    CheckConversion(F, A);
    TestConversion(F, A, m, nb_iter);

    Matrix<Complex_wp, General, RowComplexSparse> G;
    CheckConversion(G, A);
    TestConversion(G, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> H;
    CheckConversion(H, A);
    TestConversion(H, A, m, nb_iter);

    Matrix<Complex_wp, General, ColComplexSparse> J;
    CheckConversion(J, A);
    TestConversion(J, A, m, nb_iter);

  }

  {
    cout << "Testing ColSymSparse conversion" << endl;
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckConversion(B, A);
    TestConversion(B, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> C;
    CheckConversion(C, A);
    TestConversion(C, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> D;
    CheckConversion(D, A);
    TestConversion(D, A, m, nb_iter);

    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> E;
    CheckConversion(E, A);
    TestConversion(E, A, m, nb_iter);

  }
  
  return FinalizeMontjoie();
}
