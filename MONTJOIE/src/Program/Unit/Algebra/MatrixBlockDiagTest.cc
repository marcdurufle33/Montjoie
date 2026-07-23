#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-14;

template<class T>
void GetRand(T & x)
{
  x = 2.0*T(rand())/RAND_MAX - 1.0;
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
      if (IsSymmetricMatrix(A))
        A.Set(j, i, x);
    }
}

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

template<class T, class Prop1, class Storage1, class T2, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T2, Prop2, Storage2>& B, Real_wp epsilon = threshold)
{
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if (( abs(A(i, j) - B(i, j) ) > epsilon) || isnan(abs(A(i, j) - B(i, j))))
	{
          DISP(A(i, j)); DISP(B(i, j));
          //abort();
          return false;
        }
  
  return true;
}

template<class T0, class T, class Prop1, class Storage1, class Prop2, class Storage2, class Prop3, class Storage3>
bool EqualMatrix(const T0& alpha,
                 const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T, Prop2, Storage2>& B,
                 const Matrix<T, Prop3, Storage3>& C)
{
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;

  if ( (A.GetM() != C.GetM()) || (A.GetN() != C.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ( abs(C(i, j) - alpha*A(i, j) - B(i, j) ) > threshold)
	return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void PerfBlockDiagMatrix(Matrix<T, Prop, Storage, Allocator>& A,
			 int size_block, int size_mat, int nb_iter)
{
  // nombre de blocks
  int nb_blocks = size_mat/size_block;
  size_mat = nb_blocks*size_block;
  
  IVect permut;
  GenerateRandomPermutation(size_mat, permut);
  
  // initializing sparsity pattern of the matrix
  Vector<IVect> BlockRowNumber(nb_blocks);
  int nb = 0;
  for (int i = 0; i < nb_blocks; i++)
    {
      BlockRowNumber(i).Reallocate(size_block);
      for (int j = 0; j < size_block; j++)
	BlockRowNumber(i)(j) = permut(nb++);
    }
  
  A.SetPattern(BlockRowNumber);
  A.FillRand();
  Mlt(1.0/RAND_MAX, A);
  
  //int jeton_tps = glob_chrono.GetNumber();
  //glob_chrono.Start(jeton_tps);
  GetInverse(A);
  //glob_chrono.Stop(jeton_tps);
  //cout<<"Inversion took "<<glob_chrono.GetSeconds(jeton_tps)<<" s"<<endl;
  //glob_chrono.Reset(jeton_tps);
    
  // testing fastness of matrix vector product
  Vector<T> x, y;
  T alpha, beta; 
  GenerateRandomVector(x, size_mat);
  //glob_chrono.Start(jeton_tps);
  for (int k = 0; k < nb_iter; k++)
    {
      GetRand(alpha); GetRand(beta);
      GenerateRandomVector(y, size_mat);
      //DISP(Norm2(x)); DISP(Norm2(y)); DISP(alpha); DISP(beta);
      MltAdd(alpha, A, x, beta, y);      
    }
  //glob_chrono.Stop(jeton_tps);
  //cout<<"MltAdd took "<<glob_chrono.GetSeconds(jeton_tps)<<" s"<<endl;

  //glob_chrono.Reset(jeton_tps);
  //glob_chrono.Start(jeton_tps);
  for (int k = 0; k < nb_iter; k++)
    {
      GenerateRandomVector(y, size_mat);
      //DISP(Norm2(x)); DISP(Norm2(y)); DISP(alpha); DISP(beta);
      Mlt(A, y);      
    }
  //glob_chrono.Stop(jeton_tps);
  //cout<<"Mlt took "<<glob_chrono.GetSeconds(jeton_tps)<<" s"<<endl;

  DISP(Norm2(y));
}

template<class T, class Allocator>
int TestBlockDiagMatrix(Matrix<T, General, BlockDiagRow, Allocator>& A)
{  
  // for example a matrix with three blocks [4], [0, 2], [1, 3, 5]
  Vector<IVect> BlockNum(3);
  BlockNum(0).Reallocate(1); BlockNum(0)(0) = 4;
  BlockNum(1).Reallocate(2); BlockNum(1)(0) = 0; BlockNum(1)(1) = 2;
  BlockNum(2).Reallocate(3); BlockNum(2)(0) = 1; BlockNum(2)(1) = 3; BlockNum(2)(2) = 5;
  
  // A is initialized with the structure you gave, all values are set to 0
  A.SetPattern(BlockNum);
  
  // checking that pattern is correct
  if ((A.GetM() != 6) || (A.GetN() != 6))
    {
      cout << "Number of rows or columns incorrect " << endl;
      abort();
    }
  
  if (A.GetNbBlocks() != 3)
    {
      cout << "Number of blocks incorrect " << endl;
      abort();
    }
  
  if (A.GetDataSize() != 14)
    {
      cout << "Number of non-zero entries is incorrect " << endl;
      abort();
    }
  
  if ((A.GetBlockSize(0) != 1) || (A.GetBlockSize(1) != 2) || (A.GetBlockSize(2) != 3))
    {
      cout << "Size of block incorrect " << endl;
      abort();
    }
  
  if ((A.Index(0, 0) != 4) || (A.Index(1, 0) != 0) || (A.Index(1, 1) != 2) 
      || (A.Index(2, 0) != 1) || (A.Index(2, 1) != 3) || (A.Index(2, 2) != 5) )
    {
      cout << "Index incorrect " << endl;
      abort();
    }

  // now if you want to add values in this structure, you can use AddInteraction/AddInteractionRow
  A.AddInteraction(0, 0, 1.5); A.AddInteraction(0, 2, 2.5);
  A.AddInteraction(2, 0, 2.5); A.AddInteraction(2, 2, 2.0);
  A.AddInteraction(4, 4, -3.5);
  
  // several entries with AddInteractionRow
  IVect col(3); col(0) = 1; col(1) = 3; col(2) = 5;
  Vector<T> value(3); value(0) = -1.5; value(1) = 0.5; value(2) = 1.0;
  A.AddInteractionRow(1, col.GetM(), col, value);
  value(0) = 0.5; value(1) = 4.5; value(2) = -2.5;
  A.AddInteractionRow(3, col.GetM(), col, value);
  value(0) = 1.0; value(1) = -2.5; value(2) = 3.0;
  A.AddInteractionRow(5, col.GetM(), col, value);
  
  // A should be equal to |1.5   0   2.5   0   0   0  |
  //                      | 0  -1.5   0   0.5  0  1.0 |
  //                      |2.5   0   2.0   0   0   0  |
  //                      | 0   0.5   0   4.5  0 -2.5 |
  //                      | 0    0    0    0 -3.5  0  |
  //                      | 0   1.0   0  -2.5  0  3.0 |
  Matrix<T, General, RowMajor> Adense, Aref;
  Matrix<T, General, BlockDiagRow, Allocator> B(A);
  Aref.Reallocate(6, 6);
  Aref.Fill(0);
  Aref(0, 0) = 1.5; Aref(0, 2) = 2.5;
  Aref(2, 0) = 2.5; Aref(2, 2) = 2.0;
  Aref(1, 1) = -1.5; Aref(1, 3) = 0.5; Aref(1, 5) = 1.0;
  Aref(3, 1) = 0.5; Aref(3, 3) = 4.5; Aref(3, 5) = -2.5;
  Aref(5, 1) = 1.0; Aref(5, 3) = -2.5; Aref(5, 5) = 3.0;
  Aref(4, 4) = -3.5;
  
  Copy(A, Adense);
  if (!EqualMatrix(A, Aref))
    {
      cout << "Operator() incorrect " << endl;
      abort();
    }
  
  if (!EqualMatrix(A, Adense))
    {
      cout << "Conversion to dense matrix incorrect " << endl;
      abort();
    }
  
  if (!EqualMatrix(Aref, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < A.GetNbBlocks(); i++)
    {
      int p = A.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          {
            int j1 = A.Index(i, j);
            int k1 = A.Index(i, k);
            if ((abs(A(j1, k1) - A.Value(i, j, k)) > threshold)
                || isnan(abs(A(j1, k1) - A.Value(i, j, k))))
              {
                cout << "Value incorrect" << endl;
                abort();
              }
          }
    }
  
  B.Clear();
  
  if ((B.GetM() != 0) || (B.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // you can copy block-diagonal matrices
  B = A;
  
  // you can also set all non-zero entries to the same value
  B.Fill(2.0);
  // B should be equal to [2 0 2 0 0 0; 0 2 0 2 0 2; 2 0 2 0 0 0; 0 2 0 2 0 2; 0 0 0 0 2 0; 0 2 0 2 0 2]
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int p = B.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          if ((abs(B.Value(i, j, k) - 2.0) > threshold) || isnan(abs(B.Value(i, j, k) - 2.0)))
            {
              cout << "Fill incorrect" << endl;
              abort();
            }
    }

  // values of non-zero entries can be set to 0
  B.Zero();
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int p = B.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          if ((abs(B.Value(i, j, k)) > threshold) || isnan(abs(B.Value(i, j, k))))
            {
              cout << "Zero incorrect" << endl;
              abort();
            }
    }

  // or randomly set
  B.FillRand();
    
  // you can write A in ascii format (same format as sparse matrices)
  A.WriteText("mat_block.dat");
  
  // checking if it is correctly written
  Matrix<T, General, ArrayRowSparse> Acsr;
  Acsr.ReadText("mat_block.dat");
  if (!EqualMatrix(A, Acsr))
    {
      cout << "WriteText incorrect " << endl;
      abort();
    }
  
  // checking conversions
  Acsr.Clear();
  ConvertToSparse(A, Acsr);
  if (!EqualMatrix(A, Acsr))
    {
      cout << "ConvertToSparse incorrect " << endl;
      abort();
    }
  
  A.Clear();
  ConvertToBlockDiagonal(Acsr, A);
  if (!EqualMatrix(A, Acsr))
    {
      cout << "ConvertToBlockDiagonal incorrect " << endl;
      abort();
    }

  B = A;
  B.ClearRow(2);

  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetN(); j++)
      {
	T val_ref = A(i, j);
	if (i == 2)
	  SetComplexZero(val_ref);
	
	if (B(i, j) != val_ref)
	  {
	    cout << "ClearRow incorrect" << endl;
	    abort();
	  }
      }

  B = A;
  B.ClearColumn(2);

  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetN(); j++)
      {
	T val_ref = A(i, j);
	if (j == 2)
	  SetComplexZero(val_ref);
	
	if (B(i, j) != val_ref)
	  {
	    cout << "ClearColumn incorrect" << endl;
	    abort();
	  }
      }

  // you can combine matrix with Add
  // B = B + alpha*A  
  B.FillRand(); B *= 1e-9;
  Matrix<T, General, BlockDiagRow> C(B);
  Add(Real_wp(0.45), A, B);
  if (!EqualMatrix(Real_wp(0.45), A, C, B))
    {
      cout << "Add incorrect " << endl;
      abort();
    }
  
  // you can use MltAdd to perform a matrix-vector product
  Vector<T> X(6), Y(6), Yref(6);
  Vector<T> Xref;
  Y.Fill(0); X.FillRand(); X *= 1e-9; Yref.Fill(0);
  MltAdd(Real_wp(3)/Real_wp(10), A, X, Real_wp(0.0), Y);
  MltAdd(Real_wp(3)/Real_wp(10), Acsr, X, Real_wp(0.0), Yref);
  for (int i = 0; i < 6; i++)
    if ((abs(Yref(i) - Y(i)) > threshold) || isnan(abs(Yref(i) - Y(i))))
      {
        cout << "MltAdd incorrect " << endl;
        abort();
      }
  
  // you can replace A by its inverse, and then find the solution of A X = B
  Xref = X;
  GetInverse(A);
  X.Fill(0);
  MltAdd(Real_wp(10)/Real_wp(3), A, Y, Real_wp(0), X);
  // solution should be [0;1;2;3;4;5]
  for (int i = 0; i < 6; i++)
    if ((abs(X(i) - Xref(i)) > threshold) || isnan(abs(X(i) - Xref(i))))
      {
        DISP(X(i)); DISP(Xref(i)); DISP(abs(X(i) - Xref(i)));
        cout << "GetInverse incorrect " << endl;
        abort();
      }  
  
  return 0;
}

template<class T, class Allocator>
int TestBlockDiagMatrix(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
{
  // for example a matrix with three blocks [4] [0, 2], [1, 3, 5]
  Vector<IVect> BlockNum(3);
  BlockNum(0).Reallocate(1); BlockNum(0)(0) = 4;
  BlockNum(1).Reallocate(2); BlockNum(1)(0) = 0; BlockNum(1)(1) = 2;
  BlockNum(2).Reallocate(3); BlockNum(2)(0) = 1; BlockNum(2)(1) = 3; BlockNum(2)(2) = 5;
  
  // A is initialized with the structure you gave, all values are set to 0
  A.SetPattern(BlockNum);
  
  // checking that pattern is correct
  if ((A.GetM() != 6) || (A.GetN() != 6))
    {
      cout << "Number of rows or columns incorrect " << endl;
      abort();
    }
  
  if (A.GetNbBlocks() != 3)
    {
      cout << "Number of blocks incorrect " << endl;
      abort();
    }
  
  if (A.GetDataSize() != 10)
    {
      cout << "Number of non-zero entries is incorrect " << endl;
      abort();
    }
  
  if ((A.GetBlockSize(0) != 1) || (A.GetBlockSize(1) != 2) || (A.GetBlockSize(2) != 3))
    {
      cout << "Size of block incorrect " << endl;
      abort();
    }
  
  if ((A.Index(0, 0) != 4) || (A.Index(1, 0) != 0) || (A.Index(1, 1) != 2) 
      || (A.Index(2, 0) != 1) || (A.Index(2, 1) != 3) || (A.Index(2, 2) != 5) )
    {
      cout << "Index incorrect " << endl;
      abort();
    }

  // now if you want to add values in this structure, you can use AddInteraction/AddInteractionRow
  A.AddInteraction(0, 0, 4.3); A.AddInteraction(0, 2, 2.5);
  A.AddInteraction(2, 0, 2.5); A.AddInteraction(2, 2, 2.0);
  A.AddInteraction(4, 4, 3.5);
  
  // several entries with AddInteractionRow
  IVect col(3); col(0) = 1; col(1) = 3; col(2) = 5;
  Vector<T> value(3); value(0) = 2.8; value(1) = 0.5; value(2) = 1.0;
  A.AddInteractionRow(1, col.GetM(), col, value);
  value(0) = 0.5; value(1) = 4.5; value(2) = -2.5;
  A.AddInteractionRow(3, col.GetM(), col, value);
  value(0) = 1.0; value(1) = -2.5; value(2) = 3.0;
  A.AddInteractionRow(5, col.GetM(), col, value);
  
  // A should be equal to |4.3   0   2.5   0   0   0  |
  //                      | 0   2.8   0   0.5  0  1.0 |
  //                      |2.5   0   2.0   0   0   0  |
  //                      | 0   0.5   0   4.5  0 -2.5 |
  //                      | 0    0    0    0 3.5  0  |
  //                      | 0   1.0   0  -2.5  0  3.0 |
  Matrix<T, Symmetric, RowSymPacked> Adense, Aref;
  Matrix<T, Symmetric, BlockDiagRowSym, Allocator> B(A);
  Aref.Reallocate(6, 6);
  Aref.Fill(0);
  Aref(0, 0) = 4.3; Aref(0, 2) = 2.5;
  Aref(2, 0) = 2.5; Aref(2, 2) = 2.0;
  Aref(1, 1) = 2.8; Aref(1, 3) = 0.5; Aref(1, 5) = 1.0;
  Aref(3, 1) = 0.5; Aref(3, 3) = 4.5; Aref(3, 5) = -2.5;
  Aref(5, 1) = 1.0; Aref(5, 3) = -2.5; Aref(5, 5) = 3.0;
  Aref(4, 4) = 3.5;
  
  Copy(A, Adense);
  if (!EqualMatrix(Aref, A))
    {
      cout << "Operator() incorrect " << endl;
      abort();
    }

  if (!EqualMatrix(Aref, Adense))
    {
      cout << "Conversion to dense matrix incorrect " << endl;
      abort();
    }
  
  if (!EqualMatrix(Aref, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < A.GetNbBlocks(); i++)
    {
      int p = A.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          {
            int j1 = A.Index(i, j);
            int k1 = A.Index(i, k);
            if ((abs(A(j1, k1) - A.Value(i, j, k)) > threshold) || isnan(abs(A(j1, k1) - A.Value(i, j, k))))
              {
                cout << "Value incorrect" << endl;
                abort();
              }
          }
    }
  
  B.Clear();
  
  if ((B.GetM() != 0) || (B.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // you can copy block-diagonal matrices
  B = A;
  
  // you can also set all non-zero entries to the same value
  B.Fill(2.0);
  // B should be equal to [2 0 2 0 0 0; 0 2 0 2 0 2; 2 0 2 0 0 0; 0 2 0 2 0 2; 0 0 0 0 2 0; 0 2 0 2 0 2]
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int p = B.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          if ((abs(B.Value(i, j, k) - 2.0) > threshold) || isnan(abs(B.Value(i, j, k) - 2.0)))
            {
              cout << "Fill incorrect" << endl;
              abort();
            }
    }

  // values of non-zero entries can be set to 0
  B.Zero();
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int p = B.GetBlockSize(i);
      for (int j = 0; j < p; j++)
        for (int k = 0; k < p; k++)
          if ((abs(B.Value(i, j, k)) > threshold) || isnan(abs(B.Value(i, j, k))))
            {
              cout << "Zero incorrect" << endl;
              abort();
            }
    }

  // or randomly set
  B.FillRand();

  B = A;
  B.ClearRow(5);

  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetN(); j++)
      {
	T val_ref = A(i, j);
	if ((i == 5) || (j == 5))
	  SetComplexZero(val_ref);
	
	if (B(i, j) != val_ref)
	  {
	    cout << "ClearRow incorrect" << endl;
	    abort();
	  }
      }

  B = A;
  B.ClearColumn(5);

  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetN(); j++)
      {
	T val_ref = A(i, j);
	if ((i == 5) || (j == 5))
	  SetComplexZero(val_ref);
	
	if (B(i, j) != val_ref)
	  {
	    cout << "ClearColumn incorrect" << endl;
	    abort();
	  }
      }
    
  // you can write A in ascii format (same format as sparse matrices)
  A.WriteText("mat_block.dat");
  
  // checking if it is correctly written
  Matrix<T, Symmetric, ArrayRowSymSparse> Acsr;
  Acsr.ReadText("mat_block.dat");
  if (!EqualMatrix(Acsr, A))
    {
      cout << "WriteText incorrect " << endl;
      abort();
    }
  
  // checking conversions
  Acsr.Clear();
  ConvertToSparse(A, Acsr);
  if (!EqualMatrix(Acsr, A))
    {
      cout << "ConvertToSparse incorrect " << endl;
      abort();
    }
  
  A.Clear();
  ConvertToBlockDiagonal(Acsr, A);
  if (!EqualMatrix(Acsr, A))
    {
      cout << "ConvertToBlockDiagonal incorrect " << endl;
      abort();
    }
  
  // you can combine matrix with Add
  // B = B + alpha*A  
  B.FillRand(); B *= 1e-9;
  Matrix<T, Symmetric, BlockDiagRowSym> C(B);
  Add(Real_wp(0.45), A, B);
  if (!EqualMatrix(Real_wp(0.45), A, C, B))
    {
      cout << "Add incorrect " << endl;
      abort();
    }
  
  // you can use MltAdd to perform a matrix-vector product
  Vector<T> X(6), Y(6), Yref(6), Xref;
  Y.Fill(0); X.FillRand(); X *= 1e-9; Yref.Fill(0);
  MltAdd(Real_wp(3)/Real_wp(10), A, X, Real_wp(0.0), Y);
  MltAdd(Real_wp(3)/Real_wp(10), Acsr, X, Real_wp(0.0), Yref);
  for (int i = 0; i < 6; i++)
    if ((abs(Yref(i) - Y(i)) > threshold) || isnan(abs(Yref(i) - Y(i))))
      {
        cout << "MltAdd incorrect " << endl;
        abort();
      }
  
  C = A;
  
  // you can replace A by its inverse, and then find the solution of A X = B
  Xref = X;
  GetInverse(A);
  X.Fill(0);
  MltAdd(Real_wp(10)/Real_wp(3), A, Y, Real_wp(0.0), X);
  // solution should be Xref
  for (int i = 0; i < 6; i++)
    if ((abs(X(i) - Xref(i)) > threshold) || isnan(abs(X(i) - Xref(i))))
      {
        cout << "GetInverse incorrect " << endl;
        abort();
      }  
  
  // you can use Cholesky factorisation too
  GetCholesky(C);
  X = Y; Mlt(Real_wp(10)/Real_wp(3), X);
  SolveCholesky(SeldonNoTrans, C, X);
  SolveCholesky(SeldonTrans, C, X);
  for (int i = 0; i < 6; i++)
    if ((abs(X(i) - Xref(i)) > threshold) || isnan(abs(X(i) - Xref(i))))
      {
        cout << "GetCholesky/SolveCholesky incorrect " << endl;
        abort();
      }  
  
  MltCholesky(SeldonTrans, C, X);
  MltCholesky(SeldonNoTrans, C, X);
  Mlt(Real_wp(3)/Real_wp(10), X);
  for (int i = 0; i < 6; i++)
    if ((abs(X(i) - Y(i)) > threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "MltCholesky incorrect " << endl;
        abort();
      }  
  
  return 0;
}


template<class T, class Prop, class Storage, class Allocator,
         class T2, class Prop2, class Storage2, class Allocator2>
void CheckSizeMat(Matrix<T, Prop, Storage, Allocator>& A,
                  Matrix<T2, Prop2, Storage2, Allocator2>& B)
{
  int size_max = 0;
  
  for (int i = 0; i < A.GetM(); i++)
    {
      if (!IsSymmetricMatrix(A))
        if (A.GetRowSize(i) != B.GetRowSize(i))
          {
            cout << "GetRowSize incorrect" << endl;
            abort();
          }
      
      size_max = max(size_max, A.GetRowSize(i));
    }
  
  if (size_max != A.GetSizeMaxBlock())
    {
      cout << "GetSizeMaxBlock incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Storage, class Allocator,
         class T2, class Prop2, class Allocator2>
void CheckSizeMat(Matrix<T, Prop, Storage, Allocator>& A,
                  Matrix<T2, Prop2, ArrayRowComplexSparse, Allocator2>& B)
{}

template<class T, class Prop, class Storage, class Allocator,
         class T2, class Prop2, class Allocator2>
void CheckSizeMat(Matrix<T, Prop, Storage, Allocator>& A,
                  Matrix<T2, Prop2, ArrayRowSymComplexSparse, Allocator2>& B)
{}


template<class T, class Prop, class Storage, class Allocator>
void CheckCholesky(Matrix<T, Prop, Storage, Allocator>& A)
{}

template<class Allocator>
void CheckCholesky(Matrix<Real_wp, Symmetric, BlockDiagRowSym, Allocator>& A)
{
  int m = A.GetM();
  Matrix<Real_wp, Symmetric, RowSymPacked> Adense, Bdense(m, m);
  Copy(A, Adense);

  Matrix<Real_wp, Symmetric, BlockDiagRowSym, Allocator> B(A);
  Mlt(A, A, B);
  Mlt(Adense, Adense, Bdense);
  
  //A.WriteText("A.dat");
  //Adense.WriteText("Adense.dat");
  //B.WriteText("B.dat");
  //Bdense.WriteText("Bdense.dat");
  Real_wp norm_ref = MaxAbs(Bdense);
  if (!EqualMatrix(B, Bdense, norm_ref*threshold))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Copy(B, Bdense);
  GetCholesky(B);
  GetCholesky(Bdense);
  
  Vector<Real_wp> X, Y;
  GenerateRandomVector(X, m);
  Y = X;
  
  SolveCholesky(SeldonNoTrans, B, X);
  SolveCholesky(SeldonTrans, B, X);
  SolveCholesky(SeldonNoTrans, Bdense, Y);
  SolveCholesky(SeldonTrans, Bdense, Y);
  for (int i = 0; i < m; i++)
    if ((abs(X(i) - Y(i)) > norm_ref*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "SolveCholesky incorrect" << endl;
        abort();
      }

  MltCholesky(SeldonTrans, B, X);
  MltCholesky(SeldonNoTrans, B, X);
  MltCholesky(SeldonTrans, Bdense, Y);
  MltCholesky(SeldonNoTrans, Bdense, Y);
  for (int i = 0; i < m; i++)
    if ((abs(X(i) - Y(i)) > norm_ref*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "MltCholesky incorrect" << endl;
        abort();
      }  
}

template<class T, class Prop, class Storage, class Allocator,
         class T2, class Prop2, class Storage2, class Allocator2>
void CheckBlockDiagonalMatrix(Matrix<T, Prop, Storage, Allocator>& Atest,
                              Matrix<T2, Prop2, Storage2, Allocator2>& Btest)
{
  Matrix<T, Prop, Storage, Allocator> A;
  
  // random pattern
  int m = 500; IVect permut;
  GenerateRandomPermutation(m, permut);
  
  IVect size_block;
  int nb = 0;
  while (nb < m)
    {
      int k = rand()%25 + 1;
      if (nb+k > m)
        k = m-nb;
      
      size_block.PushBack(k);
      nb += k;
    }
  
  //DISP(size_block);
  Vector<IVect> pattern;
  pattern.Reallocate(size_block.GetM());
  nb = 0;
  for (int i = 0; i < size_block.GetM(); i++)
    {
      pattern(i).Reallocate(size_block(i));
      for (int j = 0; j < size_block(i); j++)
        {
          pattern(i)(j) = permut(nb);
          nb++;
        }
    }
  
  //DISP(pattern);
  A.SetPattern(pattern);

  // testing conversion between block-diagonal matrices and sparse matrices
  T coef; SetComplexReal(1.0/RAND_MAX, coef); 
  A.FillRand(); A *= coef;
  Matrix<T2, Prop2, Storage2, Allocator2> B;
  Copy(A, B);
  
  if (!EqualMatrix(A, B))
    {
      cout << "Copy incorrect" <<endl;
      abort();
    }
  
  Matrix<T, Prop, Storage, Allocator> Abis;
  Copy(B, Abis);
  
  if (!EqualMatrix(Abis, B))
    {
      cout << "Copy incorrect" <<endl;
      abort();
    }

  // testing GetRowSize and GetSizeMaxBlock
  CheckSizeMat(A, B);
  
  // testing function Mlt
  Vector<T> X, Y(m), Xblock;
  GenerateRandomVector(X, m);  
  Mlt(B, X, Y);
  
  Xblock = X;
  Mlt(A, Xblock);
  for (int i = 0; i < m; i++)
    if ((abs(Xblock(i) - Y(i)) > threshold) || isnan(abs(Xblock(i) - Y(i))))
      {
        cout << "Mlt incorrect " << endl;
        abort();
      }
  
  // testing function MltAdd for transpose, and conjugate transpose
  T alpha, beta; GetRand(alpha); GetRand(beta);
  GenerateRandomVector(Y, m);  
  Vector<T> Yblock(Y); Xblock = X;
  MltAdd(alpha, SeldonNoTrans, B, X, beta, Y);
  MltAdd(alpha, SeldonNoTrans, A, Xblock, beta, Yblock);
  for (int i = 0; i < m; i++)
    if ((abs(Yblock(i) - Y(i)) > threshold) || isnan(abs(Yblock(i) - Y(i))))
      {
        cout << "MltAdd incorrect " << endl;
        abort();
      }
  
  Yblock = Y;
  MltAdd(alpha, SeldonTrans, B, X, beta, Y);
  MltAdd(alpha, SeldonTrans, A, Xblock, beta, Yblock);
  for (int i = 0; i < m; i++)
    if ((abs(Yblock(i) - Y(i)) > threshold) || isnan(abs(Yblock(i) - Y(i))))
      {
        cout << "MltAdd incorrect (transpose)" << endl;
        abort();
      }

  Yblock = Y;
  MltAdd(alpha, SeldonConjTrans, B, X, beta, Y);
  MltAdd(alpha, SeldonConjTrans, A, Xblock, beta, Yblock);
  for (int i = 0; i < m; i++)
    if ((abs(Yblock(i) - Y(i)) > threshold) || isnan(abs(Yblock(i) - Y(i))))
      {
        cout << "MltAdd incorrect (transpose conjugate)" << endl;
        abort();
      }
  
  // testing Mlt for a matrix
  if (realpart(alpha) == Real_wp(0))
    {
      Real_wp alpha_r; GetRand(alpha_r);
      alpha += alpha_r;
    }
  
  Mlt(realpart(alpha), A);
  Mlt(realpart(alpha), B);
  if (!EqualMatrix(A, B))
    {
      cout << "Mlt incorrect" <<endl;
      abort();
    }
  
  // testing GetInverse
  Abis = A;
  GetInverse(A);

  GenerateRandomVector(X, m);  
  Y = X;
  Mlt(Abis, X);
  Mlt(A, X);
  for (int i = 0; i < m; i++)
    if ((abs(X(i) - Y(i)) > 1000.0*threshold) || isnan(abs(X(i) - Y(i))))
      {
        DISP(i); DISP(X(i)); DISP(Y(i));
        cout << "GetInverse incorrect" << endl;
        abort();
      }
  
  // testing GetCholesky
  CheckCholesky(A);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 100*epsilon_machine;
  
  int n = toInteger(-log10(epsilon_machine));
  cout.precision(n);
  
  /* srand(time(NULL));
  if (argc < 3)
    {
      cout << "Entrez trois arguments : taille bloc, taille matrice et nombre iterations" << endl;
      abort();
    }
  
  int size_block = atoi(argv[1]);
  int size_mat = atoi(argv[2]);
  int nb_iter = atoi(argv[3]);
  Matrix<double, Symmetric, BlockDiagRowSym> A;
  PerfBlockDiagMatrix(A, size_block, size_mat, nb_iter); */
  
  bool overall_success = true;

  {
    Matrix<Real_wp, General, BlockDiagRow> A;
    int success = TestBlockDiagMatrix(A);
    if (success != 0)
      {
        overall_success = false;
        cout << "Problem on BlockDiagRow" << endl;
      }
    
    Matrix<Real_wp, General, ArrayRowSparse> B;
    CheckBlockDiagonalMatrix(A, B);
  }
    
  {
    Matrix<Real_wp, Symmetric, BlockDiagRowSym> A;
    int success = TestBlockDiagMatrix(A);
    if (success != 0)
      {
        overall_success = false;
        cout << "Problem on BlockDiagRowSym" << endl;
      }

    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> B;
    CheckBlockDiagonalMatrix(A, B);
  }


  {
    Matrix<Complex_wp, General, BlockDiagRow> A;
    Matrix<Complex_wp, General, ArrayRowComplexSparse> B;
    CheckBlockDiagonalMatrix(A, B);
  }
    
  {
    Matrix<Complex_wp, Symmetric, BlockDiagRowSym> A;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckBlockDiagonalMatrix(A, B);
  }
  
  if (overall_success)
    {
      cout << "All tests passed successfully" << endl;
      return 0;
    }
  
  return FinalizeMontjoie();
}
