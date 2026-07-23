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
      GetRandNumber(x);
      A.Set(i, j, x);
    }
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - y(i)) > eps)
      return false;
  
  return true;
}


template<class T, class T2>
bool EqualMatrix(const Matrix<T, General, ColMajor>& x,
		 const Matrix<T2, General, ColMajor>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;

  if (x.GetN() != y.GetN())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    for (int j = 0; j < x.GetN(); j++)
      if ( (abs(x(i, j) - y(i, j)) > eps) || isnan(abs(x(i, j) - y(i, j))) )
	return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void CheckDistributedSolver(Matrix<T, Prop, Storage, Allocator>& A)
{
  T zero, one;
  SetComplexZero(zero);
  SetComplexOne(one);
  
  int n = 400, nnz = 5000;
  //int n = 7, nnz = 20;
  GenerateRandomMatrix(A, n, n, nnz);

  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 0;
      for (int j = 0; j < n; j++)
	sum += abs(A(i, j));
      
      A.Set(i, i, 1.0 + sum);
    }

  // matrix with same pattern but different values
  Matrix<T, Prop, Storage, Allocator> B;
  GenerateRandomMatrix(B, n, n, nnz);

  for (int i = 0; i < n; i++)
    {
      Real_wp sum = 0;
      for (int j = 0; j < n; j++)
	sum += abs(B(i, j));
      
      B.Set(i, i, 1.0 + sum);
    }

  Matrix<T, Prop, Storage, Allocator> B2(B), Bref(B);
  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetRowSize(i); j++)
      GetRandNumber(B2.Value(i, j));

  Matrix<T, Prop, Storage, Allocator> Bref2(B2);
    
  typedef typename ClassComplexType<T>::Tcplx Tcplx;
  Vector<T> x, y, b(n), bt(n);
  Vector<Tcplx> xc(n), bc(n), bt_c(n), yc;

  int nb_rhs = 20;
  Matrix<T, General, ColMajor> x_mat(n, nb_rhs), b_mat(n, nb_rhs), bt_mat(n, nb_rhs);
  
  for (int i = 0; i < nb_rhs; i ++)
    {
      GenerateRandomVector(x, n);
      Mlt(A, x, b);
      Mlt(SeldonTrans, A, x, bt);
      SetCol(b, i, b_mat);
      SetCol(bt, i, bt_mat);
      SetCol(x, i, x_mat);
    }

  Matrix<T, General, ColMajor> y_mat(x_mat);
  
  GenerateRandomVector(x, n);

  b.Fill(zero); bt.Fill(zero);
  Mlt(A, x, b);
  Mlt(SeldonTrans, A, x, bt);
  y = x;
  
  GenerateRandomVector(xc, n);
  
  Mlt(A, xc, bc);
  Mlt(SeldonTrans, A, xc, bt_c);
  yc = xc;

  Vector<T> xb, bb(n), bb2(n), xsolb;
  GenerateRandomVector(xb, n);

  bb.Zero(); bb2.Zero();
  Mlt(B, xb, bb);
  Mlt(B2, xb, bb2);

  Matrix<T, Prop, typename GenericStorage<Storage>::StorageCSR> Bs, Brefs, Bs2, Brefs2;
  Copy(Bref, Brefs);
  Copy(Bref2, Brefs2);
  
  {
    SparseDistributedSolver<T> mat_lu;
    
    for (int t = 0; t < 7; t++)
      if ((SparseDirectSolver<T>::IsAvailableSolver(t)) 
	  && (t != SparseDirectSolver<T>::ILUT) && (t != SparseDirectSolver<T>::SELDON_SOLVER))
	{
	  cout << "Testing Solver " << t << endl;
	  mat_lu.SelectDirectSolver(t);	  
	  
	  mat_lu.ShowMessages();

	  mat_lu.Factorize(A, true);
	  
	  int ierr = 0;
	  int info = mat_lu.GetInfoFactorization(ierr);
	  if (info != mat_lu.FACTO_OK)
	    {
	      cout << "Factorize incorrect" << endl;
	      abort();
	    }
	  
	  mat_lu.HideMessages();
	  
	  x = b;
	  mat_lu.Solve(x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }    
	  
	  x = bt;
	  mat_lu.TransSolve(x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  x = b;
	  mat_lu.Solve(SeldonNoTrans, x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  
	  x = bt;
	  mat_lu.Solve(SeldonTrans, x);
	  
	  if (!EqualVector(x, y))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  	  
	  // testing resolution of multiple right hand-sides
	  x_mat = b_mat;
	  mat_lu.Solve(x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }

	  x_mat = bt_mat;
	  mat_lu.TransSolve(x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }
	  
	  x_mat = b_mat;
	  mat_lu.Solve(SeldonNoTrans, x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }
	  
	  x_mat = bt_mat;
	  mat_lu.Solve(SeldonTrans, x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }

	  // testing solution with a complex vector
	  xc = bc;
	  mat_lu.Solve(xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }    
	  
	  xc = bt_c;
	  mat_lu.TransSolve(xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  xc = bc;
	  mat_lu.Solve(SeldonNoTrans, xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  
	  xc = bt_c;
	  mat_lu.Solve(SeldonTrans, xc);
	  
	  if (!EqualVector(xc, yc))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }


	  // testing with scaling
	  Matrix<T, Prop, Storage, Allocator> B(A);
	  mat_lu.Factorize(B, false, true);
	  
	  info = mat_lu.GetInfoFactorization(ierr);
	  if (info != mat_lu.FACTO_OK)
	    {
	      cout << "Factorize incorrect" << endl;
	      abort();
	    }
	  
	  x = b;
	  mat_lu.Solve(x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }    
	  
	  x = bt;
	  mat_lu.TransSolve(x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  x = b;
	  mat_lu.Solve(SeldonNoTrans, x);
	  
	  if (!EqualVector(x, y, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  
	  x = bt;
	  mat_lu.Solve(SeldonTrans, x);
	  
	  if (!EqualVector(x, y))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  	  
	  // testing resolution of multiple right hand-sides
	  x_mat = b_mat;
	  mat_lu.Solve(x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }

	  x_mat = bt_mat;
	  mat_lu.TransSolve(x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }
	  
	  x_mat = b_mat;
	  mat_lu.Solve(SeldonNoTrans, x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }
	  
	  x_mat = bt_mat;
	  mat_lu.Solve(SeldonTrans, x_mat);
	  
	  if (!EqualMatrix(x_mat, y_mat))
	    {
	      cout << "Solve of SparseDirectSolver incorrect" << endl;
	      abort();
	    }
	 
	  // testing solution with a complex vector
	  xc = bc;
	  mat_lu.Solve(xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }    
	  
	  xc = bt_c;
	  mat_lu.TransSolve(xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  xc = bc;
	  mat_lu.Solve(SeldonNoTrans, xc);
	  
	  if (!EqualVector(xc, yc, 10.0*threshold))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
	  
	  xc = bt_c;
	  mat_lu.Solve(SeldonTrans, xc);
	  
	  if (!EqualVector(xc, yc))
	    {
	      cout << "Solve of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  // factorisation is cleared
	  mat_lu.Clear();

#ifdef SELDON_WITH_MUMPS
	  // testing PerformAnalysis/PerformFactorization
	  B = Bref;
	  mat_lu.PerformAnalysis(B);
	  mat_lu.PerformFactorization(B);

	  xsolb = bb;
	  mat_lu.Solve(xsolb);

	  if (!EqualVector(xsolb, xb, 10.0*threshold))
	    {
	      cout << "PerformAnalysis/PerformFactorization of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  B2 = Bref2;
	  mat_lu.PerformFactorization(B2);

	  xsolb = bb2;
	  mat_lu.Solve(xsolb);

	  if (!EqualVector(xsolb, xb, 10.0*threshold))
	    {
	      cout << "PerformAnalysis/PerformFactorization of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  Bs = Brefs;
	  mat_lu.PerformAnalysis(Bs);
	  mat_lu.PerformFactorization(Bs);

	  xsolb = bb;
	  mat_lu.Solve(xsolb);

	  if (!EqualVector(xsolb, xb, 10.0*threshold))
	    {
	      cout << "PerformAnalysis/PerformFactorization of SparseDirectSolver incorrect " << endl;
	      abort();
	    }

	  Bs2 = Brefs2;
	  mat_lu.PerformFactorization(Bs2);

	  xsolb = bb2;
	  mat_lu.Solve(xsolb);

	  if (!EqualVector(xsolb, xb, 10.0*threshold))
	    {
	      cout << "PerformAnalysis/PerformFactorization of SparseDirectSolver incorrect " << endl;
	      abort();
	    }
#endif
	}
  }

}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e6*epsilon_machine;
  
  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckDistributedSolver(A);
  }
  
  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckDistributedSolver(A);
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckDistributedSolver(A);
  }
  
  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckDistributedSolver(A);
  }
  
  cout << "All tests passed successfully" << endl;
  return FinalizeMontjoie();
}

