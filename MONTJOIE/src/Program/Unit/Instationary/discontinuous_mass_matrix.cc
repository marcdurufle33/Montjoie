#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL

#include "Instationary/MontjoieUnsteady.hxx"

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
  x = complex<T>(rand(), rand())/RAND_MAX;
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
    {
      DISP(x.GetM()); DISP(y.GetM());
      return false;
    }

  Real_wp yref = y.GetNormInf();
  if ((Norm2(x) == 0) && (Norm2(y) == 0))
    return true;
  
  if ((Norm2(x) <= eps) || isnan(yref))
    {
      DISP(Norm2(x)); DISP(Norm2(y));
      return false;
    }

  for (int i = 0; i < x.GetM(); i++)
    if (isnan(abs(x(i) - y(i))) || (abs(x(i) - y(i))/yref > eps))
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }
  
  return true;
}

void CheckBlockMass()
{
  int n = 99;
  
  // diagonal mass matrix
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> Mh(n, n);
    Mh.Zero();
    VectReal_wp diagonal(n);
    for (int i = 0; i < n; i++)
      {
	GetRandNumber(diagonal(i));
	Mh(i, i) = diagonal(i);
      }

    DiscontinuousDiagonalMassMatrix block(diagonal);
  
    Vector<Real_wp> x0, x, y(n);
    y.Zero();
    GenerateRandomVector(x0, n);
  
    x = x0;
    block.MltMass(x);
    Mlt(Mh, x0, y);
    
    if (!EqualVector(y, x))
      {
	DISP(x); DISP(y);
	cout << "MltMass incorrect" << endl;
	abort();
      }
  
    block.SolveMass(x);
    
    if (!EqualVector(x0, x))
      {
	DISP(x); DISP(x0);
	cout << "SolveMass incorrect" << endl;
	abort();
      }
  
    x = y;
    block.SolveCholesky(SeldonNoTrans, x);
    block.SolveCholesky(SeldonTrans, x);
    
    if (!EqualVector(x0, x))
      {
	DISP(x); DISP(x0);
	cout << "SolveCholesky incorrect" << endl;
	abort();
      }
  }

  
  // block-diagonal mass matrix
  {
    VarRandomGenerator var;
    IVect permut(n); permut.Fill();
    var.GenerateRandomPermutation(n, permut);
    //DISP(permut);
    
    Matrix<Real_wp, Symmetric, RowSymPacked> Mh(n, n);
    Mh.Zero();
    Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> > block_diagonal(n/3);

    Matrix3_3 B; Matrix3_3sym A;
    for (int i = 0; i < n/3; i++)
      {
	int nx = permut(3*i), ny = permut(3*i+1), nz = permut(3*i+2);
	GetRandNumber(B(0, 0)); GetRandNumber(B(0, 1)); GetRandNumber(B(0, 2));
	GetRandNumber(B(1, 0)); GetRandNumber(B(1, 1)); GetRandNumber(B(1, 2));
	GetRandNumber(B(2, 0)); GetRandNumber(B(2, 1)); GetRandNumber(B(2, 2));
	MltTrans(B, B, A);
	block_diagonal(i) = A;
	Mh(nx, nx) = A(0, 0); 
	Mh(ny, ny) = A(1, 1);
	Mh(nz, nz) = A(2, 2);
	Mh(nx, ny) = A(0, 1);
	Mh(nx, nz) = A(0, 2);
	Mh(ny, nz) = A(1, 2);
      }
    
    DiscontinuousBlockDiagMassMatrix<3> block(permut, block_diagonal);
    
    Vector<Real_wp> x0, x, y(n);
    y.Zero();
    GenerateRandomVector(x0, n);
  
    x = x0;
    block.MltMass(x);
    Mlt(Mh, x0, y);
    
    if (!EqualVector(y, x))
      {
	DISP(x); DISP(y);
	cout << "MltMass incorrect" << endl;
	abort();
      }
  
    block.SolveMass(x);
    
    if (!EqualVector(x0, x))
      {
	DISP(x); DISP(x0);
	cout << "SolveMass incorrect" << endl;
	abort();
      }
  
    x = y;
    block.SolveCholesky(SeldonNoTrans, x);
    block.SolveCholesky(SeldonTrans, x);
    
    if (!EqualVector(x0, x))
      {
	DISP(x); DISP(x0);
	cout << "SolveCholesky incorrect" << endl;
	abort();
      }
  }

  /*
  // dense mass matrix
  Matrix<Real_wp> B(n, n), C(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(B(i, j));
        C(i, j) = 0.0;
      }
  
  MltAdd(1.0, SeldonNoTrans, B, SeldonTrans, B, 0.0, C);
  
  for (int i = 0; i < n; i++)
    {
      for (int j = i; j < n; j++)
        Mh(i, j) = C(i, j)/n;
      
      Mh(i, i) += 1.0;
    }

  block.type_matrix = block.DENSE;
  block.dense_Mh = Mh;
  
  block.Factorize();
  
  x = x0;
  block.MltMass(x);
  Mlt(Mh, x0, y);
  
  if (!EqualVector(y, x))
    {
      DISP(x); DISP(y);
      cout << "MltMass incorrect" << endl;
      abort();
    }
  
  block.SolveMass(x);
  
  if (!EqualVector(x0, x))
    {
      DISP(x); DISP(x0);
      cout << "SolveMass incorrect" << endl;
      abort();
    }
  
  x = y;
  block.SolveCholesky(SeldonNoTrans, x);
  block.SolveCholesky(SeldonTrans, x);
  
  if (!EqualVector(x0, x))
    {
      DISP(x); DISP(x0);
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  // sparse mass matrix
  Matrix<Real_wp, Symmetric, ArrayRowSymSparse> Sh;
  GenerateRandomMatrix(Sh, n, n, 3*n);
  
  VectReal_wp diag_Sh;
  GetRowSum(diag_Sh, Sh);
  for (int i = 0; i < n; i++)
    Sh.Get(i, i) = 1.0 + diag_Sh(i);
  
  Copy(Sh, Mh);
  
  block.dense_Mh.Clear();  
  block.type_matrix = block.SPARSE;
  Copy(Sh, block.sparse_Mh);
  
  block.Factorize();
  
  x = x0;
  block.MltMass(x);
  Mlt(Mh, x0, y);
  
  if (!EqualVector(y, x))
    {
      DISP(x); DISP(y);
      cout << "MltMass incorrect" << endl;
      abort();
    }
  
  block.SolveMass(x);
  
  if (!EqualVector(x0, x))
    {
      DISP(x); DISP(x0);
      cout << "SolveMass incorrect" << endl;
      abort();
    }
  
  x = y;
  block.SolveCholesky(SeldonNoTrans, x);
  block.SolveCholesky(SeldonTrans, x);
  
  if (!EqualVector(x0, x))
    {
      DISP(x); DISP(x0);
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }
  */
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e4*epsilon_machine;

  CheckBlockMass();

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
