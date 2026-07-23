#ifndef MONTJOIE_FILE_TINY_BLOCK_SOLVER_1D_CXX

namespace Seldon
{

  //! Default constructor
  template<class T, int r>
  TinyBlockSolver1D<T, r>::TinyBlockSolver1D()
  {
    nb_elem = 0;
  }


  //! the matrix is factorized
  template<class T, int r>
  void TinyBlockSolver1D<T, r>
  ::Factorize(const Matrix<T, General, BandedCol>& mat, bool hg)
  {
    vec_a12.Clear(); vec_a21.Clear(); vec_inv_a22.Clear();
    if (hg)
      {
	// homogeneous case : all blocks are identical
	// we retrieve a12, a21 and a22
	for(int i = 0; i < r; i++)
	  {
	    this->a12(0, i) = mat(0, i+1);
	    this->a12(1, i) = mat(r+1, i+1);
	    
	    this->a21(i, 0) = mat(i+1, 0);
	    this->a21(i, 1) = mat(i+1, r+1);
	    
	    for(int j = 0; j < r; j++)
	      this->inv_a22(i, j) = mat(i+1, j+1);
	  }
	
	// a22 is inverteed
	GetInverse(this->inv_a22);
	
	// initializing mat_tridiag with a11
	this->nb_elem = mat.GetM() / (r + 1);
	int N = this->nb_elem+1;
	this->mat_tridiag.Reallocate(N, N);
	this->mat_tridiag.Get(0, 0) = mat(0, 0);
	for (int i = 1; i < N; i++)
	  {
	    this->mat_tridiag.Get(i, i) = mat(i*(r+1), i*(r+1));
	    this->mat_tridiag.Get(i, i-1) = mat(i*(r+1), (i-1)*(r+1));
	    this->mat_tridiag.Get(i-1, i) = mat((i-1)*(r+1), i*(r+1));
	  }
	
	// we compute schur = a_12 a_22^{-1} a_21
	TinyMatrix<T, General, r, 2> prodTmp;
	TinyMatrix<T, General, 2, 2> schur;
	Mlt(this->inv_a22, this->a21, prodTmp);
	Mlt(this->a12, prodTmp, schur);
	
	// assembling a_11 - a_12 a_22^{-1} a_21
	for(int i = 0; i < this->nb_elem; i++)
	  {
	    this->mat_tridiag.Get(i,i) -= schur(0,0);
	    this->mat_tridiag.Get(i,i+1) -= schur(0,1);
	    this->mat_tridiag.Get(i+1,i) -= schur(1,0);
	    this->mat_tridiag.Get(i+1,i+1) -= schur(1,1);
	  }
	
	// tridiagonal matrix is factorized
	this->mat_tridiag.Factorize();

	// a_12 is replaced by a_12 a_22^{-1}
	TinyMatrix<T, General, 2, r> B;
	B = this->a12;
	Mlt(B, this->inv_a22, this->a12);

	return;
      }

    // inhomogeneous case
    this->nb_elem = mat.GetM() / (r+1);

    int N =  this->nb_elem + 1;
    this->mat_tridiag.Reallocate(N, N);
    this->mat_tridiag.Zero();

    // initializing mat_tridiag with a11
    this->mat_tridiag.Get(0, 0) = mat(0, 0);
    for (int i = 1; i < N; i++)
      {
	this->mat_tridiag.Get(i, i) = mat(i*(r+1), i*(r+1));
	this->mat_tridiag.Get(i, i-1) = mat(i*(r+1), (i-1)*(r+1));
	this->mat_tridiag.Get(i-1, i) = mat((i-1)*(r+1), i*(r+1));
      }
    
    // allocation of Schur blocks
    TinyMatrix<T, General, r, 2> prodTmp;
    TinyMatrix<T, General, 2, 2> schur;
    TinyMatrix<T, General, 2, r> B;
    vec_a12.Reallocate(N-1);
    vec_a21.Reallocate(N-1);
    vec_inv_a22.Reallocate(N-1);
    int offset = 0;
    for (int b = 0; b < N-1; b++)
      {
	// we retrieve a12, a21 and a22
	for(int i = 0; i < r; i++)
	  {
	    this->vec_a12(b)(0, i) = mat(offset, offset+i+1);
	    this->vec_a12(b)(1, i) = mat(offset+r+1, offset+i+1);
	    
	    this->vec_a21(b)(i, 0) = mat(offset+i+1, offset);
	    this->vec_a21(b)(i, 1) = mat(offset+i+1, offset+r+1);
	    
	    for(int j = 0; j < r; j++)
	      this->vec_inv_a22(b)(i, j) = mat(offset+i+1, offset+j+1);
	  }
	
	// a22 is inverteed
	GetInverse(this->vec_inv_a22(b));

	// we compute schur = a_12 a_22^{-1} a_21
	Mlt(this->vec_inv_a22(b), this->vec_a21(b), prodTmp);
	Mlt(this->vec_a12(b), prodTmp, schur);
	
	// assembling a_11 - a_12 a_22^{-1} a_21
	this->mat_tridiag.Get(b, b) -= schur(0,0);
	this->mat_tridiag.Get(b, b+1) -= schur(0,1);
	this->mat_tridiag.Get(b+1, b) -= schur(1,0);
	this->mat_tridiag.Get(b+1, b+1) -= schur(1,1);

	// a_12 is replaced by a_12 a_22^{-1}
	B = this->vec_a12(b);
	Mlt(B, this->vec_inv_a22(b), this->vec_a12(b));

	// next element
	offset += r+1;
      }

    // tridiagonal matrix is factorized
    this->mat_tridiag.Factorize();
  }


  //! solves A X = B
  template<class T, int r> template<class T0>
  void TinyBlockSolver1D<T, r>::Solve(Vector<T0>& X)
  {
    if (this->nb_elem <= 0)
      return;

    bool hg = (this->vec_a12.GetM() <= 0);    
    int N = this->nb_elem+1;
    Vector<T0> Xtri(N);
    for (int i = 0; i < N; i++)
      Xtri(i) = X(i*(r+1));
    
    TinyVector<T0, r> G, prodY;
    TinyVector<T0, 2> prodG;
    int offset = 1;
    for(int i = 0; i < this->nb_elem; i++)
      {
	TinyVectorLoop<r>::ExtractVector(X, offset, G);

	if (hg)	  
	  Mlt(this->a12, G, prodG);
	else
	  Mlt(this->vec_a12(i), G, prodG);
	
	Xtri(i) -= prodG(0);
	Xtri(i+1) -= prodG(1);
	offset += r+1;
      }

    this->mat_tridiag.Solve(Xtri);

    for (int i = 0; i < N; i++)
      X(i*(r+1)) = Xtri(i);

    T one; SetComplexOne(one);
    offset = 1;
    for(int i = 0; i < this->nb_elem; i++)
      {
	TinyVectorLoop<r>::ExtractVector(X, offset, G);
	
	// we compute G - a_12 X
	prodG(0) = Xtri(i);
	prodG(1) = Xtri(i+1);

	if (hg)
	  MltAdd(-one, this->a21, prodG, G);
	else
	  MltAdd(-one, this->vec_a21(i), prodG, G);
	
	// then we put a_22^{-1} G in X
	if (hg)
	  Mlt(this->inv_a22, G, prodY);
	else
	  Mlt(this->vec_inv_a22(i), G, prodY);
	
	TinyVectorLoop<r>::ExtractVector(prodY, offset, X);
	offset += r+1;
      }
  }
  
}

#define MONTJOIE_FILE_TINY_BLOCK_SOLVER_1D_CXX
#endif
