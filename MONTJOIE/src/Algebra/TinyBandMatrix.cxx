#ifndef SELDON_FILE_TINY_BAND_MATRIX_CXX

#include "TinyBandMatrix.hxx"

namespace Seldon
{
  
  /******************
   * TinyBandMatrix *
   ******************/
  
  
  //! Matrix is erased
  template<class T, int d>    
  void TinyBandMatrix<T, d>::Clear()
  {
    diag_.Clear();
    this->m_ = 0;
    this->n_ = 0;
  }
  
  
  //! all elements of the matrix are set to 0
  template<class T, int d>    
  void TinyBandMatrix<T, d>::Zero()
  {
    diag_.Zero();
  }
  
  
  //! no message is displayed during factorization
  template<class T, int d>    
  void TinyBandMatrix<T, d>::HideMessages()
  {
  }
  
  
  //! changes the size of the matrix while resetting previous entries
  template<class T, int d>    
  void TinyBandMatrix<T, d>::Reallocate(int m, int n)
  {
    this->m_ = m;
    this->n_ = n;
    diag_.Reallocate(n);
    FillZero(diag_);
  }
  
  
  //! adding several terms val in a_{i,j}
  template<class T, int d>    
  void TinyBandMatrix<T, d>::AddInteractionRow(int i, int n, const Vector<int>& num,
                                               const Vector<T>& val, bool sorted)
  {
    for (int k = 0; k < n; k++)
      AddInteraction(i, num(k), val(k));
  }
  
  
  //! multiplication by a scalar
  template<class T, int d>
  TinyBandMatrix<T, d>& TinyBandMatrix<T, d>::operator *=(const T& alpha)
  {
    diag_ *= alpha;
    return *this;
  }
  
  
  //! Filling non-zero entries with a same value
  template<class T, int d>  template<class T0>
  void TinyBandMatrix<T, d>::Fill(const T0& x)
  {
    T x_; SetComplexReal(x, x_);
    TinyVector<T, 2*d+1> col;
    col.Fill(x_);
    this->diag_.Fill(col);
  }


  //! Filling non-zero entries with random values
  template<class T, int d>
  void TinyBandMatrix<T, d>::FillRand()
  {
    for (int i = 0; i < this->diag_.GetM(); i++)
      this->diag_(i).FillRand();
  }
  
  
  //! Sets the matrix to the identity matrix
  template<class T, int d>    
  void TinyBandMatrix<T, d>::SetIdentity()
  {
    T one; SetComplexOne(one);
    TinyVector<T, 2*d+1> x;
    x.Zero(); x(d) = one;
    for (int i = 0; i < this->diag_.GetM(); i++)
      this->diag_(i) = x;
  }
  
  
  //! clears a row of the matrix
  template<class T, int d>    
  void TinyBandMatrix<T, d>::ClearRow(int i)
  {
    T zero; SetComplexZero(zero);
    int n = diag_.GetM();
    for (int j = max(0, i-d); j < min(i+d+1, n); j++)
      this->Val(i, j) = zero;
  }
  
  
  //! conversion from a CSR matrix
  template<class T, int d>    
  void TinyBandMatrix<T, d>::Copy(const Matrix<T, General, ArrayRowSparse>& A)
  {
    int n = A.GetM();
    Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = A.GetRowSize(i);
	for (int k = 0; k < size_row; k++)
	  {
	    int j = A.Index(i, k);
	    T val = A.Value(i, k);
	    AddInteraction(i, j, val);
	  }
      }
  }
  
  
  //! LU factorization
  template<class T, int d>
  void TinyBandMatrix<T, d>::Factorize()
  {
    int n = diag_.GetM();
    if (n <= (2*d))
      {
        cout << "Matrix too small" << endl;
        abort();
      }
       
    T one; SetComplexOne(one);
    for (int j = 0; j < n-d; j++)
      {
	this->diag_(j)(d) = one/this->diag_(j)(d);
	
	// we eliminate column j
	TinyBandMatrixLoop<d>::EliminateColumn(*this, j);
      }
    
    // for last columns, basic loops
    T pivot;
    for (int j = n-d; j < n; j++)
      {
	diag_(j)(d) = one/diag_(j)(d);
	for (int p = 0; p < n-1-j; p++)
	  {
	    pivot = diag_(j)(d-1-p)*diag_(j)(d);
	    diag_(j)(d-1-p) = pivot;
	    for (int k = j+1; k < n; k++)
	      {
		int m = k-j;
		diag_(k)(d-1-p+m) -= pivot*diag_(k)(d+m);
	      }
	  }
      }
    
  }
  
  
  //! adds alpha*A to current matrix
  template<class T, int d> template<class T0>
  void TinyBandMatrix<T, d>::Add_(const T0& alpha, const TinyBandMatrix<T, d>& A)
  {
    for (int i = 0; i < this->diag_.GetM(); i++)
      Add(alpha, A.diag_(i), this->diag_(i));
  }
  
  
  //! operation y = y + alpha*A*x
  template<class T, int d> template<class T0, class T1>
  void TinyBandMatrix<T, d>
  ::MltAdd(const T0& alpha, const SeldonTranspose& trans, const Vector<T1>& x, Vector<T1>& y) const
  {
    int n = this->diag_.GetM();
    T1 val;
    if (trans.NoTrans())
      {
        if (n <= (2*d+1))
          {
            for (int j = 0; j < n; j++)
              {
                val = alpha*x(j);
                for (int i = max(j-d, 0); i < min(j+d+1, n); i++)
                  y(i) += diag_(j)(j-i+d)*val;
              }
            return;
          }
        
        // first columns
        for (int j = 0; j < d; j++)
          {
            val = alpha*x(j);
            for (int i = 0; i < j+d+1; i++)
              y(i) += diag_(j)(j-i+d)*val;
          }
        
        // main part of the matrix
        for (int j = d; j < n-d; j++)
          TinyBandMatrixLoop<2*d+1>::MltAdd(j+d+1, diag_(j), alpha*x(j), y);
        
        // last columns
        for (int j = n-d; j < n; j++)
          {
            val = alpha*x(j);
            for (int i = j-d; i < n; i++)
              y(i) += diag_(j)(j-i+d)*val;
          }
      }
    else if (trans.Trans())
      {
        if (n <= (2*d+1))
          {
            for (int j = 0; j < n; j++)
              {
                SetComplexZero(val);
                for (int i = max(j-d, 0); i < min(j+d+1, n); i++)
                  val += diag_(j)(j-i+d)*x(i);
                
                y(j) += alpha*val;
              }
            
            return;
          }
        
        // first columns
        for (int j = 0; j < d; j++)
          {
            SetComplexZero(val);
            for (int i = 0; i < j+d+1; i++)
              val += diag_(j)(j-i+d)*x(i);
            
            y(j) += alpha*val;
          }
        
        // main part of the matrix
        for (int j = d; j < n-d; j++)
          {
            SetComplexZero(val);
            TinyBandMatrixLoop<2*d+1>::MltTransAdd(j+d+1, diag_(j), x, val);
            y(j) += alpha*val;
          }
        
        // last columns
        for (int j = n-d; j < n; j++)
          {
            SetComplexZero(val);
            for (int i = j-d; i < n; i++)
              val += diag_(j)(j-i+d)*x(i);
            
            y(j) += alpha*val;
          }
      }
    else
      {
        if (n <= (2*d+1))
          {
            for (int j = 0; j < n; j++)
              {
                SetComplexZero(val);
                for (int i = max(j-d, 0); i < min(j+d+1, n); i++)
                  val += conjugate(diag_(j)(j-i+d))*x(i);
                
                y(j) += alpha*val;
              }
            
            return;
          }
        
        // first columns
        for (int j = 0; j < d; j++)
          {
            SetComplexZero(val);
            for (int i = 0; i < j+d+1; i++)
              val += conjugate(diag_(j)(j-i+d))*x(i);
            
            y(j) += alpha*val;
          }
        
        // main part of the matrix
        for (int j = d; j < n-d; j++)
          {
            SetComplexZero(val);
            TinyBandMatrixLoop<2*d+1>::MltConjTransAdd(j+d+1, diag_(j), x, val);
            y(j) += alpha*val;
          }
        
        // last columns
        for (int j = n-d; j < n; j++)
          {
            SetComplexZero(val);
            for (int i = j-d; i < n; i++)
              val += conjugate(diag_(j)(j-i+d))*x(i);
            
            y(j) += alpha*val;
          }
      }
  }
  
  
  //! LU resolution
  template<class T, int d> template<class T0>
  void TinyBandMatrix<T, d>::Solve(Vector<T0>& x)
  {
    int n = diag_.GetM();
    // solving L x = x
    
    // basic loops for first rows
    for (int j = 1; j < d; j++)
      for (int k = 0; k < j; k++)
	x(j) -= diag_(k)(k-j+d)*x(k);
    
    // then optimized loops
    for (int j = d; j < n; j++)
      TinyBandMatrixLoop<d>::SolveLower(*this, j, x);
    
    // solving U x = x
    
    // basic loops for last rows
    for (int j = n-1; j >= n-d; j--)
      {
	for (int k = j+1; k < n; k++)
	  x(j) -= diag_(k)(k-j+d)*x(k);
	
	// dividing by coefficient of diagonal
	x(j) *= diag_(j)(d);
      }
    
    // then optimized loops
    for (int j = n-d-1; j >= 0; j--)
      {
	TinyBandMatrixLoop<d>::SolveUpper(*this, j, x);
	
	// dividing by coefficient of diagonal
	x(j) *= diag_(j)(d);
      }
    
  }
  
  
  //! matrix is written in binary format
  template<class T, int d>
  void TinyBandMatrix<T, d>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyBandMatrix::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();    
  }
  
  
  //! matrix is written in binary format
  template<class T, int d>
  void TinyBandMatrix<T, d>::Write(ostream& FileStream) const
  {
    abort();
    //this->diag_.Write(FileStream);
  }
  
  
  //! matrix is written in ascii format
  template<class T, int d>
  void TinyBandMatrix<T, d>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyBandMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! matrix is written in ascii format
  template<class T, int d>
  void TinyBandMatrix<T, d>::WriteText(ostream& FileStream) const
  {
    for (int j = 0; j < this->diag_.GetM(); j++)
      {
	for (int p = 0; p <= 2*d; p++)
	  {
	    int i = j + d - p;
	    if ((i >= 0)&&(i < this->diag_.GetM()))
	      FileStream << i+1 << " " << j+1 << " " << this->diag_(j)(p) << '\n';
	  }
      }
  }
  
  
  //! computes row_sum_i = \sum_j |a_{i,j}|
  template<class T, int d>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& row_sum, const TinyBandMatrix<T, d>& A)
  {
    row_sum.Reallocate(A.GetM());
    row_sum.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int  j = max(i-d,0); j < min(i+d+1,A.GetM()); j++)
        row_sum(i) += abs(A(i, j));    
  }
  
  
  //! scales matrix by coefficients scale (A is replaced by scale A)
  template<class T, int d>
  void ScaleLeftMatrix(TinyBandMatrix<T, d>& A, const Vector<typename ClassComplexType<T>::Treal>& scale)
  {
    if (A.GetM() <= (2*d+1))
      {
        for (int i = 0; i < A.GetM(); i++)
          for (int  j = max(i-d,0); j < min(i+d+1,A.GetM()); j++)
            A.Get(i, j) *= scale(i);
        
        return;
      }
    
    // first rows
    for (int i = 0; i < d; i++)
      for (int  j = 0; j <= i+d; j++)
        A.Get(i, j) *= scale(i);
    
    // main part of the matrix
    for (int i = d; i < A.GetM()-d; i++)
      TinyBandMatrixLoop<2*d+1>::MltRow(i, scale(i), A);

    // last rows
    for (int i = A.GetM()-d; i < A.GetM(); i++)
      for (int  j = i-d; j < A.GetM(); j++)
        A.Get(i, j) *= scale(i);
  }
  

  //! LU factorisation
  template<class T, int d>
  void GetLU(TinyBandMatrix<T, d>& A,
             TinyBandMatrix<T, d>& mat_lu, bool keep_matrix)
  {
    mat_lu = A;
    if (!keep_matrix)
      A.Clear();
    
    mat_lu.Factorize();
  }
  
  
  //! LU factorisation
  template<class T, int d>
  void GetLU(TinyBandMatrix<T, d>& A)
  {
    A.Factorize();
  }
  
  
  //! conversion from ArrayRowSparse to band matrix
  template<class T, int d>
  void Copy(const Matrix<T, General, ArrayRowSparse>& A, TinyBandMatrix<T, d>& B)
  {
    B.Copy(A);
  }
  
  
  //! resolution of A x = b, once GetLU has been called
  template<class T, int d>
  void SolveLU(TinyBandMatrix<T, d>& mat_lu, Vector<T>& x)
  {
    mat_lu.Solve(x);
  }
  
  
  //! B = B + alpha*A
  template<class T0, class T1, int d>
  void Add(const T0& alpha, const TinyBandMatrix<T1, d>& A, TinyBandMatrix<T1, d>& B)           
  {
    B.Add_(alpha, A);
  }
  
  
  //! y = beta*y + alpha*A*x
  template<class T0, class T, int d, class T2>
  void MltAdd(const T0& alpha, const TinyBandMatrix<T, d>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y)
  {
    T0 zero; SetComplexZero(zero);
    T2 czero; SetComplexZero(czero);
    T0 one; SetComplexOne(one);
    if (beta == zero)
      y.Fill(czero);
    else if (beta != one)
      Mlt(beta, y);
    
    A.MltAdd(alpha, SeldonNoTrans, x, y);
  }


  //! y = beta*y + alpha*A*x (or A^T or A^H)
  template<class T0, class T, int d, class T2>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans, const TinyBandMatrix<T, d>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y)
  {
    T0 zero; SetComplexZero(zero);
    T2 czero; SetComplexZero(czero);
    T0 one; SetComplexOne(one);
    if (beta == zero)
      y.Fill(czero);
    else if (beta != one)
      Mlt(beta, y);
    
    A.MltAdd(alpha, trans, x, y);
  }


  //! y = A*x
  template<class T, int d, class T1>
  void Mlt(const SeldonTranspose& trans, const TinyBandMatrix<T, d>& A,
	   const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero; SetComplexZero(zero);
    T1 one; SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, trans, x, y);
  }

  
  //! y = A*x
  template<class T, int d, class T1>
  void Mlt(const TinyBandMatrix<T, d>& A,
	   const Vector<T1>& x, Vector<T1>& y)
  {
    T1 zero; SetComplexZero(zero);
    T1 one; SetComplexOne(one);
    y.Fill(zero);
    A.MltAdd(one, SeldonNoTrans, x, y);
  }

  
  //! A = A*alpha
  template<class T, int d>
  void Mlt(const T& alpha, TinyBandMatrix<T, d>& A)
  {
    A *= alpha;
  }
  
  
  //! displays matrix
  template<class T, int d>
  ostream& operator<<(ostream& out, const TinyBandMatrix<T, d>& A)
  {
    A.WriteText(out);
    return out;
  }

  
  /*******************
   * TinyArrowMatrix *
   *******************/

  
  //! matrix is cleared
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::Clear()
  {
    this->diag_.Clear();
    last_row_.Clear();
    last_col_.Clear();
  }
  
  
  //! all elements of the matrix are set to 0
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::Zero()
  {
    FillZero(this->diag_);
    FillZero(last_row_);
    FillZero(last_col_);
    FillZero(last_block_);
  }
  
  
  //! changes the size of the matrix while resetting previous entries
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::Reallocate(int nb, int n)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (n <= m)
      {
        cout << "Too small matrix" << endl;
        abort();
      }
#endif

    TinyBandMatrix<T, d>::Reallocate(n-m, n-m);
    last_row_.Reallocate(n-m);
    last_col_.Reallocate(n-m);
    FillZero(last_row_);
    FillZero(last_col_);
    last_block_.Zero();
  }
  
  
  //! adding a term val in a_{i,j}
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::AddInteraction(int i, int j, const T& val)
  {
    if (i >= last_row_.GetM())
      {
	if (j >= last_row_.GetM())
	  {
	    // last block
	    int i2 = i - last_row_.GetM();
	    int j2 = j - last_row_.GetM();
	    last_block_(i2, j2) += val;
	  }
	else
	  {
	    int i2 = i - last_row_.GetM();
	    last_row_(j)(i2) += val;
	  }
      }
    else if (j >= last_col_.GetM())
      {
        int j2 = j - last_col_.GetM();
	last_col_(i)(j2) += val;
      }
    else
      {
	// (i, j) in the band
	TinyBandMatrix<T, d>::AddInteraction(i, j, val);
      }
  }

  
  //! adding terms val(k) in a_{i,col(k)}
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T>& val, bool sorted)
  {
    for (int k = 0; k < nb; k++)
      AddInteraction(i, col(k), val(k));
  }
  
  
  //! returns A(i, j)
  template<class T, int d, int m>    
  T TinyArrowMatrix<T, d, m>::operator()(int i, int j) const
  {
    if (i >= last_row_.GetM())
      {
	if (j >= last_row_.GetM())
	  {
	    // last block
	    int i2 = i - last_row_.GetM();
	    int j2 = j - last_row_.GetM();
	    return last_block_(i2, j2);
	  }
	else
	  {
	    int i2 = i - last_row_.GetM();
	    return last_row_(j)(i2);
	  }
      }
    else if (j >= last_col_.GetM())
      {
        int j2 = j - last_col_.GetM();
	return last_col_(i)(j2);
      }
    
    // (i, j) in the band
    return static_cast<const TinyBandMatrix<T, d>& >(*this)(i, j);
    
  }
  
  
  //! multiplication by a scalar
  template<class T, int d, int m>
  TinyArrowMatrix<T, d, m>& TinyArrowMatrix<T, d, m>::operator *=(const T& alpha)
  {
    this->diag_ *= alpha;
    last_row_ *= alpha;
    last_col_ *= alpha;
    last_block_ *= alpha;
    return *this;
  }

  
  //! returns A(i, j)
  template<class T, int d, int m>
  T& TinyArrowMatrix<T, d, m>::Get(int i, int j)
  {
    if (i >= last_row_.GetM())
      {
	if (j >= last_row_.GetM())
	  {
	    // last block
	    int i2 = i - last_row_.GetM();
	    int j2 = j - last_row_.GetM();
	    return last_block_(i2, j2);
	  }
	else
	  {
	    int i2 = i - last_row_.GetM();
	    return last_row_(j)(i2);
	  }
      }
    else if (j >= last_col_.GetM())
      {
        int j2 = j - last_col_.GetM();
	return last_col_(i)(j2);
      }
    
    // (i, j) in the band
    return static_cast<TinyBandMatrix<T, d>& >(*this).Get(i, j);
  }


  //! returns A(i, j)
  template<class T, int d, int m>
  const T& TinyArrowMatrix<T, d, m>::Get(int i, int j) const
  {
    if (i >= last_row_.GetM())
      {
	if (j >= last_row_.GetM())
	  {
	    // last block
	    int i2 = i - last_row_.GetM();
	    int j2 = j - last_row_.GetM();
	    return last_block_(i2, j2);
	  }
	else
	  {
	    int i2 = i - last_row_.GetM();
	    return last_row_(j)(i2);
	  }
      }
    else if (j >= last_col_.GetM())
      {
        int j2 = j - last_col_.GetM();
	return last_col_(i)(j2);
      }
    
    // (i, j) in the band
    return TinyBandMatrix<T, d>::Get(i, j);
  }
  

  //! returns A(i, j)
  template<class T, int d, int m>
  T& TinyArrowMatrix<T, d, m>::Val(int i, int j)
  {
    return Get(i, j);
  }
  

  //! returns A(i, j)
  template<class T, int d, int m>
  const T& TinyArrowMatrix<T, d, m>::Val(int i, int j) const
  {
    return Get(i, j);
  }
  

  //! returns A(i, j)
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! returns A(i, j)
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::SetEntry(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }

  
  //! multiplication by a scalar
  template<class T, int d, int m> template<class T0>
  const TinyArrowMatrix<T, d, m>& TinyArrowMatrix<T, d, m>
  ::operator *=(const T0& x)
  {
    for (int i = 0; i < this->diag_.GetM(); i++)
      this->diag_(i) *= x;

    for (int i = 0; i < this->last_row_.GetM(); i++)
      this->last_row_(i) *= x;

    for (int i = 0; i < this->last_col_.GetM(); i++)
      this->last_col_(i) *= x;
    
    this->last_block_ *= x;
    
    return *this;
  }
  
  
  //! clears a row of the matrix
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::ClearRow(int i)
  {
    TinyBandMatrix<T, d>::ClearRow(i);
    
    T zero; SetComplexZero(zero);
    if (i >= last_row_.GetM())
      {
        int i2 = i - last_row_.GetM();
        for (int j = 0; j < last_row_.GetM(); j++)
          last_row_(j)(i2) = zero;
        
        for (int j = 0; j < m; j++)
          last_block_(i2, j) = zero;
      }
    else
      {
        for (int j = 0; j < m; j++)
          last_col_(i)(j) = zero;
      }
  }
  
  
  //! sets the matrix to the identity matrix
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::SetIdentity()
  {
    TinyBandMatrix<T, d>::SetIdentity();
    
    for (int i = 0; i < this->last_row_.GetM(); i++)
      this->last_row_(i).Zero();

    for (int i = 0; i < this->last_col_.GetM(); i++)
      this->last_col_(i).Zero();
    
    this->last_block_.SetIdentity();
  }
  
  
  //! Filling non-zero entries with a same value
  template<class T, int d, int m>  template<class T0>
  void TinyArrowMatrix<T, d, m>::Fill(const T0& x)
  {
    T x_; SetComplexReal(x, x_);
    TinyVector<T, 2*d+1> col;
    col.Fill(x_);
    this->diag_.Fill(col);
    
    TinyVector<T, m> elt;
    elt.Fill(x_);
    this->last_row_.Fill(elt);
    this->last_col_.Fill(elt);
    this->last_block_.Fill(x_);
  }


  //! Filling non-zero entries with random values
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::FillRand()
  {
    for (int i = 0; i < this->diag_.GetM(); i++)
      this->diag_(i).FillRand();

    for (int i = 0; i < this->last_row_.GetM(); i++)
      this->last_row_(i).FillRand();

    for (int i = 0; i < this->last_col_.GetM(); i++)
      this->last_col_(i).FillRand();
    
    this->last_block_.FillRand();
  }

  
  //! conversion from a CSR matrix
  template<class T, int d, int m>    
  void TinyArrowMatrix<T, d, m>::Copy(const Matrix<T, General, ArrayRowSparse>& A)
  {
    int n = A.GetM();
    Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = A.GetRowSize(i);
	for (int k = 0; k < size_row; k++)
	  {
	    int j = A.Index(i, k);
	    T val = A.Value(i, k);
	    AddInteraction(i, j, val);
	  }
      }
  }
  
  
  //! LU factorization
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::Factorize()
  {
    // factorization of the band-matrix
    TinyBandMatrix<T, d>::Factorize();
    
    T pivot;
    T one; SetComplexOne(one);
    // we apply combinations done for band-matrix to last rows and columns
    for (int j = 0; j < this->diag_.GetM()-d; j++)
      {
	// we eliminate column j
	TinyArrowMatrixLoop<d>::EliminateLastColumn(*this, j);
	
	// we eliminate last rows
	TinyArrowMatrixLoop<m>::EliminateLastRow(*this, j);
      }
    
    // for last columns, basic loops
    for (int j = this->diag_.GetM()-d; j < this->diag_.GetM(); j++)
      {
	// 
	for (int p = 0; p < this->diag_.GetM()-1-j; p++)
	  {
	    pivot = this->diag_(j)(d-1-p);
	    // loop over columns
	    for (int k = 0; k < m; k++)
	      last_col_(j+1+p)(k) -= pivot*last_col_(j)(k);
	  }
	
	for (int k = 0; k < m; k++)
	  {
	    pivot = last_row_(j)(k)*this->diag_(j)(d);
	    last_row_(j)(k) = pivot;
	    for (int p = 0; p < this->diag_.GetM()-1-j; p++)
	      last_row_(j+1+p)(k) -= pivot*this->diag_(j+1+p)(d+1+p);
	    
	    for (int p = 0; p < m; p++)
	      last_block_(k, p) -= pivot*last_col_(j)(p);
	  }	    
      }

    // now treating last_block_ with basic loops
    for (int i = 0; i < m; i++)
      {
	last_block_(i, i) = one/last_block_(i, i);
	for (int j = i+1; j < m; j++)
	  {
	    pivot = last_block_(j, i)*last_block_(i, i);
	    last_block_(j, i) = pivot;
	    for (int k = i+1; k < m; k++)
	      last_block_(j, k) -= pivot*last_block_(i, k);
	  }
      }
    
  }
  
  
  //! adds alpha*A to the current matrix
  template<class T, int d, int m> template<class T0>
  void TinyArrowMatrix<T, d, m>::Add_(const T0& alpha, const TinyArrowMatrix<T, d, m>& A)
  {
    for (int i = 0; i < this->diag_.GetM(); i++)
      Add(alpha, A.diag_(i), this->diag_(i));
    
    for (int i = 0; i < this->last_row_.GetM(); i++)
      Add(alpha, A.last_row_(i), this->last_row_(i));
    
    for (int i = 0; i < this->last_col_.GetM(); i++)
      Add(alpha, A.last_col_(i), this->last_col_(i));
    
    Add(alpha, A.last_block_, this->last_block_);
  }
  
  
  //! adds alpha*A to the current matrix
  template<class T, int d, int m> template<class T0, int d2, int m2>
  void TinyArrowMatrix<T, d, m>::Add_(const T0& alpha, const TinyArrowMatrix<T, d2, m2>& A)
  {
    // diagonals are added
    for (int i = 0; i <= 2*d2; i++)
      {
        int j = i + d - d2;
        for (int k = 0; k < this->diag_.GetM(); k++)
          this->diag_(k)(j) += alpha*A.diag_(k)(i);
      }
    
    // then last rows and columns
    for (int i = 0; i < m2; i++)
      {        
        int j = i + m - m2;
        for (int k = 0; k < last_row_.GetM(); k++)
          last_row_(k)(j) += alpha*A.last_row_(k)(i);
        
        for (int k = 0; k < last_col_.GetM(); k++)
          last_col_(k)(j) += alpha*A.last_col_(k)(i);
      }
    
    // for last block, we use the method AddInteraction
    int n = A.GetM();
    for (int i = n-m; i < n; i++)
      for (int j = n-m; j < n; j++)
        AddInteraction(i, j, alpha*A(i, j));
    
    // and small blocks in last_row/last_rows
    for (int i = 0; i < m-m2; i++)
      for (int j = 0; j < d2; j++)
        {
          int iglob = n-m2-1-i;
          int jglob = n-m-1-j;
          if ((iglob >= 0) && (jglob >= 0))
            {
              AddInteraction(iglob, jglob, alpha*A(iglob, jglob));
              AddInteraction(jglob, iglob, alpha*A(jglob, iglob));
            }
        }    
  }
  
  
  //! operation y = y + alpha*A*x
  template<class T, int d, int m> template<class T0, class T1>
  void TinyArrowMatrix<T, d, m>::
  MltAdd(const T0& alpha, const SeldonTranspose& trans, const Vector<T1>& x, Vector<T1>& y) const
  {
    // band part
    TinyBandMatrix<T, d>::MltAdd(alpha, trans, x, y);
    
    // last columns and rows now
    int n = this->diag_.GetM();
    T1 val, zero;
    SetComplexZero(zero);

    if (trans.NoTrans())
      {
        // y = y + alpha*A*x
        for (int j = 0; j < n; j++)
          {
            val = zero;
            TinyArrowMatrixLoop<m>::MltAdd(n, last_col_(j), x, val);
            y(j) += alpha*val;
            
            TinyArrowMatrixLoop<m>::MltTransAdd(n, last_row_(j), alpha*x(j), y);
          }
        
        // last block
        for (int p = 0; p < m; p++)
          for (int q = 0; q < m; q++)
            y(n+p) += alpha*last_block_(p, q)*x(n+q);
        
      }
    else if (trans.Trans())
      {
        // y = y + alpha*A^T*x
        for (int j = 0; j < n; j++)
          {
            val = zero;
            TinyArrowMatrixLoop<m>::MltAdd(n, last_row_(j), x, val);
            y(j) += alpha*val;
            
            TinyArrowMatrixLoop<m>::MltTransAdd(n, last_col_(j), alpha*x(j), y);
          }
        
        // last block
        for (int p = 0; p < m; p++)
          for (int q = 0; q < m; q++)
            y(n+q) += alpha*last_block_(p, q)*x(n+p);
      }
    else
      {
        // y = y + alpha*A^T*x
        for (int j = 0; j < n; j++)
          {
            val = zero;
            TinyArrowMatrixLoop<m>::MltConjAdd(n, last_row_(j), x, val);
            y(j) += alpha*val;
            
            TinyArrowMatrixLoop<m>::MltConjTransAdd(n, last_col_(j), alpha*x(j), y);
          }
        
        // last block
        for (int p = 0; p < m; p++)
          for (int q = 0; q < m; q++)
            y(n+q) += alpha*conjugate(last_block_(p, q))*x(n+p);
      }
  }
  
  
  //! LU resolution
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::Solve(Vector<T>& x)
  {
    int n = this->diag_.GetM();
    // solving L x = x
    
    // basic loops for first rows
    for (int j = 1; j < d; j++)
      for (int k = 0; k < j; k++)
	x(j) -= this->diag_(k)(k-j+d)*x(k);
    
    // then optimized loops
    for (int j = d; j < n; j++)
	TinyArrowMatrixLoop<d>::SolveLower(*this, j, x);
    
    // then treatment of last_row and last_block
    for (int j = 0; j < n; j++)
      TinyArrowMatrixLoop<m>::SolveRow(*this, j, n, x);
    
    for (int j = 0; j < m; j++)
      for (int k = 0; k < j; k++)
	x(n+j) -= x(n+k)*last_block_(j,k);
    
    
    // solving U x = x
    for (int j = m-1; j >= 0; j--)
      {
	for (int k = j+1; k < m; k++)
	  x(n+j) -= x(n+k)*last_block_(j,k);
	
	x(n+j) *= last_block_(j, j);
      }
    
    // basic loops for last rows
    for (int j = n-1; j >= n-d; j--)
      {
	for (int k = 0; k < m; k++)
	  x(j) -= last_col_(j)(k)*x(n+k);
	
	for (int k = j+1; k < n; k++)
	  x(j) -= this->diag_(k)(k-j+d)*x(k);
	
	// dividing by coefficient of diagonal
	x(j) *= this->diag_(j)(d);
      }
    
    // then optimized loops
    for (int j = n-d-1; j >= 0; j--)
      {
	TinyArrowMatrixLoop<m>::SolveColumn(*this, j, n, x);
	TinyArrowMatrixLoop<d>::SolveUpper(*this, j, x);
	
	// dividing by coefficient of diagonal
	x(j) *= this->diag_(j)(d);
      }
    
  }
  
  
  //! matrix is written in binary format
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyArrowMatrix::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();    
  }
  
  
  //! matrix is written in binary format
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::Write(ostream& FileStream) const
  {
    abort();
    //this->diag_.Write(FileStream);
    //this->last_col_.Write(FileStream);
    //this->last_row_.Write(FileStream);
    // this->last_block_.Write(FileStream);
  }
  
  
  //! matrix is written in ascii format
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyBandMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! matrix is written in ascii format
  template<class T, int d, int m>
  void TinyArrowMatrix<T, d, m>::WriteText(ostream& FileStream) const
  {
    TinyBandMatrix<T, d>::WriteText(FileStream);
    int n = last_col_.GetM();
    for (int j = 0; j < n; j++)
      for (int p = 0; p < m; p++)
	FileStream << j+1 << " " << n+p+1 << " " << last_col_(j)(p) << '\n';

    for (int j = 0; j < n; j++)
      for (int p = 0; p < m; p++)
	FileStream << n+p+1 << " " << j+1 << " " << last_row_(j)(p) << '\n';

    for (int j = 0; j < m; j++)
      for (int p = 0; p < m; p++)
	FileStream << n+j+1 << " " << n+p+1 << " " << last_block_(j, p) << '\n';
    
  }
  
  
  //! computes row_sum_i = \sum_j |a_{i,j}|
  template<class T, int d, int m>
  void GetRowSum(Vector<typename ClassComplexType<T>::Treal>& row_sum, const TinyArrowMatrix<T, d, m>& A)
  {
    row_sum.Reallocate(A.GetM());
    row_sum.Fill(0);
    int n = A.GetM();
    for (int i = 0; i < n-m; i++)
      {
        for (int  j = max(i-d,0); j < min(i+d+1,n-m); j++)
          row_sum(i) += abs(A(i, j));    
        
        for (int j = 0; j < m; j++)
          row_sum(i) += abs(A(i, j+n-m));
      }
    
    for (int i = n-m; i < n; i++)
      for (int j = 0; j < n; j++)
        row_sum(i) += abs(A(i, j));
  }
  
  
  //! scales matrix by coefficients scale (A is replaced by scale A)
  template<class T, int d, int m>
  void ScaleLeftMatrix(TinyArrowMatrix<T, d, m>& A, const Vector<typename ClassComplexType<T>::Treal>& scale)
  {
    int n = A.GetM();
    for (int i = 0; i < n-m; i++)
      {
        for (int  j = max(i-d,0); j < min(i+d+1,n-m); j++)
          A.Get(i, j) *= scale(i);    
        
        for (int j = 0; j < m; j++)
          A.Get(i, j+n-m) *= scale(i);
      }
    
    for (int i = n-m; i < n; i++)
      for (int j = 0; j < n; j++)
        A.Get(i, j) *= scale(i);
  }

  
  //! LU factorisation
  template<class T, int d, int m>
  void GetLU(TinyArrowMatrix<T, d, m>& A,
             TinyArrowMatrix<T, d, m>& mat_lu, bool keep_matrix)
  {
    mat_lu = A;
    if (!keep_matrix)
      A.Clear();
    
    mat_lu.Factorize();
  }
  

  //! LU factorisation
  template<class T, int d, int m>
  void GetLU(TinyArrowMatrix<T, d, m>& A)
  {
    A.Factorize();
  }
  
  
  //! resolution of A x = b, once GetLU has already been called
  template<class T, int d, int m>
  void SolveLU(TinyArrowMatrix<T, d, m>& mat_lu, Vector<T>& x)
  {
    mat_lu.Solve(x);
  }
  
  
  //! conversion from ArrayRowSparse to arrow matrix
  template<class T, int d, int m>
  void Copy(const Matrix<T, General, ArrayRowSparse>& A, TinyArrowMatrix<T, d, m>& B)
  {
    B.Copy(A);
  }

  
  //! B = B + alpha*A
  template<class T0, class T1, int d, int m, int d2, int m2>
  void Add(const T0& alpha, const TinyArrowMatrix<T1, d, m>& A,
           TinyArrowMatrix<T1, d2, m2>& B)
  {
    B.Add_(alpha, A);
  }
  
  
  //! y = beta*y + alpha*A*x
  template<class T0, class T, int d, int m, class T2>
  void MltAdd(const T0& alpha, const TinyArrowMatrix<T, d, m>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y)
  {
    T0 zero; SetComplexZero(zero);
    T2 czero; SetComplexZero(czero);
    T0 one; SetComplexOne(one);
    if (beta == zero)
      y.Fill(czero);
    else if (beta != one)
      Mlt(beta, y);
    
    A.MltAdd(alpha, SeldonNoTrans, x, y);
  }


  //! y = beta*y + alpha*A*x
  template<class T0, class T, int d, int m, class T2>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans, const TinyArrowMatrix<T, d, m>& A,
	      const Vector<T2>& x, const T0& beta, Vector<T2>& y)
  {
    T0 zero; SetComplexZero(zero);
    T2 czero; SetComplexZero(czero);
    T0 one; SetComplexOne(one);
    if (beta == zero)
      y.Fill(czero);
    else if (beta != one)
      Mlt(beta, y);
    
    A.MltAdd(alpha, trans, x, y);
  }


  //! y = A*x
  template<class T, int d, int m, class T2>
  void Mlt(const SeldonTranspose& trans, const TinyArrowMatrix<T, d, m>& A,
	   const Vector<T2>& x, Vector<T2>& y)
  {
    T2 zero; SetComplexZero(zero);
    T2 one; SetComplexOne(one);

    y.Fill(zero);
    A.MltAdd(one, trans, x, y);
  }


  //! y = A*x
  template<class T, int d, int m, class T2>
  void Mlt(const TinyArrowMatrix<T, d, m>& A,
	   const Vector<T2>& x, Vector<T2>& y)
  {
    T2 zero; SetComplexZero(zero);
    T2 one; SetComplexOne(one);

    y.Fill(zero);
    A.MltAdd(one, SeldonNoTrans, x, y);
  }
  

  //! A = A*alpha
  template<class T, int d, int m>
  void Mlt(const T& alpha, TinyArrowMatrix<T, d, m>& A)
  {
    A *= alpha;
  }

  
  //! displays matrix
  template<class T, int d, int m>
  ostream& operator<<(ostream& out, const TinyArrowMatrix<T, d, m>& A)
  {
    A.WriteText(out);
    return out;
  }


  //! copies a tiny arrow matrix into a current one
  /*!
    d2 and m2 have to be lower than d and m
   */
  template<class T, int d, int m> template<int d2, int m2>    
  const TinyArrowMatrix<T, d, m>& TinyArrowMatrix<T, d, m>::
  operator =(const TinyArrowMatrix<T, d2, m2>& A)
  {
    if ((d2 > d) || (m2 > m))
      {
        cout << "The profile of the matrix A is too large to be contained in this matrix" << endl;
        abort();
      }
    
    // initialisation of the matrix
    Reallocate(A.GetM(), A.GetN());
    Zero();
    
    // diagonals are copied
    for (int i = 0; i <= 2*d2; i++)
      {
        int j = i + d - d2;
        for (int k = 0; k < this->diag_.GetM(); k++)
          this->diag_(k)(j) = A.diag_(k)(i);
      }
    
    // then last rows and columns
    for (int i = 0; i < m2; i++)
      {        
        int j = i + m - m2;
        for (int k = 0; k < last_row_.GetM(); k++)
          last_row_(k)(j) = A.last_row_(k)(i);

        for (int k = 0; k < last_col_.GetM(); k++)
          last_col_(k)(j) = A.last_col_(k)(i);
      }
    
    // for last block, we use the method Set
    int n = A.GetM();
    for (int i = n-m; i < n; i++)
      for (int j = n-m; j < n; j++)
        Set(i, j, A(i, j));
    
    // and small blocks in last_row/last_rows
    for (int i = 0; i < m-m2; i++)
      for (int j = 0; j < d2; j++)
        {
          int iglob = n-m2-1-i;
          int jglob = n-m-1-j;
          if ((iglob >= 0) && (jglob >= 0))
            {
              Set(iglob, jglob, A(iglob, jglob));
              Set(jglob, iglob, A(jglob, iglob));
            }
        }

    return *this;
  }
    
}

#define SELDON_FILE_TINY_BAND_MATRIX_CXX
#endif

