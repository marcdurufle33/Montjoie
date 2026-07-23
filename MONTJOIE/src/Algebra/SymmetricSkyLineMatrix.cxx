#ifndef SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_CXX

#include "SymmetricSkyLineMatrix.hxx"

namespace Seldon
{
  
  //! setting indices and values
  template <class T, class Prop, class Storage, class Allocator>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::
  SetData(IVect& ind_, Vector<T>& val_)
  {
    ind.Clear();
    val.Clear();
    
    ind.SetData(ind_.GetM(), ind_.GetData());
    val.SetData(val_.GetM(), val_.GetData());
    
    this->m_ = ind_.GetM()-1;
    this->n_ = this->m_;

    ind_.Nullify();
    val_.Nullify();
  }
  
  
  //! sets all non-zero entries to a given value
  template <class T, class Prop, class Storage, class Allocator> template<class T0>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    val.Fill(x);
  }
  
  
  //! sets all non-zero entries to random values
  template <class T, class Prop, class Storage, class Allocator>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::FillRand()
  {
    val.FillRand();
  }
  
  
  //! sets all non-zero entries to 0
  template <class T, class Prop, class Storage, class Allocator>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::Zero()
  {
    val.Zero();
  }
  
  
  //! Element access
  template <class T, class Prop, class Storage, class Allocator>
  T SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::operator() (int i, int j) const
  {
    T zero; SetComplexZero(zero);
    if (j <= i)
      {
	int jn = i+1-(ind(i+1)-ind(i));
	if (j >= jn)
          return val(ind(i)+j-jn);
        else
	  return zero;
      }
    else
      {
	int in = j+1-(ind(j+1)-ind(j));
	if (i >= in)
	  return val(ind(j)+i-in);
	else
	  return zero;
      }
    
    return zero;
  }
  
  
  //! writes the matrix in Matlab format
  template <class T, class Prop, class Storage, class Allocator>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::
  WriteText(const string& file_name) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("SymmetricSkyLineMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! writes the matrix in Matlab format
  template <class T, class Prop, class Storage, class Allocator>
  void SymmetricSkyLineMatrix<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream) const
  {
    for (int i = 0; i < ind.GetM()-1; i++)
      {
        int size_col = ind(i+1) - ind(i);
        int irow = i - size_col + 1;
        for (int j = ind(i); j < ind(i+1); j++)
          {
            FileStream << irow+1 << ' ' << i+1 << ' ' << val(j) << '\n';
            if (irow != i)
              FileStream << i+1 << ' ' << irow+1 << ' ' << val(j) << '\n';
            
            irow++;
          }
      }
  }
  
  
  //! displays the matrix
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, SymColSkyLine, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }
  
  
  //! Conversion from any matrix to skyline matrix
  template<class T, class Prop, class Storage,
           class Allocator, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
		  Matrix<T, Symmetric, SymColSkyLine, Allocator2>& B)
  {
    T zero; SetComplexZero(zero);
    Vector<T> val; Vector<int> ind;
    int nb, jn;

    nb = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
	jn = 0;
	while ((jn <= i) && (A(i, jn) == zero))
	  jn++;
        
        if (i >= jn)
          nb += i - jn + 1;
        
        if (nb < 0)
          {
            cout << "Matrix too large to be converted" << endl;
            abort();
          }
      }
    
    val.Reallocate(nb);
    val.Fill(zero);
    ind.Reallocate(A.GetM()+1);
    ind.Fill(0);
    
    nb = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
	jn = 0;
	while ((jn <= i) && (A(i, jn) == zero))
	  jn++;
	
	ind(i+1) = ind(i) + i - jn +1;

	for (int j = jn; j <= i; j++)
	  val(nb++) = A(i, j);
      }
    
    B.SetData(ind, val);
  }
  
  
  //! Cholesky factorisation
  template <class T, class Allocator>
  void GetCholesky(Matrix<T, Symmetric, SymColSkyLine, Allocator>& A)
  {
    T* val = A.GetData();
    int* ind = A.GetInd();
    
    T S;
    int Li, Lj, Lk, N = A.GetM();
    // fisrt diagonal element 
    val[ind[1]-1] = sqrt(val[ind[1]-1]);
    
    // fisrt column
    for (int i = 1; i < N; i++)
      {
	Li = ind[i+1]-ind[i];
	if (Li == i+1)
	  val[ind[i]] /= val[ind[1]-1];
      }
    
    // other columns
    for (int k = 1; k < N-1; k++)
      {
	Lk = ind[k+1]-ind[k];
	// diagonal term
	SetComplexZero(S);
        int off_k = ind[k] - (k-Lk) - 1;
	for (int i = k-Lk+1; i < k; i++)
          S += val[off_k + i]*val[off_k + i];
	
	val[ind[k+1] - 1] =  sqrt(val[ind[k+1]-1] - S);
	
	// Extra-diagonal terms
 	for (int j = k+1; j < N; j++)
	  {
	    Lj = ind[j+1]-ind[j];
	    if (Lj > j-k)
	      {
		SetComplexZero(S);
                int off_j = ind[j] - (j-Lj) - 1;
		for (int i = max(k-Lk+1, j-Lj+1); i < k; i++)
                  S += val[off_k + i]*val[off_j + i];
		
		val[off_j + k] = (val[off_j + k]-S)/val[ind[k+1]-1];
	      }
	  }
      }
    
    // last diagonal term
    Lk = ind[N]-ind[N-1];
    SetComplexZero(S);
    int off_n = ind[N-1] - (N-Lk);
    for (int i = N-Lk; i < N-1; i++)
      S += val[off_n + i]*val[off_n + i];
    
    val[ind[N]-1] =  sqrt(val[ind[N]-1] - S);
  }
  
  
  //! overwrites x by solution of L y = x or L^T y = x
  template<class T, class Allocator, class Alloc2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const Matrix<T, Symmetric, SymColSkyLine, Allocator>& A,
                     Vector<T, VectFull, Alloc2>& x)
  {
    T* val = A.GetData();
    int* ind = A.GetInd();

    int Lk, N = A.GetM();
    if (trans.NoTrans())
      {
        // descente
        for (int k = 0; k < N; k++)
          {
            Lk = ind[k+1]-ind[k];
            int offset = ind[k] - (k-Lk) - 1;
            for (int i = k-Lk+1; i < k; i++)
              x(k) -= x(i)*val[offset + i];
            
            x(k) /= val[ind[k+1]-1];
          }
      }
    else
      {    
        // remontee
        for (int k = N-1; k >= 0; k--)
          {
            x(k) /= val[ind[k+1]-1];
            Lk = ind[k+1]-ind[k];
            int offset = ind[k] - (k-Lk) - 1;
            for (int i = k-Lk+1; i < k; i++)
              x(i) -= x(k)*val[offset + i];
          }
      }
  }

  
  //! overwrites x with L x or L^T x
  template<class T, class Allocator, class Alloc2>
  void MltCholesky(const SeldonTranspose& trans,
                   const Matrix<T, Symmetric, SymColSkyLine, Allocator>& A,
                   Vector<T, VectFull, Alloc2>& x)
  {
    T* val = A.GetData();
    int* ind = A.GetInd();

    int Lk, N = A.GetM();
    if (trans.NoTrans())
      {
        // descente
        for (int k = N-1; k >= 0; k--)
          {
            Lk = ind[k+1]-ind[k];
            int offset = ind[k] - (k-Lk) - 1;
            x(k) *= val[ind[k+1]-1];
            for (int i = k-Lk+1; i < k; i++)
              x(k) += x(i)*val[offset + i];
          }
      }
    else
      {    
        // remontee
        for (int k = 0; k < N; k++)
          {            
            Lk = ind[k+1]-ind[k];
            int offset = ind[k] - (k-Lk) - 1;
            for (int i = k-Lk+1; i < k; i++)
              x(i) += x(k)*val[offset + i];
            
            x(k) *= val[ind[k+1]-1];
          }
      }
  }
  
  
  //! LDL^T factorisation
  template <class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, SymColSkyLine, Allocator>& A)
  {
    T* val = A.GetData();
    int* ind = A.GetInd();
    
    T S;
    int Li, Lj, Lk, N = A.GetM();
    T one; 
    SetComplexOne(one);
    
    Vector<T> invDiag(N);
    // fisrt diagonal element 
    invDiag(0) = one/val[ind[1]-1];
    
    // fisrt column
    for (int i = 1; i < N; i++)
      {
	Li = ind[i+1]-ind[i];
	if (Li == i+1)
	  val[ind[i]] *= invDiag(0);
      }
    
    // other columns
    for (int k = 1; k < N-1; k++)
      {
	Lk = ind[k+1]-ind[k];
	// diagonal term
	SetComplexZero(S);
        int off_k = ind[k] - (k-Lk) - 1;
	for (int i = k-Lk+1; i < k; i++)
          S += val[off_k + i]*val[off_k + i]*val[ind[i+1]-1];
        
	val[ind[k+1] - 1] -=  S;
        invDiag(k) = one/val[ind[k+1]-1];
	
	// Extra-diagonal terms
 	for (int j = k+1; j < N; j++)
	  {
	    Lj = ind[j+1]-ind[j];
	    if (Lj > j-k)
	      {
		SetComplexZero(S);
                int off_j = ind[j] - (j-Lj) - 1;
		for (int i = max(k-Lk+1, j-Lj+1); i < k; i++)
                  S += val[off_k + i]*val[off_j + i]*val[ind[i+1]-1];
		
		val[off_j + k] = (val[off_j + k]-S)*invDiag(k);
	      }
	  }
      }
    
    // last diagonal term
    Lk = ind[N]-ind[N-1];
    SetComplexZero(S);
    int off_n = ind[N-1] - (N-Lk);
    for (int i = N-Lk; i < N-1; i++)
      S += val[off_n+i]*val[off_n+i]*val[ind[i+1]-1];;
    
    val[ind[N]-1] -= S;
    invDiag(N-1) = one/val[ind[N]-1];
    
    // putting inverse of diagonal on diagonal of A
    for (int i = 0; i < N; i++)
      val[ind[i+1]-1] = invDiag(i);
    
  }
  
  
  //! overwrites x by solution of LDL^T y = x
  template<class T, class Allocator, class T2>
  void SolveLuVector(const SeldonTranspose& trans,
		     const Matrix<T, Symmetric, SymColSkyLine, Allocator>& A,
		     Vector<T2>& x)
  {
    T* val = A.GetData();
    int* ind = A.GetInd();

    int Lk, N = A.GetM();
    // solving by L
    for (int k = 0; k < N; k++)
      {
        Lk = ind[k+1]-ind[k];
        int offset = ind[k] - (k-Lk) - 1;
        for (int i = k-Lk+1; i < k; i++)
          x(k) -= x(i)*val[offset + i];
      }
    
    // solving by D
    for (int k = 0; k < N; k++)
      x(k) *= val[ind[k+1]-1];
    
    // solving by L^T
    for (int k = N-1; k >= 0; k--)
      {
        Lk = ind[k+1]-ind[k];
        int offset = ind[k] - (k-Lk) - 1;
        for (int i = k-Lk+1; i < k; i++)
          x(i) -= x(k)*val[offset + i];
      }
  }
  

  //! y = beta*y + alpha*A*x for symmetric skyline matrix
  template<class T1, class Allocator1, class T2,
           class Allocator2, class T4, class Allocator4>
  void MltVector(const Matrix<T1, Symmetric, SymColSkyLine, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& x,
		 Vector<T4, VectFull, Allocator4>& y)
  {
    T4 czero; SetComplexZero(czero);
    y.Fill(czero);
    
    int n = A.GetM();
    int* ind = A.GetInd();
    T1* val = A.GetData();
    
    for (int i = 0; i < n; i++)
      {
	int size_col = ind[i+1] - ind[i];
	int irow = i - size_col + 1;
	for (int j = ind[i]; j < ind[i+1]-1; j++)
	  {
	    y(i) += val[j]*x(irow);
	    y(irow) += val[j]*x(i);
	    irow++;
	  }
	
	y(i) += val[ind[i+1]-1]*x(i);
      }
  }

  
  //! y = beta*y + alpha*A*x for symmetric skyline matrix
  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, Symmetric, SymColSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y)
  {
    T3 zero; SetComplexZero(zero);
    T4 czero; SetComplexZero(czero);
    T4 one; SetComplexOne(one);
    if (beta == zero)
      y.Fill(czero);
    else
      Mlt(beta, y);
    
    int n = A.GetM();
    int* ind = A.GetInd();
    T1* val = A.GetData();
    
    if (alpha == one)
      {
        for (int i = 0; i < n; i++)
          {
            int size_col = ind[i+1] - ind[i];
            int irow = i - size_col + 1;
            for (int j = ind[i]; j < ind[i+1]-1; j++)
              {
                y(i) += val[j]*x(irow);
                y(irow) += val[j]*x(i);
                irow++;
              }
            
            y(i) += val[ind[i+1]-1]*x(i);
          }
      }
    else
      {
        for (int i = 0; i < n; i++)
          {
            int size_col = ind[i+1] - ind[i];
            int irow = i - size_col + 1;
            for (int j = ind[i]; j < ind[i+1]-1; j++)
              {
                y(i) += alpha*val[j]*x(irow);
                y(irow) += alpha*val[j]*x(i);
                irow++;
              }
            
            y(i) += alpha*val[ind[i+1]-1]*x(i);
          }
      }
  }

}

#define SELDON_FILE_SYMMETRIC_SKYLINE_MATRIX_CXX
#endif
