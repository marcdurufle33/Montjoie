#ifndef MONTJOIE_FILE_GENERAL_SKYLINE_MATRIX_CXX

#include "GeneralSkyLineMatrix.hxx"

namespace Seldon
{
  
  //! initialization of a m x n matrix
  template<class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    val.Clear();
    first_col.Clear();
    
    val.Reallocate(i);
    first_col.Reallocate(i);
    first_col.Fill(-1);
    this->m_ = i;
    this->n_ = j;
  }
  
  
  //! adds x to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  AddInteraction(int i, int j, const T& x)
  {
    int size_row = val(i).GetM();
    if (first_col(i) == -1)
      {
        // new row
        ReallocateRow(i, 1, j);
        val(i)(0) = x;
      }
    else if (j < first_col(i))
      {
        T zero; SetComplexZero(zero);
        Vector<T> row(size_row + first_col(i) - j);
        row(0) = x;
        for (int k = 1; k < first_col(i)-j; k++)
          row(k) = zero;
        
        for (int k = first_col(i)-j; k < row.GetM(); k++)
          row(k) = val(i)(k - first_col(i) + j);
        
        val(i).SetData(row.GetM(), row.GetData());
        first_col(i) = j;
        row.Nullify();
      }
    else if (j >= first_col(i)+size_row)
      {
        T zero; SetComplexZero(zero);
        val(i).Resize(j-first_col(i) + 1);
        for (int k = size_row; k < val(i).GetM()-1; k++)
          val(i)(k) = zero;
        
        val(i)(j-first_col(i)) = x;
      }
    else
      val(i)(j-first_col(i)) += x;
  }
  

  //! adds x(k) to A(i, col(k)) for each k between 0 and nb-1
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T>& x, bool sorted)
  {
    if (nb <= 0)
      return;
    
    int jmin = col(0), jmax = col(0);
    for (int j = 1; j < nb; j++)
      {
        jmin = min(jmin, col(j));
        jmax = max(jmax, col(j));
      }
    
    // resizing row if needed
    int size_row = val(i).GetM();
    T zero; SetComplexZero(zero);
    if (first_col(i) == -1)
      {
        // new row
        ReallocateRow(i, jmax-jmin+1, jmin);
        val(i).Fill(zero);
      }
    else if ((jmin < first_col(i)) || (jmax >= first_col(i)+size_row))
      {
        int jbeg = first_col(i);
        int jend = first_col(i) + size_row - 1;
        int j0 = min(jmin, jbeg);
        Vector<T> row(max(jend, jmax) + 1 - j0);
        row.Fill(zero);
             
        for (int j = jbeg; j <= jend; j++)
          row(j-j0) = val(i)(j-jbeg);
        
        val(i).SetData(row.GetM(), row.GetData());
        first_col(i) = j0;
        row.Nullify();
      }
    
    // adding values in x
    for (int j = 0; j < nb; j++)
      val(i)(col(j) - first_col(i)) += x(j);
      
  }
  

  //! Returns the size used by the object in bytes.
  template <class T, class Prop, class Storage, class Allocator>
  size_t GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = val.GetMemorySize() + first_col.GetMemorySize();
    for (int i = 0; i < val.GetM(); i++)
      taille += val(i).GetMemorySize();
    
    return taille;
  }
  
    
  //! returns number of non-zero entries stored
  template <class T, class Prop, class Storage, class Allocator>
  int GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    int nnz = 0;
    for (int i = 0; i < val.GetM(); i++)
      nnz += val(i).GetM();
    
    return nnz;
  }
  
  
  //! sets all non-zero entries to a given value
  template <class T, class Prop, class Storage, class Allocator> template<class T0>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).Fill(x);
  }
  
  
  //! sets all non-zero entries to random values
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::FillRand()
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).FillRand();
  }
  
  
  //! sets all non-zero entries to 0
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).Zero();
  }
  
  
  //! Element access
  template <class T, class Prop, class Storage, class Allocator>
  const T GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::operator() (int i, int j) const
  {
    int size_row = val(i).GetM();
    if ((j >= first_col(i)) && (j < first_col(i)+size_row))
      {
        return val(i)(j - first_col(i));
      }
    
    return T(0);
  }
  
  
  //! writes the matrix in Matlab format
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  WriteText(const string& file_name) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("GeneralSkyLineMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif
    
    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! writes the matrix in Matlab format
  template <class T, class Prop, class Storage, class Allocator>
  void GeneralSkyLineMatrix<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream) const
  {
    for (int i = 0; i < this->m_; i++)
      {
        int icol = first_col(i);
        for (int j = 0; j < val(i).GetM(); j++)
          {
            FileStream << i+1 << ' ' << icol+1 << ' ' << val(i)(j) << '\n';
            icol++;
          }
      }
  }
  
  
  //! displays the matrix
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, RowSkyLine, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }
  
  
  //! Conversion from any matrix to skyline matrix
  template<class T, class Prop, class Storage,
           class Allocator, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
		  Matrix<T, General, RowSkyLine, Allocator2>& B)
  {
    int jbeg, jend;
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
        B.Clear();
        return;
      }
    
    B.Reallocate(m, n);
    for (int i = 0; i < A.GetM(); i++)
      {
	jbeg = 0;
	while ((jbeg < n) && (A(i,jbeg)==0))
	  jbeg++;
        
        jend = n-1;
	while ((jend >= 0) && (A(i,jend)==0))
	  jend--;
        
        B.ReallocateRow(i, jend-jbeg+1, jbeg);
        for (int j = jbeg; j <= jend; j++)
          B.Value(i, j-jbeg) = A(i, j);
      }
    
  }
  
  
  //! LU factorisation
  template <class T, class Allocator>
  void GetLU(Matrix<T, General, RowSkyLine, Allocator>& A)
  {    
    int n = A.GetM();
    Vector<T> Val(n);
    
    for (int i = 0; i < n; i++)
      {
        int pos = A.GetFirstColNumber(i);
        int size_Li = A.GetRowSize(i);
        if ((i < pos) || (i >= pos+size_Li))
          {
            cout << "No diagonal entry on row " << i << endl;
            abort();
          }
      }
    
    T pivot, one, zero;
    SetComplexZero(zero);
    SetComplexOne(one);
    Val.Fill(zero);
    // loop over rows
    for (int i = 0; i < n; i++)
      {
        int size_row = A.GetRowSize(i);
        int ibeg = A.GetFirstColNumber(i);
        int iend = ibeg + size_row - 1;
        int imax = iend;
        for (int k = 0; k < size_row; k++)
          Val(k) = A.Value(i, k);
        
        // loop over elements of L (on sub-diagonals)
        for (int j = ibeg; j < i; j++)
          {
            int jbeg = A.GetFirstColNumber(j);
            int size_jrow = A.GetRowSize(j);
            int jend = jbeg + size_jrow - 1;
            
            // pivot
            pivot = Val(j-ibeg)*A.Value(j, j-jbeg);
            
            // Li <- Li - pivot*Lj
            if (jend > imax)
              {
                for (int k = imax-ibeg+1; k <= jend-ibeg; k++)
                  Val(k) = zero;
                
                imax = jend;
              }
            
            for (int k = j+1; k <= jend; k++)
              Val(k-ibeg) -= pivot*A.Value(j, k-jbeg);
            
            // storing pivot
            Val(j - ibeg) = pivot;
          }
        
        if (imax != iend)
          A.ReallocateRow(i, imax-ibeg+1, ibeg);
        
        for (int k = 0; k <= imax-ibeg; k++)
          A.Value(i, k) = Val(k);
        
        // inverting diagonal
        pivot = Val(i - ibeg);
        if (pivot == zero)
          {
            cout << "non-invertible diagonal coefficient" << endl;
            abort();
          }
        
        A.Value(i, i - ibeg) = one/pivot;
      }
  }
  
  
  //! overwrites x by solution of L U y = x or (L U)^T y = x
  template<class T, class Allocator, class T2>
  void SolveLuVector(const SeldonTranspose& trans,
		     const Matrix<T, General, RowSkyLine, Allocator>& A,
		     Vector<T2>& x)
  {
    int n = A.GetM();

    if (trans.NoTrans())
      {
        // solving by L
        for (int i = 0; i < n; i++)
          {
            int ibeg = A.GetFirstColNumber(i);
            for (int j = ibeg; j < i; j++)
              x(i) -= A.Value(i, j - ibeg)*x(j);
          }
        
        // solving by U
        for (int i = n-1; i >= 0; i--)
          {
            int ibeg = A.GetFirstColNumber(i);
            int iend = ibeg + A.GetRowSize(i) - 1;            
            for (int j = i+1; j <= iend; j++)
              x(i) -= A.Value(i, j - ibeg)*x(j);
            
            x(i) *= A.Value(i, i-ibeg);
          }
      }
    else
      {
        // solving by U^T
        for (int i = 0; i < n; i++)
          {
            int ibeg = A.GetFirstColNumber(i);
            int iend = ibeg + A.GetRowSize(i) - 1;            
            x(i) *= A.Value(i, i-ibeg);
            for (int j = i+1; j <= iend; j++)
              x(j) -= A.Value(i, j - ibeg)*x(i);
            
          }

        // solving by L^T
        for (int i = n-1; i >= 0; i--)
          {
            int ibeg = A.GetFirstColNumber(i);
            for (int j = ibeg; j < i; j++)
              x(j) -= A.Value(i, j - ibeg)*x(i);
          }
      }
  }

  
  //! overwrites x by solution of L U y = x
  template<class T, class Allocator, class T2>
  void SolveLuVector(const Matrix<T, General, RowSkyLine, Allocator>& A,
		     Vector<T2>& x)
  {
    SolveLuVector(SeldonNoTrans, A, x);
  }
  

  //! y = A*x for skyline matrix
  template<class T1, class Allocator1, class T2,
           class Allocator2, class T4, class Allocator4>
  void MltVector(const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& x,
		 Vector<T4, VectFull, Allocator4>& y)
  {
    T4 zero; SetComplexZero(zero);    
    int n = A.GetM();
    T4 vloc;
    
    for (int i = 0; i < n; i++)
      {
	int irow = A.GetFirstColNumber(i);
	vloc = zero;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    vloc += A.Value(i, j)*x(irow);
	    irow++;
	  }
	
	y(i) = vloc; 
      }
  }

  
  //! y = beta*y + alpha*A*x for skyline matrix
  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y)
  {
    T3 zero; SetComplexZero(zero);
    if (beta == zero)
      y.Fill(zero);
    else
      Mlt(beta, y);
    
    int n = A.GetM();
    T4 vloc;
    
    T0 one; SetComplexOne(one);
    if (alpha == one)
      {
        for (int i = 0; i < n; i++)
          {
            int irow = A.GetFirstColNumber(i);
            vloc = zero;
            for (int j = 0; j < A.GetRowSize(i); j++)
              {
                vloc += A.Value(i, j)*x(irow);
                irow++;
              }
            
            y(i) += vloc; 
          }
      }
    else
      {
        for (int i = 0; i < n; i++)
          {
            int irow = A.GetFirstColNumber(i);
            vloc = zero;
            for (int j = 0; j < A.GetRowSize(i); j++)
              {
                vloc += A.Value(i, j)*x(irow);
                irow++;
              }
            
            y(i) += alpha*vloc; 
          }
      }
  }
  
  
  //! y = beta*y + alpha*A^T*x for skyline matrix
  template<class T0, class T1, class Allocator1, class T2,
           class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& trans,
		    const Matrix<T1, General, RowSkyLine, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& x,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& y)
  {
    if (trans.NoTrans())
      {
	MltAddVector(alpha, A, x, beta, y);
	return;
      }
    
    T3 zero; SetComplexZero(zero);
    if (beta == zero)
      y.Fill(zero);
    else
      Mlt(beta, y);
    
    int n = A.GetM();
    
    T0 one; SetComplexOne(one);
    if (trans.Trans())
      {
	if (alpha == one)
	  {
	    for (int i = 0; i < n; i++)
	      {
		int irow = A.GetFirstColNumber(i);
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    y(irow) += A.Value(i, j)*x(i);
		    irow++;
		  }
	      }
	  }
	else
	  {
	    for (int i = 0; i < n; i++)
	      {
		int irow = A.GetFirstColNumber(i);
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    y(irow) += alpha*A.Value(i, j)*x(i);
		    irow++;
		  }
	      }
	  }
      }
    else
      {
	if (alpha == one)
	  {
	    for (int i = 0; i < n; i++)
	      {
		int irow = A.GetFirstColNumber(i);
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    y(irow) += conjugate(A.Value(i, j))*x(i);
		    irow++;
		  }
	      }
	  }
	else
	  {
	    for (int i = 0; i < n; i++)
	      {
		int irow = A.GetFirstColNumber(i);
		for (int j = 0; j < A.GetRowSize(i); j++)
		  {
		    y(irow) += alpha*conjugate(A.Value(i, j))*x(i);
		    irow++;
		  }
	      }
	  }
      }
  }

}
  
#define MONTJOIE_FILE_GENERAL_SKYLINE_MATRIX_CXX
#endif
