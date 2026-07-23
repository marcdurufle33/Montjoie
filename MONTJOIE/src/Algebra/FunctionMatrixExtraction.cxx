#ifndef SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_CXX

#include "Algebra/FunctionMatrixExtraction.hxx"

namespace Seldon
{

#ifdef SELDON_WITH_MPI
  void ExtractDistributedInfo(const DistributedMatrixIntegerArray& A,
			      const Vector<int>& num, const Vector<int>& index,
			      DistributedMatrixIntegerArray& B)
  {
    int nb_proc; MPI_Comm_size(A.comm, &nb_proc);
    if (nb_proc == 1)
      {
	B.nloc = num.GetM();
	B.nglob = num.GetM();
	B.nb_unknowns_scal = 1;
	B.nodl_scalar = num.GetM();
	B.comm = A.comm;
	
	return;
      }
	
    Vector<bool> IsOverlappedRow(A.nloc); Vector<int> proc_overlap(A.nloc);
    IsOverlappedRow.Fill(false);
    proc_overlap.Fill(-1);
    for (int i = 0; i < A.OverlapRowNumbers.GetM(); i++)
      {
	IsOverlappedRow(A.OverlapRowNumbers(i)) = true;
	proc_overlap(A.OverlapRowNumbers(i)) = A.OverlapProcNumbers(i);
      }
    
    // local number of rows for B
    B.nloc = num.GetM();

    // we count the number of original rows (owned by the processor)
    int ndof_orig = 0;
    B.GlobalRowNumbers.Reallocate(B.nloc);
    B.GlobalRowNumbers.Fill(-1);
    for (int i = 0; i < num.GetM(); i++)
      if (!IsOverlappedRow(num(i)))
	{
	  B.GlobalRowNumbers(i) = ndof_orig;
	  ndof_orig++;
	}

    // the root processor receives the number of rows and original rows for each proc
    Vector<int> nb_rows_per_proc;
    int rank; MPI_Comm_rank(A.comm, &rank);
    if (rank == 0)
      nb_rows_per_proc.Reallocate(2*nb_proc);
    
    int nb_row[2];
    nb_row[0] = num.GetM(); nb_row[1] = ndof_orig;
    MPI_Gather(nb_row, 2, MPI_INTEGER, nb_rows_per_proc.GetData(), 2, MPI_INTEGER, 0, A.comm);

    // the root processor counts the number of global rows
    int nglob = 0;
    Vector<int> new_proc_num, offset_proc;
    new_proc_num.Reallocate(nb_proc+1);
    if (rank == 0)
      {
	for (int i = 0; i < nb_proc; i++)
	  nglob += nb_rows_per_proc(2*i+1);

	// new_proc_num contains the new processor number (in the MPI communicator of B)
	new_proc_num.Fill(-1);
	new_proc_num(nb_proc) = nglob;
	int np = 0;
	for (int i = 0; i < nb_proc; i++)
	  if (nb_rows_per_proc(2*i) > 0)
	    {
	      new_proc_num(i) = np;
	      np++;
	    }

	// offset_proc contains the offset to obtain global row numbers for each proc
	offset_proc.Reallocate(nb_proc);
	offset_proc(0) = 0;
	for (int i = 0; i < nb_proc-1; i++)
	  offset_proc(i+1) = offset_proc(i) + nb_rows_per_proc(2*i+1);
      }

    // offset_proc is scattered to other processors
    int offset_row;
    MPI_Scatter(offset_proc.GetData(), 1, MPI_INTEGER, &offset_row, 1, MPI_INTEGER, 0, A.comm);

    // global row numbers are updated
    for (int i = 0; i < num.GetM(); i++)
      if (!IsOverlappedRow(num(i)))
	B.GlobalRowNumbers(i) += offset_row;

    // new_proc_num is broadcasted to all processors
    MPI_Bcast(new_proc_num.GetData(), nb_proc+1, MPI_INTEGER, 0, A.comm);

    // global number of rows
    B.nglob = new_proc_num(nb_proc);

    // we create the communicator for B
    int rank_B = new_proc_num(rank);
    int color = MPI_UNDEFINED;
    if (rank_B >= 0)
      color = 0;
    
    MPI_Comm_split(A.comm, color, rank_B, &B.comm);

    // setting nodl_scalar and nb_unknowns_scal
    B.nb_unknowns_scal = A.nb_unknowns_scal;
    B.nodl_scalar = 0;
    for (int i = 0; i < num.GetM(); i++)
      if (num(i) < A.nodl_scalar)
	B.nodl_scalar++;

    // overlapped rows for B
    int noverlap = 0;
    for (int i = 0; i < num.GetM(); i++)
      if (IsOverlappedRow(num(i)))
	noverlap++;

    if (num.GetM() == 0)
      {
	B.OverlapRowNumbers.Clear();
	B.OverlapProcNumbers.Clear();
	B.ProcSharingRows.Clear();
	B.SharingRowNumbers.Clear();
      }
    else
      {
	B.OverlapRowNumbers.Reallocate(noverlap);
	B.OverlapProcNumbers.Reallocate(noverlap);
	noverlap = 0;
	for (int i = 0; i < num.GetM(); i++)
	  if (IsOverlappedRow(num(i)))
	    {
	      B.OverlapRowNumbers(noverlap) = i;
	      B.OverlapProcNumbers(noverlap) = new_proc_num(proc_overlap(num(i)));
	      noverlap++;
	    }
	
	// shared rows of B
	int nb_shared_proc = 0;
	Vector<int> nb_shared_rows(A.ProcSharingRows.GetM());
	nb_shared_rows.Zero();
	
	for (int i = 0; i < A.ProcSharingRows.GetM(); i++)
	  {	
	    for (int j = 0; j < A.SharingRowNumbers(i).GetM(); j++)
	      {
		int jloc = index(A.SharingRowNumbers(i)(j));
		if (jloc >= 0)
		  nb_shared_rows(i)++;
	      }
	    
	    if (nb_shared_rows(i) > 0)
	      nb_shared_proc++;
	  }
	
	B.ProcSharingRows.Reallocate(nb_shared_proc);
	B.SharingRowNumbers.Reallocate(nb_shared_proc);
	
	nb_shared_proc = 0;
	for (int i = 0; i < A.ProcSharingRows.GetM(); i++)
	  {
	    int ib = nb_shared_proc;
	    if (nb_shared_rows(i) > 0)
	      {
		B.ProcSharingRows(ib) = new_proc_num(A.ProcSharingRows(i));
		B.SharingRowNumbers(ib).Reallocate(nb_shared_rows(i));
		nb_shared_rows(i) = 0;
		for (int j = 0; j < A.SharingRowNumbers(i).GetM(); j++)
		  {
		    int jloc = index(A.SharingRowNumbers(i)(j));
		    if (jloc >= 0)
		      {
			B.SharingRowNumbers(ib)(nb_shared_rows(i)) = jloc;
			nb_shared_rows(i)++;
		      }
		  }
		
		nb_shared_proc++;
	      }
	  }
	
	// for GlobalRowNumbers, we recover the dof numbers for shared dofs
	AssembleVector(B.GlobalRowNumbers, MPI_MAX, B.ProcSharingRows,
		       B.SharingRowNumbers, B.comm, B.nodl_scalar, B.nb_unknowns_scal, 17);
      }
  }

  
  template<class T>
  void ExtractDistributedInfo(const DistributedMatrix_Base<T>& A,
			      const Vector<int>& num, const Vector<int>& index,
			      DistributedMatrixIntegerArray& B)
  {
    DistributedMatrixIntegerArray Ainfo;

    int nb_proc; MPI_Comm_size(A.GetCommunicator(), &nb_proc);
    if (nb_proc > 1)
      Ainfo.SetData(A.GetLocalM(), A.GetGlobalM(), A.GetNodlScalar(), A.GetNbScalarUnknowns(),
		    A.GetCommunicator(), A.GetOverlapRowNumber(), A.GetOverlapProcNumber(),
		    A.GetGlobalRowNumber(), A.GetProcessorSharingRows(), A.GetSharingRowNumbers());
    else
      Ainfo.comm = A.GetCommunicator();
    
    ExtractDistributedInfo(Ainfo, num, index, B);
    
    if (nb_proc > 1)
      Ainfo.Nullify();
  }
#endif

  
  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSymSparse, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    int m = col_number.GetM();
    V.Reallocate(m);
    // counting size of columns
    IVect size_col(m);
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    size_col.Fill(0); int num, k;
    for (int i = 0; i < A.GetM(); i++)
      {
	if (index(i) != -1)
	  for (int j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i,j) != i)
	      size_col(index(i))++;
	
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    num = index(A.Index(i,j));
	    if (num != -1)
	      size_col(num)++;
	  }
      }
    
    // now we fill columns
    for (int i = 0; i < m; i++)
      {
	V.Index(i) = col_number(i);
	V.Value(i).Reallocate(size_col(i));
      }
    
    size_col.Zero();
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index(A.Index(i, j)) != -1)
	    {
	      num = index(A.Index(i, j));
	      k = size_col(num);
	      V.Value(num).Index(k) = i;
	      V.Value(num).Value(k) = A.Value(i, j);
	      size_col(num)++;
	    }
	
	if (index(i) != -1)
	  for (int j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) != i)
	      {
		num = index(i);
		k = size_col(num);
		V.Value(num).Index(k) = A.Index(i,j);
		V.Value(num).Value(k) = A.Value(i, j);
		size_col(num)++;
	      }
      }
  }


  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSparse, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    int m = col_number.GetM();
    V.Reallocate(m);
    // counting size of columns
    IVect size_col(m);
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    size_col.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	if (index(A.Index(i, j)) != -1)
	  size_col(index(A.Index(i, j)))++;
    
    // now we fill columns
    for (int i = 0; i < m; i++)
      {
	V.Index(i) = col_number(i);
	V.Value(i).Reallocate(size_col(i));
      }
    
    size_col.Zero();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	if (index(A.Index(i, j)) != -1)
	  {
	    int num = index(A.Index(i, j));
	    int k = size_col(num);
	    V.Value(num).Index(k) = i;
	    V.Value(num).Value(k) = A.Value(i, j);
	    size_col(num)++;
	  }
  }
  

  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    int m = col_number.GetM();
    V.Reallocate(m);
    // counting size of columns
    IVect size_col(m);
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    size_col.Fill(0); int num, k;
    for (int i = 0; i < A.GetM(); i++)
      {
	if (index(i) != -1)
	  for (int j = 0; j < A.GetRealRowSize(i); j++)
	    if (A.IndexReal(i,j) != i)
	      size_col(index(i))++;
	
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    num = index(A.IndexReal(i,j));
	    if (num != -1)
	      size_col(num)++;
	  }

	if (index(i) != -1)
	  for (int j = 0; j < A.GetImagRowSize(i); j++)
	    if (A.IndexImag(i,j) != i)
	      size_col(index(i))++;
	
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  {
	    num = index(A.IndexImag(i,j));
	    if (num != -1)
	      size_col(num)++;
	  }
      }
    
    // now we fill columns
    for (int i = 0; i < m; i++)
      {
	V.Index(i) = col_number(i);
	V.Value(i).Reallocate(size_col(i));
      }
    
    size_col.Zero();
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  if (index(A.IndexReal(i, j)) != -1)
	    {
	      num = index(A.IndexReal(i, j));
	      k = size_col(num);
	      V.Value(num).Index(k) = i;
	      V.Value(num).Value(k) = T1(A.ValueReal(i, j), 0);
	      size_col(num)++;
	    }
	
	if (index(i) != -1)
	  for (int j = 0; j < A.GetRealRowSize(i); j++)
	    if (A.IndexReal(i, j) != i)
	      {
		num = index(i);
		k = size_col(num);
		V.Value(num).Index(k) = A.IndexReal(i,j);
		V.Value(num).Value(k) = T1(A.ValueReal(i, j), 0);
		size_col(num)++;
	      }

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  if (index(A.IndexImag(i, j)) != -1)
	    {
	      num = index(A.IndexImag(i, j));
	      k = size_col(num);
	      V.Value(num).Index(k) = i;
	      V.Value(num).Value(k) = T1(0, A.ValueImag(i, j));
	      size_col(num)++;
	    }
	
	if (index(i) != -1)
	  for (int j = 0; j < A.GetImagRowSize(i); j++)
	    if (A.IndexImag(i, j) != i)
	      {
		num = index(i);
		k = size_col(num);
		V.Value(num).Index(k) = A.IndexImag(i,j);
		V.Value(num).Value(k) = T1(0, A.ValueImag(i, j));
		size_col(num)++;
	      }
      }
    
    // assembling vector to have sorted row numbers
    for (int i = 0; i < m; i++)
      V.Value(i).Assemble();
  }


  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */  
  template<class T1, class Prop, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T1, Prop, ArrayRowComplexSparse, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    int m = col_number.GetM();
    V.Reallocate(m);
    // counting size of columns
    IVect size_col(m);
    IVect index(A.GetM()); index.Fill(-1);
    for (int i = 0; i < m; i++)
      index(col_number(i)) = i;
    
    size_col.Fill(0);
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  if (index(A.IndexReal(i, j)) != -1)
	    size_col(index(A.IndexReal(i, j)))++;

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  if (index(A.IndexImag(i, j)) != -1)
	    size_col(index(A.IndexImag(i, j)))++;
      }
    
    // now we fill columns
    for (int i = 0; i < m; i++)
      {
	V.Index(i) = col_number(i);
	V.Value(i).Reallocate(size_col(i));
      }
    
    size_col.Zero();
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  if (index(A.IndexReal(i, j)) != -1)
	    {
	      int num = index(A.IndexReal(i, j));
	      int k = size_col(num);
	      V.Value(num).Index(k) = i;
	      V.Value(num).Value(k) = T1(A.ValueReal(i, j), 0);
	      size_col(num)++;
	    }

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  if (index(A.IndexImag(i, j)) != -1)
	    {
	      int num = index(A.IndexImag(i, j));
	      int k = size_col(num);
	      V.Value(num).Index(k) = i;
	      V.Value(num).Value(k) = T1(0, A.ValueImag(i, j));
	      size_col(num)++;
	    }
      }
    
    // assembling vector to have sorted row numbers
    for (int i = 0; i < m; i++)
      V.Value(i).Assemble();
  }
  

  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */  
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, BlockDiagRowSym, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    abort();
  }
  
  
  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, BlockDiagRow, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    abort();
  }
  
  
  //! extracts several columns of a sparse matrix
  /*!
    \param[in] A sparse matrix where columns are extracted
    \param[in] col_number numbers of the columns to be extracted
    \param[out] V extracted sparse columns    
   */
  template<class T, class Prop, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void GetCol(const Matrix<T, Prop, DiagonalRow, Allocator>& A, const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>, VectSparse, Allocator3>& V)
  {
    int m = col_number.GetM();
    if (m <= 0)
      {
        V.Clear();
        return;
      }
    
    V.Reallocate(m);
    for (int ic = 0; ic < m; ic++)
      {
        int i = col_number(ic);
        V.Index(0) = i;
        V.Value(0).Reallocate(1);
        V.Value(0).Index(0) = i;
        V.Value(0).Value(0) = A(i, i);
      }
  }
  
  
  //! clears several columns of a sparse matrix
  /*!
    \param[in] col_number numbers of the columns to be cleared
    \param[inout] A sparse matrix where columns are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, DiagonalRow, Allocator>& A)
  {
    for (int ic = 0; ic < col_number.GetM(); ic++)
      {
        int i = col_number(ic);
        A.Get(i, i) = 0.0;
      }
  }
   

  //! clears several rows of a sparse matrix
  /*!
    \param[in] col_number numbers of the rows to be cleared
    \param[inout] A sparse matrix where rows are erased
   */
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, DiagonalRow, Allocator>& A)
  {
    EraseCol(col_number, A);
  }
  

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BandedCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM();
    int n = A.GetN();
    Vector<bool> SelectedCol(n);
    SelectedCol.Fill(false);
    for (int i = 0; i < col_number.GetM(); i++)
      SelectedCol(col_number(i)) = true;
    
    T1 zero; SetComplexZero(zero);
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
	  if (SelectedCol(j))
	    A.Set(i, j, zero);
        }
  }
  
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& row_number,
		Matrix<T1, Prop, BandedCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    //int m = A.GetM();
    int n = A.GetN();
    T1 zero; SetComplexZero(zero);
    for (int i1 = 0; i1 < row_number.GetM(); i1++)
      {
	int i = row_number(i1);
	for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
	  {
	    int j = i + k;
	    A.Set(i, j, zero);
	  }
      }
  }


  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRow, Allocator>& A)
  {
    for (int ic = 0; ic < col_number.GetM(); ic++)
      {
        int i = col_number(ic);
	A.ClearColumn(i);
      }
  }

  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRow, Allocator>& A)
  {
    for (int ic = 0; ic < col_number.GetM(); ic++)
      {
        int i = col_number(ic);
	A.ClearRow(i);
      }
  }
  

  template<class T1, class Prop, class Allocator>
  void EraseCol(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRowSym, Allocator>& A)
  {
    for (int ic = 0; ic < col_number.GetM(); ic++)
      {
        int i = col_number(ic);
	A.ClearColumn(i);
      }
  }
  
  
  template<class T1, class Prop, class Allocator>
  void EraseRow(const IVect& col_number,
		Matrix<T1, Prop, BlockDiagRowSym, Allocator>& A)
  {
    for (int ic = 0; ic < col_number.GetM(); ic++)
      {
        int i = col_number(ic);
	A.ClearRow(i);
      }
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale_left vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[in] mat given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowSum(Vector<T>& diagonal_scale_left,
		 const Matrix<Complexe, General, ArrowCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM() - A.GetNbLastRow();
    int n = A.GetN() - A.GetNbLastCol(); 
    diagonal_scale_left.Reallocate(A.GetM());
    diagonal_scale_left.Fill(0);
    
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          diagonal_scale_left(i) += abs(A(i, j));
        }
    
    for (int i = 0; i < A.GetNbLastRow(); i++)
      for (int j = 0; j < A.GetN(); j++)
        diagonal_scale_left(m+i) += abs(A(m+i, j));
    
    for (int j = 0; j < A.GetNbLastCol(); j++)
      for (int i = 0; i < m; i++)
        diagonal_scale_left(i) += abs(A(i, n+j));    
  }
  

  //! For each column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] diagonal_scale vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] mat given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetColSum(Vector<T>& diagonal_scale,
		 const Matrix<Complexe, General, ArrowCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM() - A.GetNbLastRow();
    int n = A.GetN() - A.GetNbLastCol(); 
    diagonal_scale.Reallocate(A.GetN());
    diagonal_scale.Fill(0);
    
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          diagonal_scale(j) += abs(A(i, j));
        }
    
    for (int i = 0; i < A.GetNbLastRow(); i++)
      for (int j = 0; j < A.GetN(); j++)
        diagonal_scale(j) += abs(A(m+i, j));
    
    for (int j = 0; j < A.GetNbLastCol(); j++)
      for (int i = 0; i < m; i++)
        diagonal_scale(n+j) += abs(A(i, n+j));    
  }

  
  //! For each row and column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row, Vector<T>& sum_col,
                    const Matrix<Complexe, General, BandedCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM();
    int n = A.GetN();
    sum_row.Reallocate(A.GetM());
    sum_row.Fill(0);
    sum_col.Reallocate(A.GetN());
    sum_col.Fill(0);
    
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          sum_row(i) += abs(A(i, j));
          sum_col(j) += abs(A(i, j));
        }
  }


  //! For each row and column of the matrix, computation of the sum of absolute values
  /*!
    \param[out] sum_row vector containing the sum of
                            the magnitudes of non-zero entries of each row
    \param[out] sum_col vector containing the sum of
                            the magnitudes of non-zero entries of each column
    \param[in] A given matrix
   */
  template<class T, class Complexe, class Allocator>
  void GetRowColSum(Vector<T>& sum_row, Vector<T>& sum_col,
                    const Matrix<Complexe, General, ArrowCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM() - A.GetNbLastRow();
    int n = A.GetN() - A.GetNbLastCol(); 
    sum_row.Reallocate(A.GetM());
    sum_row.Fill(0);
    sum_col.Reallocate(A.GetN());
    sum_col.Fill(0);
    
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          sum_row(i) += abs(A(i, j));
          sum_col(j) += abs(A(i, j));
        }
    
    for (int i = 0; i < A.GetNbLastRow(); i++)
      for (int j = 0; j < A.GetN(); j++)
        {
          sum_row(m+i) += abs(A(m+i, j));
          sum_col(j) += abs(A(m+i, j));
        }
    
    for (int j = 0; j < A.GetNbLastCol(); j++)
      for (int i = 0; i < m; i++)
        {
          sum_row(i) += abs(A(i, n+j));
          sum_col(n+j) += abs(A(i, n+j));
        }
  }
  
  
  //! Multiplication of rows by coefficients contained in coef_row
  /*!
    Equivalent Matlab function A = coef_row*A
   */
  template<class T, class Allocator, class T0>
  void ScaleLeftMatrix(Matrix<T, General, ArrowCol, Allocator>& A,
                       const Vector<T0>& coef_row)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM() - A.GetNbLastRow();
    int n = A.GetN() - A.GetNbLastCol(); 
    
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          A.Get(i, j) *= coef_row(i);
        }
    
    for (int i = 0; i < A.GetNbLastRow(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A.Get(m+i, j) *= coef_row(m+i);
    
    for (int j = 0; j < A.GetNbLastCol(); j++)
      for (int i = 0; i < m; i++)
        A.Get(i, n+j) *= coef_row(i);
  }
  
  
  //! copy real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSymComplexSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSymSparse, Allocator2>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int nb = A.GetRealRowSize(i);
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = A.IndexReal(i, j);
            B.Value(i, j) = A.ValueReal(i, j);
          }
      }
  }
  
  
  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSymSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSymSparse, Allocator2>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int nb = A.GetRowSize(i);
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = real(A.Value(i, j));
          }
      }    
  }
  

  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSymComplexSparse, Allocator1>& A,
		Matrix<T, Prop, RowSymSparse, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetRealDataSize();
    if ((m <= 0) || (n <= 0))
      {
	B.Clear();
	return;
      }
    
    Vector<T> Val(nnz);
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);

    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* val_real = A.GetRealData();

    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_real[i];
    
    for (long i = 0; i < nnz; i++)
      {
	Val(i) = val_real[i];
	Ind(i) = ind_real[i];
      }

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSymSparse, Allocator1>& A,
		Matrix<T, Prop, RowSymSparse, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    if ((m <= 0) || (n <= 0))
      {
	B.Clear();
	return;
      }
    
    Vector<T> Val(nnz);
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);

    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    complex<T>* val_ = A.GetData();
    
    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_[i];
    
    for (long i = 0; i < nnz; i++)
      {
	Val(i) = real(val_[i]);
	Ind(i) = ind_[i];
      }
    
    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowComplexSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSparse, Allocator2>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int nb = A.GetRealRowSize(i);
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = A.IndexReal(i, j);
            B.Value(i, j) = A.ValueReal(i, j);
          }
      }
  }


  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, ArrayRowSparse, Allocator1>& A,
		Matrix<T, Prop, ArrayRowSparse, Allocator2>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int nb = A.GetRowSize(i);
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = real(A.Value(i, j));
          }
      }    
  }


  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowComplexSparse, Allocator1>& A,
		Matrix<T, Prop, RowSparse, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetRealDataSize();
    if ((m <= 0) || (n <= 0))
      {
	B.Clear();
	return;
      }
    
    Vector<T> Val(nnz);
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);

    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* val_real = A.GetRealData();

    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_real[i];
    
    for (long i = 0; i < nnz; i++)
      {
	Val(i) = val_real[i];
	Ind(i) = ind_real[i];
      }

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! copies real part of A in B
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyReal(const Matrix<complex<T>, Prop, RowSparse, Allocator1>& A,
		Matrix<T, Prop, RowSparse, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    long nnz = A.GetDataSize();
    if ((m <= 0) || (n <= 0))
      {
	B.Clear();
	return;
      }
    
    Vector<T> Val(nnz);
    Vector<long> Ptr(m+1); Vector<int> Ind(nnz);

    long* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    complex<T>* val_ = A.GetData();
    
    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_[i];
    
    for (long i = 0; i < nnz; i++)
      {
	Val(i) = real(val_[i]);
	Ind(i) = ind_[i];
      }
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! B = real(A)
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void CopyReal(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
		Matrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    cout << "not implemented" << endl;
    abort();
  }
  
  
  //! Extracts a sub-matrix from a given sparse matrix
  /*!
    \param[in] A given sparse matrix
    \param[in] m first row/column to extract
    \param[in] n last row/column+1 to extract
    \param[out] B extracted sub-matrix of size (m-n) x (m-n)
    Equivalent Matlab function B = A(m:n-1, m:n-1)
   */
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m, int n, Matrix<T, General, ArrayRowSparse, Allocator>& B)
  {
    GetSubMatrix(A, m, n, m, n, B);
  }

  
  //! Extracts a sub-matrix from a given sparse matrix
  /*!
    \param[in] A given sparse matrix
    \param[in] m1 first row to extract
    \param[in] m2 last row+1 to extract
    \param[in] n1 first column to extract
    \param[in] n2 last column+1 to extract
    \param[out] B extracted sub-matrix of size (m2-m1) x (n2-n1)
    Equivalent Matlab function B = A(m1:m2-1, n1:n2-1)
   */  
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, General, ArrayRowSparse, Allocator>& B)
  {
    int M = m2 - m1;
    int N = n2 - n1;
    B.Reallocate(M, N);
    for (int i = m1; i < m2; i++)
      {
        int size_row = A.GetRowSize(i);
        int nb = 0;
        // counting the number of elements in the row
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2))
              nb++;
          }
        
        // filling the row
        B.ReallocateRow(i-m1, nb);
        nb = 0;
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2))
              {
                B.Index(i-m1, nb) = icol - n1;
                B.Value(i-m1, nb) = A.Value(i, j);
                nb++;
              }
          }
      }
  }


  //! Extracts a sub-matrix from a given sparse matrix
  /*!
    \param[in] A given sparse matrix
    \param[in] m1 first row to extract
    \param[in] m2 last row+1 to extract
    \param[in] n1 first column to extract
    \param[in] n2 last column+1 to extract
    \param[out] B extracted sub-matrix of size (m2-m1) x (n2-n1)
    Equivalent Matlab function B = A(m1:m2-1, n1:n2-1)
   */  
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B)
  {
    int M = m2 - m1;
    int N = n2 - n1;
    B.Reallocate(M, N);
    for (int i = m1; i < m2; i++)
      {
        int size_row = A.GetRowSize(i);
        int nb = 0;
        // counting the number of elements in the row
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2) && (icol-n1 >= i-m1))
              nb++;
          }
        
        // filling the row
        B.ReallocateRow(i-m1, nb);
        nb = 0;
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2) && (icol-n1 >= i-m1))
              {
                B.Index(i-m1, nb) = icol - n1;
                B.Value(i-m1, nb) = A.Value(i, j);
                nb++;
              }
          }
      }
  }
  

  //! Extracts a sub-matrix from a given sparse matrix
  /*!
    \param[in] A given sparse matrix
    \param[in] m first row/column to extract
    \param[in] n last row/column+1 to extract
    \param[out] B extracted sub-matrix of size (m-n) x (m-n)
    Equivalent Matlab function B = A(m:n-1, m:n-1)
   */
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
                    int m, int n, Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B)
  {
    GetSubMatrix(A, m, n, m, n, B);
  }

  
  //! Extracts a sub-matrix from a given sparse matrix
  /*!
    \param[in] A given sparse matrix
    \param[in] m1 first row to extract
    \param[in] m2 last row+1 to extract
    \param[in] n1 first column to extract
    \param[in] n2 last column+1 to extract
    \param[out] B extracted sub-matrix of size (m2-m1) x (n2-n1)
    Equivalent Matlab function B = A(m1:m2-1, n1:n2-1)
   */  
  template<class T, class Allocator>
  void GetSubMatrix(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& B)
  {
    int M = m2 - m1;
    int N = n2 - n1;
    if (M != N)
      {
        cout << "A symmetric matrix must be squared " << endl;
        abort();
      }
    
    B.Reallocate(M, N);
    for (int i = m1; i < m2; i++)
      {
        int size_row = A.GetRowSize(i);
        int nb = 0;
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2) && (icol-n1 >= i-m1))
              nb++;
          }
        
        B.ReallocateRow(i-m1, nb);
        nb = 0;
        for (int j = 0; j < size_row; j++)
          {
            int icol = A.Index(i, j);
            if ((icol >= n1)&&(icol < n2)&& (icol-n1 >= i-m1))
              {
                B.Index(i-m1, nb) = icol - n1;
                B.Value(i-m1, nb) = A.Value(i, j);
                nb++;
              }
          }
      }
  }
  

  //! Extracts a sub-matrix from a given sparse matrix  
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
                    int m, int n, Matrix<T, Prop, Storage, Allocator>& B)
  {
    cout << "not implemented" << endl;
    abort();
  }
  
  
  //! Extracts a sub-matrix from a given sparse matrix
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Prop, Storage, Allocator>& B)
  {
    cout << "not implemented" << endl;
    abort();
  }
  

#ifdef SELDON_WITH_MPI
  //! extracts a sub-matrix B from matrix A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m, int n,
		    DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    MPI_Comm& comm = B.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return GetSubMatrix(static_cast<const Matrix<T, Prop,
			  Storage, Allocator>& >(A),
                          m, n, static_cast<Matrix<T, Prop,
			  Storage, Allocator>& >(B));

    cout << "GetSubMatrix not implemented for distributed matrices" << endl;
    abort();    
  }
  
  
  //! extracts a sub-matrix B from matrix A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    MPI_Comm& comm = B.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return GetSubMatrix(static_cast<const Matrix<T, Prop,
			  Storage, Allocator>& >(A),
                          m1, m2, n1, n2,
                          static_cast<Matrix<T, Prop,
			  Storage, Allocator>& >(B));

    cout << "GetSubMatrix not implemented for distributed matrices" << endl;
    abort();    
  }
#endif
  
  
  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    long* ptrA = A.GetPtr();
    int* indA = A.GetInd();
    T* dataA = A.GetData();
    int n = 0; long nnz = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  nnz += ptrA[i+1] - ptrA[i];
	  n++;
	}
    
    Vector<T> values(nnz); Vector<long> ptr(n+1); Vector<int> ind(nnz);
    n = 0; ptr(0) = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  ptr(n+1) = ptr(n) + ptrA[i+1] - ptrA[i];
	  for (long j = ptr(n); j < ptr(n+1); j++)
	    {
	      long j2 = j - ptr(n) + ptrA[i];
	      ind(j) = IndexRow(indA[j2]);
	      values(j) = dataA[j2];
	    }
	  
	  // next row
	  n++;
	}
    
    A.SetData(n, n, values, ptr, ind);
  }
  

  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSymSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    long* ptrA = A.GetPtr();
    int* indA = A.GetInd();
    T* dataA = A.GetData();
    int n = 0; long nnz = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  nnz += ptrA[i+1] - ptrA[i];
	  n++;
	}
    
    Vector<T> values(nnz); Vector<long> ptr(n+1); Vector<int> ind(nnz);
    n = 0;
    ptr(0) = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  ptr(n+1) = ptr(n) + ptrA[i+1] - ptrA[i];
	  for (long j = ptr(n); j < ptr(n+1); j++)
	    {
	      long j2 = j - ptr(n) + ptrA[i];
	      ind(j) = IndexRow(indA[j2]);
	      values(j) = dataA[j2];
	    }
	  
	  // next row
	  n++;
	}
    
    A.SetData(n, n, values, ptr, ind);
  }

  
  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    int n = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  for (int j = 0; j < A.GetRowSize(i); j++)
	    A.Index(i, j) = IndexRow(A.Index(i, j));
	  
	  if (i != n)
	    A.SwapRow(i, n);
	  
	  // next row
	  n++;
	}

    A.Resize(n, n);
  }


  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    int n = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  for (int j = 0; j < A.GetRowSize(i); j++)
	    A.Index(i, j) = IndexRow(A.Index(i, j));
	  
	  if (i != n)
	    A.SwapRow(i, n);

	  // next row
	  n++;
	}

    A.Resize(n, n);
  }


  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowComplexSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    int n = 0; long nnz_real = 0, nnz_imag = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  nnz_real += ptr_real[i+1] - ptr_real[i];
 	  nnz_imag += ptr_imag[i+1] - ptr_imag[i];
	  n++;
	}

    Vector<Treal> ValReal(nnz_real), ValImag(nnz_imag);
    Vector<long> PtrReal(n+1), PtrImag(n+1);
    Vector<int> IndReal(nnz_real), IndImag(nnz_imag);
    n = 0; PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  PtrReal(n+1) = PtrReal(n) + ptr_real[i+1] - ptr_real[i];
	  PtrImag(n+1) = PtrImag(n) + ptr_imag[i+1] - ptr_imag[i];
	  for (long j = PtrReal(n); j < PtrReal(n+1); j++)
	    {
	      long j2 = j - PtrReal(n) + ptr_real[i];
	      IndReal(j) = IndexRow(ind_real[j2]);
	      ValReal(j) = data_real[j2];
	    }

	  for (long j = PtrImag(n); j < PtrImag(n+1); j++)
	    {
	      long j2 = j - PtrImag(n) + ptr_imag[i];
	      IndImag(j) = IndexRow(ind_imag[j2]);
	      ValImag(j) = data_imag[j2];
	    }

	  // next row
	  n++;
	}
    
    A.SetData(n, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }
  

  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, RowSymComplexSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    long* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    long* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    int n = 0; long nnz_real = 0, nnz_imag = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  nnz_real += ptr_real[i+1] - ptr_real[i];
 	  nnz_imag += ptr_imag[i+1] - ptr_imag[i];
	  n++;
	}

    Vector<Treal> ValReal(nnz_real), ValImag(nnz_imag);
    Vector<long> PtrReal(n+1), PtrImag(n+1);
    Vector<int> IndReal(nnz_real), IndImag(nnz_imag);
    n = 0; PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  PtrReal(n+1) = PtrReal(n) + ptr_real[i+1] - ptr_real[i];
	  PtrImag(n+1) = PtrImag(n) + ptr_imag[i+1] - ptr_imag[i];
	  for (long j = PtrReal(n); j < PtrReal(n+1); j++)
	    {
	      long j2 = j - PtrReal(n) + ptr_real[i];
	      IndReal(j) = IndexRow(ind_real[j2]);
	      ValReal(j) = data_real[j2];
	    }

	  for (long j = PtrImag(n); j < PtrImag(n+1); j++)
	    {
	      long j2 = j - PtrImag(n) + ptr_imag[i];
	      IndImag(j) = IndexRow(ind_imag[j2]);
	      ValImag(j) = data_imag[j2];
	    }

	  // next row
	  n++;
	}
    
    A.SetData(n, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    int n = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  for (int j = 0; j < A.GetRealRowSize(i); j++)
	    A.IndexReal(i, j) = IndexRow(A.IndexReal(i, j));

	  for (int j = 0; j < A.GetImagRowSize(i); j++)
	    A.IndexImag(i, j) = IndexRow(A.IndexImag(i, j));

	  if (i != n)
	    {
	      A.SwapRealRow(i, n);
	      A.SwapImagRow(i, n);
	    }
	  
	  // next row
	  n++;
	}
    
    A.Resize(n, n);
  }


  //! removes unnecessary rows from matrix A
  //! (it is assumed that these rows do not communicate with other rows)
  template<class T, class Prop, class Allocator>
  void CompressMatrix(DistributedMatrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
		      const IVect& IndexRow)
  {
    int n = 0;
    for (int i = 0; i < IndexRow.GetM(); i++)
      if (IndexRow(i) >= 0)
	{
	  if (IndexRow(i) != n)
	    {
	      cout << "Case not handled" << endl;
	      abort();
	    }

	  for (int j = 0; j < A.GetRealRowSize(i); j++)
	    A.IndexReal(i, j) = IndexRow(A.IndexReal(i, j));

	  for (int j = 0; j < A.GetImagRowSize(i); j++)
	    A.IndexImag(i, j) = IndexRow(A.IndexImag(i, j));

	  if (i != n)
	    {
	      A.SwapRealRow(i, n);
	      A.SwapImagRow(i, n);
	    }
	  
	  // next row
	  n++;
	}
    
    A.Resize(n, n);
  }


  //! extracts a sub-matrix Asub from A by giving row and column numbers
  template<class T>
  void ExtractSubMatrix(const Matrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			const IVect& col_num, const IVect& index_col,
			Matrix<T, General, ArrayRowSparse>& Asub)
  {
    Asub.Clear();
    if ((row_num.GetM() <= 0) || (col_num.GetM() <= 0))
      return;
    
    Asub.Reallocate(row_num.GetM(), col_num.GetM());
    for (int i2 = 0; i2 < row_num.GetM(); i2++)
      {
	int i = row_num(i2);
	int size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index_col(A.Index(i, j)) >= 0)
	    size_row++;

	Asub.ReallocateRow(i2, size_row);
	size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index_col(A.Index(i, j)) >= 0)
	    {
	      Asub.Index(i2, size_row) = index_col(A.Index(i, j));
	      Asub.Value(i2, size_row) = A.Value(i, j);
	      size_row++;
	    }
      }
  }


  //! extracts a sub-matrix Asub from A by giving row and column numbers
  template<class T>
  void ExtractSubMatrix(const Matrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			const IVect& col_num, const IVect& index_col,
			Matrix<T, Symmetric, ArrayRowSymSparse>& Asub)
  {
    Asub.Clear();
    if ((row_num.GetM() <= 0) || (col_num.GetM() <= 0))
      return;
    
    Asub.Reallocate(row_num.GetM(), col_num.GetM());
    for (int i2 = 0; i2 < row_num.GetM(); i2++)
      {
	int i = row_num(i2);
	int size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if ((index_col(A.Index(i, j)) >= 0) && (index_col(A.Index(i, j)) >= i2))
	    size_row++;
	      
	Asub.ReallocateRow(i2, size_row);
	size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if ((index_col(A.Index(i, j)) >= 0) && (index_col(A.Index(i, j)) >= i2))
	    {
	      Asub.Index(i2, size_row) = index_col(A.Index(i, j));
	      Asub.Value(i2, size_row) = A.Value(i, j);
	      size_row++;
	    }
      }
  }


  template<class T>
  void ExtractSubMatrix(const DistributedMatrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			DistributedMatrixIntegerArray& info,
			DistributedMatrix<T, General, ArrayRowSparse>& Asub)
  {
#ifdef SELDON_WITH_MPI
    if ((A.GetMaxDataSizeDistantCol() > 0) || (A.GetMaxDataSizeDistantRow() > 0))
      {
	cout << "Function not implemented for distant interactions" << endl;
	abort();
      }
    
    ExtractSubMatrix(A, row_num, index_row, row_num, index_row, Asub);

    Asub.Init(info.nglob, &info.GlobalRowNumbers, &info.OverlapRowNumbers,
	      &info.OverlapProcNumbers, info.nodl_scalar, info.nb_unknowns_scal,
	      &info.ProcSharingRows, &info.SharingRowNumbers, info.comm);
    
#else
    ExtractSubMatrix(A, row_num, index_row, row_num, index_row, Asub);
#endif
  }
  

  template<class T>
  void ExtractSubMatrix(const DistributedMatrix<T, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			DistributedMatrixIntegerArray& info,
			DistributedMatrix<T, Symmetric, ArrayRowSymSparse>& Asub)
  {
#ifdef SELDON_WITH_MPI
    if ((A.GetMaxDataSizeDistantCol() > 0) || (A.GetMaxDataSizeDistantRow() > 0))
      {
	cout << "Function not implemented for distant interactions" << endl;
	abort();
      }
    
    ExtractSubMatrix(A, row_num, index_row, row_num, index_row, Asub);

    Asub.Init(info.nglob, &info.GlobalRowNumbers, &info.OverlapRowNumbers,
	      &info.OverlapProcNumbers, info.nodl_scalar, info.nb_unknowns_scal,
	      &info.ProcSharingRows, &info.SharingRowNumbers, info.comm);
    
#else
    ExtractSubMatrix(A, row_num, index_row, row_num, index_row, Asub);
#endif

  }

}

#define SELDON_FILE_FUNCTION_MATRIX_EXTRACTION_CXX
#endif

