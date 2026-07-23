// Copyright (C) 2001-2010 Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.

#ifndef SELDON_FILE_PASTIX_CXX

#include "Pastix.hxx"

namespace Seldon
{

  //! Default constructor.
  template<class T>
  MatrixPastix<T>::MatrixPastix()
  {
    comm_ = MPI_COMM_WORLD;
    pastix_data = NULL;
    spm = NULL;
    n = 0;
    for (int i = 0; i < IPARM_SIZE; i++)
      iparm[i] = 0;
    
    for (int i = 0; i < DPARM_SIZE; i++)
      dparm[i] = 0;
    
    // Factorization of a matrix on a single processor.
    distributed = false;

    // No refinement by default.
    refine_solution = false;

    // initializing parameters
    pastixInitParam(iparm, dparm);
    spm = (spmatrix_t*)malloc(sizeof(spmatrix_t));
    
    iparm[IPARM_VERBOSE] = PastixVerboseNot;
    iparm[IPARM_FREE_CSCUSER] = 1;
    
    threshold_pivot = 0.0;
    adjust_threshold_pivot = false;
    cholesky = false;
  }


  //! destructor
  template<class T>
  MatrixPastix<T>::~MatrixPastix()
  {
    Clear();
    free(spm);
  }


  template<class T>
  bool MatrixPastix<T>::UseInteger8() const  
  {
    if (sizeof(pastix_int_t) == 8)
      return true;
    
    return false;
  }


  //! Clearing factorization.
  template<class T>
  void MatrixPastix<T>::Clear()
  {
    if (n > 0)
      {
        pastixFinalize(&pastix_data);
        spmExit(spm);
	n = 0;
        pastix_data = NULL;
        distributed = false;
      }
  }


  //! no message will be displayed
  template<class T>
  void MatrixPastix<T>::HideMessages()
  {
    iparm[IPARM_VERBOSE] = PastixVerboseNot;
  }


  //! Low level of display.
  template<class T>
  void MatrixPastix<T>::ShowMessages()
  {
    iparm[IPARM_VERBOSE] = PastixVerboseNo;
  }


  //! Displaying all messages.
  template<class T>
  void MatrixPastix<T>::ShowFullHistory()
  {
    iparm[IPARM_VERBOSE] = PastixVerboseYes;
  }


  //! selects the algorithm used for reordering
  template<class T>
  void MatrixPastix<T>::SelectOrdering(int type)
  {
    iparm[IPARM_ORDERING] = type;
  }

  
  //! provides a permutation array (instead of using Scotch reordering)
  template<class T>
  void MatrixPastix<T>::SetPermutation(const IVect& permut)
  {
    iparm[IPARM_ORDERING] = PastixOrderPersonal;
    pastixOrderAlloc(&ord, permut.GetM(), 0);
    for (int i = 0; i < permut.GetM(); i++)
      {
	ord.permtab[i] = permut(i);
	ord.peritab[permut(i)] = i;
      }
  }

  
  //! sets Cholesky factorisation
  template<class T>
  void MatrixPastix<T>::SetCholeskyFacto(bool chol)
  {
    cholesky = chol;
  }
  
  
  //! you can change the threshold used for static pivoting
  template<class T>
  void MatrixPastix<T>::SetPivotThreshold(double eps)
  {
    adjust_threshold_pivot = true;
    threshold_pivot = eps;
  }
  
  
  //! You can require that solution is refined after LU resolution.
  template<class T>
  void MatrixPastix<T>::RefineSolution()
  {
    refine_solution = true;
  }


  //! You can require that solution is not refined (faster).
  template<class T>
  void MatrixPastix<T>::DoNotRefineSolution()
  {
    refine_solution = false;
  }

  
  //! Returns the size of memory used by the factorisation in bytes
  template<class T>
  size_t MatrixPastix<T>::GetMemorySize() const
  {
    size_t taille = 0;
    if (n <= 0)
      return taille;
    
    // assuming that for each term, a value and an index is needed
    taille += (sizeof(T)+sizeof(pastix_int_t))*iparm[IPARM_NNZEROS];
    return taille;
  }
  

  template<class T>
  int MatrixPastix<T>::GetInfoFactorization() const
  {
    return 0;
  }

  
  //! Returning ordering found by Scotch.
  template<class T>
  template<class T0, class Prop, class Storage, class Allocator, class Tint>
  void MatrixPastix<T>::
  FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
               Vector<Tint>& numbers, bool keep_matrix)
  {
    // We clear the previous factorization, if any.
    Clear();

    distributed = false;

    n = mat.GetN();
    if (n <= 0)
      return;

    Vector<pastix_int_t> Ptr, Ind;

    GetSymmetricPattern(mat, Ptr, Ind);
    if (!keep_matrix)
      mat.Clear();
    
    cout << "Not implemented" << endl;
    abort();
  }


  //! Factorization of unsymmetric matrix
  template<class T> template<class T0, class Storage, class Allocator>
  void MatrixPastix<T>
  ::FactorizeMatrix(Matrix<T0, General, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    Vector<pastix_int_t> Ptr, IndRow;
    Vector<T> Val;

    General prop;
    ConvertToCSC(mat, prop, Ptr, IndRow, Val, true);
    if (!keep_matrix)
      mat.Clear();

    FactorizeCSC(Ptr, IndRow, Val, false);
  }


  template<class T>
  void MatrixPastix<T>::FillSpmMatrix(Vector<pastix_int_t>& Ptr, Vector<pastix_int_t>& IndRow,
                                      Vector<T>& Val, bool sym, const Vector<int>& glob_number,
                                      MPI_Comm comm_facto)
  {
    spmInit(spm);
    
    pastix_int_t nnz = IndRow.GetM();
    spm->fmttype = SpmCSC;
    spm->dof = 1;
    if (glob_number.GetM() < n)
      spm->loc2glob = NULL;
    else
      {
        Vector<int> glob_num(glob_number);
        spm->loc2glob = glob_num.GetData();
        spm->comm = comm_facto;
        int nb_proc, rank_proc;
        MPI_Comm_size(comm_facto, &nb_proc);
        MPI_Comm_rank(comm_facto, &rank_proc);
        spm->clustnum = rank_proc;
        spm->clustnbr = nb_proc;
        glob_num.Nullify();
      }
    
    spm->baseval = 0;
    spm->n = n;
    spm->nnz = nnz;
    spm->layout = SpmColMajor;

    if (sizeof(T) == 8)
      spm->flttype = SpmDouble;
    else
      spm->flttype = SpmComplex64;

    if (sym)
      spm->mtxtype = SpmSymmetric;
    else
      spm->mtxtype = SpmGeneral;

    spm->colptr = Ptr.GetData();
    spm->rowptr = IndRow.GetData();
    spm->values = (double*)Val.GetData();

    spmUpdateComputedFields( spm );

    spmatrix_t spm2;
    int rc = spmCheckAndCorrect( spm, &spm2 );
    if (rc != 0)
      {
        cout << "Error in the matrix" << endl;
        abort();
      }        
    
    Ptr.Nullify(); IndRow.Nullify(); Val.Nullify();
  }
  

  template<class T>
  void MatrixPastix<T>
  ::FactorizeCSC(Vector<pastix_int_t>& Ptr, Vector<pastix_int_t>& IndRow,
		 Vector<T>& Val, bool sym)
  {
    distributed = false;
    n = Ptr.GetM()-1;
    if (n <= 0)
      return;

    Vector<int> glob_num; // void because of sequential matrix
    FillSpmMatrix(Ptr, IndRow, Val, sym, glob_num, MPI_COMM_SELF);
    
    if (sym)
      {
	if (cholesky)
	  iparm[IPARM_FACTORIZATION] = PastixFactLLT;
	else
	  iparm[IPARM_FACTORIZATION] = PastixFactLDLT;
      }
    else
      iparm[IPARM_FACTORIZATION] = PastixFactLU;

    // pivot threshold
    if (adjust_threshold_pivot)
      dparm[DPARM_EPSILON_MAGN_CTRL] = threshold_pivot;
    
    pastixInit(&pastix_data, MPI_COMM_SELF, iparm, dparm);

    // analyze
    pastix_task_analyze(pastix_data, spm);

    // numerical factorization
    pastix_task_numfact(pastix_data, spm);
    
    if (iparm[IPARM_VERBOSE] != PastixVerboseNot)
      cout << "Factorization successful" << endl;
  }


  //! Factorization of symmetric matrix.
  template<class T> template<class T0, class Storage, class Allocator>
  void MatrixPastix<T>::
  FactorizeMatrix(Matrix<T0, Symmetric, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    Vector<pastix_int_t> Ptr, IndRow;
    Vector<T> Val;

    Symmetric prop;
    ConvertToCSR(mat, prop, Ptr, IndRow, Val);

    FactorizeCSC(Ptr, IndRow, Val, true);
  }


  //! solving A x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA,
                              Vector<T, VectFull, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), 1);
  }


  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA, T* x_ptr, int nrhs_)
  {
    pastix_int_t nrhs = nrhs_;

    if (cholesky)
      {
        spm_coeftype_t flttype;
        if (sizeof(T) == 8)
          flttype = SpmDouble;
        else
          flttype = SpmComplex64;

        if (TransA.NoTrans())
          {
            pastix_subtask_applyorder(pastix_data, flttype, PastixDirForward,
                                      n, nrhs, (void*)x_ptr, n);

            pastix_subtask_trsm(pastix_data, flttype, PastixLeft, PastixLower,
                                PastixNoTrans, PastixNonUnit, nrhs, (void*)x_ptr, n);
          }
        else
          {
            pastix_subtask_trsm(pastix_data, flttype, PastixLeft, PastixLower,
                                PastixTrans, PastixNonUnit, nrhs, (void*)x_ptr, n);

            pastix_subtask_applyorder(pastix_data, flttype, PastixDirBackward,
                                      n, nrhs, (void*)x_ptr, n);
          }

        return;
      }
    
    if (TransA.Trans())
      iparm[IPARM_TRANSPOSE_SOLVE] = PastixTrans;
    else
      iparm[IPARM_TRANSPOSE_SOLVE] = PastixNoTrans;        
    
    T* b_ptr = NULL;
    if (refine_solution)
      {
        b_ptr = (T*) malloc(nrhs*n*sizeof(T));
        memcpy(b_ptr, x_ptr, nrhs*n*sizeof(T));
      }
    
    pastix_task_solve(pastix_data, nrhs, (void*)x_ptr, n);
    if (refine_solution)
      pastix_task_refine(pastix_data, n, nrhs, (void*)b_ptr, n, (void*)x_ptr, n);
  }
  

  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(const SeldonTranspose& TransA,
                              Matrix<T, General, ColMajor, Allocator2>& x)
  {
    Solve(TransA, x.GetData(), x.GetN());
  }


  //! Computes x = L b or x = L^T b (A = L L^T)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Mlt(const SeldonTranspose& TransA,
			    Vector<T, VectFull, Allocator2>& x)
  {
    Mlt(TransA, x.GetData(), 1);
  }


  //! Computes x = L b or x = L^T b (A = L L^T)
  template<class T>
  void MatrixPastix<T>::Mlt(const SeldonTranspose& TransA, T* x_ptr, int nrhs_)
  {
    //    pastix_int_t nrhs = nrhs_;

    if (cholesky)
      {
        cout << "Multiplication with L Not implemented in MatrixPastix" << endl;
        abort();
        /*        if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_UTRMV;
        else
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_LTRMV;
        
          iparm[IPARM_END_TASK] = API_TASK_SOLVE; */
      }
    else
      {
        cout << "Mlt is defined for Cholesky only" << endl;
        abort();
      }
  }


  //! Modifies the number of threads per node.
  template<class T>
  void MatrixPastix<T>::SetNumberOfThreadPerNode(int num_thread)
  {
    iparm[IPARM_THREAD_NBR] = num_thread;
  }


  template<class T>
  void MatrixPastix<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<long>& Ptr,
                             Vector<int>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    Vector<int> PtrInt(Ptr.GetM());
    for (int i = 0; i < Ptr.GetM(); i++)
      PtrInt(i) = Ptr(i);

    FactorizeParallel(comm_facto, PtrInt, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }

  
  template<class T>
  void MatrixPastix<T>::
  FactorizeDistributedMatrix(MPI_Comm& comm_facto, Vector<int64_t>& Ptr,
                             Vector<int64_t>& IndRow, Vector<T>& Val,
                             const Vector<int>& glob_number,
                             bool sym, bool keep_matrix)
  {
    FactorizeParallel(comm_facto, Ptr, IndRow, Val, 
                      glob_number, sym, keep_matrix);
  }
  
  
  template<class T> template<class Tint>
  void MatrixPastix<T>::  
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<Tint>&, Vector<Tint>&,
                    Vector<T>&,
                    const Vector<int>& glob_num,
                    bool sym, bool keep_matrix)
  {
    cout << "Not available for this type of integers " << endl;
    cout << "Size of Tint = " << sizeof(Tint) << endl;
    abort();
  }
  

  //! Distributed factorization (on several nodes).
  template<class T>
  void MatrixPastix<T>::
  FactorizeParallel(MPI_Comm& comm_facto,
                    Vector<pastix_int_t>& Ptr,
                    Vector<pastix_int_t>& IndRow,
                    Vector<T>& Val, const Vector<int>& glob_number,
                    bool sym, bool keep_matrix)
  {
    Clear();

    iparm[IPARM_FREE_CSCUSER] = 0;
    distributed = true; n = Ptr.GetM()-1;
    FillSpmMatrix(Ptr, IndRow, Val, sym, glob_number, comm_facto);

    if (sym)
      {
	if (cholesky)
	  iparm[IPARM_FACTORIZATION] = PastixFactLLT;
	else
	  iparm[IPARM_FACTORIZATION] = PastixFactLDLT;
      }
    else
      iparm[IPARM_FACTORIZATION] = PastixFactLU;

    // pivot threshold
    if (adjust_threshold_pivot)
      dparm[DPARM_EPSILON_MAGN_CTRL] = threshold_pivot;
    
    pastixInit(&pastix_data, comm_facto, iparm, dparm);

    // analyze
    pastix_task_analyze(pastix_data, spm);

    // numerical factorization
    pastix_task_numfact(pastix_data, spm);
  }

  
  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
                                         const SeldonTranspose& TransA,
                                         Vector<T, VectFull, Allocator2>& x,
                                         const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
					 const SeldonTranspose& TransA,
					 T* x_ptr, int nrhs_,
					 const IVect& glob_num)
  {
    pastix_int_t nrhs = nrhs_;

    if (cholesky)
      {
        cout << "Not implemented" << endl;
        abort();
        /*if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_BACKWARD_ONLY;
        else
        iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_FORWARD_ONLY; */
      }
    else
      {
        if (TransA.Trans())
          iparm[IPARM_TRANSPOSE_SOLVE] = PastixTrans;
        else
          iparm[IPARM_TRANSPOSE_SOLVE] = PastixNoTrans;        
      }

    T* b_ptr;
    if (refine_solution)
      {
        b_ptr = (T*) malloc(nrhs*n*sizeof(T));
        memcpy(b_ptr, x_ptr, nrhs*n*sizeof(T));
      }
    
    pastix_task_solve(pastix_data, nrhs, (void*)x_ptr, n);
    if (refine_solution)
      pastix_task_refine(pastix_data, n, nrhs, (void*)b_ptr, n, (void*)x_ptr, n);
  }
  

  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::SolveDistributed(MPI_Comm& comm_facto,
                                         const SeldonTranspose& TransA,
                                         Matrix<T, General, ColMajor, Allocator2>& x,
                                         const Vector<int>& glob_num)
  {
    SolveDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }

  
  //! computes L x or L^T x
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       Vector<T, VectFull, Allocator2>& x,
				       const Vector<int>& glob_num)
  {
    MltDistributed(comm_facto, TransA, x.GetData(), 1, glob_num);
  }


  //! computes L x or L^T x
  template<class T>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       T* x_ptr, int nrhs_,
				       const IVect& glob_num)
  {
    /*
    pastix_int_t nrhs = nrhs_;

    for (int k = 0; k < nrhs_; k++)
      {
	// workaround : each linear system is solved separately
	nrhs = 1;
        
        if (cholesky)
          {
            if (TransA.Trans())
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_UTRMV;
            else
              iparm[IPARM_TRANSPOSE_SOLVE] = API_SOLVE_LTRMV;
            
            iparm[IPARM_END_TASK] = API_TASK_SOLVE;
          }
        else
          {
            cout << "Mlt is defined for Cholesky only" << endl;
            abort();
          }
        
        iparm[IPARM_START_TASK] = API_TASK_SOLVE;
        //CallPastix(comm_facto, NULL, NULL, NULL, x_ptr, nrhs);
        
        CallPastix(comm_facto, NULL, NULL, NULL, x_ptr, nrhs);
        x_ptr += n;
      }
    */
  }


  //! solves A x = b or A^T x = b in parallel
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::MltDistributed(MPI_Comm& comm_facto,
				       const SeldonTranspose& TransA,
				       Matrix<T, General, ColMajor, Allocator2>& x,
				       const Vector<int>& glob_num)
  {
    MltDistributed(comm_facto, TransA, x.GetData(), x.GetN(), glob_num);
  }


  //! Factorization of a matrix of same type T as the Pastix object 
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<T>& mat_lu, bool keep_matrix, T& x)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  

  //! Factorization of a complex matrix with a real Pastix object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<T>& mat_lu, bool keep_matrix, complex<T>& x)
  {
    throw WrongArgument("GetLU(Matrix<complex<T> >& A, MatrixPastix<T>& mat_lu, bool)",
			"The LU matrix must be complex");
  }

  
  //! Factorization of a real matrix with a complex Pastix object
  template<class MatrixSparse, class T>
  void GetLU(MatrixSparse& A, MatrixPastix<complex<T> >& mat_lu, bool keep_matrix, T& x)
  {
    throw WrongArgument("GetLU(Matrix<T>& A, MatrixMumps<Pastix<T> >& mat_lu, bool)",
			"The sparse matrix must be complex");
  }
  
  
  //! Factorization of a general matrix with Pastix
  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixPastix<T>& mat_lu,
	     bool keep_matrix)
  {
    // we check if the type of non-zero entries of matrix A
    // and of the Pastix object (T) are different
    // we call one of the GetLUs written above
    // such a protection avoids to compile the factorisation of a complex
    // matrix with a real Pastix object
    typename Matrix<T0, Prop, Storage, Allocator>::entry_type x;
    GetLU(A, mat_lu, keep_matrix, x);
  }
  

  //! LU resolution with a vector whose type is the same as Pastix object
  template<class T, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  //! LU resolution with a vector whose type is the same as Pastix object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! LU resolution with a matrix whose type is the same as Pastix object
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  //! LU resolution with a matrix whose type is the same as UmfPack object
  //! Solves transpose system A^T x = b or A x = b depending on TransA
  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  //! Solves A x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(MatrixPastix<double>& mat_lu,
               Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));
  }
  

  //! Solves A x = b or A^T x = b, where A is real and x is complex
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<double>& mat_lu, Vector<complex<double>, VectFull, Allocator>& x)
  {
    Matrix<double, General, ColMajor> y(x.GetM(), 2);
    
    for (int i = 0; i < x.GetM(); i++)
      {
	y(i, 0) = real(x(i));
	y(i, 1) = imag(x(i));
      }
    
    SolveLU(TransA, mat_lu, y);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<double>(y(i, 0), y(i, 1));

  }


  //! Solves A x = b, where A is complex and x is real => Forbidden
  template<class Allocator>
  void SolveLU(MatrixPastix<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPastix<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }

  
  //! Solves A x = b or A^T x = b, where A is complex and x is real => Forbidden  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixPastix<complex<double> >& mat_lu,
               Vector<double, VectFull, Allocator>& x)
  {
    throw WrongArgument("SolveLU(MatrixPastix<complex<double> >, Vector<double>)", 
			"The result should be a complex vector");
  }


  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixPastix<T>& mat_chol, bool keep_matrix)
  {
    mat_chol.SetCholeskyFacto(true);
    //IVect permut(A.GetM()); permut.Fill();
    //mat_chol.SetPermutation(permut);
    mat_chol.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Allocator>
  void
  SolveCholesky(const SeldonTranspose& TransA,
                MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Solve(TransA, x);
  }


  template<class T, class Allocator>
  void
  MltCholesky(const SeldonTranspose& TransA,
              MatrixPastix<T>& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Mlt(TransA, x);
  }

} // end namespace

#define SELDON_FILE_PASTIX_CXX
#endif
