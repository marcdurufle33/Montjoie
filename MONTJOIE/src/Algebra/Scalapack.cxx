#ifndef SELDON_FILE_SCALAPACK_CXX

namespace Seldon
{
  //! default constructor
  BlacsHandle::BlacsHandle()
  {
    ictxt = 0;
    nprow = 0;
    npcol = 0;
    myrow = 0;
    mycol = 0;
  }
  
  
  //! returns Blacs context
  int BlacsHandle::GetContext() const
  {
    return ictxt;
  }
  
  
  //! returns the number of processors over which rows are distributed
  int BlacsHandle::GetNumOfRowProcs() const
  {
    return nprow;
  }
  
  
  //! returns the number of processors over which columns are distributed
  int BlacsHandle::GetNumOfColProcs() const
  {
    return npcol;
  }
  
  
  //! returns the processor number for rows
  int BlacsHandle::GetRowProcNum() const
  {
    return myrow;
  }
  
  
  //! returns the processor number for columns
  int BlacsHandle::GetColProcNum() const
  {
    return mycol;
  }
  
  
  //! initialisation of Blacs grid
  void BlacsHandle::Init()
  {
    int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    DistributeBlacsComputation(nb_proc,
                               ictxt, nprow, npcol, myrow, mycol);
  }

  
  //! default constructor
  DistributedMatrixDense_Base::DistributedMatrixDense_Base()
  {
    mb_ = 64; nb_ = 64;
    mglob_ = 0; nglob_ = 0;
    ictxt_ = 0;
  }    

  
  //! initialisation of a distributed matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_Pointers<T, Prop, Storage, Allocator>::
  Init(BlacsHandle& bl, int m, int n, int mb, int nb)
  {
    this->mglob_ = m;
    this->nglob_ = n;
    this->mb_ = mb;
    this->nb_ = nb;
    
    this->ictxt_ = bl.GetContext();
    this->nprow_ = bl.GetNumOfRowProcs();
    this->npcol_ = bl.GetNumOfColProcs();
    this->myrow_ = bl.GetRowProcNum();
    this->mycol_ = bl.GetColProcNum();
    
    int mloc, nloc;
    GetScalapackDescriptor(m, n, mb, nb, this->ictxt_,
                           this->nprow_, this->npcol_, this->myrow_, this->mycol_,
                           this->desc_, mloc, nloc);
    
    this->Reallocate(mloc, nloc);
  }

  
  //! performs A(i, j) = val with i, j global numbers
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_Pointers<T, Prop, Storage, Allocator>::
  SetGlobal(int i, int j, const T& val)
  {
    int iloc, jloc, iarow, jacol;
    GetLocalScalapackIndex(i, j, this->desc_, this->nprow_,
                           this->npcol_, this->myrow_, this->mycol_,
                           iloc, jloc, iarow, jacol);
    
    if ((iarow == this->myrow_) && (jacol == this->mycol_))
      this->Set(iloc, jloc, val);
  }


  //! returns A(i, j) with i, j global numbers
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_Pointers<T, Prop, Storage, Allocator>::
  GetGlobal(int i, int j, T& val)
  {
    int iloc, jloc, iarow, jacol;
    GetLocalScalapackIndex(i, j, this->desc_, this->nprow_,
                           this->npcol_, this->myrow_, this->mycol_,
                           iloc, jloc, iarow, jacol);
    
    if ((iarow == this->myrow_) && (jacol == this->mycol_))
      val = this->Get(iloc, jloc);
  }
  
  
  //! initialisation of the Blas grid
  void DistributeBlacsComputation(int nb_procs, int& ictxt,
                                  int& nprow, int& npcol,
                                  int& myrow, int& mycol)
  {
    // context is retrieved
    Cblacs_get(-1, 0, &ictxt);
    
    // we search nprow and npcol such that nb_procs = nprow * npcol
    Montjoie::FindTwoFactors(nb_procs, nprow, npcol);
    
    // then initialisation of the 2-D grid
    char format = 'R';
    Cblacs_gridinit(&ictxt, &format, nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  }

  
  //! initialisaion of the Scalapack descriptor for a dense matrix
  void GetScalapackDescriptor(int m, int n, int mb, int nb, int ictxt,
                              int nprow, int npcol, int myrow, int mycol,
                              IVect& desc, int& np, int& nq)
  {
    if ((myrow < nprow) && (mycol < npcol))
      {
        // number of rows/columns 
        int rank = 0;
        np = numroc_(&m, &mb, &myrow, &rank, &nprow);
        nq = numroc_(&n, &nb, &mycol, &rank, &npcol);
        
        // calling descinit
        int info;
        desc.Reallocate(9);
        descinit_(desc.GetData(), &m, &n, &mb, &nb, &rank, &rank, &ictxt, &np, &info);
        
        if (info != 0)
          {
            cout << "Failed to distribute matrix "<< endl;
            abort();
          }
      }
  }
  
  
  //! returns local numbers iloc, jloc from global numbers i, j
  void GetLocalScalapackIndex(int i, int j, const IVect& desc, int nprow, int npcol, 
                              int myrow, int mycol, int& iloc, int& jloc, int& iarow, int& jacol)
  {
    i++; j++;
    infog2l_(&i, &j, desc.GetData(), &nprow, &npcol, &myrow, &mycol, &iloc, &jloc, &iarow, &jacol);
    iloc--; jloc--; 
  }


  //! vector is used to fill the first column of a distributed matrix
  template<class T, class Prop, class Allocator, class Allocator2>
  void DistributeVector(Vector<T, VectFull, Allocator>& X,
                        DistributedMatrix<T, Prop, ColMajor, Allocator2>& A)
  {
    for (int i = 0; i < X.GetM(); i++)
      A.SetGlobal(i, 0, X(i));
  }

  
  //! First column of matrix A is assembled into the vector X
  template<class T, class Prop, class Allocator, class Allocator2>
  void GatherVector(DistributedMatrix<T, Prop, ColMajor, Allocator>& A,
                    Vector<T, VectFull, Allocator2>& X)
  {
    T val;
    for (int i = 0; i < X.GetM(); i++)
      {
        val = 0.0;
        A.GetGlobal(i, 0, val);
        X(i) = val;
      }
       
    Vector<T> Y(X);
    MPI_Allreduce(Y.GetData(), X.GetData(), X.GetM(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  
  //! Cholesky factorisation of a matrix A
  template<class Prop, class Allocator>
  void GetCholesky(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                   LapackInfo& info)
  {
    int n = A.GetGlobalM();
    char uplo('U'); int ia = 1;
    pdpotrf_(&uplo, &n, A.GetData(), &ia, &ia, A.GetDescriptor(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (A.GetRowProcNum() == 0)
      if (info.GetInfo() != 0)
        throw LapackError(info.GetInfo(), "GetCholesky",
                          "An error occured during the factorization.");
#endif
  }

  
  //! Cholesky resolution A y = x, result y overwrites x
  template<class Prop, class Allocator, class Allocator2>
  void SolveCholesky(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                     DistributedMatrix<double, Prop, ColMajor, Allocator2>& x,
                     LapackInfo& info)
  {
    int n = A.GetGlobalM();
    char uplo('U'); int ia = 1, ib = 1;
    int nrhs = x.GetGlobalN();
    pdpotrs_(&uplo, &n, &nrhs, A.GetData(), &ia, &ia, A.GetDescriptor(),
             x.GetData(), &ib, &ib, x.GetDescriptor(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (A.GetRowProcNum() == 0)
      if (info.GetInfo() != 0)
        throw LapackError(info.GetInfo(), "GetCholesky",
                          "An error occured during the resolution.");
#endif
  }
    

  template<class Prop, class Allocator>
  void GetEigenvalues(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                      Vector<double>& wr, Vector<double>& wi)
  {
    int n = A.GetGlobalM(), ilo = 1, ihi = n, ia = 1, lwork = -1, info;
    
    int nb = A.GetSizeBlockCol(), itmp = A.GetDescriptor()[6], jtmp = A.GetDescriptor()[7];
    int myrow = A.GetRowProcNum(), nprow = A.GetNumberOfRowProc();
    int mycol = A.GetColProcNum(), npcol = A.GetNumberOfColProc();
    int ihip = numroc_(&ihi, &nb, &myrow, &itmp, &nprow);
    int ilrow = indxg2p_(&ilo, &nb, &myrow, &itmp, &nprow);
    int ilcol = indxg2p_(&ilo, &nb, &mycol, &jtmp, &npcol);
    int ihlp = numroc_(&n, &nb, &myrow, &ilrow, &nprow);
    int inlq = numroc_(&n, &nb, &mycol, &ilcol, &npcol);
    lwork = nb*( nb + max( ihip+1, ihlp+inlq) );
    
    Vector<double> tau(n), work(lwork);
    Matrix<double, Prop, ColMajor> T(A.GetM(), A.GetN()); T.Zero();
    pdgehrd_(&n, &ilo, &ihi, A.GetData(), &ia, &ia, A.GetDescriptor(),
             tau.GetData(), work.GetData(), &lwork, &info);

    /* A.WriteText("Hess" + to_str(MPI_COMM_WORLD.Get_rank())+".dat");

    int zero_p(0);
    int size_r0 = numroc_(&n, &nb, &zero_p, &itmp, &nprow);
    int size_c0 = numroc_(&n, &nb, &zero_p, &jtmp, &npcol);

    int i0 = size_r0*myrow, j0 = size_c0*mycol;
    DISP(i0); DISP(j0);
    for (int i = i0; i < i0 + nprow; i++)
      for (int j = j0; j < j0 + npcol; j++)
        if (j < i-1)
        A(i-i0, j-j0) = 0.0; */

    for (int i = 0; i < n; i++)
      for (int j = 0; j < i-1; j++)
        A.SetGlobal(i, j, 0.0);
    
    wr.Reallocate(n); wi.Reallocate(n);
    wr.Zero(); wi.Zero();
    
    bool wantt = false, wantz = true; int ilwork = n;
    Vector<int> iwork(ilwork); work.Reallocate(lwork);
    pdlahqr_(&wantt, &wantz, &n, &ilo, &ihi, A.GetData(), A.GetDescriptor(),
             wr.GetData(), wi.GetData(), &ia, &ia, T.GetData(), A.GetDescriptor(),
             work.GetData(), &lwork, iwork.GetData(), &ilwork, &info); 
    
  }
  
}

#define SELDON_FILE_SCALAPACK_CXX
#endif

