#ifndef SELDON_FILE_SCALAPACK_HXX

extern "C"
{
  void Cblacs_get(int ConTxt, int what, int *val);
  int Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);
  void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);
  
  int numroc_(int*, int*, int*, int*, int*);
  int indxg2p_(int*, int*, int*, int*, int*);
  int descinit_(int*, int*, int*, int*, int*,
                int*, int*, int*, int*, int*);
  
  void infog2l_(int*, int*, int*, int*, int*, int*, int*,
                int*, int*, int*, int*);
  
  void pdpotrf_(char*, int*, double*, int*, int*, int*, int*);  
  void pdpotrs_(char*, int*, int*, double*, int*, int*, int*,
                double*, int*, int*, int*, int*);

  void pdgehrd_(int*, int*, int*, double*, int*, int*, int*,
                double*, double*, int*, int*);
  
  void pdlahqr_(bool*, bool*, int*, int*, int*, double*, int*,
                double*, double*, int*, int*, double*, int*,
                double*, int*, int*, int*, int*);
}

namespace Seldon
{
  //! class for initialisation of Blacs grid
  class BlacsHandle
  {
  protected :
    int ictxt; //!< Blacs context (equivalent of MPI communicator)
    int nprow; //!< number of processors over which rows are distributed
    int npcol; //!< number of processors over which columns are distributed
    int myrow; //!< processor number for rows
    int mycol; //!< processor number for columns
    
  public :
    BlacsHandle();
    
    int GetContext() const;
    
    int GetNumOfRowProcs() const;
    int GetNumOfColProcs() const;

    int GetRowProcNum() const;
    int GetColProcNum() const;

    void Init();
    
  };

  
  //! Base class for a dense distributed matrix
  class DistributedMatrixDense_Base
  {
  protected :
    int ictxt_; //!< Blacs context
    int mb_, nb_; //!< size of blocks
    int mglob_; //!< global number of rows
    int nglob_; //!< global number of columns
    IVect desc_; //!< Scalapack descriptor
    int nprow_; //!< number of processors over which rows are distributed
    int npcol_; //!< number of processors over which columns are distributed
    int myrow_; //!< processor number for rows
    int mycol_; //!< processor number for columns
    
  public :
    
    DistributedMatrixDense_Base();
    
    //! returns Scalapack descriptor
    inline int* GetDescriptor() const { return desc_.GetData(); }
    
    //! returns the global number of rows
    inline int GetGlobalM() const { return mglob_; }
    //! returns the global number of columns
    inline int GetGlobalN() const { return nglob_; }

    //! returns the processor number for rows
    inline int GetRowProcNum() const { return myrow_; }
    //! returns the processor number for columns
    inline int GetColProcNum() const { return mycol_; }

    //! returns the number of processors for rows
    inline int GetNumberOfRowProc() const { return nprow_; }
    //! returns the processor number for columns
    inline int GetNumberOfColProc() const { return npcol_; }
    
    inline int GetSizeBlockRow() const { return mb_; }
    inline int GetSizeBlockCol() const { return nb_; }

  };
  

  //! Base class for a dense general distributed matrix
  template<class T, class Prop, class Storage, class Allocator
           = typename SeldonDefaultAllocator<Storage, T>::allocator >
  class DistributedMatrix_Pointers : public Matrix<T, Prop, Storage, Allocator>,
                                     public DistributedMatrixDense_Base
  {
  protected :

  public :
    void Init(BlacsHandle& bl, int m, int n, int mb = 64, int nb = 64);
    
    void SetGlobal(int i, int j, const T& val);
      
    void GetGlobal(int i, int j, T& val);
    
  };
  
  
  //! class for a general dense distributed matrix stored by columns
  template<class T, class Prop, class Allocator>
  class DistributedMatrix<T, Prop, ColMajor, Allocator>
    : public DistributedMatrix_Pointers<T, Prop, ColMajor, Allocator>
  {
    
  };

  void DistributeBlacsComputation(int nb_procs, int& ictxt,
                                  int& nprow, int& npcol,
                                  int& myrow, int& mycol);

  void GetScalapackDescriptor(int m, int n, int mb, int nb, int ictxt,
                              int nprow, int npcol, int myrow, int mycol,
                              IVect& desc, int& np, int& nq);
  
  void GetLocalScalapackIndex(int i, int j, const IVect& desc, int nprow, int npcol, 
                              int myrow, int mycol, int& iloc, int& jloc, int& iarow, int& jacol);
  
  template<class T, class Prop, class Allocator, class Allocator2>
  void DistributeVector(Vector<T, VectFull, Allocator>& X,
                        DistributedMatrix<T, Prop, ColMajor, Allocator2>& A);

  template<class T, class Prop, class Allocator, class Allocator2>
  void GatherVector(DistributedMatrix<T, Prop, ColMajor, Allocator>& A,
                    Vector<T, VectFull, Allocator2>& X);
  
  template<class Prop, class Allocator>
  void GetCholesky(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator, class Allocator2>
  void SolveCholesky(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                     DistributedMatrix<double, Prop, ColMajor, Allocator2>& x,
                     LapackInfo& info = lapack_info);
  
  template<class Prop, class Allocator>
  void GetEigenvalues(DistributedMatrix<double, Prop, ColMajor, Allocator>& A,
                      Vector<double>& wr, Vector<double>& wi);

}

#define SELDON_FILE_SCALAPACK_HXX
#endif

