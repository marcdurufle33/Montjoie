#ifndef MONTJOIE_FILE_NON_LINEAR_LEAST_SQUARES_HXX

#ifdef SELDON_WITH_MKL
extern "C"
{
  typedef void* _TRNSP_HANDLE_t;
  
  extern int dtrnlsp_init     (_TRNSP_HANDLE_t*, int*, int*, double*, double*, int*, int*, double*);
  extern int dtrnlsp_check    (_TRNSP_HANDLE_t*, int*, int*, double*, double*, double*, int*);
  extern int dtrnlsp_solve    (_TRNSP_HANDLE_t*, double*, double*, int*);
  extern int dtrnlsp_get      (_TRNSP_HANDLE_t*, int*, int*, double*, double*);
  extern int dtrnlsp_delete   (_TRNSP_HANDLE_t*);
}
#endif

namespace Montjoie
{

  //! base class to specify a least squares problem for mkl
  template<class T, class VectorSol = Vector<T>, class VectorRhs = Vector<T>,
	   class MatrixJac = Matrix<T, General, ColMajor>,
	   class MatrixAAt = Matrix<T, Symmetric, RowSymPacked> >
  class VirtualLeastSquaresFunction
  {
  protected:
    int m_, n_;
    
  public:
    inline VirtualLeastSquaresFunction() { m_ = 0; n_ = 0; }
    
    inline int GetM() { return m_; }
    inline int GetN() { return n_; }
    
    virtual void FindInitGuess(VectorSol& param) = 0;
    virtual void EvaluateF(const VectorSol& x, VectorRhs& feval) = 0;
    virtual void EvaluateJacobian(const VectorSol& x, VectorRhs& feval, MatrixJac& fjac) = 0;
    
    virtual inline T GetNorm2_Rhs(const VectorRhs& f) { return Norm2(f); }
    virtual void InitJacobian(MatrixJac& fjac) {}
    virtual void GiveIterate(int step, VectorSol& x, const T& normRes) {}
    
  };
  
  
  template<class T, class VectorSol, class VectorRhs, class MatrixJac, class MatrixAAt>
  T SolveLeastSquaresLvm(VirtualLeastSquaresFunction<T, VectorSol, VectorRhs, MatrixJac, MatrixAAt>& fct,
			 VectorSol& xsol, Real_wp epsilon1 = 1e-12, Real_wp epsilon2 = -1,
                         unsigned nb_max_iter = 5000, int print_level = 0, Real_wp tau = 1e-3);
  
#ifdef SELDON_WITH_MKL
  template<class T>
  T SolveLeastSquaresMkl(VirtualLeastSquaresFunction<T>& fct, Vector<T>& xsol,
                         Real_wp epsilon = 1e-12, unsigned nb_max_iter = 5000);
#endif
    
}

#define MONTJOIE_FILE_NON_LINEAR_LEAST_SQUARES_HXX
#endif
