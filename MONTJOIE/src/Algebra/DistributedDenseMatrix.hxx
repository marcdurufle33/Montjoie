#ifndef SELDON_FILE_DISTRIBUTED_DENSE_MATRIX_HXX

namespace Seldon
{
  
  //! basic distributed matrix (only rows are distributed)
  template<class T, class Prop, class Allocator>
  class DistributedMatrix<T, Prop, RowMajor, Allocator>
    : public Matrix<T, Prop, RowMajor, Allocator>
  {
  protected:
    MPI_Comm comm_;
    
  public:
    DistributedMatrix();
    DistributedMatrix(int i, int j);
    
    void SetCommunicator(const MPI_Comm&);
    const MPI_Comm& GetCommunicator() const;
    
  };

  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0, class T1, class T2>
  void MltVector(const SeldonTranspose& trans,
                 const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
                 const Vector<T1>& X, Vector<T2>& Y);
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class T1, class Prop1, class Allocator1, class T2, class T3, class T4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
                    const DistributedMatrix<T1, Prop1, RowMajor, Allocator1>& A,
                    const Vector<T2>& X, const T3& beta, Vector<T4>& Y);
  
  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y,
	      bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const T0& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);

  template<class T0, class T1, class Prop1, class Allocator1, class T2, class Prop2, class Allocator2,
           class T3, class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha, const SeldonTranspose& transA,
                    const DistributedMatrix<T1, Prop1, RowMajor, Allocator1>& A,
                    const SeldonTranspose& transB,
                    const DistributedMatrix<T2, Prop2, RowMajor, Allocator2>& B,
                    const T3& beta, Matrix<T4, Prop4, RowSymPacked, Allocator4>& C);
  
  template<class T0, class Prop1, class Allocator1,
           class Prop2, class Allocator2, class Prop4, class Allocator4>
  void MltAdd(const T0& alpha, const SeldonTranspose& transA,
              const DistributedMatrix<T0, Prop1, RowMajor, Allocator1>& A,
              const SeldonTranspose& transB,
              const DistributedMatrix<T0, Prop2, RowMajor, Allocator2>& B,
              const T0& beta, Matrix<T0, Prop4, RowSymPacked, Allocator4>& C);
  
}

#define SELDON_FILE_DISTRIBUTED_DENSE_MATRIX_HXX
#endif

