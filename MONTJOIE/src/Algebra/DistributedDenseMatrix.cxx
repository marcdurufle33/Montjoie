#ifndef SELDON_FILE_DISTRIBUTED_DENSE_MATRIX_CXX

namespace Seldon
{
  
  //! default constructor
  template<class T, class Prop, class Allocator>
  DistributedMatrix<T, Prop, RowMajor, Allocator>::DistributedMatrix()
    : Matrix<T, Prop, RowMajor, Allocator>()
  {
    comm_ = MPI_COMM_SELF;
  }

  
  //! constructor with the shape of the matrix
  template<class T, class Prop, class Allocator>
  DistributedMatrix<T, Prop, RowMajor, Allocator>::DistributedMatrix(int i, int j)
    : Matrix<T, Prop, RowMajor, Allocator>(i, j)
  {
    comm_ = MPI_COMM_SELF;
  }
  
  
  //! sets the MPI communicator
  template<class T, class Prop, class Allocator>
  void DistributedMatrix<T, Prop, RowMajor, Allocator>::SetCommunicator(const MPI_Comm& comm)
  {
    comm_ = comm;
  }

  
  //! returns the MPI communicator
  template<class T, class Prop, class Allocator>
  const MPI_Comm& DistributedMatrix<T, Prop, RowMajor, Allocator>::GetCommunicator() const
  {
    return comm_;
  }

  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, Y);
  }


  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<complex<T0>, Prop0, RowMajor, Allocator0>& >(A), X, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<complex<T0>, Prop0, RowMajor, Allocator0>& >(A), X, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0, class T1, class T2>
  void MltVector(const SeldonTranspose& trans,
                 const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
                 const Vector<T1>& X, Vector<T2>& Y)
  {
    if (trans.NoTrans())
      return Mlt(static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, Y);

    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc(1); MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return Mlt(trans, static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, Y);
    
    Y.Zero();
    T2 vloc;
    for (int i = 0; i < A.GetN(); i++)
      {
        SetComplexZero(vloc);
        if (trans.ConjTrans())
          for (int j = 0; j < A.GetM(); j++)
            vloc += conjugate(A(j, i))*X(j);
        else
          for (int j = 0; j < A.GetM(); j++)
            vloc += A(j, i)*X(j);
        
        Y(i) = vloc;
      }

    Vector<T2> Y0(Y);
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, Y0, xtmp, Y, Y.GetM(), MPI_SUM);    
  }


  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    MltVector(trans, A, X, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    MltVector(trans, A, X, Y);
  }
  

  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    MltVector(trans, A, X, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }
  

  template<class T0, class T1, class Prop1, class Allocator1, class T2, class T3, class T4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
                    const DistributedMatrix<T1, Prop1, RowMajor, Allocator1>& A,
                    const Vector<T2>& X, const T3& beta, Vector<T4>& Y)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc(1); MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return MltAddVector(alpha, trans, static_cast<const Matrix<T1, Prop1, RowMajor, Allocator1>& >(A), X, beta, Y);
    
    Vector<T2> Y0(Y.GetM());

    MltVector(trans, A, X, Y0);
    
    Y = beta*Y + alpha*Y0;
  }
  

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, beta, Y);
  }


  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
              Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<complex<T0>, Prop0, RowMajor, Allocator0>& >(A), X, beta, Y);
  }


  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const T0& beta,
	      Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<complex<T0>, Prop0, RowMajor, Allocator0>& >(A), X, beta, Y);
  }
  
  
  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<T0, Prop0, RowMajor, Allocator0>& >(A), X, beta, Y);
  }


  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<complex<T0>, Prop0, RowMajor, Allocator0>& >(A), X, beta, Y);
  }
  

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y);
  }
  

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y);
  }
  

  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y);
  }


  template<class T0, class Prop0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }


  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2, class T3, class T4, class Prop4, class Allocator4>
  void MltAddMatrix(const T0& alpha, const SeldonTranspose& transA,
                    const DistributedMatrix<T1, Prop1, RowMajor, Allocator1>& A,
                    const SeldonTranspose& transB,
                    const DistributedMatrix<T2, Prop2, RowMajor, Allocator2>& B,
                    const T3& beta, Matrix<T4, Prop4, RowSymPacked, Allocator4>& C)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc(1); MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return MltAddMatrix(alpha, transA, static_cast<const Matrix<T1, Prop1, RowMajor, Allocator1>& >(A),
                          transB, static_cast<const Matrix<T2, Prop2, RowMajor, Allocator2>& >(B), beta, C);
    
    if ((!transA.Trans()) && (!transB.NoTrans()))
      {
        cout << "Not implemented" << endl;
        abort();
      }

    int m = A.GetN();
    Matrix<T4, Prop4, RowSymPacked, Allocator4> C0(m, m);
    C0.Zero();

    T2 vloc;
    for (int i = 0; i < m; i++)
      for (int k = i; k < m; k++)
        {
          SetComplexZero(vloc);
          for (int j = 0; j < A.GetM(); j++)
            vloc += A(j, i)*B(j, k);
          
          C0(i, k) = vloc;
        }
    
    int size_C1 = m*(m+1)/2;
    Vector<int64_t> xtmp;

    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexOne(one);

    if (beta == zero)
      {
        MpiAllreduce(comm, C0.GetData(), xtmp, C.GetData(), size_C1, MPI_SUM);    
        if (alpha != one)
          Mlt(alpha, C);
      }
    else
      {
        Matrix<T4, Prop4, RowSymPacked, Allocator4> C1(m, m);
        
        MpiAllreduce(comm, C0.GetData(), xtmp, C1.GetData(), size_C1, MPI_SUM);    
        
        if (beta != one)
          Mlt(beta, C);
        
        Add(alpha, C1, C);
      }
  }
  
  
  template<class T0, class Prop1, class Allocator1,
           class Prop2, class Allocator2, class Prop4, class Allocator4>
  void MltAdd(const T0& alpha, const SeldonTranspose& transA,
              const DistributedMatrix<T0, Prop1, RowMajor, Allocator1>& A,
              const SeldonTranspose& transB,
              const DistributedMatrix<T0, Prop2, RowMajor, Allocator2>& B,
              const T0& beta, Matrix<T0, Prop4, RowSymPacked, Allocator4>& C)
  {
    MltAddMatrix(alpha, transA, A, transB, B, beta, C);
  }
  
}

#define SELDON_FILE_DISTRIBUTED_DENSE_MATRIX_CXX
#endif

