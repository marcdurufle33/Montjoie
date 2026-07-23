#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_SEQUENTIAL_INLINE_CXX

namespace Seldon
{

  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>
  ::DistributedMatrix() : Matrix<T, Prop, Storage, Allocator>()
  {
  }
  
  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>
  ::DistributedMatrix(int m, int n) : Matrix<T, Prop, Storage, Allocator>(m, n)
  {
  }
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {    
    Mlt(static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {    
    Mlt(trans, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(trans, static_cast<const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble)
  {
    Mlt(trans, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A), X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y,
		     bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A),
	   X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& >(A),
	   X, beta, Y);
  }
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const T0& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& >(A),
	   X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A),
	   X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAdd(alpha, trans, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A)
	   , X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, trans, static_cast<const Matrix<complex<T0>, Prop0, Storage0, Allocator0>& >(A)
	   , X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
 		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, trans, static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A)
	   , X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    cout << "Incompatible types" << endl;
    abort();
  }

}

#define SELDON_FILE_DISTRIBUTED_MATRIX_SEQUENTIAL_INLINE_CXX
#endif
