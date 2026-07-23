#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_SEQUENTIAL_HXX

namespace Seldon
{

  class DistributedMatrixIntegerArray
  {
  };
  
  template<class T, class Prop, class Storage, class Allocator
	   = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class DistributedMatrix : public Matrix<T, Prop, Storage, Allocator>
  {
  public :
    typedef typename Matrix<T, Prop,
			    Storage, Allocator>::entry_type entry_type;
    
    DistributedMatrix();
    DistributedMatrix(int m, int n);
    
  };

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y,
	      bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const T0& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);


}

#define SELDON_FILE_DISTRIBUTED_MATRIX_SEQUENTIAL_HXX
#endif
