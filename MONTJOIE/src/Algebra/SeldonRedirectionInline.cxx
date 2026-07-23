#ifndef SELDON_FILE_SELDON_REDIRECTION_INLINE_CXX

namespace Seldon
{
  
  /*****************
   * Dense vectors *
   *****************/
  

  template<>
  inline void Vector<bool, VectFull, MallocAlloc<bool> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int, VectFull, MallocAlloc<int> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int64_t, VectFull, MallocAlloc<int64_t> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<double, VectFull, MallocAlloc<double> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<complex<double>, VectFull, MallocAlloc<complex<double> > >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }

  
  template<>
  inline void Vector<bool, VectFull, MallocAlloc<bool> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int, VectFull, MallocAlloc<int> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int64_t, VectFull, MallocAlloc<int64_t> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<double, VectFull, MallocAlloc<double> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<complex<double>, VectFull, MallocAlloc<complex<double> > >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }

  
  /******************
   * Sparse vectors *
   ******************/
  

  template<>
  inline void Vector<int, VectSparse, MallocAlloc<int> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int64_t, VectSparse, MallocAlloc<int64_t> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<double, VectSparse, MallocAlloc<double> >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<complex<double>, VectSparse, MallocAlloc<complex<double> > >::Reallocate(size_t n)
  {
    if (n < this->GetSize())
      {
	Clear();
	ReallocateVector(n);
      }
    else
      ReallocateVector(n);
  }

  
  template<>
  inline void Vector<int, VectSparse, MallocAlloc<int> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<int64_t, VectSparse, MallocAlloc<int64_t> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<double, VectSparse, MallocAlloc<double> >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }


  template<>
  inline void Vector<complex<double>, VectSparse, MallocAlloc<complex<double> > >::Resize(size_t n)
  {
    if (n < this->GetSize())
      ResizeVector(n);
    else
      ReallocateVector(n);
  }

  
  /*****************
   * Zero for mpfr *
   *****************/
  
  
#ifdef MONTJOIE_WITH_MPFR
  template<>
  inline void Vector<mpfr::mpreal>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Vector<complex<mpfr::mpreal> >::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Vector<TinyVector<mpfr::mpreal, 1 > >::Zero()
  {
    TinyVector<mpfr::mpreal, 1> zero;
    Fill(zero);
  }


  template<>
  inline void Vector<TinyVector<complex<mpfr::mpreal>, 1> >::Zero()
  {
    TinyVector<complex<mpfr::mpreal>, 1> zero;
    Fill(zero);
  }

  
  template<>
  inline void Vector<TinyVector<mpfr::mpreal, 2> >::Zero()
  {
    TinyVector<mpfr::mpreal, 2> zero;
    Fill(zero);
  }


  template<>
  inline void Vector<TinyVector<complex<mpfr::mpreal>, 2> >::Zero()
  {
    TinyVector<complex<mpfr::mpreal>, 2> zero;
    Fill(zero);
  }

  
  template<>
  inline void Vector<TinyVector<mpfr::mpreal, 3> >::Zero()
  {
    TinyVector<mpfr::mpreal, 3> zero;
    Fill(zero);
  }


  template<>
  inline void Vector<TinyVector<complex<mpfr::mpreal>, 3> >::Zero()
  {
    TinyVector<complex<mpfr::mpreal>, 3> zero;
    Fill(zero);
  }

  
  template<>
  inline void Matrix_Pointers<mpfr::mpreal, General, RowMajor>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Pointers<complex<mpfr::mpreal>, General, RowMajor>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Pointers<mpfr::mpreal, General, ColMajor>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Pointers<complex<mpfr::mpreal>, General, ColMajor>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Symmetric<mpfr::mpreal, Symmetric, RowSym>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Symmetric<complex<mpfr::mpreal>, Symmetric, RowSym>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_SymPacked<mpfr::mpreal, Symmetric, RowSymPacked>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_SymPacked<complex<mpfr::mpreal>, Symmetric, RowSymPacked>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_HermPacked<mpfr::mpreal, Hermitian, RowHermPacked>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_HermPacked<complex<mpfr::mpreal>, Hermitian, RowHermPacked>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Sparse<mpfr::mpreal, General, RowSparse>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Sparse<complex<mpfr::mpreal>, General, RowSparse>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_SymSparse<mpfr::mpreal, Symmetric, RowSymSparse>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_SymSparse<complex<mpfr::mpreal>, Symmetric, RowSymSparse>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }

  template<>
  inline void Matrix_ComplexSparse<complex<mpfr::mpreal>, General, RowComplexSparse>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_SymComplexSparse<complex<mpfr::mpreal>, Symmetric, RowSymComplexSparse>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Band<mpfr::mpreal, General, BandedCol>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Band<complex<mpfr::mpreal>, General, BandedCol>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Arrow<mpfr::mpreal, General, ArrowCol>::Zero()
  {
    mpfr::mpreal zero(0);
    Fill(zero);
  }


  template<>
  inline void Matrix_Arrow<complex<mpfr::mpreal>, General, ArrowCol>::Zero()
  {
    complex<mpfr::mpreal> zero(0, 0);
    Fill(zero);
  }
#endif
  
  
  /*************
   * MltVector *
   *************/


  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltAddVector(const Real_wp& alpha, const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltAddVector(const Complex_wp& alpha, const Vector<Complex_wp>& x,
		 const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  { 
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Complex_wp>& x, const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltVector(const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltVector(const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, RowMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Real_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltAddVector(const Real_wp& alpha, const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltAddVector(const Complex_wp& alpha, const Vector<Complex_wp>& x,
		 const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Complex_wp>& x, const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltVector(const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltVector(const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, RowMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Complex_wp, General, RowMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltAddVector(const Real_wp& alpha, const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltAddVector(const Complex_wp& alpha, const Vector<Complex_wp>& x,
		 const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Complex_wp>& x, const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltVector(const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltVector(const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Real_wp, General, ColMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Real_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltAddVector(const Real_wp& alpha, const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltAddVector(const Complex_wp& alpha, const Vector<Complex_wp>& x,
		 const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose& trans,
		 const Vector<Complex_wp>& x, const Complex_wp& beta, Vector<Complex_wp>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this),
	   x, beta, y);
  }

  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltVector(const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltVector(const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Real_wp>& x, Vector<Real_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this), x, y);
  }


  template <>
  inline void Matrix_Pointers<Complex_wp, General, ColMajor>
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Complex_wp>& x, Vector<Complex_wp>& y) const
  {
    Mlt(trans, static_cast<const Matrix<Complex_wp, General, ColMajor>& >(*this), x, y);
  }

};

#define SELDON_FILE_SELDON_REDIRECTION_INLINE_CXX
#endif

