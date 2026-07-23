#ifndef SELDON_FILE_EIGENVALUE_HXX

namespace Seldon
{

#ifndef SELDON_WITH_LAPACK
  
  /***************************************
   * C++ translations of Lapack routines *
   ***************************************/
  
  ////////////////////////////
  // General complex matrix //
  ////////////////////////////
  
  
  template<class Matrix1, class Vector1, class reel>
  void TransformMatrix_to_HessenbergForm(int low, int high, Matrix1& A,
                                         Vector1& Q, const reel& ftmp);
  
  template<class Matrix1, class Vector1, class reel>
  int ComputeEigenValuesQR_Complex(int low, int high, Matrix1& H, Vector1& W,
                                   const reel& norm_tmp);

  template<class Matrix1, class Vector1, class Matrix2, class Matrix3, class reel>
  int ComputeEigenValuesQR_Complex(int low, int high, Matrix1& H, Matrix2& Q,
                                   Vector1& W, Matrix3& Z, const reel& norm_tmp);
  
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<T>, Prop, Storage, Allocator1>& A,
		      Vector<complex<T>, VectFull, Allocator2>& L,
		      LapackInfo& info = lapack_info);
  
  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<T>, General, Storage, Allocator1>& A,
				  Vector<complex<T>, VectFull, Allocator2>& L,
				  Matrix<complex<T>, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info = lapack_info);
  
  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Allocator3, class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<T, General, Storage, Allocator1>& A,
				  Vector<T, VectFull, Allocator2>& Lreal,
				  Vector<T, VectFull, Allocator3>& Limag,
				  Matrix<T, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info = lapack_info);

  template<class T, class Prop, class Storage, class Allocator1, class Allocator2,
	   class Allocator3>
  void GetEigenvalues(Matrix<T, Prop, Storage, Allocator1>& A,
                      Vector<T, VectFull, Allocator2>& Lreal,
                      Vector<T, VectFull, Allocator3>& Limag,
		      LapackInfo& info = lapack_info);

  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<T, Prop, Storage, Allocator1>& A,
                      Vector<T, VectFull, Allocator2>& Lreal,
		      LapackInfo& info = lapack_info);

  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<T, Symmetric, Storage, Allocator1>& A,
				  Vector<T, VectFull, Allocator2>& L,
				  Matrix<T, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info = lapack_info);

  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<T>, Symmetric, Storage, Allocator1>& A,
				  Vector<complex<T>, VectFull, Allocator2>& L,
				  Matrix<complex<T>, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info = lapack_info);
  
#endif

  
  template<class Storage>
  class ClassStorageTypes
  {
  public :
    typedef RowMajor DenseStorage;
  };

  template<>
  class ClassStorageTypes<ColSym>
  {
  public :
    typedef ColMajor DenseStorage;
  };

  template<>
  class ClassStorageTypes<ColSymPacked>
  {
  public :
    typedef ColMajor DenseStorage;
  };
  
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrixReal(Matrix<T, General, Storage, Allocator1>& A, complex<T> (&f)(const complex<T>&));
    
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrix(Matrix<complex<T>, General, Storage, Allocator1>& A, complex<T> (&f)(const complex<T>&));

  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrixReal(Matrix<T, Symmetric, Storage, Allocator1>& A, complex<T> (&f)(const complex<T>&));
    
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrix(Matrix<complex<T>, Symmetric, Storage, Allocator1>& A, complex<T> (&f)(const complex<T>&));

  template<class T, class Prop, class Storage, class Allocator>
  void GetAbsoluteValue(Matrix<T, Prop, Storage, Allocator>& A, bool take_abs=true);

  template<class T, class Prop, class Storage, class Allocator>
  void GetAbsoluteValue(Matrix<complex<T>, Prop, Storage, Allocator>& A, bool take_abs=true);

  template<class T, class Prop, class Storage, class Allocator>
  void GetExponential(Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void GetExponential(Matrix<complex<T>, Prop, Storage, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void GetSquareRoot(Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void GetSquareRoot(Matrix<complex<T>, Prop, Storage, Allocator>& A);

  template<class T, int m>
  void GetAbsoluteValue(TinyMatrix<T, Symmetric, m, m>& A, bool take_abs = true);
  
  template<class T, int m>
  void GetAbsoluteValue(TinyMatrix<T, General, m, m>& A, bool take_abs=true);

}

#define SELDON_FILE_EIGENVALUE_HXX
#endif
