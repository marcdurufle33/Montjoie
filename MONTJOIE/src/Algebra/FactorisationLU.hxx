#ifndef SELDON_FILE_FACTORISATION_LU_HXX

namespace Seldon
{
  
#ifndef SELDON_WITH_LAPACK
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T, General, Storage, Allocator>& A, IVect& IPivot,
	     LapackInfo& info = lapack_info);
  
  template<class T, class Allocator,
           class Storage, class Prop, class Allocator2>
  void SolveLuVector(const SeldonTranspose& TransA,
                     const Matrix<T, General, Storage, Allocator>& A,
                     const IVect& IPivot, Vector<T, Prop, Allocator2>& b,
		     LapackInfo& info = lapack_info);
  
  template<class T, class Allocator, class Storage, class Prop, class Allocator2>
  void SolveLuVector(const Matrix<T, General, Storage, Allocator>& A,
                     const IVect& IPivot, Vector<T, Prop, Allocator2>& b,
		     LapackInfo& info = lapack_info);
  
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T, Symmetric, Storage, Allocator>& A, IVect& IPivot,
	     LapackInfo& info = lapack_info);
  
  template<class T, class Allocator, class Storage, class Prop, class Allocator2>
  void SolveLuVector(const Matrix<T, Symmetric, Storage, Allocator>& A,
                     const IVect& IPivot, Vector<T, Prop, Allocator2>& b,
		     LapackInfo& info = lapack_info);
  
  template<class T, class Allocator>
  void GetInverse(Matrix<T, Symmetric, RowSymPacked, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Storage, class Allocator>
  void GetInverse(Matrix<T, General, Storage, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, RowUpTriang, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, RowUpTriangPacked, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, ColUpTriang, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, ColUpTriangPacked, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, RowLoTriang, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, RowLoTriangPacked, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, ColLoTriang, Allocator>& A,
		  LapackInfo& info = lapack_info);

  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, ColLoTriangPacked, Allocator>& A,
		  LapackInfo& info = lapack_info);
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<complex<T>, Prop, Storage, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class T, class Prop,
           class Storage, class Allocator,
	   class T2, class Storage2, class Allocator2>
  void SolveCholesky(const SeldonTranspose& TransA,
		     const Matrix<T, Prop, Storage, Allocator>& A,
		     Vector<T2, Storage2, Allocator2>& x,
		     LapackInfo& info = lapack_info);
  
  template<class T,
           class Prop, class Storage, class Allocator,
	   class T2, class Storage2, class Allocator2>
  void MltCholesky(const SeldonTranspose& TransA,
                   const Matrix<T, Prop, Storage, Allocator>& A,
                   Vector<T2, Storage2, Allocator2>& x,
		   LapackInfo& info = lapack_info);

  template<class T, class Prop,
           class Storage, class Allocator,
	   class T2, class Storage2, class Allocator2>
  void SolveCholesky(const SeldonTranspose& TransA,
		     const Matrix<complex<T>, Prop, Storage, Allocator>& A,
		     Vector<complex<T2>, Storage2, Allocator2>& x,
		     LapackInfo& info = lapack_info);
  
  template<class T,
           class Prop, class Storage, class Allocator,
	   class T2, class Storage2, class Allocator2>
  void MltCholesky(const SeldonTranspose& TransA,
                   const Matrix<complex<T>, Prop, Storage, Allocator>& A,
                   Vector<complex<T2>, Storage2, Allocator2>& x,
		   LapackInfo& info = lapack_info);
  
  template<class T, class Prop, class Storage, class Allocator, class Allocator2>
  void GetQR(Matrix<T, Prop, Storage, Allocator>& A, Vector<T, VectFull, Allocator2>& tau,
	     LapackInfo& info = lapack_info);
  
  template<class T, class Prop, class Storage, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void SolveQR(const Matrix<T, Prop, Storage, Allocator>& A,
               const Vector<T, VectFull, Allocator2>& tau,
               Vector<T2, VectFull, Allocator3>& X,
	       LapackInfo& info = lapack_info);

  template<class IsTranspose, class T, class Prop,
	   class Storage, class Allocator,
           class T2, class Allocator2, class Allocator3>
  void MltQ_FromQR(const IsTranspose& trans,
		   const Matrix<T, Prop, Storage, Allocator>& A,
		   const Vector<T, VectFull, Allocator2>& tau,
		   Vector<T2, VectFull, Allocator3>& X,
		   LapackInfo& info = lapack_info);
  
  template<class T, class Allocator1, class Storage, class Allocator2>
  void SolveUpper(const Matrix<T, General, Storage, Allocator1>& A,
                  Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void SolveUpper(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
                  const Matrix<T, General, Storage, Allocator1>& A,
                  Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void SolveLower(const Matrix<T, General, Storage, Allocator1>& A,
                  Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void SolveLower(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
                  const Matrix<T, General, Storage, Allocator1>& A,
                  Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, ColUpTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColUpTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, ColUpTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColUpTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, RowUpTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowUpTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, RowUpTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowUpTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, ColLoTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColLoTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, ColLoTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColLoTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, RowLoTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowLoTriang, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const Matrix<T, General, RowLoTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Solve(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowLoTriangPacked, Allocator1>& A,
             Vector<T, VectFull, Allocator2>& x);

  template<class T, class Prop, class Storage, class Allocator,
           class Prop0, class Storage0, class Allocator0>
  void Solve(const SeldonSide& side, const T& alpha,
             const Matrix<T, Prop, Storage, Allocator>& A,
             Matrix<T, Prop0, Storage0, Allocator0>& B);

  template<class T, class Prop, class Storage, class Allocator,
           class Prop0, class Storage0, class Allocator0>
  void Solve(const SeldonSide& side, const T& alpha,
             const SeldonTranspose& TransA,
             const SeldonDiag& DiagA,
             const Matrix<T, Prop, Storage, Allocator>& A,
             Matrix<T, Prop0, Storage0, Allocator0>& B);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void MltUpper(const Matrix<T, General, Storage, Allocator1>& A,
		Vector<T, VectFull, Allocator2>& x);


  template<class T, class Allocator1, class Storage, class Allocator2>
  void MltUpper(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
                const Matrix<T, General, Storage, Allocator1>& A,
                Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void MltLower(const Matrix<T, General, Storage, Allocator1>& A,
		Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Storage, class Allocator2>
  void MltLower(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
                const Matrix<T, General, Storage, Allocator1>& A,
                Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, ColUpTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	   const Matrix<T, General, ColUpTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, ColUpTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColUpTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, RowUpTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowUpTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, RowUpTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowUpTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, ColLoTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColLoTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, ColLoTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, ColLoTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, RowLoTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowLoTriang, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Allocator1, class Allocator2>
  void Mlt(const Matrix<T, General, RowLoTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);
  
  template<class T, class Allocator1, class Allocator2>
  void Mlt(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
             const Matrix<T, General, RowLoTriangPacked, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x);

  template<class T, class Prop, class Storage, class Allocator,
           class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonSide& side, const T& alpha,
           const Matrix<T, Prop, Storage, Allocator>& A,
           Matrix<T, Prop0, Storage0, Allocator0>& B);

  template<class T, class Prop, class Storage, class Allocator,
           class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonSide& side, const T& alpha,
           const SeldonTranspose& TransA,
           const SeldonDiag& DiagA,
           const Matrix<T, Prop, Storage, Allocator>& A,
           Matrix<T, Prop0, Storage0, Allocator0>& B);


  template<class T, class Alloc1, class Alloc2, class Storage, class Alloc0>
  void Rank1Update(const T& alpha, const Vector<T, VectFull, Alloc1>& X,
                   const Vector<T, VectFull, Alloc2>& Y,
                   Matrix<T, General, Storage, Alloc0>& A);
  
  template<class T, class Alloc1, class Alloc2, class Storage, class Alloc0>
  void Rank1Update(const T& alpha, const Vector<T, VectFull, Alloc1>& X,
                   const SeldonConjugate& ConjY, const Vector<T, VectFull, Alloc2>& Y,
                   Matrix<T, General, Storage, Alloc0>& A);

  template<class T, class Alloc1, class Storage, class Alloc0>
  void Rank1Update(const T& alpha, const Vector<T, VectFull, Alloc1>& X,
                   Matrix<T, Symmetric, Storage, Alloc0>& A);

  template<class T, class Alloc1, class Storage, class Alloc0>
  void Rank1Update(const T& alpha, const Vector<complex<T>, VectFull, Alloc1>& X,
                   Matrix<complex<T>, Hermitian, Storage, Alloc0>& A);

  template<class T, class Alloc1, class Alloc2, class Storage, class Alloc0>
  void Rank2Update(const T& alpha, const Vector<T, VectFull, Alloc1>& X,
                   const Vector<T, VectFull, Alloc2>& Y,
                   Matrix<T, Symmetric, Storage, Alloc0>& A);

  template<class T, class Alloc1, class Alloc2, class Storage, class Alloc0>
  void Rank2Update(const T& alpha, const Vector<complex<T>, VectFull, Alloc1>& X,
                   const Vector<complex<T>, VectFull, Alloc2>& Y,
                   Matrix<complex<T>, Hermitian, Storage, Alloc0>& A);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3, class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const SeldonSide& Side, const T0& alpha,
              const Matrix<T1, Prop1, Storage1, Allocator1>& A,
              const Matrix<T2, Prop2, Storage2, Allocator2>& B,
              const T3& beta, Matrix<T4, Prop4, Storage4, Allocator4>& C);
  
#else
  template<class T, class Prop, class Allocator1,
           class Storage2, class Allocator2>
  void SolveCholesky(const SeldonTranspose& TransA,
		     const Matrix<T, Prop, RowSymPacked, Allocator1>& A,
		     Vector<complex<T>, Storage2, Allocator2>& x,
		     LapackInfo& info = lapack_info);
#endif

  template<class T, class Prop, class Storage, class Allocator,
           class Allocator2, class Allocator3>
  void SolveQR(const Matrix<T, Prop, Storage, Allocator>& A,
               const Vector<T, VectFull, Allocator2>& tau,
               Vector<complex<T>, VectFull, Allocator3>& x,
	       LapackInfo& info = lapack_info);
  
  template<class T>
  T GetHouseholderNormale(Vector<T>& U);

  template<class T>
  void GetReorthogonalization(Vector<Vector<T> >& H, int size, Vector<T>& Q);
  
  template<class T>
  void ApplyHouseholderTransformation(const Vector<T>& N, int k, Vector<T>& X);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  void GetPseudoInverse(Matrix<T0, Prop1, Storage1, Allocator1>& A, const T1& epsilon,
			LapackInfo& info = lapack_info);
  
}

#define SELDON_FILE_FACTORISATION_LU_HXX
#endif

