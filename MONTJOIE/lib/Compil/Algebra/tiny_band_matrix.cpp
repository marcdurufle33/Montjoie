#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#include "Algebra/TinyBlockSolver1D.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/TinyBandMatrix.cxx"
#include "Algebra/TinyBlockSolver1D.cxx"
#endif

namespace Seldon
{
  /* TinyBandMatrix */

  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 1>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 1>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 1>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 1>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 1>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template ostream& operator <<(ostream&, const TinyBandMatrix<Real_wp, 1>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Real_wp, 1>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Real_wp, 1>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 1>&, TinyBandMatrix<Real_wp, 1>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 1>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Real_wp, 1>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Real_wp, 1>&, TinyBandMatrix<Real_wp, 1>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 1>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 1>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Real_wp, 1>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Real_wp, 1>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyBandMatrix<Real_wp, 1>&);
  
  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 2>;
  
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 2>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 2>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 2>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 2>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Real_wp, 2>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Real_wp, 2>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 2>&, TinyBandMatrix<Real_wp, 2>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 2>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Real_wp, 2>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Real_wp, 2>&, TinyBandMatrix<Real_wp, 2>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 2>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 2>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 2>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Real_wp, 2>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Real_wp, 2>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyBandMatrix<Real_wp, 2>&);


  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 3>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 3>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 3>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 3>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 3>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Real_wp, 3>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Real_wp, 3>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 3>&, TinyBandMatrix<Real_wp, 3>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 3>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Real_wp, 3>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Real_wp, 3>&, TinyBandMatrix<Real_wp, 3>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 3>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 3>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Real_wp, 3>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Real_wp, 3>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyBandMatrix<Real_wp, 3>&);


  SELDON_EXTERN template class TinyBandMatrix<Complex_wp, 2>;

  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 2>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 2>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 2>::Fill(const Complex_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 2>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Complex_wp, 2>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Complex_wp, 2>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 2>&, TinyBandMatrix<Complex_wp, 2>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, TinyBandMatrix<Complex_wp, 2>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Complex_wp, 2>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Complex_wp, 2>&, TinyBandMatrix<Complex_wp, 2>&);
  SELDON_EXTERN template void Add(const Complex_wp&, const TinyBandMatrix<Complex_wp, 2>&, TinyBandMatrix<Complex_wp, 2>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Complex_wp, 2>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const Complex_wp&, TinyBandMatrix<Complex_wp, 2>&);


  SELDON_EXTERN template class TinyBandMatrix<Complex_wp, 3>;

  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 3>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 3>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 3>::Fill(const Complex_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 3>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Complex_wp, 3>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Complex_wp, 3>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 3>&, TinyBandMatrix<Complex_wp, 3>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, TinyBandMatrix<Complex_wp, 3>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Complex_wp, 3>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Complex_wp, 3>&, TinyBandMatrix<Complex_wp, 3>&);
  SELDON_EXTERN template void Add(const Complex_wp&, const TinyBandMatrix<Complex_wp, 3>&, TinyBandMatrix<Complex_wp, 3>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Complex_wp, 3>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const Complex_wp&, TinyBandMatrix<Complex_wp, 3>&);

  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 4>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 4>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 4>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 4>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 4>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 4>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 4>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 4>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 6>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 6>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 6>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 6>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 6>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 6>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 6>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 6>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  
  /*SELDON_EXTERN template class TinyBandMatrix<Real_wp, 8>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 8>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 8>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 8>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 8>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 8>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 8>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 8>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 9>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 9>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 9>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 9>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 9>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Real_wp, 9>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Real_wp, 9>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 9>&, TinyBandMatrix<Real_wp, 9>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 9>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 9>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Real_wp, 9>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Real_wp, 9>&, TinyBandMatrix<Real_wp, 9>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 9>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 9>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 9>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Real_wp, 9>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Real_wp, 9>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyBandMatrix<Real_wp, 9>&);


  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 10>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 10>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 10>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 10>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 10>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Real_wp, 10>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Real_wp, 10>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 10>&, TinyBandMatrix<Real_wp, 10>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Real_wp, 10>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 10>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Real_wp, 10>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Real_wp, 10>&, TinyBandMatrix<Real_wp, 10>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 10>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 10>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 10>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Real_wp, 10>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Real_wp, 10>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyBandMatrix<Real_wp, 10>&);


  SELDON_EXTERN template class TinyBandMatrix<Complex_wp, 10>;

  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 10>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 10>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 10>::Fill(const Complex_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Complex_wp, 10>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyBandMatrix<Complex_wp, 10>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyBandMatrix<Complex_wp, 10>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 10>&, TinyBandMatrix<Complex_wp, 10>&, bool);
  SELDON_EXTERN template void GetLU(TinyBandMatrix<Complex_wp, 10>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, TinyBandMatrix<Complex_wp, 10>&);
  SELDON_EXTERN template void SolveLU(TinyBandMatrix<Complex_wp, 10>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyBandMatrix<Complex_wp, 10>&, TinyBandMatrix<Complex_wp, 10>&);
  SELDON_EXTERN template void Add(const Complex_wp&, const TinyBandMatrix<Complex_wp, 10>&, TinyBandMatrix<Complex_wp, 10>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const TinyBandMatrix<Complex_wp, 10>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const Complex_wp&, TinyBandMatrix<Complex_wp, 10>&); */

  /* SELDON_EXTERN template class TinyBandMatrix<Real_wp, 12>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 12>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 12>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 12>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 12>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 12>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 12>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 12>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 14>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 14>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 14>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 14>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 14>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 14>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 14>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 14>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBandMatrix<Real_wp, 16>;

  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 16>::Fill(const int&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 16>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 16>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBandMatrix<Real_wp, 16>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyBandMatrix<Real_wp, 16>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 16>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyBandMatrix<Real_wp, 16>&, const Vector<Complex_wp>&, const Real_wp&, Vector<Complex_wp>&);
  */  
  
  /* TinyArrowMatrix */

  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 2, 1>;
  SELDON_EXTERN template void TinyArrowMatrix<Real_wp, 2, 1>::Fill(const int&);
  SELDON_EXTERN template void TinyArrowMatrix<Real_wp, 2, 1>::Fill(const Real_wp&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyArrowMatrix<Real_wp, 2, 1>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyArrowMatrix<Real_wp, 2, 1>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Real_wp, 2, 1>&, TinyArrowMatrix<Real_wp, 2, 1>&, bool);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Real_wp, 2, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 2, 1>&);
  SELDON_EXTERN template void SolveLU(TinyArrowMatrix<Real_wp, 2, 1>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyArrowMatrix<Real_wp, 2, 1>&, TinyArrowMatrix<Real_wp, 2, 1>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyArrowMatrix<Real_wp, 2, 1>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyArrowMatrix<Real_wp, 2, 1>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyArrowMatrix<Real_wp, 2, 1>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyArrowMatrix<Real_wp, 2, 1>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyArrowMatrix<Real_wp, 2, 1>&);


  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 3, 5>;
  SELDON_EXTERN template void TinyArrowMatrix<Real_wp, 3, 5>::Fill(const int&);
  SELDON_EXTERN template void TinyArrowMatrix<Real_wp, 3, 5>::Fill(const Real_wp&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyArrowMatrix<Real_wp, 3, 5>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyArrowMatrix<Real_wp, 3, 5>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Real_wp, 3, 5>&, TinyArrowMatrix<Real_wp, 3, 5>&, bool);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Real_wp, 3, 5>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 3, 5>&);
  SELDON_EXTERN template void SolveLU(TinyArrowMatrix<Real_wp, 3, 5>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const TinyArrowMatrix<Real_wp, 3, 5>&, TinyArrowMatrix<Real_wp, 3, 5>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const TinyArrowMatrix<Real_wp, 3, 5>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const TinyArrowMatrix<Real_wp, 3, 5>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyArrowMatrix<Real_wp, 3, 5>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const TinyArrowMatrix<Real_wp, 3, 5>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void Mlt(const Real_wp&, TinyArrowMatrix<Real_wp, 3, 5>&);


  SELDON_EXTERN template class TinyArrowMatrix<Complex_wp, 2, 1>;
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 2, 1>::Fill(const int&);
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 2, 1>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 2, 1>::Fill(const Complex_wp&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyArrowMatrix<Complex_wp, 2, 1>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyArrowMatrix<Complex_wp, 2, 1>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Complex_wp, 2, 1>&, TinyArrowMatrix<Complex_wp, 2, 1>&, bool);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Complex_wp, 2, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Complex_wp, 2, 1>&);
  SELDON_EXTERN template void SolveLU(TinyArrowMatrix<Complex_wp, 2, 1>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Add(const Complex_wp&, const TinyArrowMatrix<Complex_wp, 2, 1>&, TinyArrowMatrix<Complex_wp, 2, 1>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const TinyArrowMatrix<Complex_wp, 2, 1>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const SeldonTranspose&, const TinyArrowMatrix<Complex_wp, 2, 1>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyArrowMatrix<Complex_wp, 2, 1>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const TinyArrowMatrix<Complex_wp, 2, 1>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const Complex_wp&, TinyArrowMatrix<Complex_wp, 2, 1>&);

  SELDON_EXTERN template class TinyArrowMatrix<Complex_wp, 3, 5>;
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 3, 5>::Fill(const int&);
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 3, 5>::Fill(const Real_wp&);
  SELDON_EXTERN template void TinyArrowMatrix<Complex_wp, 3, 5>::Fill(const Complex_wp&);
  
  SELDON_EXTERN template void GetRowSum(Vector<Real_wp>&, const TinyArrowMatrix<Complex_wp, 3, 5>&);
  SELDON_EXTERN template void ScaleLeftMatrix(TinyArrowMatrix<Complex_wp, 3, 5>&, const Vector<Real_wp>&);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Complex_wp, 3, 5>&, TinyArrowMatrix<Complex_wp, 3, 5>&, bool);
  SELDON_EXTERN template void GetLU(TinyArrowMatrix<Complex_wp, 3, 5>&);
  SELDON_EXTERN template void Copy(const Matrix<Complex_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Complex_wp, 3, 5>&);
  SELDON_EXTERN template void SolveLU(TinyArrowMatrix<Complex_wp, 3, 5>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Add(const Complex_wp&, const TinyArrowMatrix<Complex_wp, 3, 5>&, TinyArrowMatrix<Complex_wp, 3, 5>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const TinyArrowMatrix<Complex_wp, 3, 5>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void MltAdd(const Complex_wp&, const SeldonTranspose&, const TinyArrowMatrix<Complex_wp, 3, 5>&, const Vector<Complex_wp>&, const Complex_wp&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const TinyArrowMatrix<Complex_wp, 3, 5>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const TinyArrowMatrix<Complex_wp, 3, 5>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void Mlt(const Complex_wp&, TinyArrowMatrix<Complex_wp, 3, 5>&);

  /* TinyBlockSolver1D */
  
  SELDON_EXTERN template class TinyBlockSolver1D<Real_wp, 1>;
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 1>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 1>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Real_wp, 3>;
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 3>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 3>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBlockSolver1D<Real_wp, 5>;
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 5>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 5>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBlockSolver1D<Real_wp, 7>;
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 7>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 7>::Solve(Vector<Complex_wp>&);

  SELDON_EXTERN template class TinyBlockSolver1D<Real_wp, 9>;
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 9>::Solve(Vector<Real_wp>&);
  SELDON_EXTERN template void TinyBlockSolver1D<Real_wp, 9>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Complex_wp, 1>;
  SELDON_EXTERN template void TinyBlockSolver1D<Complex_wp, 1>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Complex_wp, 3>;
  SELDON_EXTERN template void TinyBlockSolver1D<Complex_wp, 3>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Complex_wp, 5>;
  SELDON_EXTERN template void TinyBlockSolver1D<Complex_wp, 5>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Complex_wp, 7>;
  SELDON_EXTERN template void TinyBlockSolver1D<Complex_wp, 7>::Solve(Vector<Complex_wp>&);
  
  SELDON_EXTERN template class TinyBlockSolver1D<Complex_wp, 9>;
  SELDON_EXTERN template void TinyBlockSolver1D<Complex_wp, 9>::Solve(Vector<Complex_wp>&);
  
    
  /* For Piano */
  
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 4>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 6>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyBandMatrix<Real_wp, 8>&);
  
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 1, 1>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 3, 1>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 4, 1>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 6, 1>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 8, 1>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 10, 1>;

  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 1, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 3, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 4, 1>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 6, 1>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 8, 1>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 10, 1>&);

  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 1, 2>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 2, 2>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 3, 2>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 4, 2>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 6, 2>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 8, 2>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 10, 2>;

  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 1, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 2, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 3, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 4, 2>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 6, 2>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 8, 2>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 10, 2>&);

  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 1, 3>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 2, 3>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 3, 3>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 4, 3>;
  SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 6, 3>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 8, 3>;
  //SELDON_EXTERN template class TinyArrowMatrix<Real_wp, 10, 3>;

  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 1, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 2, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 3, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 4, 3>&);
  SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 6, 3>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 8, 3>&);
  //SELDON_EXTERN template void Copy(const Matrix<Real_wp, General, ArrayRowSparse>&, TinyArrowMatrix<Real_wp, 10, 3>&);
  
}
