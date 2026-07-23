#include "MontjoieFlag.hxx"

#include "Algebra/MontjoieAlgebraHeader.hxx"
#include "Algebra/MontjoieAlgebraInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/TinyBandMatrix.cxx"
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
  
}
