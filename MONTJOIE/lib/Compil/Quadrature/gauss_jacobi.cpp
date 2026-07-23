#include "MontjoieFlag.hxx"

#include "Quadrature/MontjoieQuadratureHeader.hxx"
#include "Quadrature/MontjoieQuadratureInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/GaussJacobi.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void GetJacobiPolynomial(Matrix<Real_wp>&, int, const Real_wp&, const Real_wp&);
  
  SELDON_EXTERN template void EvaluateJacobiPolynomial(const Matrix<Real_wp>&, int, const Real_wp&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void EvaluateJacobiPolynomial(const Matrix<Real_wp>&, int, const Real_wp&, Vector<Real_wp>&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void ComputeGaussLegendre(Vector<Real_wp>&, Vector<Real_wp>&, int);
  
  SELDON_EXTERN template void ComputeGaussLobatto(Vector<Real_wp>&, Vector<Real_wp>&, int);
  
  SELDON_EXTERN template void ComputeGaussJacobi(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&, const Real_wp&);

  SELDON_EXTERN template void ComputeGaussRadauJacobi(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&, const Real_wp&, bool);

  SELDON_EXTERN template void ComputeGaussLobattoJacobi(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&, const Real_wp&);

  SELDON_EXTERN template void ComputeGaussChebyshev(Vector<Real_wp>&, Vector<Real_wp>&, int);

  SELDON_EXTERN template void ComputeGaussLogarithmic(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&);

  SELDON_EXTERN template void ComputeGaussLaguerre(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&);

  SELDON_EXTERN template void ComputeGaussHermite(Vector<Real_wp>&, Vector<Real_wp>&, int, const Real_wp&);
  
}

