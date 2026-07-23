#include "MontjoieFlag.hxx"

#include "Share/MontjoieCommonHeader.hxx"
#include "Share/MontjoieCommonInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Share/UnivariatePolynomial.cxx"
#include "Share/MultivariatePolynomial.cxx"
#endif

namespace Montjoie
{
  
  // univariate polynomials
  SELDON_EXTERN template class UnivariatePolynomial<Real_wp>;
  SELDON_EXTERN template class UnivariatePolynomial<Complex_wp>;
  SELDON_EXTERN template ostream& operator <<(ostream& out, const UnivariatePolynomial<Real_wp>&);
  SELDON_EXTERN template ostream& operator <<(ostream& out, const UnivariatePolynomial<Complex_wp>&);
  
  SELDON_EXTERN template void MltAdd(const Real_wp&, const UnivariatePolynomial<Real_wp>&,
				     const UnivariatePolynomial<Real_wp>&, const Real_wp&, 
				     UnivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void MltAdd(const Complex_wp&, const UnivariatePolynomial<Complex_wp>&,
                                     const UnivariatePolynomial<Complex_wp>&, const Complex_wp&, 
                                     UnivariatePolynomial<Complex_wp>&);

  SELDON_EXTERN template UnivariatePolynomial<Real_wp> Pow(const UnivariatePolynomial<Real_wp>&, int);
  SELDON_EXTERN template void ComputeLegendrePolynome(int, int, Vector<UnivariatePolynomial<Real_wp> >&);

  // multivariate polynomials
  SELDON_EXTERN template class MultivariatePolynomial<Real_wp>;
  SELDON_EXTERN template Real_wp MultivariatePolynomial<Real_wp>::Evaluate(const TinyVector<Real_wp, 2>&) const;
  SELDON_EXTERN template Real_wp MultivariatePolynomial<Real_wp>::Evaluate(const TinyVector<Real_wp, 3>&) const;
  
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> operator *(const Real_wp& alpha, const MultivariatePolynomial<Real_wp>& P);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> operator +(const Real_wp& alpha, const MultivariatePolynomial<Real_wp>& P);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> operator -(const Real_wp& alpha, const MultivariatePolynomial<Real_wp>& P);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> operator -(const MultivariatePolynomial<Real_wp>& P);

  SELDON_EXTERN template MultivariatePolynomial<Real_wp> MultivariatePolynomial<Real_wp>::operator *(const Real_wp&);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> MultivariatePolynomial<Real_wp>::operator -(const Real_wp&);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> MultivariatePolynomial<Real_wp>::operator +(const Real_wp&);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp>& MultivariatePolynomial<Real_wp>::operator *=(const Real_wp&);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp>& MultivariatePolynomial<Real_wp>::operator -=(const Real_wp&);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp>& MultivariatePolynomial<Real_wp>::operator +=(const Real_wp&);

  SELDON_EXTERN template MultivariatePolynomial<Real_wp> Pow(const MultivariatePolynomial<Real_wp>&, int);
  SELDON_EXTERN template MultivariatePolynomial<Real_wp> Derivate(const MultivariatePolynomial<Real_wp>&, int);

  SELDON_EXTERN template void GetCurlPolynomial(const TinyVector<MultivariatePolynomial<Real_wp>, 3>&, TinyVector<MultivariatePolynomial<Real_wp>, 3>&);

  SELDON_EXTERN template void GetDivPolynomial(const TinyVector<MultivariatePolynomial<Real_wp>, 3>&, MultivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void DerivatePolynomial(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, int);
  SELDON_EXTERN template void IntegratePolynomial(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, int);
  
  SELDON_EXTERN template void Mlt(const Real_wp&, MultivariatePolynomial<Real_wp>&);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const MultivariatePolynomial<Real_wp>&, const MultivariatePolynomial<Real_wp>&, const Real_wp&, MultivariatePolynomial<Real_wp>&);
  
  SELDON_EXTERN template Real_wp IntegrateOnSimplex(const MultivariatePolynomial<Real_wp>&);
  SELDON_EXTERN template void Add(const Real_wp&, const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&);
  
  SELDON_EXTERN template void ScalePolynomVariable(MultivariatePolynomial<Real_wp>&, const Real_wp&);
  SELDON_EXTERN template void CompressPolynom(MultivariatePolynomial<Real_wp>&, const Real_wp&);

  SELDON_EXTERN template void OppositeX(MultivariatePolynomial<Real_wp>&);
  SELDON_EXTERN template void OppositeY(MultivariatePolynomial<Real_wp>&);
  SELDON_EXTERN template void PermuteXY(MultivariatePolynomial<Real_wp>&);
  
  SELDON_EXTERN template void GenerateSymPol(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void GenerateSym4_Vertices(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void GenerateSym4_Edges(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void GenerateSym8(const MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&, MultivariatePolynomial<Real_wp>&);

  SELDON_EXTERN template void SolvePolynomialEquation(const UnivariatePolynomial<Real_wp>&, Vector<Real_wp>&, Vector<Real_wp>&);

    SELDON_EXTERN template void SolvePolynomialEquation(const UnivariatePolynomial<Complex_wp>&, Vector<Complex_wp>&);
  
  // Legendre polynomials
  SELDON_EXTERN template class LegendrePolynomial<Real_wp>;
  SELDON_EXTERN template class AssociatedLegendrePolynomial<Real_wp>;
  
}
