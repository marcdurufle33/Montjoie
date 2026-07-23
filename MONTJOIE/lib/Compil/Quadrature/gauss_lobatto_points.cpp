#include "MontjoieFlag.hxx"

#include "Quadrature/MontjoieQuadratureHeader.hxx"
#include "Quadrature/MontjoieQuadratureInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/GaussLobattoPoints.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template void 
  ComputeFourierCoef(const int, Real_wp&, Vector<Real_wp>&,
		     Vector<Real_wp>&, Vector<Real_wp>&);
  
  SELDON_EXTERN template void
  ComputeLegendrePol_and_Derivative(int, Real_wp, const Real_wp&, const Vector<Real_wp>&,
				    const Vector<Real_wp>&, const Vector<Real_wp>&,
				    Real_wp&, Real_wp&, Real_wp&);
  
  SELDON_EXTERN template void
  ComputeGaussFormulas(int, Vector<Real_wp>&, Vector<Real_wp>&, bool);

  SELDON_EXTERN template void
  ComputeGaussBlendedFormulas(int nb_points, Vector<Real_wp>& points,
			      Vector<Real_wp>& weights, const Real_wp& tau);  
  
}

