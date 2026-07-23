#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/SingularIntegration2D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension1>;
  
  SELDON_EXTERN template void SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension1>::
  ComputeElemMatrix(const VectR2&, const VectR2&, bool, const VectReal_wp&, const VectReal_wp&,
		    const VectR2&, const VectR2&, const VectR2&, const VectR2&,
		    const VectR2&, const VectR2&, const VectR2&, const VectR2&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const EuclidianDistanceClass_Base<Dimension2>&, Matrix<Real_wp>&,
		    Matrix<Real_wp>&, Matrix<Real_wp>&);
  
  SELDON_EXTERN template void SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension1>::
  ComputeStiffnessMatrix(const VectR2&, bool, const VectReal_wp&, const VectR2&,
			 const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
			 const Complex_wp&, const Real_wp&, Matrix<Complex_wp>& mat_elem);
  
}
