#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Quadrature/SingularIntegration3D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension2>;
  
  SELDON_EXTERN template void SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension2>::
  ComputeElemMatrix(const VectR3&, const VectR3&, bool, const VectReal_wp&, const VectReal_wp&,
		    const VectR3&, const VectR3&, const VectR3&, const VectR3&,
		    const VectR3&, const VectR3&, const VectR3&, const VectR3&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
		    const EuclidianDistanceClass_Base<Dimension3>&, Matrix<Real_wp>&,
		    Matrix<Real_wp>&, Matrix<Real_wp>&);
  
  SELDON_EXTERN template void SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension2>::
  ComputeStiffnessMatrix(const VectR3&, bool, const VectReal_wp&, const VectR3&,
			 const Vector<TinyMatrix<Real_wp, General, 3, 2> >&,
			 const Complex_wp&, const Real_wp&, Matrix<Complex_wp>& mat_elem);
  
}
